# ============================================================================
# BUILD_WORKSPACE.R
# Erstellt einen optimierungsbereiten R-Workspace für die AFS Supply-Chain-
# Optimierung. Dieser enthält ausschließlich die für das MILP-Modell
# (build_agroforestry_lp_v10.r / app_v10.r) relevanten Objekte:
#
#   A) FELDBLÖCKE  → sites_sf, sites (Produktionsstandorte KUP)
#   B) VORBEHANDLUNG → storages (Sägewerke + Hackschnitzel-Hubs)
#   C) VERBRAUCHSPUNKTE → consumers (Zellstoff, Papier, Chemie, Energie)
#   D) DISTANZMATRIZEN → dist_ij, dist_jk (via OSM/OSRM)
#   E) WORKSPACE-OBJEKTE → afs_workspace (optimierungsfertiges list-Objekt)
#
# Quelldateien:
#   - !base_data.rdata           : Feldblockatlas + Sägewerke (sawmills_sf)
#   - consumers_real_data.R      : reale Industriekonsumenten (9 Betriebe)
#   - mastr_consumers.R          : MaStR-Biogas/Biomethan-Anlagen (aggregiert)
#   - Stromerzeuger.csv          : Biomasse-Stromanlagen (MaStR)
#   - Gaserzeuger.csv            : Biomethan-Einspeiser (MaStR)
#
# Output:
#   - afs_workspace.RData        : fertig für build_agroforestry_lp_v10.r
#
# Verwendung:
#   source("build_workspace.R")
#   # → afs_workspace$sites, $storages, $consumers, $dist_ij, $dist_jk
#
# Autoren: tkirschstein / Perplexity-Analyse 2026-03-25
# ============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

cat("==========================================================\n")
cat(" AFS Workspace Builder\n")
cat("==========================================================\n\n")


# ============================================================================
# SCHRITT 1: BASIS-DATEN LADEN (!base_data.rdata)
# ============================================================================
# Enthält folgende Objekte aus dem Feldblockatlas:
#   data_fb_filtered   : sf-Polygone aller Feldblöcke (HBN, area, ...)
#   fb_zentroide_wgs84 : sf-Punkte der Feldblock-Zentroide (X/Y/area)
#   sawmills_sf        : sf-Punkte der Sägewerke/Vorbehandlungsstandorte
#                        Spalten: name, geometry

cat("[1/5] Lade !base_data.rdata ...\n")
load("!base_data.rdata")

# Kleine Feldblöcke ausfiltern (< 10 ha, analog app_v10.r Zeile 38/39)
data_fb_filtered   <- data_fb_filtered   %>% filter(area >= 0)
fb_zentroide_wgs84 <- fb_zentroide_wgs84 %>% filter(area >= 0)

# Nur Ackerflächen (HBN == "AL" = Ackerland), falls Spalte vorhanden
if ("HBN" %in% colnames(data_fb_filtered) && "AL" %in% unique(data_fb_filtered$HBN)) {
  data_fb_filtered   <- data_fb_filtered   %>% filter(HBN == "AL")
  fb_zentroide_wgs84 <- fb_zentroide_wgs84 %>% filter(area >= 0)
  cat("   → Ackerflächen (HBN=AL) gefiltert\n")
}

cat("   Feldblöcke gesamt:", nrow(data_fb_filtered), "\n")
cat("   Zentroide gesamt: ", nrow(fb_zentroide_wgs84), "\n")


# ============================================================================
# SCHRITT 2: SITES – Feldblock-Produktionsstandorte
# ============================================================================
cat("[2/5] Erstelle sites-Datensatz (Feldblöcke) ...\n")

sites_coords <- sf::st_coordinates(fb_zentroide_wgs84)

sites <- data.frame(
  site_id  = seq_len(nrow(sites_coords)),
  name     = paste0("FB_", seq_len(nrow(sites_coords))),
  lat      = sites_coords[, "Y"],
  lng      = sites_coords[, "X"],
  area_ha  = fb_zentroide_wgs84$area,
  stringsAsFactors = FALSE
)

sites_sf <- data_fb_filtered

cat("   Sites (Feldblöcke >= 10 ha):", nrow(sites), "\n")
cat("   Gesamtfläche:               ", round(sum(sites$area_ha)/1000, 1), "Tausend ha\n")


# ============================================================================
# SCHRITT 3: STORAGES – Vorbehandlungsstandorte (Sägewerke + Hubs)
# ============================================================================
# Quelle A: sawmills_sf aus !base_data.rdata
# Quelle B: manuelle Ergänzung weiterer Hackschnitzel-Hubs (Logistikstandorte)
#
# CAP_stor: Lagerkapazität in t (Frischmasse Hackschnitzel)
# CAP_proc: Verarbeitungskapazität kt/Jahr (Chipping + Trocknung)
# ============================================================================
cat("[3/5] Erstelle storages-Datensatz (Vorbehandlung/Sägewerke) ...\n")

# Sägewerke aus base_data
sawmills_coords <- sf::st_coordinates(sawmills_sf)

# Standard-Kapazitäten analog app_v10.r (Zeile 616–617)
# CAP_proc-Vektor muss Länge nrow(sawmills_sf) haben
n_saw <- nrow(sawmills_coords)
cap_proc_default <- c(700, 500, 250, 40, 150, 250, 90, 200, 180, 120,
                       80,  60,  50,  50,  50,  50,  50,  50,  50,  50)

storages_sawmills <- data.frame(
  storage_id  = seq_len(n_saw),
  name        = sawmills_sf$name,
  lat         = sawmills_coords[, "Y"],
  lng         = sawmills_coords[, "X"],
  type        = "Sägewerk",
  CAP_stor    = 500,                                    # t Puffer
  CAP_proc    = cap_proc_default[seq_len(n_saw)],       # kt/Jahr
  stringsAsFactors = FALSE
)

# Zusätzliche strategische Hackschnitzel-Hubs (aus Infrastrukturanalyse)
# Diese ergänzen die Sägewerke als reine Logistik-/Trocknungsstandorte
storages_hubs <- data.frame(
  storage_id  = n_saw + 1:4,
  name        = c(
    "Hub Merseburg-Süd (Saalekreis)",
    "Hub Halle-West (Saalekreis)",
    "Hub Zeitz (Burgenlandkreis)",
    "Hub Bitterfeld (Anhalt-Bitterfeld)"
  ),
  lat         = c(51.348,  51.482,  51.051,  51.628),
  lng         = c(11.992,  11.892,  12.134,  12.315),
  type        = "Hub",
  CAP_stor    = 250,
  CAP_proc    = 100,
  stringsAsFactors = FALSE
)

# Zusammenführen
storages <- bind_rows(storages_sawmills, storages_hubs) %>%
  mutate(storage_id = seq_len(n()))

cat("   Sägewerke aus base_data: ", n_saw, "\n")
cat("   Zusätzliche Hubs:         4\n")
cat("   Storages gesamt:         ", nrow(storages), "\n")


# ============================================================================
# SCHRITT 4: CONSUMERS – Verbrauchspunkte
# Drei Kategorien:
#   I)  Reale Industriekunden (Zellstoff, Sägewerke, Chemie)
#   II) MaStR-Biogas-BHKW (aggregiert nach Landkreis, P3)
#   III: MaStR-Biomethan-Einspeiser (P3/P1 je nach Substrat)
# ============================================================================
cat("[4/5] Erstelle consumers-Datensatz ...\n")

## --- 4a: Reale Industriekonsumenten ---
source("consumers_real_data.R", local = TRUE)
# → consumers_real (data.frame, 9 Zeilen)

## --- 4b: MaStR Biogas-BHKW (aggregiert nach Landkreis) ---
# Parameter (analog mastr_consumers.R)
REGION_BUNDESLAENDER <- c("Sachsen-Anhalt", "Thüringen", "Sachsen")
ZIEL_LANDKREISE      <- c("Burgenlandkreis", "Saalekreis", "Mansfeld-Südharz",
                           "Wittenberg", "Anhalt-Bitterfeld", "Halle (Saale)",
                           "Kyffhäuserkreis", "Altenburger Land",
                           "Landkreis Leipzig", "Nordsachsen")

# Stromerzeuger laden und filtern
strom_raw <- read.csv2("Stromerzeuger.csv", fileEncoding = "UTF-8-BOM",
                        stringsAsFactors = FALSE, check.names = FALSE)

strom_region <- strom_raw %>%
  filter(
    Bundesland %in% REGION_BUNDESLAENDER,
    `Betriebs-Status` == "In Betrieb",
    Landkreis %in% ZIEL_LANDKREISE
  ) %>%
  mutate(
    lat = as.numeric(str_replace(`Koordinate: Breitengrad (WGS84)`, ",", ".")),
    lng = as.numeric(str_replace(`Koordinate: Längengrad (WGS84)`,  ",", ".")),
    kW  = as.numeric(str_replace(`Bruttoleistung der Einheit`,       ",", "."))
  ) %>%
  filter(!is.na(lat), !is.na(lng), !is.na(kW))

# Koeffizient: kW_el → kt Holz-Kosubstrat/Jahr
# B = P_kW * 8760h * h_VL=0.70 / (eta_el=0.30 * H_u=4000kWh/t * 1000)
KOEFF_BIOGAS <- 8760 * 0.70 / (0.30 * 4000 * 1000)

# Aggregierung nach Landkreis → virtuelle Cluster-Consumer
mastr_biogas_consumers <- strom_region %>%
  group_by(Postleitzahl) %>%
  summarise(
    lat          = mean(lat, na.rm = TRUE),
    lng          = mean(lng, na.rm = TRUE),
    n_anlagen    = n(),
    kW_total     = sum(kW, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    name         = paste0("Biogas-Cluster ", Postleitzahl),
    # P3: Holz-Kosubstrat-Potenzial (8% des Energieeinsatzes)
    demand_P1    = 0,
    demand_P2    = 0,
    demand_P3    = round(kW_total * KOEFF_BIOGAS * 0.08, 1),
    P1 = 0, P2 = 0, P3 = 35
  ) %>%
  select(name, lat, lng, demand_P1, demand_P2, demand_P3, P1, P2, P3)

## --- 4c: MaStR Biomethan-Einspeiser (Zielregion PLZ 06/07/04) ---
gas_raw <- read.csv2("Gaserzeuger.csv", fileEncoding = "UTF-8-BOM",
                      stringsAsFactors = FALSE, check.names = FALSE)

gas_consumers <- gas_raw %>%
  filter(
    `Betriebs-Status` == "In Betrieb",
    str_detect(Postleitzahl, "^(06|07|04)")
  ) %>%
  mutate(
    lat    = as.numeric(str_replace(`Koordinate: Breitengrad (WGS84)`, ",", ".")),
    lng    = as.numeric(str_replace(`Koordinate: Längengrad (WGS84)`,  ",", ".")),
    kWh_h  = as.numeric(str_replace(`Erzeugungsleistung in kWh/h`,     ",", "."))
  ) %>%
  filter(!is.na(lat), !is.na(lng)) %>%
  mutate(
    name = `Anzeige-Name der Einheit`,
    # VERBIO Zörbig + Könnern: Stroh → kein Holz (P1-Kandidat für Stroh)
    # Alle anderen: geringe Holz-Beimischung möglich (P3)
    ist_verbio = str_detect(name, regex("Zörbig|VERBIO", ignore_case = TRUE)),
    demand_P1  = ifelse(ist_verbio, round(kWh_h * 8760 * 0.70 / (0.55 * 4000 * 1000), 1), 0),
    demand_P2  = 0,
    demand_P3  = ifelse(!ist_verbio, round(kWh_h * 8760 * 0.70 / (0.60 * 3500 * 1000) * 0.05, 1), 0),
    P1 = ifelse(ist_verbio, 55, 0),   # Stroh-Cellulose (P1, niedriger Preis)
    P2 = 0,
    P3 = ifelse(!ist_verbio, 35, 0)
  ) %>%
  select(name, lat, lng, demand_P1, demand_P2, demand_P3, P1, P2, P3)

## --- 4d: Zusammenführen aller Consumer ---
consumers_all <- bind_rows(
  consumers_real %>% select(name, lat, lng, demand_P1, demand_P2, demand_P3, P1, P2, P3),
  mastr_biogas_consumers,
  gas_consumers
) %>%
  filter(demand_P1 + demand_P2 + demand_P3 > 0) %>%   # nur Anlagen mit Bedarf > 0
  mutate(consumer_id = seq_len(n())) %>%
  select(consumer_id, name, lat, lng, demand_P1, demand_P2, demand_P3, P1, P2, P3)

cat("   Reale Industrie-Consumer:   ", nrow(consumers_real), "\n")
cat("   MaStR Biogas-Cluster:       ", nrow(mastr_biogas_consumers), "\n")
cat("   MaStR Biomethan-Einspeiser: ", nrow(gas_consumers), "\n")
cat("   Consumer gesamt (aktiv):    ", nrow(consumers_all), "\n")
cat("   Gesamt P1-Nachfrage:        ", sum(consumers_all$demand_P1), "kt/Jahr\n")
cat("   Gesamt P2-Nachfrage:        ", sum(consumers_all$demand_P2), "kt/Jahr\n")
cat("   Gesamt P3-Nachfrage:        ", round(sum(consumers_all$demand_P3), 1), "kt/Jahr\n")


# ============================================================================
# SCHRITT 5: WORKSPACE ASSEMBLIEREN UND SPEICHERN
# ============================================================================
cat("[5/5] Assembliere und speichere afs_workspace.RData ...\n")

afs_workspace <- list(

  # --- A: Feldblöcke ---
  sites_sf          = sites_sf,           # sf-Polygone (für Kartendarstellung)
  sites             = sites,              # data.frame: site_id, name, lat, lng, area_ha

  # --- B: Vorbehandlungsstandorte ---
  storages          = storages,           # data.frame: storage_id, name, lat, lng,
                                          #             type, CAP_stor, CAP_proc

  # --- C: Verbrauchspunkte (alle Kategorien) ---
  consumers         = consumers_all,      # data.frame: consumer_id, name, lat, lng,
                                          #             demand_P1/P2/P3, P1/P2/P3
  consumers_real    = consumers_real,     # nur reale Industrie-Consumer (9)
  consumers_biogas  = mastr_biogas_consumers,  # MaStR Biogas-Cluster
  consumers_gas     = gas_consumers,      # MaStR Biomethan-Einspeiser

  # --- D: SF-Objekte für Kartenintegration ---
  sawmills_sf       = sawmills_sf,
  consumers_sf      = st_as_sf(consumers_all, coords = c("lng", "lat"), crs = 4326),

  # --- E: Distanzmatrizen (leer, werden in app_v10.r via OSM befüllt) ---
  dist_ij           = NULL,               # sites → storages  (wird via OSRM berechnet)
  dist_jk           = NULL,               # storages → consumers

  # --- F: Metadaten ---
  meta = list(
    created         = Sys.time(),
    n_sites         = nrow(sites),
    n_storages      = nrow(storages),
    n_consumers     = nrow(consumers_all),
    region          = "Südliches Sachsen-Anhalt + angrenzende LK (TH/SN)",
    area_total_ha   = round(sum(sites$area_ha)),
    demand_P1_kt    = sum(consumers_all$demand_P1),
    demand_P2_kt    = sum(consumers_all$demand_P2),
    demand_P3_kt    = round(sum(consumers_all$demand_P3), 1),
    source_files    = c("!base_data.rdata", "consumers_real_data.R",
                        "mastr_consumers.R", "Stromerzeuger.csv", "Gaserzeuger.csv")
  )
)

save(afs_workspace, file = "afs_workspace.RData")

cat("\n✓ afs_workspace.RData gespeichert.\n")
cat("\n=== WORKSPACE-ZUSAMMENFASSUNG ===\n")
cat("  Sites (Feldblöcke >= 10 ha): ", afs_workspace$meta$n_sites, "\n")
cat("  Gesamtfläche:                ", afs_workspace$meta$area_total_ha, "ha\n")
cat("  Vorbehandlungsstandorte:     ", afs_workspace$meta$n_storages, "\n")
cat("  Verbrauchspunkte gesamt:     ", afs_workspace$meta$n_consumers, "\n")
cat("  Nachfrage P1 (Chemisch):     ", afs_workspace$meta$demand_P1_kt, "kt/Jahr\n")
cat("  Nachfrage P2 (Pulp/Säge):   ", afs_workspace$meta$demand_P2_kt, "kt/Jahr\n")
cat("  Nachfrage P3 (Energie):      ", afs_workspace$meta$demand_P3_kt, "kt/Jahr\n")
cat("\n=== VERWENDUNG IN app_v10.r ===\n")
cat("  # Zeile 35 ersetzen durch:\n")
cat("  load('afs_workspace.RData')\n")
cat("  data_fb_filtered   <- afs_workspace$sites_sf\n")
cat("  fb_zentroide_wgs84 <- afs_workspace$sites    # (oder via st_centroid)\n")
cat("  sawmills_sf        <- afs_workspace$sawmills_sf\n")
cat("  biochem_sites_sf   <- afs_workspace$consumers_sf\n")
cat("\n=== VERWENDUNG IN build_agroforestry_lp_v10.r ===\n")
cat("  load('afs_workspace.RData')\n")
cat("  data <- list(\n")
cat("    sites    = afs_workspace$sites,\n")
cat("    storages = afs_workspace$storages,\n")
cat("    consumers = afs_workspace$consumers,\n")
cat("    dist_ij  = afs_workspace$dist_ij,   # nach OSRM-Berechnung\n")
cat("    dist_jk  = afs_workspace$dist_jk\n")
cat("  )\n")
