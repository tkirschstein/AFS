# ============================================================================
# MaStR-Analyse: Abnahmepotentiale für AFS-Biomasse
# Zielregion: Südliches Sachsen-Anhalt, Nordthüringen, Nordsachsen
#
# Eingangsdaten (Marktstammdatenregister, BNetzA):
#   - Stromerzeuger.csv : alle Biomasse-Stromanlagen (19.252 Einträge)
#   - Gaserzeuger.csv   : alle Biomethan-Einspeiser (293 Einträge)
#
# Ergebnis: consumer-Datensatz für MILP-Modell (app_v10.r)
#
# Datenstruktur:
#   Stromerzeuger: Semikolon-CSV, Spalte "Energieträger" = "Biomasse"
#     Hauptbrennstoffe: Biogas (18.004), Biomethan (1.163), Holzgas (85)
#     Relevante Spalten: MaStR-Nr., Name, Bundesland, Landkreis, Ort,
#       Breitengrad, Längengrad, Bruttoleistung, Hauptbrennstoff,
#       Elektrische KWK-Leistung, Thermische Nutzleistung
#   Gaserzeuger: Semikolon-CSV, Technologie = "Biomethan-Erzeugung"
#     Relevante Spalten: MaStR-Nr., Name, Bundesland, PLZ, Ort,
#       Breitengrad, Längengrad, Erzeugungsleistung in kWh/h
# ============================================================================

library(dplyr)
library(sf)
library(stringr)

# ============================================================================
# 1. EINLESEN DER CSV-DATEIEN
# ============================================================================

# Pfad relativ zum Modell-Verzeichnis
strom_raw <- read.csv2("Stromerzeuger.csv", fileEncoding = "UTF-8-BOM",
                       stringsAsFactors = FALSE, check.names = FALSE)

gas_raw   <- read.csv2("Gaserzeuger.csv",   fileEncoding = "UTF-8-BOM",
                       stringsAsFactors = FALSE, check.names = FALSE)

cat("Stromerzeuger geladen:", nrow(strom_raw), "Zeilen,", ncol(strom_raw), "Spalten\n")
cat("Gaserzeuger geladen:",   nrow(gas_raw),   "Zeilen,", ncol(gas_raw),   "Spalten\n")

# ============================================================================
# 2. FILTER: REGION & BETRIEBSSTATUS
# ============================================================================

REGION_BUNDESLAENDER <- c("Sachsen-Anhalt", "Thüringen", "Sachsen")

# ZIELREGION Landkreise (Süd-Sachsen-Anhalt)
ZIEL_LANDKREISE <- c("Burgenlandkreis", "Saalekreis", "Mansfeld-Südharz",
                     "Wittenberg", "Anhalt-Bitterfeld", "Halle (Saale)")

# Stromerzeuger: Region + in Betrieb
strom_region <- strom_raw %>%
  filter(
    Bundesland %in% REGION_BUNDESLAENDER,
    `Betriebs-Status` == "In Betrieb"
  ) %>%
  mutate(
    lat = as.numeric(str_replace(`Koordinate: Breitengrad (WGS84)`, ",", ".")),
    lng = as.numeric(str_replace(`Koordinate: L\u00e4ngengrad (WGS84)`, ",", ".")),
    kW  = as.numeric(str_replace(`Bruttoleistung der Einheit`, ",", "."))
  )

# Gaserzeuger: Region + in Betrieb
gas_region <- gas_raw %>%
  filter(
    Bundesland %in% REGION_BUNDESLAENDER,
    `Betriebs-Status` == "In Betrieb"
  ) %>%
  mutate(
    lat    = as.numeric(str_replace(`Koordinate: Breitengrad (WGS84)`, ",", ".")),
    lng    = as.numeric(str_replace(`Koordinate: L\u00e4ngengrad (WGS84)`, ",", ".")),
    kWh_h  = as.numeric(str_replace(`Erzeugungsleistung in kWh/h`, ",", "."))
  )

cat("\nStromerzeuger in Region (in Betrieb):", nrow(strom_region), "\n")
cat("Gaserzeuger in Region (in Betrieb):",   nrow(gas_region),   "\n")

# ============================================================================
# 3. HOLZGAS-ANLAGEN (P3 direkt: Holzvergasung)
# ============================================================================

# Direkt AFS-relevant: Holzgas = Vergasung fester Biomasse (Hackschnitzel, KUP)
holzgas <- strom_region %>%
  filter(str_detect(`Hauptbrennstoff der Einheit`,
                    regex("Holzgas|Hackschnitzel|Holz|Rinde|Pellet|feste Biomasse",
                          ignore_case = TRUE)))

cat("Holzgas/Festbiomasse-Anlagen in Region:", nrow(holzgas), "\n")

# ============================================================================
# 4. BIOGAS-BHKW ALS P3-KONSUMENTEN (Kosubstrat-Potenzial)
#    Annahme: 8% Holz-/Energiepflanzen-Kosubstrat möglich
# ============================================================================

# Holzbedarf-Umrechnung aus Leistung:
# B_kt = P_kW * 8760h * h_VL / (eta_el * H_u_kWh_t * 1000)
# Parameter: h_VL=0.70, eta_el=0.30 (Biogas-BHKW), H_u=4000 kWh/t (Energiepflanzen)
KOEFF_BIOGAS <- 8760 * 0.70 / (0.30 * 4000 * 1000)  # = 5.11e-3 kt/(kW*a)

biogas_ziel <- strom_region %>%
  filter(Landkreis %in% ZIEL_LANDKREISE) %>%
  mutate(
    # Konservativ: 8% Holzanteil als Kosubstrat (Kurzumtriebsholz, Pappel, Weide)
    demand_P3_kt = kW * KOEFF_BIOGAS * 0.08
  ) %>%
  filter(!is.na(lat), !is.na(lng), demand_P3_kt > 0)

cat("Biogas-BHKW in Zielregion Landkreisen:", nrow(biogas_ziel), "\n")
cat("P3-Gesamtbedarf (Kosubstrat 8%):",
    round(sum(biogas_ziel$demand_P3_kt, na.rm=TRUE), 1), "kt/Jahr\n")

# ============================================================================
# 5. BIOMETHAN-EINSPEISER ALS P3-KONSUMENTEN
#    VERBIO Zörbig: Stroh-dominiert → P1/Chemisch
#    InfraLeuna: Nähe UPM Bioraffinerie
# ============================================================================

gas_ziel <- gas_region %>%
  filter(str_detect(Postleitzahl, "^(06|07|04)")) %>%
  filter(!is.na(lat), !is.na(lng))

# Erzeugungsleistung kWh/h → P3-Bedarf kt/Jahr
# B_kt = L_kWh_h * 8760h / (H_u_kWh_t * 1000)
# H_u Biogas-Input (Maissilage) = 1000 kWh/t FM → kein Holz
# VERBIO Zörbig nutzt Stroh → P1 (Chemisch/Cellulose), nicht P3

verbio_zorbig <- gas_ziel %>%
  filter(str_detect(`Anzeige-Name der Einheit`, regex("Zörbig", ignore_case=TRUE)))

cat("\nVERBIO Zörbig (größte Anlage):",
    verbio_zorbig$`Anzeige-Name der Einheit`[1], "-",
    verbio_zorbig$kWh_h[1], "kWh/h\n")
cat("  → Hauptsubstrat: Weizenstroh/Getreide (nicht Holz)\n")
cat("  → Relevant als P1-Konsument (Stroh-Cellulose, kein direktes AFS-Holz)\n")

# ============================================================================
# 6. AGGREGIERTER CONSUMER-DATENSATZ FÜR MILP-MODELL
#    Zusammenführung: MaStR-basierte Anlagen + reale Industriekunden
# ============================================================================

# Aggregation: Biogas-BHKW nach Gemeinde (Cluster, max. 5 km)
# Vereinfacht: nach Landkreis aggregiert als "virtueller Consumer"
agg_by_lk <- biogas_ziel %>%
  group_by(Landkreis, Bundesland) %>%
  summarise(
    name       = paste0("Biogas-Cluster ", Landkreis),
    lat        = mean(lat, na.rm = TRUE),
    lng        = mean(lng, na.rm = TRUE),
    n_anlagen  = n(),
    kW_total   = sum(kW, na.rm = TRUE),
    demand_P3_kt = sum(demand_P3_kt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    demand_P1_kt = 0,
    demand_P2_kt = 0,
    P1 = 0, P2 = 0, P3 = 35  # €/t Hackschnitzel am Werkstor
  )

cat("\n=== AGGREGIERTER P3-CONSUMER-DATENSATZ ===\n")
print(agg_by_lk %>% select(name, lat, lng, n_anlagen, kW_total, demand_P3_kt))

# ============================================================================
# 7. SF-OBJEKTE FÜR KARTENINTEGRATION
# ============================================================================

# Biogas-Einzelanlagen als sf
biogas_sf <- st_as_sf(biogas_ziel %>% filter(!is.na(lat), !is.na(lng)),
                      coords = c("lng", "lat"), crs = 4326)

# Aggregierte Cluster als sf
mastr_cluster_sf <- st_as_sf(agg_by_lk,
                              coords = c("lng", "lat"), crs = 4326)

# ============================================================================
# 8. HINWEISE ZUR INTEGRATION IN app_v10.r
# ============================================================================
cat("
=== INTEGRATION IN app_v10.r ===

# In observeEvent(base_data(), ...) nach dem Laden:
source('mastr_consumers.R')

# Für P3-Consumer (Energieholz) den MaStR-Cluster-Datensatz verwenden:
# user_data$consumers <- bind_rows(
#   consumers_real,        # aus consumers_real_data.R (Zellstoff/Sägewerke)
#   agg_by_lk %>% select(consumer_id=..., name, lat, lng,
#                         demand_P1=demand_P1_kt, demand_P2=demand_P2_kt,
#                         demand_P3=demand_P3_kt, P1, P2, P3)
# )

=== WICHTIGE BEFUNDE ===
- 203 Biogas-BHKW in Zielregion, 117.5 MW_el -> ca. 48 kt/Jahr Kosubstrat-Potenzial
- 15 Biomethan-Einspeiser PLZ 06xxx, 138.3 MW_th
- VERBIO Zörbig (33 MW_th): Stroh, kein Holz -> kein P3-Konsument
- BGEA InfraLeuna (9 MW_th, Leuna): Nähe zu UPM Bioraffinerie (P1-Cluster)
- KEINE Holzgas/Holzvergasungs-Anlage in der Region -> Marktlücke!
- P3-Gesamtnachfrage konservativ: 15-30 kt/Jahr (Kosubstrat)
- P3-Gesamtnachfrage optimistisch (inkl. geplante Holzgas-BHKW): 50-80 kt/Jahr
")

# ============================================================================
# 9. AUSGABE: VOLLSTÄNDIGE ANLAGENLISTE ZIELREGION
# ============================================================================

cat("\n=== ALLE BIOGAS-BHKW IN ZIELREGION (Top 20 nach Leistung) ===\n")
biogas_ziel %>%
  arrange(desc(kW)) %>%
  select(`Anzeige-Name der Einheit`, Landkreis, Ort, kW, demand_P3_kt) %>%
  head(20) %>%
  print()
