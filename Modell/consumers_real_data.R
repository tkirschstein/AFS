# ============================================================================
# REAL CONSUMER DATA: Holzverarbeitende Betriebe
# Untersuchungsgebiet: Südliches Sachsen-Anhalt, Nordthüringen, Nordsachsen
#
# Quellen:
#   - Mercer Stendal (Arneburg): https://renewable-carbon.eu/news/sachsen-anhalt-modernstes-zellstoffwerk-europas-soll-wirtschaft-ankurbeln/
#   - UPM Biochemicals (Leuna): https://www.upmbiochemicals.com/de/
#   - Mercer Torgau (Sägewerk + Pellets): https://www.frankenbrennstoffe.de/marken/mercer-torgau-gmbh-co.-kg/
#   - Mercer Rosenthal (Blankenstein): https://de.wikipedia.org/wiki/Zellstoff-_und_Papierfabrik_Rosenthal
#   - Sägewerk Nedlitz: https://www.saegewerk-nedlitz.de
#   - Sägewerk Selle (Holzthaleben): https://saegewerk-selle.de
#
# Spalten (entsprechend app_v10.r / generate_instance_v8_final):
#   consumer_id  : eindeutige ID
#   name         : Unternehmensname
#   lat / lng    : WGS84-Koordinaten (Werksstandort)
#   demand_P1    : Jahresbedarf Produkt 1 – Chemisch/Zellstoff (kt/Jahr)
#   demand_P2    : Jahresbedarf Produkt 2 – Papierholz/Pulp (kt/Jahr)
#   demand_P3    : Jahresbedarf Produkt 3 – Energieholz/Pellets (kt/Jahr)
#   P1 / P2 / P3 : Erlöspreise je Produkt (€/t)
#
# HINWEIS zu Kapazitäten:
#   Alle Kapazitäten sind als kt Holz-Input angegeben (nicht Output).
#   Mercer Stendal verarbeitet ~3 Mio. Festmeter ≈ 1.500 kt/Jahr Nadelholz.
#   UPM Biochemicals benötigt ca. 500 kt/Jahr Buchenholz (Laubholz → P1).
#   Mercer Rosenthal: ~1.000 kt/Jahr Nadelholz.
#   Mercer Torgau: Sägewerk ~300 kt/Jahr + Pellets ~50 kt/Jahr.
#   Kleinere Sägewerke: 10–50 kt/Jahr.
#
# HINWEIS zu Preisen (€/t Holzrohstoff, Einkaufspreise am Werkstor):
#   P1 Chemisch/Zellstoff: 80–120 €/t (Faserstoff-Qualität)
#   P2 Pulp/Sägerundholz:  60–90 €/t  (Langholz, sägefähig)
#   P3 Energieholz/Pellets: 30–50 €/t (Hackgut, Restholz)
# ============================================================================

library(sf)
library(dplyr)

consumers_real <- data.frame(
  consumer_id = 1:9,
  name = c(
    "Mercer Stendal GmbH (Arneburg)",          # 1 – Großzellstoff Nadel
    "UPM Biochemicals (Leuna)",                 # 2 – Bioraffinerie Laub
    "Mercer Torgau – Sägewerk & Pellets",       # 3 – Sägewerk + Pellets
    "Zellstoff Rosenthal / Mercer (Blankenstein)", # 4 – Zellstoff Nadel
    "Sägewerk Nedlitz GmbH",                   # 5 – kl. Sägewerk
    "Sägewerk Selle (Holzthaleben)",            # 6 – kl. Sägewerk Thüringen
    "Gebrüder Machemehl (Gernrode/Harz)",      # 7 – Sägewerk + Zimmerei
    "Holzhof Helbra (Mansfeld-Südharz)",        # 8 – Holzhandel/Verarbeitung
    "IMHOLZ / Klöpfer Holzhandel (Leipzig)"    # 9 – Holzgroßhandel
  ),
  lat = c(
    52.7817,   # Arneburg, Stendal
    51.3314,   # Leuna, Saalekreis
    51.5633,   # Torgau, Nordsachsen
    50.5014,   # Blankenstein/Rennsteig, Thüringen
    51.9650,   # Nedlitz b. Burg, S-Anhalt
    51.3200,   # Holzthaleben, Kyffhäuserkreis
    51.7240,   # Gernrode, Harz
    51.5278,   # Helbra, Mansfeld-Südharz
    51.3397    # Leipzig-Plagwitz
  ),
  lng = c(
    11.9417,   # Arneburg
    12.0017,   # Leuna
    12.9883,   # Torgau
    11.6364,   # Blankenstein
    11.8700,   # Nedlitz
    10.9200,   # Holzthaleben
    11.1410,   # Gernrode
    11.5294,   # Helbra
    12.3331    # Leipzig
  ),
  # Nachfrage in kt/Jahr (Holzinput-Äquivalent)
  demand_P1 = c(1500, 500,   0,  1000,   0,  0,  0,  0,  10),  # Chemisch/Zellstoff
  demand_P2 = c(   0,   0, 300,     0,  30, 15, 10,  5,  50),  # Pulp/Sägerundholz
  demand_P3 = c(   0,   0,  50,     0,   5,  5,  3,  8,  20),  # Energieholz/Pellets
  # Einkaufspreise (€/t Holzinput am Werkstor)
  P1 = c(100,  120,   0,   95,   0,   0,   0,   0,  80),
  P2 = c(  0,    0,  75,    0,  70,  65,  60,  55,  70),
  P3 = c(  0,    0,  40,    0,  35,  35,  30,  30,  38),
  stringsAsFactors = FALSE
)

# Als sf-Objekt (für Kartenintegration, analog zu biochem_sites_sf in app_v10.r)
biochem_sites_real_sf <- st_as_sf(
  consumers_real,
  coords = c("lng", "lat"),
  crs    = 4326
)

# Zusammenfassung ausgeben
cat("=== Reale Consumer-Daten geladen ===\n")
cat("Anzahl Betriebe:", nrow(consumers_real), "\n\n")
print(consumers_real[, c("consumer_id", "name", "demand_P1", "demand_P2", "demand_P3")])
cat("\nGesamt-Nachfrage P1 (Chemisch/Zellstoff): ", sum(consumers_real$demand_P1), "kt/Jahr\n")
cat("Gesamt-Nachfrage P2 (Pulp/Sägerundholz):  ", sum(consumers_real$demand_P2), "kt/Jahr\n")
cat("Gesamt-Nachfrage P3 (Energie/Pellets):    ", sum(consumers_real$demand_P3), "kt/Jahr\n")

# ============================================================================
# VERWENDUNG IN app_v10.r:
#
#   source("consumers_real_data.R")
#
#   # In observeEvent(base_data(), ...) ersetzen:
#   user_data$consumers <- consumers_real
#
# VERWENDUNG IN paper.Rmd:
#
#   source("../Modell/consumers_real_data.R")
#   # consumers_real und biochem_sites_real_sf stehen dann zur Verfügung
# ============================================================================
