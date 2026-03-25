# ============================================================================
# MAP_WORKSPACE.R
# Leaflet-Karte des AFS-Workspaces:
#   - Feldblock-Polygone (KUP-Potenzialflächen)
#   - Feldblock-Zentroide (Produktionsstandorte)
#   - Vorbehandlungsstandorte / Umschlagsplätze (Storages)
#   - Verbrauchspunkte (Consumers: Zellstoff, Papier, Chemie, Energie)
#
# Voraussetzung:
#   source("build_workspace.R")  # oder: load("afs_workspace.RData")
#
# Verwendung:
#   source("map_workspace.R")
#   afs_map   # interaktive Leaflet-Karte
#
# Konsistent mit app_v10.r:
#   Farben/Radius/Popups analog zum Shiny-Map-Output
#
# Paketabhängigkeiten:
#   leaflet          -- CRAN (stabil)
#   leaflet.extras2  -- CRAN (stabil, Nachfolger von leaflet.extras)
#   sf, dplyr, htmltools -- CRAN
#
# HINWEIS: leaflet.extras wurde im August 2024 von CRAN entfernt.
#   Workaround A (empfohlen): leaflet.extras2 aus CRAN
#     install.packages("leaflet.extras2")
#   Workaround B: Direkt von GitHub
#     remotes::install_github("trafficonese/leaflet.extras")
# ============================================================================

suppressPackageStartupMessages({
  library(leaflet)
  library(sf)
  library(dplyr)
  library(htmltools)
})

# Optionales Laden von leaflet.extras2 (CRAN) oder leaflet.extras (GitHub)
# -- leaflet.extras2 ist für diese Karte ausreichend (addSearchFeatures etc.)
if (requireNamespace("leaflet.extras2", quietly = TRUE)) {
  library(leaflet.extras2)
} else if (requireNamespace("leaflet.extras", quietly = TRUE)) {
  library(leaflet.extras)
} else {
  message(
    "[map_workspace] Hinweis: weder leaflet.extras2 noch leaflet.extras verfügbar.\n",
    "  Karte wird ohne Suchfunktion gerendert. Installieren mit:\n",
    "  install.packages('leaflet.extras2')  # empfohlen (CRAN)\n",
    "  # ODER: remotes::install_github('trafficonese/leaflet.extras')"
  )
}

# ============================================================================
# 0. DATEN LADEN (falls noch nicht im Environment)
# ============================================================================

if (!exists("afs_workspace")) {
  if (file.exists("afs_workspace.RData")) {
    cat("Lade afs_workspace.RData ...\n")
    load("afs_workspace.RData")
  } else {
    cat("afs_workspace nicht gefunden - starte build_workspace.R ...\n")
    source("build_workspace.R")
  }
}

# Kurzreferenzen
sites      <- afs_workspace$sites
sites_sf   <- afs_workspace$sites_sf
storages   <- afs_workspace$storages
consumers  <- afs_workspace$consumers
meta       <- afs_workspace$meta

# ============================================================================
# 1. FARBKODIERUNG (analog app_v10.r + erweitert für Consumer-Kategorien)
# ============================================================================

COL_SITES_BORDER  <- "#2d6a2d"   # darkgreen
COL_SITES_FILL    <- "#74c476"   # lightgreen
COL_STORAGE       <- "#e65c00"   # darkorange
COL_STORAGE_FILL  <- "#fd8d3c"   # orange
COL_CONSUMER_P1   <- "#6a0dad"   # darkviolet  = Chemisch/Zellstoff
COL_CONSUMER_P2   <- "#08519c"   # darkblue    = Pulp/Sägewerk
COL_CONSUMER_P3   <- "#a50026"   # darkred     = Energie/Biogas

# Consumer-Farbe nach dominantem Produkttyp
get_consumer_color <- function(p1, p2, p3) {
  dplyr::case_when(
    p1 >= p2 & p1 >= p3 & p1 > 0 ~ COL_CONSUMER_P1,
    p2 >= p1 & p2 >= p3 & p2 > 0 ~ COL_CONSUMER_P2,
    p3 > 0                        ~ COL_CONSUMER_P3,
    TRUE                          ~ "grey60"
  )
}

consumers <- consumers %>%
  mutate(
    color        = get_consumer_color(demand_P1, demand_P2, demand_P3),
    total_demand = demand_P1 + demand_P2 + demand_P3,
    radius       = pmin(18, pmax(5, 5 + log1p(total_demand) * 1.8)),
    kategorie    = dplyr::case_when(
      demand_P1 >= demand_P2 & demand_P1 >= demand_P3 & demand_P1 > 0 ~ "P1 Chemisch/Zellstoff",
      demand_P2 >= demand_P1 & demand_P2 >= demand_P3 & demand_P2 > 0 ~ "P2 Pulp/Sägewerk",
      demand_P3 > 0 ~ "P3 Energie/Biogas",
      TRUE ~ "Unbekannt"
    )
  )

storages <- storages %>%
  mutate(
    color  = dplyr::if_else(type == "Hub", "#c05000", COL_STORAGE),
    fill   = dplyr::if_else(type == "Hub", "#fdae6b", COL_STORAGE_FILL),
    radius = pmin(14, pmax(5, 4 + log1p(CAP_proc) * 1.2))
  )

# ============================================================================
# 2. POPUP-TEXTE
# ============================================================================

site_polygon_popup <- paste0(
  "<b>Feldblock</b><br>",
  "Fläche: ", round(sites_sf$area, 1), " ha"
)

site_point_popup <- paste0(
  "<b>", sites$name, "</b><br>",
  "Fläche: ", round(sites$area_ha, 1), " ha"
)

storage_popup <- paste0(
  "<b>", storages$name, "</b><br>",
  "Typ: ", storages$type, "<br>",
  "Lagerkapazität: ",       storages$CAP_stor, " t<br>",
  "Verarbeitungskapazität: ", storages$CAP_proc, " kt/Jahr"
)

consumer_popup <- paste0(
  "<b>", consumers$name, "</b><br>",
  "Kategorie: ",    consumers$kategorie, "<br>",
  "P1 Chemisch: ",  consumers$demand_P1, " kt/Jahr (", consumers$P1, " €/t)<br>",
  "P2 Pulp/Säge: ", consumers$demand_P2, " kt/Jahr (", consumers$P2, " €/t)<br>",
  "P3 Energie: ",   consumers$demand_P3, " kt/Jahr (", consumers$P3, " €/t)<br>",
  "<b>Gesamt: ",    round(consumers$total_demand, 1), " kt/Jahr</b>"
)

# ============================================================================
# 3. LEAFLET-KARTE ERSTELLEN
# ============================================================================

cat("Erstelle Leaflet-Karte ...\n")
cat("  Feldblöcke: ", nrow(sites_sf),  "Polygone\n")
cat("  Zentroide:  ", nrow(sites),     "Punkte\n")
cat("  Storages:   ", nrow(storages),  "Standorte\n")
cat("  Consumers:  ", nrow(consumers), "Verbrauchsorte\n")

afs_map <- leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%

  # --- Basiskarten ---
  addProviderTiles("CartoDB.Positron",  group = "CartoDB (hell)") %>%
  addProviderTiles("OpenStreetMap.DE",  group = "OpenStreetMap DE") %>%
  addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%

  setView(lng = 11.85, lat = 51.55, zoom = 8) %>%
  addScaleBar(position = "bottomleft") %>%

  # -------------------------------------------------------
  # LAYER 1: Feldblock-Polygone
  # -------------------------------------------------------
  addPolygons(
    data             = sites_sf,
    color            = COL_SITES_BORDER,
    weight           = 0.5,
    fillColor        = COL_SITES_FILL,
    fillOpacity      = 0.35,
    popup            = site_polygon_popup,
    label            = ~paste0(round(area, 0), " ha"),
    labelOptions     = labelOptions(noHide = FALSE, textsize = "11px"),
    group            = "Feldblöcke (Polygone)",
    highlightOptions = highlightOptions(
      weight = 2, color = "#005a00", fillOpacity = 0.6, bringToFront = FALSE
    )
  ) %>%

  # -------------------------------------------------------
  # LAYER 2: Zentroide
  # -------------------------------------------------------
  addCircleMarkers(
    data        = sites,
    lng         = ~lng, lat = ~lat,
    radius      = 3,
    color       = COL_SITES_BORDER,
    fillColor   = COL_SITES_FILL,
    weight      = 1.5, opacity = 0.9, fillOpacity = 0.7,
    popup       = site_point_popup,
    group       = "Produktionsstandorte (Zentroide)"
  ) %>%

  # -------------------------------------------------------
  # LAYER 3: Vorbehandlung / Umschlag
  # -------------------------------------------------------
  addCircleMarkers(
    data         = storages,
    lng          = ~lng, lat = ~lat,
    radius       = ~radius,
    color        = ~color, fillColor = ~fill,
    weight       = 2.5, opacity = 1.0, fillOpacity = 0.85,
    popup        = storage_popup,
    label        = ~name,
    labelOptions = labelOptions(noHide = FALSE, textsize = "12px", direction = "top"),
    group        = "Vorbehandlung / Umschlag"
  ) %>%

  # -------------------------------------------------------
  # LAYER 4a: P1 Chemisch/Zellstoff
  # -------------------------------------------------------
  addCircleMarkers(
    data        = consumers %>% filter(kategorie == "P1 Chemisch/Zellstoff"),
    lng = ~lng, lat = ~lat, radius = ~radius,
    color = COL_CONSUMER_P1, fillColor = COL_CONSUMER_P1,
    weight = 2.5, opacity = 1.0, fillOpacity = 0.75,
    popup = ~paste0(
      "<b>", name, "</b><br>Kategorie: ", kategorie, "<br>",
      "P1: ", demand_P1, " kt/Jahr (", P1, " €/t)<br>",
      "<b>Gesamt: ", round(total_demand, 1), " kt/Jahr</b>"
    ),
    label = ~name,
    labelOptions = labelOptions(noHide = FALSE, textsize = "12px", direction = "top"),
    group = "P1 Chemisch / Zellstoff"
  ) %>%

  # -------------------------------------------------------
  # LAYER 4b: P2 Pulp/Sägewerk
  # -------------------------------------------------------
  addCircleMarkers(
    data        = consumers %>% filter(kategorie == "P2 Pulp/Sägewerk"),
    lng = ~lng, lat = ~lat, radius = ~radius,
    color = COL_CONSUMER_P2, fillColor = COL_CONSUMER_P2,
    weight = 2.5, opacity = 1.0, fillOpacity = 0.75,
    popup = ~paste0(
      "<b>", name, "</b><br>Kategorie: ", kategorie, "<br>",
      "P2: ", demand_P2, " kt/Jahr (", P2, " €/t)<br>",
      "<b>Gesamt: ", round(total_demand, 1), " kt/Jahr</b>"
    ),
    label = ~name,
    labelOptions = labelOptions(noHide = FALSE, textsize = "12px", direction = "top"),
    group = "P2 Pulp / Sägewerk"
  ) %>%

  # -------------------------------------------------------
  # LAYER 4c: P3 Energie/Biogas
  # -------------------------------------------------------
  addCircleMarkers(
    data        = consumers %>% filter(kategorie == "P3 Energie/Biogas"),
    lng = ~lng, lat = ~lat, radius = ~radius,
    color = COL_CONSUMER_P3, fillColor = COL_CONSUMER_P3,
    weight = 2.5, opacity = 1.0, fillOpacity = 0.75,
    popup = ~paste0(
      "<b>", name, "</b><br>Kategorie: ", kategorie, "<br>",
      "P3: ", round(demand_P3, 1), " kt/Jahr (", P3, " €/t)<br>",
      "<b>Gesamt: ", round(total_demand, 1), " kt/Jahr</b>"
    ),
    label = ~name,
    labelOptions = labelOptions(noHide = FALSE, textsize = "12px", direction = "top"),
    group = "P3 Energie / Biogas-BHKW"
  ) %>%

  # -------------------------------------------------------
  # LEGENDE
  # -------------------------------------------------------
  addLegend(
    position = "bottomright",
    colors   = c(COL_SITES_FILL, COL_STORAGE, "#fdae6b",
                 COL_CONSUMER_P1, COL_CONSUMER_P2, COL_CONSUMER_P3),
    labels   = c(
      "Feldblöcke / KUP-Potenzialflächen",
      "Vorbehandlung (Sägewerk)",
      "Vorbehandlung (Logistik-Hub)",
      "P1 Chemisch / Zellstoff",
      "P2 Pulp / Sägewerk",
      "P3 Energie / Biogas-BHKW"
    ),
    opacity = 0.85,
    title   = HTML("<b>AFS Supply Chain</b><br><small>Süd-Sachsen-Anhalt</small>")
  ) %>%

  # -------------------------------------------------------
  # LAYER-CONTROL
  # -------------------------------------------------------
  addLayersControl(
    baseGroups    = c("CartoDB (hell)", "OpenStreetMap DE", "Satellite"),
    overlayGroups = c(
      "Feldblöcke (Polygone)",
      "Produktionsstandorte (Zentroide)",
      "Vorbehandlung / Umschlag",
      "P1 Chemisch / Zellstoff",
      "P2 Pulp / Sägewerk",
      "P3 Energie / Biogas-BHKW"
    ),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%

  # Polygone standardmäßig ausgeblendet (Performance)
  hideGroup("Feldblöcke (Polygone)")

# ============================================================================
# 4. ZUSAMMENFASSUNG
# ============================================================================

cat("\n=== KARTENINHALT ===\n")
cat("Region:        ", meta$region, "\n")
cat("Feldblöcke:    ", meta$n_sites, "Flächen |",
    format(meta$area_total_ha, big.mark = "."), "ha Gesamt\n")
cat("Vorbehandlung: ", meta$n_storages, "Standorte\n")
cat("Consumer:      ", meta$n_consumers, "Verbrauchspunkte\n")
cat("  P1 Chemisch: ", meta$demand_P1_kt, "kt/Jahr\n")
cat("  P2 Pulp:     ", meta$demand_P2_kt, "kt/Jahr\n")
cat("  P3 Energie:  ", meta$demand_P3_kt, "kt/Jahr\n")
cat("\nObjekt: 'afs_map' (interaktive Leaflet-Karte)\n")
cat("Aufruf:  afs_map\n")
cat("Export:  htmlwidgets::saveWidget(afs_map, 'afs_map.html')\n")

afs_map
