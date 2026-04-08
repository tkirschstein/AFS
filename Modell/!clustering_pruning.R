# ==============================================================================
# TWO-STAGE SITE CLUSTERING FOR AFS MILP
# Stage 1: DBSCAN — prune isolated / non-viable peripheral sites
# Stage 2: HAC (Ward.D2) with max-radius cut — form spatial super-sites
#
# Input:  afs_workspace$sites  (site_id, lat, lng, area_ha)
# Output: sites_clustered      (cluster_id, lat, lng, area_ha, n_sites)
#         replaces afs_workspace$sites before dist_ij is recomputed
#
# Authors: tkirschstein / SmartAgroforst 2026
# ==============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(dbscan)   # install.packages("dbscan")
  library(cluster)  # for silhouette diagnostics
  library(ggplot2)
  library(patchwork)     # install.packages("patchwork")
  library(leaflet)
  library(leaflet.extras2)
  library(htmltools)
  library(RColorBrewer)
  library(scales)
})

# ── Parameters (adjust to study region) ──────────────────────────────────────
DBSCAN_EPS_KM    <- 5    # neighbourhood radius for DBSCAN (km); sites beyond
# this distance from any neighbour = outlier
DBSCAN_MIN_PTS   <- 3     # minimum cluster size in Stage 1
HAC_MAX_RADIUS   <- 10   # maximum within-cluster geodesic radius (km) for
# Stage 2 cut; all members ≤ this distance to medoid
# ─────────────────────────────────────────────────────────────────────────────

load("../Modell/afs_workspace.RData")

sites_fb <- afs_workspace$sites_sf   # columns: site_id, lat, lng, area_ha

sites <- st_centroid(sites_fb)

centroids_coords <- st_coordinates(sites)
sites$lng <- centroids_coords[, "X"]
sites$lat <- centroids_coords[, "Y"]

sites <- sites %>%
  st_drop_geometry() %>%
  mutate(site_id = seq_len(nrow(sites))) %>%
  select(site_id, lat, lng, area)

cat("Sites before clustering:", nrow(sites), "\n")


# ==============================================================================
# HELPER: great-circle distance matrix (km) from lat/lng data.frame
# ==============================================================================
geodist_km <- function(df) {
  pts <- sf::st_as_sf(df, coords = c("lng", "lat"), crs = 4326)
  m   <- sf::st_distance(pts, pts) / 1000          # metres → km
  units::drop_units(m)
}


# ==============================================================================
# STAGE 1: DBSCAN — identify and remove spatial outliers
# ==============================================================================
# Convert degrees to approximate Euclidean km for fast DBSCAN scan
# (valid for moderate-extent regions, e.g. Sachsen-Anhalt ~100×200 km)
sites_xy <- sites %>%
  mutate(
    x_km = lng * cos(mean(lat) * pi / 180) * 111.32,
    y_km = lat * 111.32
  )

db <- dbscan::dbscan(
  x       = as.matrix(sites_xy[, c("x_km", "y_km")]),
  eps     = DBSCAN_EPS_KM,
  minPts  = DBSCAN_MIN_PTS
)

sites_xy$dbscan_cluster <- db$cluster          # 0 = noise/outlier

n_outliers <- sum(db$cluster == 0)
n_core     <- sum(db$cluster  > 0)

cat(sprintf("[Stage 1] DBSCAN eps=%.0f km, minPts=%d\n",
            DBSCAN_EPS_KM, DBSCAN_MIN_PTS))
cat(sprintf("  Outlier sites removed : %d  (%.1f%%)\n",
            n_outliers, 100 * n_outliers / nrow(sites)))
cat(sprintf("  Core sites retained  : %d\n", n_core))

# Retain only non-outlier sites
sites_core <- sites_xy %>%
  filter(dbscan_cluster > 0) %>%
  select(site_id, lat, lng, area, x_km, y_km)


# ==============================================================================
# STAGE 2: HAC with Ward.D2 linkage, cut at HAC_MAX_RADIUS
# ==============================================================================
# Geodesic distance matrix on retained sites
cat("[Stage 2] Computing geodesic distance matrix ...\n")
D_geo <- geodist_km(sites_core)                 # n_core × n_core matrix

# Ward.D2 hierarchical clustering
hc <- hclust(as.dist(D_geo), method = "ward.D2")

# Iteratively find cut height such that max within-cluster distance ≤ threshold
# We use the cophenetic approach: cut at height h and check cluster compactness
find_cut_k <- function(hc, D, max_radius_km) {
  # Binary search over number of clusters k
  n <- nrow(D)
  k_lo <- 1L; k_hi <- n
  while (k_lo < k_hi) {
    k_mid <- (k_lo + k_hi) %/% 2L
    labels <- cutree(hc, k = k_mid)
    # Maximum within-cluster geodesic half-diameter
    max_r <- max(sapply(unique(labels), function(cl) {
      idx <- which(labels == cl)
      if (length(idx) < 2) return(0)
      max(D[idx, idx]) / 2          # radius = half the max pairwise distance
    }))
    if (max_r <= max_radius_km) k_hi <- k_mid else k_lo <- k_mid + 1L
  }
  k_lo
}

K_opt <- find_cut_k(hc, D_geo, max_radius_km = HAC_MAX_RADIUS)
cat(sprintf("  HAC cut: K = %d super-sites (max radius ≤ %.0f km)\n",
            K_opt, HAC_MAX_RADIUS))

sites_core$hac_cluster <- cutree(hc, k = K_opt)

n_clusters <- max(sites_core$hac_cluster)
# ==============================================================================
# AGGREGATE: build super-site data.frame for MILP input
# Each super-site:
#   lat/lng  → area-weighted centroid
#   area_ha  → sum of member site areas
#   n_sites  → cardinality of cluster
# ==============================================================================
sites_clustered <- sites_core %>%
  group_by(cluster_id = hac_cluster) %>%
  summarise(
    lat      = weighted.mean(lat, area),   # area-weighted centroid
    lng      = weighted.mean(lng, area),
    area_ha  = sum(area),
    n_sites  = n(),
    .groups  = "drop"
  ) %>%
  mutate(
    site_id = seq_len(n()),                   # re-index for MILP
    name    = paste0("SC_", site_id)          # super-site label
  ) %>%
  select(site_id, name, lat, lng, area_ha, n_sites, cluster_id)

cat(sprintf("\n=== CLUSTERING SUMMARY ===\n"))
cat(sprintf("  Original sites           : %d\n",  nrow(sites)))
cat(sprintf("  Stage 1 outliers removed : %d\n",  n_outliers))
cat(sprintf("  Stage 2 super-sites      : %d\n",  nrow(sites_clustered)))
cat(sprintf("  Reduction factor         : %.1fx\n", nrow(sites) / nrow(sites_clustered)))
cat(sprintf("  Mean cluster size        : %.1f sites\n",
            mean(sites_clustered$n_sites)))
cat(sprintf("  Total area preserved     : %.0f ha (%.1f%%)\n",
            sum(sites_clustered$area_ha),
            100 * sum(sites_clustered$area_ha) / sum(sites$area_ha)))


# ==============================================================================
# PLOT_CLUSTERING.R
# Visualisation of the two-stage DBSCAN + HAC site clustering result
#
# Produces:
#   (A) Static ggplot2 — side-by-side before/after comparison (PNG export)
#   (B) Interactive Leaflet — layered map with original + super-sites + storages
#
# Depends on: output of the two-stage clustering script (sites_core,
#             sites_clustered, db, afs_workspace)
# Colour conventions consistent with map_workspace.R / app_v10.r
# ==============================================================================

# ── Colour constants (consistent with map_workspace.R) ────────────────────────
COL_ORIG_CORE    <- "#74c476"    # light green  — retained original sites
COL_ORIG_OUTLIER <- "#fc9272"    # salmon red   — DBSCAN outliers (removed)
COL_SUPER_FILL   <- "#005a32"    # dark green   — super-site centroids
COL_SUPER_BORDER <- "#00441b"
COL_STORAGE      <- "#e65c00"    # orange       — storages (unchanged)
COL_VORONOI      <- "#add8e6"    # light blue   — cluster Voronoi fill

# Reference objects (produced by two-stage clustering script)
sites_orig     <- afs_workspace$sites      # original ~2000 sites
# → save before overwriting:
# afs_workspace_backup <- afs_workspace
sites_outlier  <- sites_xy %>% filter(dbscan_cluster == 0)
sites_inlier   <- sites_xy %>% filter(dbscan_cluster  > 0) %>%
  left_join(sites_core %>% select(site_id, hac_cluster), by = "site_id")

storages       <- afs_workspace$storages
super          <- sites_clustered


# ==============================================================================
# POLYGON-BASED VISUALISATION — replaces addCircleMarkers for site layers
#
# Requires: afs_workspace$feldblocks_sf  — sf object with columns
#             site_id  (join key to sites_xy / sites_core)
#             geometry (POLYGON or MULTIPOLYGON, CRS = 4326)
#
# If feldblocks_sf is in a different CRS, transform first:
#   afs_workspace$feldblocks_sf <- st_transform(afs_workspace$feldblocks_sf, 4326)
# ==============================================================================

# ── Join cluster labels onto the polygon layer ──────────────────────────────
fb_sf <- afs_workspace$sites_sf %>%
  mutate(site_id = seq_len(nrow(sites))) %>%
  rename(area_ha = area) %>%
  st_transform(4326) %>%
  left_join(
    sites_xy %>% select(site_id, dbscan_cluster),
    by = "site_id"
  ) %>%
  left_join(
    sites_core %>% select(site_id, hac_cluster),
    by = "site_id"
  ) %>%
  mutate(
    layer_type = case_when(
      is.na(dbscan_cluster) | dbscan_cluster == 0 ~ "outlier",
      TRUE                                         ~ "inlier"
    ),
    hac_cluster = if_else(layer_type == "inlier", hac_cluster, NA_integer_)
  )

# Separate layers for convenience
fb_outlier <- fb_sf %>% filter(layer_type == "outlier")
fb_inlier  <- fb_sf %>% filter(layer_type == "inlier")

# Polygon-level popups
popup_fb_orig <- paste0(
  "<b>", fb_sf$name, "</b><br>",
  "Fläche: ", round(fb_sf$area_ha, 1), " ha"
)

popup_fb_outlier <- paste0(
  "<b>Outlier (entfernt)</b><br>",
  "site_id: ", fb_outlier$site_id, "<br>",
  "Fläche: ", round(fb_outlier$area_ha, 1), " ha"
)

popup_fb_inlier <- paste0(
  "<b>", fb_inlier$name, "</b><br>",
  "HAC-Cluster: ", fb_inlier$hac_cluster, "<br>",
  "Fläche: ", round(fb_inlier$area_ha, 1), " ha"
)

# ── Leaflet colour palette (same Spectral ramp as before) ───────────────────
pal_cluster_poly <- colorFactor(
  palette = colorRampPalette(brewer.pal(11, "Spectral"))(n_clusters),
  domain  = as.factor(fb_inlier$hac_cluster),
  na.color = "transparent"
)

# ==============================================================================
# BUILD MAP — polygon version
# ==============================================================================
cluster_map <- leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
  
  # --- Base tiles ---
  addProviderTiles("CartoDB.Positron",  group = "CartoDB (hell)") %>%
  addProviderTiles("OpenStreetMap.DE",  group = "OpenStreetMap DE") %>%
  addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
  setView(lng = 11.85, lat = 51.55, zoom = 8) %>%
  addScaleBar(position = "bottomleft") %>%
  
  # ── Layer 1: All original Feldblock polygons ────────────────────────────────
  addPolygons(
    data        = fb_sf,
    fillColor   = COL_ORIG_CORE,
    fillOpacity = 0.45,
    color       = "#2d6a2d",
    weight      = 0.6,
    opacity     = 0.8,
    smoothFactor = 0.5,
    popup       = popup_fb_orig,
    group       = "Originalstandorte (alle)"
  ) %>%
  
  # ── Layer 2: DBSCAN outlier polygons ──────────────────────────────────────
  addPolygons(
    data        = fb_outlier,
    fillColor   = COL_ORIG_OUTLIER,
    fillOpacity = 0.65,
    color       = "#cb181d",
    weight      = 1.2,
    opacity     = 1.0,
    smoothFactor = 0.5,
    popup       = popup_fb_outlier,
    group       = "Stage 1: DBSCAN-Outlier (entfernt)"
  ) %>%
  
  # ── Layer 3: Cluster-member polygons coloured by HAC cluster ──────────────
  addPolygons(
    data        = fb_inlier,
    fillColor   = ~pal_cluster_poly(as.factor(hac_cluster)),
    fillOpacity = 0.55,
    color       = ~pal_cluster_poly(as.factor(hac_cluster)),
    weight      = 0.8,
    opacity     = 0.9,
    smoothFactor = 0.5,
    popup       = popup_fb_inlier,
    highlightOptions = highlightOptions(
      weight      = 2.5,
      color       = "#333333",
      fillOpacity = 0.85,
      bringToFront = TRUE
    ),
    group = "Stage 2: Cluster-Mitglieder (Polygone)"
  ) %>%
  
  # ── Layer 4: Super-site centroids (area-weighted, unchanged) ──────────────
  addCircleMarkers(
    data        = super,
    lng = ~lng, lat = ~lat,
    radius      = ~pmin(20, pmax(6, 4 + log1p(area_ha) * 1.5)),
    color       = COL_SUPER_BORDER,
    fillColor   = COL_SUPER_FILL,
    weight = 2.0, opacity = 1.0, fillOpacity = 0.9,
    popup       = popup_fb_inlier,
    label       = ~paste0("SC_", site_id, " (", n_sites, " sites, ",
                          round(area_ha, 0), " ha)"),
    labelOptions = labelOptions(noHide = FALSE, direction = "top",
                                textsize = "11px"),
    group       = "Super-Sites (aggregiert)"
  ) %>%
  
  # ── Layer 5: Storages (unchanged) ─────────────────────────────────────────
  addCircleMarkers(
    data        = storages,
    lng = ~lng, lat = ~lat,
    radius      = 9,
    color       = "#c05000", fillColor = "#fd8d3c",
    weight = 2.5, opacity = 1.0, fillOpacity = 0.9,
    #popup       = popup_storage,
    label       = ~name,
    labelOptions = labelOptions(noHide = FALSE, direction = "top",
                                textsize = "12px"),
    group       = "Vorbehandlung / Umschlag"
  ) %>%
  
  # ── Legend ─────────────────────────────────────────────────────────────────
  addLegend(
    position = "bottomright",
    colors   = c(COL_ORIG_CORE, COL_ORIG_OUTLIER,
                 COL_SUPER_FILL, COL_STORAGE),
    labels   = c(
      paste0("Feldblöcke original (n = ", nrow(sites_orig), ")"),
      paste0("DBSCAN-Outlier entfernt (n = ", n_outliers, ")"),
      paste0("Super-Sites, Zentroid (K = ", nrow(super), ")"),
      "Vorbehandlung / Umschlag"
    ),
    opacity  = 0.9,
    title    = HTML("<b>Clustering-Ergebnis</b><br>
                    <small>DBSCAN + HAC Ward.D2</small>")
  ) %>%
  
  # ── Layer control ──────────────────────────────────────────────────────────
  addLayersControl(
    baseGroups    = c("CartoDB (hell)", "OpenStreetMap DE", "Satellite"),
    overlayGroups = c(
      "Originalstandorte (alle)",
      "Stage 1: DBSCAN-Outlier (entfernt)",
      "Stage 2: Cluster-Mitglieder (Polygone)",
      "Super-Sites (aggregiert)",
      "Vorbehandlung / Umschlag"
    ),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  hideGroup("Originalstandorte (alle)") %>%
  hideGroup("Stage 1: DBSCAN-Outlier (entfernt)")

cluster_map

# ==============================================================================
# UPDATE WORKSPACE — replace sites, clear stale distance matrix
# afs_workspace$dist_ij must be recomputed via OSRM after this step
# ==============================================================================


source("../Modell/!helper_func.r")
source("../Modell/!helper_instance_builder_v8a.R")
source("../Modell/build_agroforestry_lp_v10.R")


afs_workspace$sites_clustered <- sites_clustered %>%
  select(site_id, name, lat, lng, area_ha, n_sites)

afs_workspace$site_cluster_assig <- sites_core %>%
  select(site_id, lat, lng, area, hac_cluster)


# add storage cost column to storages (if column c_stor is not existing)
if("c_stor" %in% colnames(afs_workspace$storages)) {
  afs_workspace$storages <- afs_workspace$storages %>%
    mutate(c_stor = 10)  # example: 10 €/t and year storage cost
}

# calculate new distance matrix with OSRM for the clustered sites

dist_ij_clustered <- calculate_distance_matrix_osm(
  starts       = afs_workspace$sites_clustered %>% select(lat, lng),
  destinations = afs_workspace$storages %>% select(lat, lng),
  max_entries  = 100
) 


dist_ij <- calculate_distance_matrix_osm(
  starts       = afs_workspace$sites %>% select(lat, lng),
  destinations = afs_workspace$storages %>% select(lat, lng),
  max_entries  = 100
) 

dist_jk <- calculate_distance_matrix_osm(
  starts       = afs_workspace$storages %>% select(lat, lng),
  destinations = afs_workspace$consumers %>% select(lat, lng),
  max_entries  = 100
) 

afs_workspace$dist_ij_clust <- dist_ij_clustered$distance_matrix_km
afs_workspace$dist_ij <- dist_ij$distance_matrix_km
afs_workspace$dist_jk <- dist_jk$distance_matrix_km

afs_workspace$meta$n_sites          <- nrow(afs_workspace$sites_clustered)
afs_workspace$meta$clustering <- list(
  method       = "DBSCAN + HAC Ward.D2",
  dbscan_eps   = DBSCAN_EPS_KM,
  dbscan_minPts= DBSCAN_MIN_PTS,
  hac_radius   = HAC_MAX_RADIUS,
  K            = nrow(sites_clustered),
  n_outliers   = n_outliers,
  date         = Sys.time()
)

save(afs_workspace, file = "afs_workspace_red.RData")
cat("\n✓ afs_workspace.RData updated with clustered sites.\n")
cat("  → Recompute dist_ij via OSRM before running build_agroforestry_lp_v10.r\n")


