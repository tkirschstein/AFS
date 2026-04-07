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
})

# ── Parameters (adjust to study region) ──────────────────────────────────────
DBSCAN_EPS_KM    <- 15    # neighbourhood radius for DBSCAN (km); sites beyond
# this distance from any neighbour = outlier
DBSCAN_MIN_PTS   <- 3     # minimum cluster size in Stage 1
HAC_MAX_RADIUS   <- 5    # maximum within-cluster geodesic radius (km) for
# Stage 2 cut; all members ≤ this distance to medoid
# ─────────────────────────────────────────────────────────────────────────────

load("../Modell/afs_workspace.RData")

sites <- afs_workspace$sites   # columns: site_id, lat, lng, area_ha

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
  select(site_id, lat, lng, area_ha, x_km, y_km)


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
    lat      = weighted.mean(lat, area_ha),   # area-weighted centroid
    lng      = weighted.mean(lng, area_ha),
    area_ha  = sum(area_ha),
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

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(sf)
  library(patchwork)     # install.packages("patchwork")
  library(leaflet)
  library(leaflet.extras2)
  library(htmltools)
  library(RColorBrewer)
  library(scales)
})

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
# (A) STATIC GGPLOT2 — side-by-side before/after
# ==============================================================================

# Cluster colour palette (one colour per HAC cluster)
n_clusters    <- nrow(super)
cluster_pal   <- colorRampPalette(
  brewer.pal(min(n_clusters, 11), "Spectral"))(n_clusters)
names(cluster_pal) <- as.character(seq_len(n_clusters))

sites_inlier <- sites_inlier %>%
  mutate(cluster_col = cluster_pal[as.character(hac_cluster)])

# --- Panel A: original sites (before clustering) ---
p_before <- ggplot() +
  geom_point(
    data = sites_orig,
    aes(x = lng, y = lat),
    colour  = COL_ORIG_CORE, fill = COL_ORIG_CORE,
    shape   = 21, size = 0.6, alpha = 0.5, stroke = 0
  ) +
  geom_point(
    data = storages,
    aes(x = lng, y = lat),
    colour = COL_STORAGE, shape = 17, size = 2.5
  ) +
  coord_quickmap() +
  labs(
    title    = "Before: Original Sites",
    subtitle = paste0(nrow(sites_orig), " Feldblöcke"),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid       = element_line(colour = "grey92"),
    plot.title       = element_text(face = "bold"),
    axis.text        = element_text(size = 7),
    legend.position  = "none"
  )

# --- Panel B: after clustering —  original points coloured by cluster,
#              super-site centroids overlaid, outliers shown in red ---
p_after <- ggplot() +
  # DBSCAN outliers
  geom_point(
    data  = sites_outlier,
    aes(x = lng, y = lat),
    colour = COL_ORIG_OUTLIER, shape = 4,
    size = 1.0, alpha = 0.6, stroke = 0.4
  ) +
  # Cluster members coloured by HAC cluster
  geom_point(
    data = sites_inlier,
    aes(x = lng, y = lat, colour = as.factor(hac_cluster)),
    shape = 16, size = 0.7, alpha = 0.6
  ) +
  # Super-site centroids
  geom_point(
    data = super,
    aes(x = lng, y = lat, size = area_ha),
    colour = COL_SUPER_BORDER, fill = COL_SUPER_FILL,
    shape = 21, stroke = 1.2, alpha = 0.9
  ) +
  # Super-site labels (cluster ID)
  geom_text(
    data = super %>% filter(n_sites >= 5),
    aes(x = lng, y = lat, label = cluster_id),
    size = 2, colour = "white", fontface = "bold",
    vjust = 0.5, hjust = 0.5
  ) +
  # Storages
  geom_point(
    data = storages,
    aes(x = lng, y = lat),
    colour = COL_STORAGE, shape = 17, size = 2.5
  ) +
  scale_colour_manual(values = cluster_pal) +
  scale_size_continuous(
    name   = "Area (ha)",
    range  = c(2, 8),
    breaks = c(50, 200, 500, 1000)
  ) +
  coord_quickmap() +
  labs(
    title    = "After: Clustered Super-Sites",
    subtitle = paste0(
      nrow(super),     " super-sites  |  ",
      n_outliers,      " outliers removed  |  ",
      "radius ≤ ", HAC_MAX_RADIUS, " km"
    ),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid      = element_line(colour = "grey92"),
    plot.title      = element_text(face = "bold"),
    axis.text       = element_text(size = 7),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(colour = "none")   # suppress cluster colour legend (too many levels)

# --- Combine with patchwork and add annotation ---
p_combined <- p_before + p_after +
  plot_annotation(
    title   = "AFS Site Clustering: DBSCAN + HAC Ward.D2",
    caption = paste0(
      "Stage 1 DBSCAN (eps=", DBSCAN_EPS_KM, " km, minPts=", DBSCAN_MIN_PTS,
      "):  ", n_outliers, " outliers removed.  ",
      "Stage 2 HAC (max radius ", HAC_MAX_RADIUS, " km):  ",
      nrow(sites_orig), " → ", nrow(super), " sites  ",
      "(reduction factor ", round(nrow(sites_orig) / nrow(super), 1), "x)."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", size = 13),
      plot.caption = element_text(size = 7, colour = "grey50")
    )
  )

ggsave(
  filename = "plot_clustering_comparison.png",
  plot     = p_combined,
  width    = 14, height = 7, dpi = 180, bg = "white"
)
cat("✓ Saved: plot_clustering_comparison.png\n")


# ==============================================================================
# (B) INTERACTIVE LEAFLET — layers for original, outliers, super-sites
# ==============================================================================

# Leaflet colour palette: one colour per cluster (max 256 distinct)
pal_cluster <- colorFactor(
  palette = colorRampPalette(brewer.pal(11, "Spectral"))(n_clusters),
  domain  = as.factor(sites_inlier$hac_cluster)
)

# Popup helpers
popup_orig <- paste0(
  "<b>", sites_orig$name, "</b><br>",
  "Fläche: ", round(sites_orig$area_ha, 1), " ha"
)

popup_super <- paste0(
  "<b>Super-Site SC_", super$site_id, "</b><br>",
  "Cluster-ID: ", super$cluster_id, "<br>",
  "Mitglieder: ", super$n_sites, " Feldblöcke<br>",
  "Gesamtfläche: ", round(super$area_ha, 0), " ha<br>",
  "Zentroid: ", round(super$lat, 4), "° N, ", round(super$lng, 4), "° E"
)

popup_outlier <- paste0(
  "<b>Outlier (entfernt)</b><br>",
  "site_id: ", sites_outlier$site_id, "<br>",
  "Fläche: ", round(sites_outlier$area_ha, 1), " ha"
)

popup_storage <- paste0(
  "<b>", storages$name, "</b><br>",
  "Typ: ", storages$type, "<br>",
  "CAP_proc: ", storages$CAP_proc, " kt/Jahr"
)

cluster_map <- leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
  
  # --- Base tiles ---
  addProviderTiles("CartoDB.Positron",  group = "CartoDB (hell)") %>%
  addProviderTiles("OpenStreetMap.DE",  group = "OpenStreetMap DE") %>%
  addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
  setView(lng = 11.85, lat = 51.55, zoom = 8) %>%
  addScaleBar(position = "bottomleft") %>%
  
  # ── Layer 1: Original sites (all) ──────────────────────────────────────────
  addCircleMarkers(
    data        = sites_orig,
    lng = ~lng, lat = ~lat,
    radius      = 3,
    color       = "#2d6a2d", fillColor = COL_ORIG_CORE,
    weight = 0.8, opacity = 0.9, fillOpacity = 0.5,
    popup       = popup_orig,
    group       = "Originalstandorte (alle)"
  ) %>%
  
  # ── Layer 2: DBSCAN outliers (removed sites) ───────────────────────────────
  addCircleMarkers(
    data        = sites_outlier,
    lng = ~lng, lat = ~lat,
    radius      = 4,
    color       = "#cb181d", fillColor = COL_ORIG_OUTLIER,
    weight = 1.5, opacity = 1.0, fillOpacity = 0.7,
    popup       = popup_outlier,
    group       = "Stage 1: DBSCAN-Outlier (entfernt)"
  ) %>%
  
  # ── Layer 3: Cluster members coloured by HAC cluster ──────────────────────
  addCircleMarkers(
    data        = sites_inlier,
    lng = ~lng, lat = ~lat,
    radius      = 3,
    color       = ~pal_cluster(as.factor(hac_cluster)),
    fillColor   = ~pal_cluster(as.factor(hac_cluster)),
    weight = 0.5, opacity = 0.8, fillOpacity = 0.6,
    popup = ~paste0(
      "<b>", name, "</b><br>",
      "HAC-Cluster: ", hac_cluster, "<br>",
      "Fläche: ", round(area_ha, 1), " ha"
    ),
    group = "Stage 2: Cluster-Mitglieder"
  ) %>%
  
  # ── Layer 4: Super-site centroids (area-weighted) ─────────────────────────
  addCircleMarkers(
    data        = super,
    lng = ~lng, lat = ~lat,
    radius      = ~pmin(20, pmax(6, 4 + log1p(area_ha) * 1.5)),
    color       = COL_SUPER_BORDER,
    fillColor   = COL_SUPER_FILL,
    weight = 2.0, opacity = 1.0, fillOpacity = 0.85,
    popup       = popup_super,
    label       = ~paste0("SC_", site_id, " (", n_sites, " sites, ",
                          round(area_ha, 0), " ha)"),
    labelOptions= labelOptions(noHide = FALSE, direction = "top",
                               textsize = "11px"),
    group       = "Super-Sites (aggregiert)"
  ) %>%
  
  # ── Layer 5: Storages ──────────────────────────────────────────────────────
  addCircleMarkers(
    data        = storages,
    lng = ~lng, lat = ~lat,
    radius      = 9,
    color       = "#c05000", fillColor = "#fd8d3c",
    weight = 2.5, opacity = 1.0, fillOpacity = 0.9,
    popup       = popup_storage,
    label       = ~name,
    labelOptions= labelOptions(noHide = FALSE, direction = "top",
                               textsize = "12px"),
    group       = "Vorbehandlung / Umschlag"
  ) %>%
  
  # ── Legend ─────────────────────────────────────────────────────────────────
  addLegend(
    position = "bottomright",
    colors   = c(COL_ORIG_CORE, COL_ORIG_OUTLIER,
                 COL_SUPER_FILL, COL_STORAGE),
    labels   = c(
      paste0("Original-Standorte (n = ", nrow(sites_orig), ")"),
      paste0("DBSCAN-Outlier entfernt (n = ", n_outliers, ")"),
      paste0("Super-Sites, fläcengewichteter Zentroid (K = ", nrow(super), ")"),
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
      "Stage 2: Cluster-Mitglieder",
      "Super-Sites (aggregiert)",
      "Vorbehandlung / Umschlag"
    ),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  hideGroup("Originalstandorte (alle)") %>%   # hidden by default (performance)
  hideGroup("Stage 1: DBSCAN-Outlier (entfernt)")

# Export as self-contained HTML
htmlwidgets::saveWidget(
  widget   = cluster_map,
  file     = "plot_clustering_map.html",
  selfcontained = TRUE
)
cat("✓ Saved: plot_clustering_map.html\n")
cat("  → Open in browser or htmltools::browsable(cluster_map) in RStudio\n")

# Print in RStudio viewer
cluster_map


# ==============================================================================
# UPDATE WORKSPACE — replace sites, clear stale distance matrix
# afs_workspace$dist_ij must be recomputed via OSRM after this step
# ==============================================================================
afs_workspace$sites   <- sites_clustered
afs_workspace$dist_ij <- NULL               # invalidated; recompute with OSRM

afs_workspace$meta$n_sites          <- nrow(sites_clustered)
afs_workspace$meta$clustering <- list(
  method       = "DBSCAN + HAC Ward.D2",
  dbscan_eps   = DBSCAN_EPS_KM,
  dbscan_minPts= DBSCAN_MIN_PTS,
  hac_radius   = HAC_MAX_RADIUS,
  K            = nrow(sites_clustered),
  n_outliers   = n_outliers,
  date         = Sys.time()
)

save(afs_workspace, file = "afs_workspace.RData")
cat("\n✓ afs_workspace.RData updated with clustered sites.\n")
cat("  → Recompute dist_ij via OSRM before running build_agroforestry_lp_v10.r\n")


# ==============================================================================
# OPTIONAL: Silhouette diagnostic plot
# ==============================================================================
if (requireNamespace("cluster", quietly = TRUE) && K_opt > 1) {
  sil <- cluster::silhouette(
    x    = sites_core$hac_cluster,
    dist = as.dist(D_geo)
  )
  avg_sil <- mean(sil[, "sil_width"])
  cat(sprintf("  Average silhouette width : %.3f  (>0.5 = good structure)\n", avg_sil))
}


