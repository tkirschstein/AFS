# function to calculate shortest route between two location via openstreetmap

#' @title Calculate Shortest Route
#' @description This function calculates the shortest route between two locations using OpenStreetMap data.
#' @param start A vector of two numeric values representing the latitude and longitude of the starting point.
#' @param end A vector of two numeric values representing the latitude and longitude of the destination point.
#' @param mode A character string indicating the mode of transportation. Options include "driving", "walking", "cycling", and "public_transport".
#' @return A list containing the route information, including the distance and duration of the route.
#' @examples
#' #' start <- c(37.7749, -122.4194)  # San Francisco
#' #' end <- c(34.0522, -118.2437)    # Los Angeles
#' #' #' route_info <- calculate_shortest_route(start, end, mode = "driving")
#' #' @export




calculate_route_info <- function(src_lon, src_lat, dst_lon, dst_lat) {
  # Start- und Zielkoordinaten als Vektor (lon, lat)
  start <- c(src_lon, src_lat)
  ziel <- c(dst_lon, dst_lat)
  
  # Route berechnen (Profil: Auto)
  route <- osrmRoute(src = start, dst = ziel, overview = "full") 
  #, returnclass = "sf"
  
  # Extrahiere Distanz (km) und Dauer (min)
  distance_km <- route$distance
  duration_min <- route$duration
  
  # Rückgabe als Liste mit Route (sf-Objekt), Distanz und Dauer
  return(list(route = route, distance_km = distance_km, duration_min = duration_min))
}


#result <- calculate_route_info(11.966, 51.482, 12.387, 51.343)

# Beispielkoordinaten:
# Quelle: Halle (Saale)
# Ziel: Leipzig
# route_sf <- get_shortest_path(
#   src_lon = 11.966,
#   src_lat = 51.482,
#   dst_lon = 12.387,
#   dst_lat = 51.343
# )


# m <- m %>% addPolylines(
#   data = result$route,
#   color = 'blue',
#   weight = 4,
#   opacity = 0.7,
#   group = 'Route',
#   popup = paste0("Distanz: ", round(result$distance_km, 2), " km<br>",
#                  "Fahrzeit: ", round(result$duration_min, 1), " min")
# )
# m
# 
# 
# # Beispielaufruf der Funktion
# # feldblock koordinaten
# fb_coords <- st_coordinates(fb_zentroide_wgs84)
# fb_coords <- data.frame(
#   lat = fb_coords[, "Y"],
#   lon = fb_coords[, "X"]
# )
# 
# # wald-zentroide koordinaten
# wald_coords <- st_coordinates(wald_zentroide_wgs84)
# wald_coords <- data.frame(
#   lat = wald_coords[, "Y"],
#   lon = wald_coords[, "X"]
# )


calculate_distance_matrix <- function(starts, destinations, max_entries = 100) {
  # Überprüfe Eingabeformat
  if (!all(c('lat', 'lon') %in% names(starts))) 
    stop('Starts müssen Spalten lat und lon enthalten')
  if (!all(c('lat', 'lon') %in% names(destinations))) 
    stop('Ziele müssen Spalten lat und lon enthalten')
  
  n_starts <- nrow(starts)
  n_dests <- nrow(destinations)
  
  # Initialisiere Ergebnis-Matrizen
  duration_matrix_min <- dist_matrix_km <- matrix(NA, nrow = n_starts, ncol = n_dests)
  
  # Startpunkte in Blöcke unterteilen
  start_blocks <- split(
    1:n_starts, 
    ceiling(seq_len(n_starts) / max_entries)
  )
  
  # Zielpunkte in Blöcke unterteilen
  dest_blocks <- split(
    1:n_dests, 
    ceiling(seq_len(n_dests) / max_entries)
  )
  
  # Iteriere über alle Block-Kombinationen
  for (s_block in start_blocks) {
    for (d_block in dest_blocks) {
      # Extrahiere aktuelle Blockdaten
      starts_block <- starts[s_block, , drop = FALSE]
      dests_block <- destinations[d_block, , drop = FALSE]
      
      # Konvertiere zu sf-Objekten
      starts_sf <- st_as_sf(starts_block, coords = c("lon", "lat"), crs = 4326)
      dests_sf <- st_as_sf(dests_block, coords = c("lon", "lat"), crs = 4326)
      
      # API-Anfrage mit Fehlerbehandlung
      tmp.res <- tryCatch({
        osrmTable(
          src = starts_sf,
          dst = dests_sf,
          measure = c("distance", "duration")
        )
      }, error = function(e) {
        message("Fehler bei Block: ", paste(s_block, collapse=","), " zu ", paste(d_block, collapse=","))
        return(NULL)
      })
      
      if (!is.null(tmp.res)) {
        # Setze Ergebnisse in Gesamtmatrix ein
        dist_matrix_km[s_block, d_block] <- tmp.res$distances / 1000
        duration_matrix_min[s_block, d_block] <- tmp.res$durations
      }
    }
  }
  
  # Benennung der Zeilen und Spalten
  rownames(dist_matrix_km) <- rownames(duration_matrix_min) <- paste0('start_', 1:n_starts)
  colnames(dist_matrix_km) <- colnames(duration_matrix_min) <- paste0('dest_', 1:n_dests)
  
  return(list(
    distance_matrix_km = dist_matrix_km,
    duration_matrix_min = duration_matrix_min
  ))
}

#dist.wald.sawmills <- calculate_distance_matrix(wald_coords, sawmills)
#dist.fb.sawmills <- calculate_distance_matrix(fb_coords, sawmills)
#dist.fb.biochem <- calculate_distance_matrix(fb_coords, biochem_sites)
#dist.wald.biochem <- calculate_distance_matrix(wald_coords, biochem_sites)
#dist.sawmills.biochem <- calculate_distance_matrix(sawmills, biochem_sites)



# Calculate distance matrix between two point sets
calculate_distance_matrix <- function(from_points, to_points) {
  if (!inherits(from_points, "sf")) {
    coords_from <- as.matrix(from_points[, c("lng", "lat")])
  } else {
    coords_from <- sf::st_coordinates(from_points)
  }
  
  if (!inherits(to_points, "sf")) {
    coords_to <- as.matrix(to_points[, c("lng", "lat")])
  } else {
    coords_to <- sf::st_coordinates(to_points)
  }
  
  # Simple Euclidean distance (replace with geospatial distance if needed)
  as.matrix(dist(rbind(coords_from, coords_to)))[
    1:nrow(coords_from), 
    (nrow(coords_from) + 1):(nrow(coords_from) + nrow(coords_to))
  ]
}

