# ============================================================================
# CALCULATE DISTANCE MATRIX (simple version)
# ============================================================================

calculate_distance_matrix <- function(from_points, to_points) {
  if (!inherits(from_points, "sf")) {
    coords_from <- as.matrix(from_points[, c("lng", "lat")])
  } else {
    coords_from <- sf::st_coordinates(from_points)[, c("X", "Y")]
  }
  
  if (!inherits(to_points, "sf")) {
    coords_to <- as.matrix(to_points[, c("lng", "lat")])
  } else {
    coords_to <- sf::st_coordinates(to_points)[, c("X", "Y")]
  }
  
  # Euclidean distance
  as.matrix(dist(rbind(coords_from, coords_to)))[
    1:nrow(coords_from), 
    (nrow(coords_from) + 1):(nrow(coords_from) + nrow(coords_to))
  ]
}

# ============================================================================
# CALCULATE DISTANCE MATRIX (OSM version)
# ============================================================================

calculate_distance_matrix_osm <- function(starts, destinations, max_entries = 100) {
  library(osrm)
  
  # Überprüfe Eingabeformat
  if (!all(c('lat', 'lng') %in% names(starts))) 
    stop('Starts müssen Spalten lat und lng enthalten')
  if (!all(c('lat', 'lng') %in% names(destinations))) 
    stop('Ziele müssen Spalten lat und lng enthalten')
  
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
      starts_sf <- st_as_sf(starts_block, coords = c("lng", "lat"), crs = 4326)
      dests_sf <- st_as_sf(dests_block, coords = c("lng", "lat"), crs = 4326)
      
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


# ============================================================================
# OPTIMIZATION INSTANCE BUILDER
# ============================================================================

build_optimization_instance <- function(data, params) {
  
  # check params completeness: n_periods, max_age,   min_age,  c_tr_raw,  c_tr_pre
  required_params <- c("n_periods", "max_age", "min_age", "c_tr_raw", "c_tr_pre")
  missing_params <- setdiff(required_params, names(params))
  if (length(missing_params) > 0) {
    stop(paste("Missing parameters:", paste(missing_params, collapse = ", "
    )))
  }
  
  # Extract parameters
  n_periods <- params$n_periods
  max_age <- params$max_age
  min_age <- params$min_age
  c_tr_raw <- params$c_tr_raw
  c_tr_pre <- params$c_tr_pre
  
  # check data completeness: sites, storages, consumers, dist_ij, dist_jk, yields_by_age
  required_data <- c("sites", "storages", "consumers", "dist_ij", "dist_jk","yields_by_age")
  missing_data <- setdiff(required_data, names(data))
  if (length(missing_data) > 0) {
    stop(paste("Missing data components:", paste(missing_data, collapse = ", ")))
  }
  
  # extract data objects
  d_ij <- data$dist_ij
  d_jk <- data$dist_jk
  sites <- data$sites
  storages <- data$storages
  consumers <- data$consumers
  yields_by_age <- data$yields_by_age
  
  
  # check matching sizes and extract numbers of entities
  n_sites <- nrow(sites)
  n_storages <- nrow(storages)
  n_consumers <- nrow(consumers)
  n_products <- length(unique(yields_by_age$product))
  
  if (nrow(d_ij) != n_sites || ncol(d_ij) != n_storages) {
    stop("Dimension mismatch: d_ij should be n_sites x n_storages")
  }
  if (nrow(d_jk) != n_storages || ncol(d_jk) != n_consumers) {
    stop("Dimension mismatch: d_jk should be n_storages x n_consumers")
  }
  if (!all(c("product", "age", "yield_ha") %in% colnames(yields_by_age))) {
    stop("yields_by_age must have columns: product, age, yield_ha")
  }
  
  if (max(yields_by_age$age) > max_age) {
    # set yields to 0 for ages > max_age
    warning(paste("yields_by_age contains ages > max_age; setting yields to 0 for those ages"))
    yields_by_age <- yields_by_age %>%  mutate(yield_ha = ifelse(age > max_age, 0, yield_ha))
  }
  
   if (max(yields_by_age$product) > n_products) {
    warning(paste("yields_by_age contains products > n_products; trimming to", n_products))
    yields_by_age <- yields_by_age %>% filter(product <= n_products)
  }  
  
  # check columns in sites: Should contain site_id, lat, lng, area_ha, C_est, C_harv, C_main, C_opp
  required_site_cols <- c("site_id", "lat", "lng", "area_ha", "C_est", "C_harv", "C_main", "C_opp")
  if (!all(required_site_cols %in% colnames(sites))) {
    stop(paste("Sites data frame must contain columns:", paste(required_site_cols, collapse = ", ")))
  }
  
  # check columns in storages: Should contain storage_id, CAP_stor, CAP_proc, c_stor
  required_storage_cols <- c("storage_id", "CAP_stor", "CAP_proc", "c_stor")
  if (!all(required_storage_cols %in% colnames(storages))) {
    stop(paste("Storages data frame must contain columns:", paste(required_storage_cols, collapse = ",")))
  }
  
  # check columns in consumers: Should contain consumer_id, demand_P1, demand_P2, demand_P3, P1, P2, P3
  required_consumer_cols <- c("consumer_id", "demand_P1", "demand_P2", "demand_P3", "P1", "P2", "P3")
  if (!all(required_consumer_cols %in% colnames(consumers))) {
    stop(paste("Consumers data frame must contain columns:", paste(required_consumer_cols, collapse=","
    )))
  }    
  
  # ids to consecutive numbers (if not already)
  sites <- sites %>% mutate(site_id = as.integer(factor(site_id)))
  storages <- storages %>% mutate(storage_id = as.integer(factor(storage_id)))
  consumers <- consumers %>% mutate(consumer_id = as.integer(factor(consumer_id)))
  # distances should be numeric matrices with appropriate row/col names
  if (!is.matrix(d_ij) || !is.numeric(d_ij)) {
    stop("d_ij must be a numeric matrix")
  }
  if (!is.matrix(d_jk) || !is.numeric(d_jk)) {
    stop("d_jk must be a numeric matrix")
  }
  # set row and col names of matrices to integers as ids
  rownames(d_ij) <- as.character(1:n_sites)
  colnames(d_ij) <- as.character(1:n_storages)
  rownames(d_jk) <- as.character(1:n_storages)
  colnames(d_jk) <- as.character(1:n_consumers)
  
  
  # build time-exanded demand matrix based on consumers data (demand_P1 demand_P2 demand_P3) with columns consumer_id, product, period and demand
  demand <- expand.grid(
    consumer_id = consumers$consumer_id,
    period = 1:n_periods
  ) %>%
    left_join(
      consumers %>% 
        select(consumer_id, demand_P1, demand_P2, demand_P3) %>%
        pivot_longer(
          cols = starts_with("demand_P"),
          names_to = "product_col",
          values_to = "D_max"
        ) %>%
        mutate(product = as.numeric(substring(product_col, nchar(product_col), nchar(product_col)))),
       by = "consumer_id") %>%
    select(consumer_id, product, period, D_max) 
    #mutate(D_max = D_max / n_periods)  # distribute demand evenly across periods
  
  # ========================================================================
  # PRICES: Extract from consumer table (P1, P2, P3) if available
  # ========================================================================
  
  # Check if consumer table has price columns
  price_cols <- c("P1", "P2", "P3")
  has_prices <- all(price_cols %in% colnames(consumers))
  
  consumer_prices <- NULL
  if (has_prices) {
    # Create price table: consumer_id x product x price
    consumer_prices <- consumers %>%
      select(consumer_id, all_of(price_cols)) %>%
      pivot_longer(
        cols = all_of(price_cols),
        names_to = "price_col",
        values_to = "price"
      ) %>%
      mutate(
        product = as.numeric(gsub("P", "", price_col))
      ) %>%
      select(consumer_id, product, price)
  } else {
    # Fallback: uniform prices by product (optional, can be set later)
    consumer_prices <- expand.grid(
      consumer_id = consumers$consumer_id,
      product = 1:n_products
    ) %>%
      mutate(price = 100)  # Default uniform price
  }
  
  
  # Build instance list
  list(
    n_sites = n_sites,
    n_storages = n_storages,
    n_consumers = n_consumers,
    n_periods = n_periods,
    n_products = n_products,
    max_age = max_age,
    min_age = min_age,
    sites = sites,
    storages = storages,
    consumers = consumers,
    consumer_prices = consumer_prices,
    yields_by_age = yields_by_age,
    demand = demand,
    d_ij = d_ij,
    d_jk = d_jk,
    c_tr_raw = c_tr_raw,
    c_tr_pre = c_tr_pre
  )
}


# ============================================================================
# Solve function
# ============================================================================

solve_agroforestry_sparse <- function(instance,
                                      time_limit = 600,
                                      mip_gap = 0.05,
                                      verbose = TRUE, 
                                      solver = "highs") {
  
  
  # Build LP
  built <- build_AFS_milp(instance)
  
  # Solve with Gurobi via ROI
  
  if(solver == "gurobi"){
    result <- ROI_solve(
      built$model,
      solver = "gurobi", 
      control = list(
        TimeLimit = time_limit,
        MIPGap = mip_gap,
        OutputFlag = ifelse(verbose, 1, 0)
      )
    )
  }
  if(solver == "highs"){
    result <- ROI_solve(
      built$model,
      solver = "highs",
      control = list(
        time_limit  = as.numeric(time_limit),
        mip_rel_gap  = as.numeric(mip_gap),
        log_to_console  = ifelse(verbose, TRUE, FALSE)
      )
    )
  }
  
  
  
  # convert results
  solution <- lapply(built$var_maps, function(x) {
    tmp <- data.frame(x, 
                      value = result$solution[x$col])
    tmp %>% select(-col)
  })
  
  # Return complete result
  list(
    result = result,
    objective = result$objval,
    solution = solution,
    status = result$status,
    var_maps = built$var_maps
    #instance_info = built$instance_info,
    #yield_matrix = built$yield_matrix
  )
}


#######################################################################
# Extract solutions
#######################################################################
extract_result <- function(opt_result){
  
  z_solution <- opt_result$solution$z
  
  # filter for z variables with value > 0.5 and drop col column
  z_solution_filtered <- z_solution %>%
    filter(value > 0.5) %>%
    arrange(i, s, t)
  
  # sites with AFS
  
  sites_est <- unique(z_solution_filtered$i)
  
  # storage quantities
  S_solution <- opt_result$solution$S
  
  S_solution_filtered <- S_solution %>% 
    filter(value > 0.01) %>%
    arrange(j, t, p)
  
  # Flows
  xij <- opt_result$solution$Xij
  
  xij_filtered <- xij %>% 
    filter(value > 0.01) %>% 
    arrange(i,j,t,p)
  
  xjk <- opt_result$solution$Xjk
  
  xjk_filtered <- xjk %>% 
    filter(value > 0.01) %>% 
    arrange(j,k,t,p)
  
  # Demand fulfillment
  if("Dkpt" %in% names(opt_result$solution)){
    Dkpt_solution <- opt_result$solution$Dkpt
    Dkpt_solution <- Dkpt_solution %>% 
      filter(value > 0.01) %>% 
      arrange(k,p,t)
  }else{
    Dkpt_solution <- NULL
  }
  
  
  
  # return results
  list(
    sites_est = sites_est, 
    z = z_solution_filtered,
    S = S_solution_filtered,
    Xij = xij_filtered,
    Xjk = xjk_filtered,
    Dkpt_solution = Dkpt_solution
  )
  
  
}