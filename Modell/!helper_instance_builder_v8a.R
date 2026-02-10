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
  # Extract parameters
  n_periods <- params$n_periods
  max_age <- params$max_age
  revenue_scale <- params$revenue_scale
  c_tr_raw <- params$c_tr_raw
  c_tr_pre <- params$c_tr_pre
  
  # Filter sites based on selection
  selected_sites <- data$sites[data$sites$site_id %in% params$selected_sites, ]
  
  # Calculate distances
  d_ij <- data$dist_ij
  
  d_jk <- data$dist_jk
  
  # Define product yields (simplified: 1 product with age-dependent yields)
  # Format: product, age, yield_ha
  yields_by_age <- data.frame(
    product = rep(1, n_periods),
    age = 1:n_periods,
    yield_ha = c(0, 3, 8, 14, 20, 25, 28, 30, 31, 31, 30, 28, 25, 23, 20)[1:n_periods]
  )
  
  # Define demand (simple constant demand)
  demand <- expand.grid(
    consumer_id = data$consumers$consumer_id,
    product = 1,
    period = 1:n_periods
  ) %>%
    left_join(data$consumers[, c("consumer_id", "demand")], by = "consumer_id") %>%
    mutate(D_max = demand / n_periods) %>%
    select(consumer_id, product, period, D_max)
  
  # Create revenue scale vector
  R_p <- rep(revenue_scale, 1)  # 1 product
  
  # Build instance list
  list(
    n_sites = nrow(selected_sites),
    n_storages = nrow(data$storages),
    n_consumers = nrow(data$consumers),
    n_periods = n_periods,
    n_products = 1,
    max_age = max_age,
    sites = selected_sites,
    storages = data$storages,
    consumers = data$consumers,
    yields_by_age = yields_by_age,
    demand = demand,
    d_ij = d_ij,
    d_jk = d_jk,
    R_p = R_p,
    c_tr_raw = c_tr_raw,
    c_tr_pre = c_tr_pre
  )
}



## ============================================================
## Instance generators for V8 and V9
## ============================================================

generate_instance_v8 <- function(data, params, seed = 123) {
  set.seed(seed)
  
  n_periods <- params$n_periods
  max_age <- params$max_age
  revenue_scale <- params$revenue_scale
  c_tr_raw <- params$c_tr_raw
  c_tr_pre <- params$c_tr_pre
  
  # Filter sites based on selection
  #selected_sites <- data$sites[data$sites$site_id %in% params$selected_sites, ]
  
  n_sites    <- nrow(data$sites)
  n_storages <- nrow(data$storages)
  n_consumers <-  nrow(data$consumers)
  n_periods  <- params$n_periods
  n_products <- params$n_products
  max_age    = params$max_age
  
  
  sites <- data$sites
  
  storages <- data$storages
    
  consumers <- data$consumers
  
  # Distances site→storage
  d_ij <- data$dist_ij
  
  
  # Distances storage→consumer
  d_jk <- data$dist_jk
  
  # Product revenues (p=1 chemical, 2 pulp, 3 energy)
  R_p <- c(
    runif(1, 150, 250),  # chem
    runif(1, 100, 180),  # pulp
    runif(1, 60, 120)    # energy
  ) * revenue_scale
  
  
  # Yield table: η_{p,a}, a = 1..max_age
  yields_by_age <- expand_grid(
    product = 1:n_products,
    age = 1:max_age
  ) %>%
    mutate(
      # crude shape: low at age 1, peak around 5–8, then decline
      base = case_when(
        age <= 3 ~ runif(n(), 0, 0),
        age <= 6 ~ runif(n(), 2, 4),
        age <= 9 ~ runif(n(), 10, 15),
        age <= max_age ~ runif(n(), 8, 10),
        age > max_age ~ 0,
        TRUE ~ runif(n(), 4, 10)
      ),
      # scale by product
      yield_ha = base * case_when(
        product == 1 ~ runif(n(), 0.8, 1.2),
        product == 2 ~ runif(n(), 1.0, 1.5),
        TRUE         ~ runif(n(), 1.5, 2.0)
      )
    )
  
  # Demand: max demand per consumer, product, period
  demand <- expand_grid(
    consumer_id = consumers$consumer_id,
    product = 1:n_products,
    period = 1:n_periods
  ) %>%
    mutate(
      D_max = runif(n(), 200, 800)   # t/yr
    )
  
  list(
    n_sites    = n_sites,
    n_storages = n_storages,
    n_consumers = n_consumers,
    n_periods  = n_periods,
    n_products = n_products,
    max_age    = max_age,
    sites      = sites,
    storages   = storages,
    consumers  = consumers,
    d_ij       = d_ij,
    d_jk       = d_jk,
    R_p        = R_p,
    c_tr_raw   = c_tr_raw,
    c_tr_pre   = c_tr_pre,
    yields_by_age = yields_by_age,
    demand     = demand
  )
}


generate_instance_v8_extended <- function(data, params, yields_by_age ) {
  
  # Extract parameters
  n_periods <- params$n_periods
  max_age <- params$max_age
  c_tr_raw <- params$c_tr_raw
  c_tr_pre <- params$c_tr_pre
  n_products <- params$n_products
  
  # Extract data components
  selected_sites <- data$sites
  storages <- data$storages
  consumers <- data$consumers
  
  n_sites <- nrow(selected_sites)
  n_storages <- nrow(storages)
  n_consumers <- nrow(consumers)
  
  # Get distances
  d_ij <- data$dist_ij  # site to storage
  d_jk <- data$dist_jk  # storage to consumer
  
  # ========================================================================
  # YIELDS TABLE: Use provided yields_by_age or default to internal generation
  # ========================================================================
  
  if (is.null(yields_by_age)) {
    # Fallback: generate default yields (original behavior)
    yields_by_age <- expand.grid(
      product = 1:n_products,
      age = 1:max_age
    ) %>%
      mutate(
        yield_ha = case_when(
          age <= 3 ~ runif(n(), 0, 2),
          age <= 6 ~ runif(n(), 2, 4),
          age <= 9 ~ runif(n(), 10, 15),
          age <= max_age ~ runif(n(), 8, 10),
          TRUE ~ runif(n(), 4, 10)
        )
      ) %>%
      mutate(
        yield_ha = case_when(
          product == 1 ~ yield_ha * runif(n(), 0.8, 1.2),
          product == 2 ~ yield_ha * runif(n(), 1.0, 1.5),
          TRUE ~ yield_ha * runif(n(), 1.5, 2.0)
        )
      )
  } else {
    # Use provided yields_by_age table
    # Ensure it has the right columns: product, age, yield_ha
    if (!all(c("product", "age", "yield_ha") %in% colnames(yields_by_age))) {
      stop("yields_by_age must have columns: product, age, yield_ha")
    }
    
    # Verify consistency with max_age and n_products
    if (max(yields_by_age$age) > max_age) {
      warning(paste("yields_by_age contains ages > max_age; trimming to", max_age))
      yields_by_age <- yields_by_age %>% filter(age <= max_age)
    }
    
    if (max(yields_by_age$product) > n_products) {
      warning(paste("yields_by_age contains products > n_products; trimming to", n_products))
      yields_by_age <- yields_by_age %>% filter(product <= n_products)
    }
  }
  
  # ========================================================================
  # DEMAND: Expand grid for all consumer x product x period combinations
  # ========================================================================
  
  demand <- expand.grid(
    consumer_id = consumers$consumer_id,
    product = 1:n_products,
    period = 1:n_periods
  ) %>%
    left_join(consumers %>% select(consumer_id, demand), by = "consumer_id") %>%
    mutate(D_max = demand) %>%
    select(consumer_id, product, period, D_max)
  
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
  
  # ========================================================================
  # REVENUE VECTOR: Average prices by product (or use defaults)
  # ========================================================================
  
  if (has_prices) {
    # Use average consumer prices for revenue vector
    R_p <- consumer_prices %>%
      group_by(product) %>%
      summarise(avg_price = mean(price, na.rm = TRUE), .groups = "drop") %>%
      arrange(product) %>%
      pull(avg_price)
  } else {
    # Default: set generic prices
    R_p <- c(150, 100, 60)  # Chemical, Pulp, Energy
  }
  
  # ========================================================================
  # BUILD INSTANCE LIST
  # ========================================================================
  
  instance <- list(
    n_sites = n_sites,
    n_storages = n_storages,
    n_consumers = n_consumers,
    n_periods = n_periods,
    n_products = n_products,
    max_age = max_age,
    
    # Entity data
    sites = selected_sites,
    storages = storages,
    consumers = consumers,
    
    # Distance matrices
    d_ij = d_ij,
    d_jk = d_jk,
    
    # Yields and demands
    yields_by_age = yields_by_age,
    demand = demand,
    
    # Prices
    consumer_prices = consumer_prices,
    R_p = R_p,
    
    # Transport costs
    c_tr_raw = c_tr_raw,
    c_tr_pre = c_tr_pre
  )
  
  return(instance)
}


# consumer and product specific demands and prices
generate_instance_v8_final <- function(data, params, yields_by_age ) {
  
  # Extract parameters
  n_periods <- params$n_periods
  max_age <- params$max_age
  min_age <- params$min_age
  c_tr_raw <- params$c_tr_raw
  c_tr_pre <- params$c_tr_pre
  c_opp <- params$c_opp
  n_products <- params$n_products
  
  # Extract data components
  selected_sites <- data$sites
  storages <- data$storages
  consumers <- data$consumers
  
  n_sites <- nrow(selected_sites)
  n_storages <- nrow(storages)
  n_consumers <- nrow(consumers)
  
  # Get distances
  d_ij <- data$dist_ij  # site to storage
  d_jk <- data$dist_jk  # storage to consumer
  
  # ========================================================================
  # YIELDS TABLE: Use provided yields_by_age or default to internal generation
  # ========================================================================
  
  if (is.null(yields_by_age)) {
    # Fallback: generate default yields (original behavior)
    yields_by_age <- expand.grid(
      product = 1:n_products,
      age = 1:max_age
    ) %>%
      mutate(
        yield_ha = case_when(
          age <= 3 ~ runif(n(), 0, 2),
          age <= 6 ~ runif(n(), 2, 4),
          age <= 9 ~ runif(n(), 10, 15),
          age <= max_age ~ runif(n(), 8, 10),
          TRUE ~ runif(n(), 4, 10)
        )
      ) %>%
      mutate(
        yield_ha = case_when(
          product == 1 ~ yield_ha * runif(n(), 0.8, 1.2),
          product == 2 ~ yield_ha * runif(n(), 1.0, 1.5),
          TRUE ~ yield_ha * runif(n(), 1.5, 2.0)
        )
      )
  } else {
    # Use provided yields_by_age table
    # Ensure it has the right columns: product, age, yield_ha
    if (!all(c("product", "age", "yield_ha") %in% colnames(yields_by_age))) {
      stop("yields_by_age must have columns: product, age, yield_ha")
    }
    
    # Verify consistency with max_age and n_products
    if (max(yields_by_age$age) > max_age) {
      warning(paste("yields_by_age contains ages > max_age; trimming to", max_age))
      yields_by_age <- yields_by_age %>% filter(age <= max_age)
    }
    
    if (max(yields_by_age$product) > n_products) {
      warning(paste("yields_by_age contains products > n_products; trimming to", n_products))
      yields_by_age <- yields_by_age %>% filter(product <= n_products)
    }
  }
  
  # ========================================================================
  # DEMAND: Expand grid for all consumer x product x period combinations
  # ========================================================================
  demand_cols <- c("demand_P1", "demand_P2", "demand_P3")
  
  tmp.demand <- consumers %>%
    select(consumer_id, all_of(demand_cols)) %>%
    pivot_longer(
      cols = all_of(demand_cols),
      names_to = "demand_col",
      values_to = "D_max"
    ) %>%
    mutate( 
      product = as.numeric(substring(demand_col, first =nchar(demand_col),  last = nchar(demand_col))) )
  
  demand <- expand.grid(
    consumer_id = consumers$consumer_id,
    product = 1:n_products,
    period = 1:n_periods
  ) %>%
    left_join(
      tmp.demand  , by = c("consumer_id","product")
    )  %>%
    select(consumer_id, product, period, D_max)
  
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
  
  # ========================================================================
  # REVENUE VECTOR: Average prices by product (or use defaults)
  # ========================================================================
  
  if (has_prices) {
    # Use average consumer prices for revenue vector
    R_p <- consumer_prices %>%
      group_by(product) %>%
      summarise(avg_price = mean(price, na.rm = TRUE), .groups = "drop") %>%
      arrange(product) %>%
      pull(avg_price)
  } else {
    # Default: set generic prices
    R_p <- c(150, 100, 60)  # Chemical, Pulp, Energy
  }
  
  # ========================================================================
  # BUILD INSTANCE LIST
  # ========================================================================
  
  instance <- list(
    n_sites = n_sites,
    n_storages = n_storages,
    n_consumers = n_consumers,
    n_periods = n_periods,
    n_products = n_products,
    max_age = max_age,
    min_age = min_age,
    
    # Entity data
    sites = selected_sites,
    storages = storages,
    consumers = consumers,
    
    # Distance matrices
    d_ij = d_ij,
    d_jk = d_jk,
    
    # Yields and demands
    yields_by_age = yields_by_age,
    demand = demand,
    
    # Prices
    consumer_prices = consumer_prices,
    R_p = R_p,
    
    # Transport costs
    c_tr_raw = c_tr_raw,
    c_tr_pre = c_tr_pre,
    
    # opportunity cost (optional, can be set later)
     c_opp = c_opp
  )
  
  return(instance)
}

