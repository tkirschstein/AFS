# OP with consumer-specific revenues per product - OPTIMIZED VERSION (MILP v10)
# Implements speed improvements: data.table, pre-allocation, vectorization, indexing
# Matches MILP_v10: X_jkppt conversion variables, opportunity cost, updated cascade
build_agroforestry_lp_sparse_v10 <- function(instance) {
  
  #------------------------------------------------------------------------
  # INITIALIZATION
  #------------------------------------------------------------------------
  cat("Building sparse LP for agroforestry SCD problem...\n")
  
  # Extract instance parameters
  ns  <- instance$n_sites
  nj  <- instance$n_storages
  nk  <- instance$n_consumers
  Tm  <- instance$n_periods
  P   <- instance$n_products
  Amax <- instance$max_age
  Amin <- instance$min_age
  Copp <- instance$c_opp
  
  # Index sets
  I   <- 1:ns
  J   <- 1:nj
  K   <- 1:nk
  Tset <- 1:Tm
  Pset <- 1:P
  T_ext <- 0:(Tm + 1)
  
  # Precomputed data
  area_ha <- setNames(instance$sites$area_ha, I)
  
  # Yield matrix: compute once [product x age]
  yield_matrix <- matrix(0, nrow = P, ncol = Tm)
  for (p in Pset) {
    for (age in 1:Tm) {
      rows <- instance$yields_by_age[
        instance$yields_by_age$product == p &
          instance$yields_by_age$age == age,
      ]
      if (nrow(rows) > 0) {
        yield_matrix[p, age] <- rows$yield_ha[1]
      }
    }
  }
  
  yield_max <- apply(yield_matrix,1,max)
  
  cat(sprintf("  Dimensions: %d sites, %d storages, %d consumers, %d periods, %d products\n",
              ns, nj, nk, Tm, P))
  
  #------------------------------------------------------------------------
  # STEP 1: CREATE VARIABLE INDEXING SCHEME
  #------------------------------------------------------------------------
  cat("Step 1: Indexing variables...\n")
  
  # Build arcs: all pairs (s,t) where s < t in extended time domain
  T_ext_pairs <- which(outer(T_ext, T_ext, "<"), arr.ind = TRUE)
  colnames(T_ext_pairs) <- c("s_idx", "t_idx")
  
  # z[i, s, t] for arc variables
  z_tuples <- expand.grid(i = I, row = seq_len(nrow(T_ext_pairs)))
  z_tuples$s <- T_ext[T_ext_pairs[z_tuples$row, "s_idx"]]
  z_tuples$t <- T_ext[T_ext_pairs[z_tuples$row, "t_idx"]]
  z_tuples$row <- NULL
  z_tuples$ub <- 1
  # set ub to 0 for  s > t
  z_tuples$ub[z_tuples$s > z_tuples$t] <- 0
  # set ub to 0 for  (t-s) < Amin or (t-s) > Amax
  z_tuples$ub[z_tuples$s > 0 & ((z_tuples$t - z_tuples$s) < Amin | (z_tuples$t - z_tuples$s) > Amax )] <- 0
  
  # drop variable with ub=0  (invalid arcs)
  z_tuples <- z_tuples[z_tuples$ub > 0, ]
  z_tuples$col <- seq_len(nrow(z_tuples))
  n_z <- nrow(z_tuples)
  
  # Y[i, p, t] - harvest quantity upper bounds
  # periods with potential harvests
  Tharv <- Tset[Tset > max(1,Amin)]
  
  Y_tuples <- expand.grid(i = I, p = Pset, t = Tharv)
  # set ub using yields and area
  Y_tuples$ub <- area_ha[as.character(Y_tuples$i)] * yield_max[Y_tuples$p]
  
  Y_tuples$col <- n_z + seq_len(nrow(Y_tuples))
  
  n_Y <- nrow(Y_tuples)
  
  # X_ij[i, j, p, t] - flows site to hub
  Xij_tuples <- expand.grid(i = I, j = J, p = Pset, t = Tharv)
  Xij_tuples$col <- max(Y_tuples$col) + seq_len(nrow(Xij_tuples))
  # set ub to inf
  Xij_tuples$ub <- Inf
  n_Xij <- nrow(Xij_tuples)
  
  # S_jpt[j, p, t] - inventory at storage hubs
  S_tuples <- expand.grid(j = J, p = Pset, t = Tset)
  S_tuples$col <- max(Xij_tuples$col) + seq_len(nrow(S_tuples))
  S_tuples$ub <- instance$storages$CAP_stor[S_tuples$j]
  n_S <- nrow(S_tuples)
  
  # X_jk[j, k, p, p', t] - flows hub to consumer
  Xjk_tuples <- expand.grid(j = J, k = K, p = Pset, pp = Pset, t = Tharv)
  # set ub to inf
  Xjk_tuples$ub <- Inf
  # set ub to 0 for pp < p (cascading demand)
  Xjk_tuples$ub[Xjk_tuples$pp < Xjk_tuples$p] <- 0
  
  # drop variables with ub=0 (invalid flows)
  Xjk_tuples <- Xjk_tuples[Xjk_tuples$ub > 0, ]
  
  Xjk_tuples$col <- max(S_tuples$col) + seq_len(nrow(Xjk_tuples))
  n_Xjk <- nrow(Xjk_tuples)
  
  n_vars <- max(Xjk_tuples$col)
  
  cat(sprintf("Total variables: %d (z:%d Y:%d Xij:%d S:%d Xjk:%d)\n",
              n_vars, n_z, n_Y, n_Xij, n_S, n_Xjk))
  
  #------------------------------------------------------------------------
  # STEP 2: BUILD OBJECTIVE VECTOR
  #------------------------------------------------------------------------
  cat("Step 2: Building objective vector...\n")
  
  c_vec <- numeric(n_vars)
  
  
  #consumer prices
  consumer_prices <- instance$consumer_prices
  # rename columns
  colnames(consumer_prices) <- c("k", "pp", "price")
  
  # find prices for each (k,p) pair
  revenue.vec <- Xjk_tuples %>%
    left_join(consumer_prices, by = c("k", "pp")) %>%
    mutate(col_idx = col) %>%
    select(col_idx, price) 
  
  # Revenue from shipments to consumers: R_p[k,p] * X_jk[j,k,p,t]
  #c_vec[Xjk_tuples$col] <- instance$R_p[Xjk_tuples$p]
  c_vec[revenue.vec$col_idx] <- revenue.vec$price
  
  # Establishment cost: -C_est[i] * z[i,0,t]
  z_est <- z_tuples[z_tuples$s == 0 & z_tuples$t %in% Tset, ]
  c_vec[z_est$col] <- c_vec[z_est$col] - instance$sites$C_est[z_est$i]
  
  # opportunity cost: -C_opp[i] * area_ha[i] * sum_t z[i,0,t] * (T-t) for t in Tset
  z_opp <- z_tuples[z_tuples$s == 0 & z_tuples$t %in% Tset, ]
  opp_cost <- Copp * area_ha[as.character(z_opp$i)] * (Tm - z_opp$t)
  c_vec[z_opp$col] <- c_vec[z_opp$col] - opp_cost
  
  # Harvest cost: -C_harv[i] * area_ha[i] * z[i,s,t] for s >= 1
  z_harv <- z_tuples[z_tuples$s >= 1 & z_tuples$t %in% Tset & z_tuples$t > z_tuples$s, ]
  harv_cost <- instance$sites$C_harv[z_harv$i] * area_ha[as.character(z_harv$i)]
  c_vec[z_harv$col] <- c_vec[z_harv$col] - harv_cost
  
  # Transport site->hub: -c_tr_raw * d_ij[i,j] * X_ij
  d_ij_vals <- numeric(nrow(Xij_tuples))
  for (k in seq_len(nrow(Xij_tuples))) {
    d_ij_vals[k] <- instance$d_ij[Xij_tuples$i[k], Xij_tuples$j[k]]
  }
  c_vec[Xij_tuples$col] <- c_vec[Xij_tuples$col] - instance$c_tr_raw * d_ij_vals
  
  # Transport hub->consumer: -c_tr_pre * d_jk[j,k] * X_jk
  d_jk_vals <- numeric(nrow(Xjk_tuples))
  for (k in seq_len(nrow(Xjk_tuples))) {
    d_jk_vals[k] <- instance$d_jk[Xjk_tuples$j[k], Xjk_tuples$k[k]]
  }
  c_vec[Xjk_tuples$col] <- c_vec[Xjk_tuples$col] - instance$c_tr_pre * d_jk_vals
  
  # Storage cost: -c_stor[j] * S_jpt
  c_vec[S_tuples$col] <- c_vec[S_tuples$col] - instance$storages$c_stor[S_tuples$j]
  
  cat(sprintf("  Objective: %d nonzeros\n", sum(c_vec != 0)))
  
  #------------------------------------------------------------------------
  # STEP 3: BUILD CONSTRAINTS (stored as structured list)
  #------------------------------------------------------------------------
  cat("Step 3: Building constraints...\n")
  
  const_count <- 0
  
  row_idx <- numeric()
  col_idx <- numeric()
  val     <- numeric()
  rhs_vec <- numeric()
  sense_vec <- character()
  
  
  # ====== C1: Path establishment (at most one path per site) ======
  # sum_t z_{i,0,t} <= 1 for i in I
  for (i in I) {
    z_est_i <- z_tuples[z_tuples$i == i & z_tuples$s == 0 & z_tuples$t %in% Tset, ]
    if (nrow(z_est_i) > 0) {
      
      const_count <- const_count + 1L
      
      row_idx <- c(row_idx, rep(const_count , nrow(z_est_i)))
      col_idx <- c(col_idx, z_est_i$col)
      val     <- c(val, rep(1, nrow(z_est_i)))
      sense_vec <- c(sense_vec, "<=")
      rhs_vec <- c(rhs_vec, 1)
      
      
    }
  }
  cat("  C1: Path establishment\n")
  
  # ====== C2: Path connectivity ======
  # For each site i and time t: sum_{s<t} z_{i,s,t} = sum_{u>t} z_{i,t,u}
  for (i in I) {
    for (t in Tset) {
      z_in <- z_tuples[z_tuples$i == i & z_tuples$s < t & z_tuples$t == t, ]
      z_out <- z_tuples[z_tuples$i == i & z_tuples$s == t & z_tuples$t > t, ]
      
      if (nrow(z_in) > 0 || nrow(z_out) > 0) {
        col_all <- c(z_in$col, z_out$col)
        coef_all <- c(rep(1, nrow(z_in)), rep(-1, nrow(z_out)))
        
        const_count <- const_count + 1L
        row_idx <- c(row_idx, rep(const_count , length(col_all)))
        col_idx <- c(col_idx, col_all)
        val     <- c(val, coef_all)
        sense_vec <- c(sense_vec, "==")
        rhs_vec <- c(rhs_vec, 0)
        
      }
    }
  }
  cat("  C2: Path connectivity\n")
  
  # ====== C3a: Biomass yield constraint ======
  # Y_{i,p,t} <= sum_{s=1}^{t-1} yield_matrix[p, t-s] * area_ha[i] * z_{i,s,t}
  for (i in I) {
    for (p in Pset) {
      for (t in Tharv) {
          # t >= 2
          Y_col <- Y_tuples$col[Y_tuples$i == i & Y_tuples$p == p & Y_tuples$t == t]
          z_arcs <- z_tuples[z_tuples$i == i & z_tuples$t == t & 
                               z_tuples$s >= 1 & z_tuples$s < t, ]
          
          if (length(Y_col) > 0 && nrow(z_arcs) > 0) {
            col_all <- c(Y_col, z_arcs$col)
            age_vec <- t - z_arcs$s
            coef_all <- c(1, -yield_matrix[p, age_vec] * area_ha[i])
            
            const_count <- const_count + 1L
            row_idx <- c(row_idx, rep(const_count , length(col_all)))
            col_idx <- c(col_idx, col_all)
            val     <- c(val, coef_all)
            sense_vec <- c(sense_vec, "<=")
            rhs_vec <- c(rhs_vec, 0)
            
          }
        }
    }
  }
  cat("  C3a: Biomass yield\n")
  
  # ====== C4: Site->Hub flow limited by harvest ======
  # sum_j X_{ijpt} <= Y_{ipt}
  for (i in I) {
    for (p in Pset) {
      for (t in Tharv) {
        Y_col <- Y_tuples$col[Y_tuples$i == i & Y_tuples$p == p & Y_tuples$t == t]
        X_cols <- Xij_tuples$col[Xij_tuples$i == i & Xij_tuples$p == p & 
                                   Xij_tuples$t == t]
        
        if (length(Y_col) > 0 && length(X_cols) > 0) {
          col_all <- c(X_cols, Y_col)
          coef_all <- c(rep(1, length(X_cols)), -1)
          
          const_count <- const_count + 1L
          row_idx <- c(row_idx, rep(const_count , length(col_all)))
          col_idx <- c(col_idx, col_all)
          val     <- c(val, coef_all)
          sense_vec <- c(sense_vec, "<=")
          rhs_vec <- c(rhs_vec, 0)
          
        }
      }
    }
  }
  cat("  C4: Flow constraints\n")
  
  # ====== C6: Inventory balance at storage hubs ======
  # S_{j,p,t} = S_{j,p,t-1} + sum_i X_{ijpt} - sum_k sum_{p'>=p} X_{jkpp't}
  for (j in J) {
    for (p in Pset) {
      # t = 1
      S_col_1 <- S_tuples$col[S_tuples$j == j & S_tuples$p == p & S_tuples$t == 1]
      X_in_1 <- Xij_tuples$col[Xij_tuples$j == j & Xij_tuples$p == p & 
                                 Xij_tuples$t == 1]
      X_out_1 <- Xjk_tuples$col[Xjk_tuples$j == j & Xjk_tuples$p == p & Xjk_tuples$pp >= p & 
                                  Xjk_tuples$t == 1]
      
      if (length(S_col_1) > 0) {
        col_all <- c(S_col_1, X_in_1, X_out_1)
        coef_all <- c(1, rep(-1, length(X_in_1)), rep(1, length(X_out_1)))
        
        
        const_count <- const_count + 1L
        row_idx <- c(row_idx, rep(const_count , length(col_all)))
        col_idx <- c(col_idx, col_all)
        val     <- c(val, coef_all)
        sense_vec <- c(sense_vec, "==")
        rhs_vec <- c(rhs_vec, 0)
        
      }
      
      # t >= 2
      for (t in Tset[Tset >= 2]) {
        S_col_t <- S_tuples$col[S_tuples$j == j & S_tuples$p == p & 
                                  S_tuples$t == t]
        S_col_tm1 <- S_tuples$col[S_tuples$j == j & S_tuples$p == p & 
                                    S_tuples$t == t-1]
        X_in_t <- Xij_tuples$col[Xij_tuples$j == j & Xij_tuples$p == p & 
                                   Xij_tuples$t == t]
        X_out_t <- Xjk_tuples$col[Xjk_tuples$j == j & Xjk_tuples$p == p & Xjk_tuples$pp >= p &
                                    Xjk_tuples$t == t]
        
        if (length(S_col_t) > 0) {
          col_all <- c(S_col_t, S_col_tm1, X_in_t, X_out_t)
          coef_all <- c(1, -1, rep(-1, length(X_in_t)), rep(+1, length(X_out_t)))
          
          const_count <- const_count + 1L
          row_idx <- c(row_idx, rep(const_count , length(col_all)))
          col_idx <- c(col_idx, col_all)
          val     <- c(val, coef_all)
          sense_vec <- c(sense_vec, "==")
          rhs_vec <- c(rhs_vec, 0)
          
        }
      }
    }
  }
  cat("  C6: Inventory balance\n")
  
  # ====== C7: Storage capacity ======
  # sum_p S_{jpt} <= CAP_stor_j
  for (j in J) {
    for (t in Tset) {
      S_cols <- S_tuples$col[S_tuples$j == j & S_tuples$t == t]
      if (length(S_cols) > 0) {
        
        const_count <- const_count + 1L
        row_idx <- c(row_idx, rep(const_count , length(S_cols)))
        col_idx <- c(col_idx, S_cols)
        val     <- c(val, rep(1, length(S_cols)))
        sense_vec <- c(sense_vec, "<=")
        rhs_vec <- c(rhs_vec, instance$storages$CAP_stor[j])
        
        
      }
    }
  }
  cat("  C7: Storage capacity\n")
  
  # ====== C8: Processing capacity ======
  # sum_{i,p} X_{ijpt} <= CAP_proc_j
  for (j in J) {
    for (t in Tset) {
      X_cols <- Xij_tuples$col[Xij_tuples$j == j & Xij_tuples$t == t]
      if (length(X_cols) > 0) {
        
        const_count <- const_count + 1L
        row_idx <- c(row_idx, rep(const_count , length(X_cols)))
        col_idx <- c(col_idx, X_cols)
        val     <- c(val, rep(1, length(X_cols)))
        sense_vec <- c(sense_vec, "<=")
        rhs_vec <- c(rhs_vec, instance$storages$CAP_proc[j])
        
      }
    }
  }
  cat("  C8: Processing capacity\n")
  
  # ====== C9: Demand with cascade ======
  # sum_j sum_{q <= p} X_{jkqpt} <= D_max(k,p,t)
  for (k in K) {
    for (p in Pset) {
      for (t in Tharv) {
        # Find demand
        dem_row <- instance$demand[
          instance$demand$consumer_id == k &
            instance$demand$product == p &
            instance$demand$period == t,
        ]
        if (nrow(dem_row) == 0) next
        D_max_kpt <- dem_row$D_max[1]
        
        # Get X_jk variables with product q <= p
        X_cols <- Xjk_tuples$col[
          Xjk_tuples$k == k &
            Xjk_tuples$t == t &
            Xjk_tuples$p <= p &
            Xjk_tuples$pp == p
        ]
        
        if (length(X_cols) > 0) {
          
          const_count <- const_count + 1L
          row_idx <- c(row_idx, rep(const_count , length(X_cols)))
          col_idx <- c(col_idx, X_cols)
          val     <- c(val, rep(1, length(X_cols)))
          sense_vec <- c(sense_vec, "<=")
          rhs_vec <- c(rhs_vec, D_max_kpt)
          
        }
      }
    }
  }
  cat("  C9: Demand satisfaction\n")
  
  # ====== C10: Drop blind flows ======
  # sum_t sum_j sum_{q <= p} X_{jkpqt} == 0
  
  
  
  cat(sprintf("  Total constraints: %d\n", const_count))
  
  
  #------------------------------------------------------------------------
  # STEP 4: BUILD ROI OP OBJECT
  #------------------------------------------------------------------------
  cat("Step 4: Building ROI OP object...\n")
  
  A <- simple_triplet_matrix(
    i = row_idx,
    j = col_idx,
    v = val,
    nrow = const_count,
    ncol = n_vars
  )
  
  model <- OP(
    objective   = L_objective(c_vec),
    constraints = L_constraint(
      L = A,
      dir = sense_vec,
      rhs = rhs_vec
    ),
    maximum    = TRUE,
    bounds     = V_bound(
      li = seq_len(n_vars),
      ui = seq_len(n_vars),
      lb = rep(0, n_vars),
      ub = c(z_tuples$ub, Y_tuples$ub, Xij_tuples$ub, S_tuples$ub, Xjk_tuples$ub),
      #ub = rep(Inf, n_vars),
      nobj = n_vars
    ),
    types      = c(rep("B", n_z),
                   rep("C", n_Y + n_Xij + n_Xjk + n_S))
  )
  
  cat(sprintf("  Sparse matrix: %d x %d, %d nonzeros (%.2f%% density)\n",
              nrow(A), ncol(A), nnzero(A), 
              100*nnzero(A)/(nrow(A)*ncol(A))))
  
  
  
  cat("OP model built successfully!\n\n")
  
  # Return everything needed for solving
  list(
    model = model,
    var_maps = list(
      z = z_tuples,
      Y = Y_tuples,
      Xij = Xij_tuples,
      Xjk = Xjk_tuples,
      S = S_tuples
    ),
    instance_info = list(
      n_vars = n_vars,
      n_constrs = const_count,
      n_z = n_z,
      n_Y = n_Y,
      n_Xij = n_Xij,
      n_Xjk = n_Xjk,
      n_S = n_S
    )
  )
}


