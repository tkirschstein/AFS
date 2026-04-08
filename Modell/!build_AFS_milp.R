build_agroforestry_lp_sparse_v10_optimized <- function(instance) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install with: install.packages('data.table')")
  }
  library(data.table)
  library(slam)
  library(ROI)
  
  cat("Building OPTIMIZED sparse LP for agroforestry SCD problem...\n")
  
  ## -------------------------------------------------------------------
  ## INITIALIZATION
  ## -------------------------------------------------------------------
  ns   <- instance$n_sites
  nj   <- instance$n_storages
  nk   <- instance$n_consumers
  Tm   <- instance$n_periods
  P    <- instance$n_products
  Amax <- instance$max_age
  Amin <- instance$min_age
  Copp <- instance$c_opp
  
  I        <- 1:ns
  J        <- 1:nj
  K        <- 1:nk
  Tset     <- 1:Tm
  Pset     <- 1:P
  T_ext    <- 0:(Tm + 1)
  area_ha  <- setNames(instance$sites$area_ha, I)
  
  # Yield matrix [product x age]
  yield_matrix <- matrix(0, nrow = P, ncol = Tm)
  for (pprod in Pset) {
    for (age_idx in 1:Tm) {
      rows <- instance$yields_by_age[
        instance$yields_by_age$product == pprod &
          instance$yields_by_age$age    == age_idx,
      ]
      if (nrow(rows) > 0) {
        yield_matrix[pprod, age_idx] <- rows$yield_ha[1]
      }
    }
  }
  yield_max <- apply(yield_matrix, 1, max)
  
  cat(sprintf("  Dimensions: %d sites, %d storages, %d consumers, %d periods, %d products\n",
              ns, nj, nk, Tm, P))
  
  ## -------------------------------------------------------------------
  ## STEP 1: VARIABLE INDEXING (data.table + filtering)
  ## -------------------------------------------------------------------
  cat("Step 1: Indexing variables...\n")
  
  # Arcs (s,t)
  T_ext_pairs <- which(outer(T_ext, T_ext, "<"), arr.ind = TRUE)
  colnames(T_ext_pairs) <- c("s_idx", "t_idx")
  
  # z[ii,s,t]
  z_tuples <- as.data.table(expand.grid(ii = I, row = seq_len(nrow(T_ext_pairs))))
  z_tuples[, `:=`(
    s   = T_ext[T_ext_pairs[row, "s_idx"]],
    t   = T_ext[T_ext_pairs[row, "t_idx"]],
    ub  = 1
  )]
  z_tuples[, row := NULL]
  
  # infeasible arcs: s>t or age outside [Amin,Amax] (for s>0)
  z_tuples[s > t, ub := 0]
  z_tuples[s > 0 & ((t - s) < Amin | (t - s) > Amax), ub := 0]
  
  # drop invalid arcs
  z_tuples <- z_tuples[ub > 0]
  z_tuples[, col := .I]
  n_z <- nrow(z_tuples)
  
  # Harvest periods
  Tharv <- Tset[Tset > max(1, Amin)]
  
  # Xij[ii,jj,pprod,tt]
  Xij_tuples <- as.data.table(expand.grid(ii = I, jj = J, pprod = Pset, tt = Tharv))
  Xij_tuples[, `:=`(ub = Inf,
                    col = .I + n_z + n_Y)]
  n_Xij <- nrow(Xij_tuples)
  
  # S[jj,pprod,tt]
  S_tuples <- as.data.table(expand.grid(jj = J, pprod = Pset, tt = Tharv))
  S_tuples[, `:=`(
    ub  = instance$storages$CAP_stor[jj],
    col = .I + n_z + n_Y + n_Xij
  )]
  n_S <- nrow(S_tuples)
  
  # Xjk[jj,kk,pprod,pp,tt]
  Xjk_tuples <- as.data.table(
    expand.grid(jj = J, kk = K, pprod = Pset, pp = Pset, tt = Tharv)
  )
  Xjk_tuples[, ub := Inf]
  Xjk_tuples[pp < pprod, ub := 0]
  Xjk_tuples <- Xjk_tuples[ub > 0]
  Xjk_tuples[, col := .I + n_z + n_Y + n_Xij + n_S]
  n_Xjk <- nrow(Xjk_tuples)
  
  n_vars <- max(Xjk_tuples$col)
  
  cat(sprintf("  Total variables: %d (z:%d Xij:%d S:%d Xjk:%d)\n",
              n_vars, n_z, n_Xij, n_S, n_Xjk))
  
  # keys for fast lookup
  setkey(z_tuples,  ii, s, t)
  setkey(Xij_tuples,ii, jj, pprod, tt)
  setkey(S_tuples,  jj, pprod, tt)
  setkey(Xjk_tuples,jj, kk, pprod, pp, tt)
  
  ## -------------------------------------------------------------------
  ## STEP 2: OBJECTIVE (vectorized)
  ## -------------------------------------------------------------------
  cat("Step 2: Building objective vector...\n")
  
  c_vec <- numeric(n_vars)
  
  # Revenue: consumer_prices (k, pp, price)
  consumer_prices <- as.data.table(instance$consumer_prices)
  setnames(consumer_prices, c("kk", "pp", "price"))
  setkey(consumer_prices, kk, pp)
  
  Xjk_tuples[consumer_prices, price := i.price, on = .(kk, pp)]
  c_vec[Xjk_tuples$col] <- Xjk_tuples$price
  
  # Establishment cost: -C_est * z[ii,0,t]
  z_est <- z_tuples[s == 0 & t %in% Tset]
  if (nrow(z_est) > 0) {
    c_vec[z_est$col] <- c_vec[z_est$col] - instance$sites$C_est[z_est$ii]
  }
  
  # Opportunity cost: -C_opp * area_ha[ii] * (Tmax - t) * z[ii,0,t]
  z_opp <- z_tuples[s == 0 & t %in% Tset]
  if (nrow(z_opp) > 0) {
    opp_cost <- Copp * area_ha[as.character(z_opp$ii)] * (Tm - z_opp$t)
    c_vec[z_opp$col] <- c_vec[z_opp$col] - opp_cost
  }
  
  # Harvest cost: -C_harv * area_ha[ii] * z[ii,s,t] for s>=1
  z_harv <- z_tuples[s >= 1 & t %in% Tset & t > s]
  if (nrow(z_harv) > 0) {
    harv_cost <- instance$sites$C_harv[z_harv$ii] * area_ha[as.character(z_harv$ii)]
    c_vec[z_harv$col] <- c_vec[z_harv$col] - harv_cost
  }
  
  # Transport site->hub: vectorized
  d_ij_vals <- instance$d_ij[cbind(Xij_tuples$ii, Xij_tuples$jj)]
  c_vec[Xij_tuples$col] <- c_vec[Xij_tuples$col] - instance$c_tr_raw * d_ij_vals
  
  # Transport hub->consumer: vectorized
  d_jk_vals <- instance$d_jk[cbind(Xjk_tuples$jj, Xjk_tuples$kk)]
  c_vec[Xjk_tuples$col] <- c_vec[Xjk_tuples$col] - instance$c_tr_pre * d_jk_vals
  
  # Storage cost
  c_vec[S_tuples$col] <- c_vec[S_tuples$col] - instance$storages$c_stor[S_tuples$jj]
  
  cat(sprintf("  Objective: %d nonzeros\n", sum(c_vec != 0)))
  
  ## -------------------------------------------------------------------
  ## STEP 3: CONSTRAINTS (pre-allocated)
  ## -------------------------------------------------------------------
  cat("Step 3: Building constraints (pre-allocated)...\n")
  
  est_constraints <- ns + ns * Tm + ns * P * length(Tharv) * 2 +
    ns * P * length(Tharv) + nj * P * Tm * 2 + nj * Tm * 2 + nk * P * length(Tharv)
  est_nonzeros <- est_constraints * 25
  
  row_idx   <- integer(est_nonzeros)
  col_idx   <- integer(est_nonzeros)
  val       <- numeric(est_nonzeros)
  rhs_vec   <- numeric(est_constraints)
  sense_vec <- character(est_constraints)
  
  const_count <- 0L
  nnz_ptr     <- 0L
  
  add_constraint <- function(cols, coeffs, sense, rhs) {
    n_new <- length(cols)
    if (n_new == 0) return()
    const_count <<- const_count + 1L
    
    if (nnz_ptr + n_new > length(row_idx)) {
      new_size <- length(row_idx) * 2L
      length(row_idx) <<- new_size
      length(col_idx) <<- new_size
      length(val)     <<- new_size
    }
    
    idx <- (nnz_ptr + 1L):(nnz_ptr + n_new)
    row_idx[idx] <<- const_count
    col_idx[idx] <<- cols
    val[idx]     <<- coeffs
    nnz_ptr <<- nnz_ptr + n_new
    
    if (const_count > length(rhs_vec)) {
      new_size <- length(rhs_vec) * 2L
      length(rhs_vec)   <<- new_size
      length(sense_vec) <<- new_size
    }
    sense_vec[const_count] <<- sense
    rhs_vec[const_count]   <<- rhs
  }
  
  ## C1: Path establishment
  for (ii_idx in I) {
    z_est_i <- z_tuples[ii== ii_idx & s == 0]
    if (nrow(z_est_i) > 0) {
      add_constraint(z_est_i$col, rep(1, nrow(z_est_i)), "<=", 1)
    }
  }
  cat("  C1: Path establishment\n")
  
  ## C2: Path connectivity
  for (ii_idx in I) {
    zi <- z_tuples[.(ii_idx)]
    if (nrow(zi) == 0) next
    for (tt_idx in Tset) {
      z_in  <- zi[s < tt_idx & t == tt_idx]
      z_out <- zi[s == tt_idx & t > tt_idx]
      if (nrow(z_in) > 0 || nrow(z_out) > 0) {
        col_all  <- c(z_in$col, z_out$col)
        coef_all <- c(rep(1, nrow(z_in)), rep(-1, nrow(z_out)))
        add_constraint(col_all, coef_all, "==", 0)
      }
    }
  }
  cat("  C2: Path connectivity\n")
  
  ## C3: Biomass yield 
  for (ii_idx in I) {
    area_i <- area_ha[ii_idx]
    for (pp_idx in Pset) {
      for (tt_idx in Tharv) {
        X_row <- Xij_tuples[ii ==  ii_idx & pprod == pp_idx & tt == tt_idx]
        if (nrow(X_row) == 0) next
        X_col  <- X_row$col
        z_arcs <- z_tuples[ii == ii_idx & t == tt_idx & s >= 1 & s < tt_idx]
        if (nrow(z_arcs) > 0) {
          age_vec  <- tt_idx - z_arcs$s
          col_all  <- c(X_col, z_arcs$col)
          coef_all <- c(1, -yield_matrix[pprod, age_vec] * area_i)
          add_constraint(col_all, coef_all, "<=", 0)
        }
      }
    }
  }
  cat("  C3: Biomass yield\n")
  
  ## C6: Inventory balance
  for (jj_idx in J) {
    for (p in Pset) {
      # t = Amin+1 (first harvest period)
      S_row_1 <- S_tuples[.(jj_idx, p, Amin+1), nomatch = 0]
      if (nrow(S_row_1) > 0) {
        S_col_1 <- S_row_1$col
        X_in_1  <- Xij_tuples[jj == jj_idx & pprod == p & tt == Amin+1]$col
        X_out_1 <- Xjk_tuples[jj == jj_idx & pprod == p & pp >= p & tt == Amin+1]$col
        col_all  <- c(S_col_1, X_in_1, X_out_1)
        coef_all <- c(1, rep(-1, length(X_in_1)), rep(1, length(X_out_1)))
        add_constraint(col_all, coef_all, "==", 0)
      }
      # t >= Amin + 1
      for (tt_idx in Tharv[Tharv >= Amin+1]) {
        S_row_t   <- S_tuples[.(jj_idx, p, tt_idx), nomatch = 0]
        if (nrow(S_row_t) == 0) next
        S_col_t   <- S_row_t$col
        S_col_tm1 <- S_tuples[.(jj_idx, p, tt_idx - 1), nomatch = 0]$col
        X_in_t    <- Xij_tuples[jj == jj_idx & pprod == p & tt == tt_idx]$col
        X_out_t   <- Xjk_tuples[jj == jj_idx & pprod == p & pp >= p & tt == tt_idx]$col
        col_all  <- c(S_col_t, S_col_tm1, X_in_t, X_out_t)
        coef_all <- c(1, -1, rep(-1, length(X_in_t)), rep(1, length(X_out_t)))
        add_constraint(col_all, coef_all, "==", 0)
      }
    }
  }
  cat("  C6: Inventory balance\n")
  
  ## C7: Storage capacity
  for (jj_idx in J) {
    for (tt_idx in Tharv) {
      S_cols <- S_tuples[jj == jj_idx & tt == tt_idx]$col
      if (length(S_cols) > 0) {
        add_constraint(S_cols, rep(1, length(S_cols)), "<=", instance$storages$CAP_stor[jj_idx])
      }
    }
  }
  cat("  C7: Storage capacity\n")
  
  ## C8: Processing capacity
  for (jj_idx in J) {
    for (tt_idx in Tharv) {
      X_cols <- Xij_tuples[jj == jj_idx & tt == tt_idx]$col
      if (length(X_cols) > 0) {
        add_constraint(X_cols, rep(1, length(X_cols)),
                       "<=", instance$storages$CAP_proc[jj_idx])
      }
    }
  }
  cat("  C8: Processing capacity\n")
  
  ## C9: Demand with cascade
  demand_dt <- as.data.table(instance$demand)
  setnames(demand_dt, c("consumer_id", "product", "period", "D_max"),
           c("kk", "pp", "tt", "D_max"))
  setkey(demand_dt, kk, pp, tt)
  
  for (kk_idx in K) {
    for (p in Pset) {  # demanded product
      for (tt_idx in Tharv) {
        dem_row <- demand_dt[.(kk_idx, pprod, tt_idx), nomatch = 0]
        if (nrow(dem_row) == 0) next
        D_max_kpt <- dem_row$D_max[1]
        Xk <- Xjk_tuples[kk == kk_idx & tt == tt_idx & pp == p & pprod <= p]$col
        if (length(Xk) > 0) {
          add_constraint(Xk, rep(1, length(Xk)), "<=", D_max_kpt)
        }
      }
    }
  }
  cat("  C9: Demand satisfaction\n")
  
  # Trim preallocated
  row_idx   <- row_idx[1:nnz_ptr]
  col_idx   <- col_idx[1:nnz_ptr]
  val       <- val[1:nnz_ptr]
  rhs_vec   <- rhs_vec[1:const_count]
  sense_vec <- sense_vec[1:const_count]
  
  cat(sprintf("  Total constraints: %d\n", const_count))
  
  ## -------------------------------------------------------------------
  ## STEP 4: ROI MODEL
  ## -------------------------------------------------------------------
  cat("Step 4: Building ROI OP object...\n")
  
  A <- simple_triplet_matrix(
    i = row_idx,
    j = col_idx,
    v = val,
    nrow = const_count,
    ncol = n_vars
  )
  
  ub_vec <- c(z_tuples$ub, Y_tuples$ub, Xij_tuples$ub, S_tuples$ub, Xjk_tuples$ub)
  
  model <- OP(
    objective   = L_objective(c_vec),
    constraints = L_constraint(
      L   = A,
      dir = sense_vec,
      rhs = rhs_vec
    ),
    maximum  = TRUE,
    bounds   = V_bound(
      li = seq_len(n_vars),
      ui = seq_len(n_vars),
      lb = rep(0, n_vars),
      ub = ub_vec,
      nobj = n_vars
    ),
    types    = c(rep("B", n_z),
                 rep("C", n_Y + n_Xij + n_S + n_Xjk))
  )
  
  cat(sprintf("  Sparse matrix: %d x %d, %d nonzeros (%.2f%% density)\n",
              nrow(A), ncol(A), nnzero(A),
              100 * nnzero(A) / (nrow(A) * ncol(A))))
  
  list(
    model = model,
    var_maps = list(
      z   = z_tuples,
      Xij = Xij_tuples,
      Xjk = Xjk_tuples,
      S   = S_tuples
    ),
    instance_info = list(
      n_vars   = n_vars,
      n_constrs = const_count,
      n_z      = n_z,
      n_Xij    = n_Xij,
      n_Xjk    = n_Xjk,
      n_S      = n_S
    )
  )
}
