library(ompr)
library(gurobi)
library(ompr.roi)
library(ROI.plugin.glpk)  # or gurobi/xpress plugin you use
library(ROI.plugin.gurobi)  # or gurobi/xpress plugin you use


solve_v8a_ompr_gurobi <- function(instance, 
                           time_limit = 600,
                           mip_gap = 0.05,
                           verbose = TRUE) {
  
  ns <- instance$n_sites
  nj <- instance$n_storages
  nk <- instance$n_consumers
  Tm <- instance$n_periods
  P  <- instance$n_products
  Amax <- instance$max_age
  
  sites    <- instance$sites
  storages <- instance$storages
  demand   <- instance$demand
  yields   <- instance$yields_by_age
  d_ij     <- instance$d_ij
  d_jk     <- instance$d_jk
  R_p      <- instance$R_p   # scale revenues
  c_tr_raw <- instance$c_tr_raw
  c_tr_pre <- instance$c_tr_pre
  
  
  I <- 1:ns
  J <- 1:nj
  K <- 1:nk
  Tset <- 1:Tm
  Pset <- 1:P
  
  # define extended time set for arcs: 0..Tmax+1
  T_ext <- 0:(Tm + 1)
  
  # helper lookups
  area_ha <- setNames(sites$area_ha, sites$site_id)
  
  # yield lookup: returns yield_ha for (p, age), or 0 if missing
  get_yield <- function(p, age) {
    if (age < 1L || age > Tm) return(0)
    rows <- instance$yields_by_age[
      instance$yields_by_age$product == p &
        instance$yields_by_age$age == age,
    ]
    if (nrow(rows) == 0) 0 else rows$yield_ha[1]
  }
  
  # create yield matrix p x age for faster access
  yield_matrix <- matrix(0, nrow = P, ncol = Tm)
  for (p in Pset) {
    for (age in 1:Tm) {
      yield_matrix[p, age] <- get_yield(p, age)
    }
  }
  
  
  
  
  model <- MIPModel() %>%
    # 1. Decision variables ----------------------------------------------------
  
  # Path / arc variables z_{i,s,t}, s,t in {0..Tmax+1}, with s < t and t in 1..Tmax+1
  add_variable(z[i, s, t],
               i = I, s = T_ext, t = T_ext, s < t,
               type = "binary",
               lb = 0, ub = 1
  ) %>%
    
    # Harvest quantity upper bound Y_{i,p,t}
    add_variable(Y[i, p, t],
                 i = I, p = Pset, t = Tset,
                 type = "continuous", lb = 0
    ) %>%
    
    # Flows site -> hub
    add_variable(X_ij[i, j, p, t],
                 i = I, j = J, p = Pset, t = Tset,
                 type = "continuous", lb = 0
    ) %>%
    
    # Flows hub -> consumer
    add_variable(X_jk[j, k, p, t],
                 j = J, k = K, p = Pset, t = Tset,
                 type = "continuous", lb = 0
    ) %>%
    
    # Storage inventory
    add_variable(S_jpt[j, p, t],
                 j = J, p = Pset, t = Tset,
                 type = "continuous", lb = 0
    )
  
  # 2. Objective function (same structure as v8 but using z_i0t for establishment) ----
  # Z_i0t = 1 means the path from 0 to t is chosen -> establishment at that site.
  model <- model %>%
    set_objective(
      # Revenue from shipments
      sum_expr(R_p[p] * X_jk[j, k, p, t],
        j = J, k = K, p = Pset, t = Tset
      )
      -
        # Establishment cost: use z_{i,0,t} as indicator
        sum_expr(
          instance$sites$C_est[i] * z[i, 0, t],
          i = I, t = Tset
        )
      -
        # Harvest cost: proportional to total harvest over p; use sum_p Y_{i,p,t}
        sum_expr(
          instance$sites$C_harv[i] * area_ha[i] *
            sum_expr(z[i, s, t], s = 1:(t-1)),
          i = I, t = 2:Tm
        )
      -
        # Transport site->hub
        sum_expr(
          c_tr_raw * instance$d_ij[i, j] * X_ij[i, j, p, t],
          i = I, j = J, p = Pset, t = Tset
        )
      -
        # Transport hub->consumer
        sum_expr(
          c_tr_pre * instance$d_jk[j, k] * X_jk[j, k, p, t],
          j = J, k = K, p = Pset, t = Tset
        )
      -
        # Storage cost
        sum_expr(
          instance$storages$c_stor[j] * S_jpt[j, p, t],
          j = J, p = Pset, t = Tset
        ),
      sense = "max"
    )
  
  # 3. Constraints -------------------------------------------------------------
  
  # C1: path establishment (at most one path per site; z_{i,0,Tmax+1}=0 means no AFS)
  # sum_t z_{i,0,t} <= 1
  model <- model %>%
    add_constraint(
     sum_expr(z[i, 0, t], t = Tset) <= 1,
      i = I
    )
  
  # Optional: forbid arcs ending at terminal T+1 except as path terminator
  # Typically, you enforce degree constraints in C2 instead.
  
  # C2: path connectivity
  # For each site i and each intermediate time t in 1..Tmax:
  #   if t is used in the path, it must have exactly one incoming and one outgoing arc.
  # Here implemented as:
  # sum_{s < t} z_{i,s,t} = sum_{u > t} z_{i,t,u}, and both <= 1.
  model <- model %>%
    add_constraint(
      sum_expr(z[i, s, t], s = T_ext[T_ext < t]) ==
        sum_expr(z[i, t, u], u = T_ext[T_ext > t]),
      i = I, t = Tset
    )  
    # add_constraint(
    #   sum_expr(z[i, s, t], s = T_ext[T_ext < t]) <= 1,
    #   i = I, t = Tset
    # ) %>%
    # add_constraint(
    #   sum_expr(z[i, t, u], u = T_ext[T_ext > t]) <= 1,
    #   i = I, t = Tset
    # )
  
  # C3a: biomass yield calculation (same as v8, but using z path)
  # Y_{i,p,t} <= sum_{s=1}^{t-1} (eta_{p, t-s} * AREA_i * z_{i,s,t})
  model <- model %>%
    add_constraint(Y[i, p, t] <= sum_expr(yield_matrix[p, t-s] * area_ha[i] * z[i, s, t], s = 1:(t-1)) , 
                 i = I, p = Pset, t = 2:Tm )
  
  # no harvest quantities in t=1
  model <- model %>%
    add_constraint(Y[i, p, 1] == 0,
                   i = I, p = Pset)
  
  # C4 (v8a numbering): site→hub flow limited by harvest
  # sum_j X_ijpt <= Y_ipt
  model <- model %>%
    add_constraint(
      sum_expr(X_ij[i, j, p, t], j = J) <= Y[i, p, t],
      i = I, p = Pset, t = Tset
    )
  
  # Inventory balance at storage hubs (C6)
  # S_{j,p,t} = S_{j,p,t-1} + sum_i X_ijpt - sum_k X_jkpt
  # with S_{j,p,0} = 0
  for (j in J) {
    for (p in Pset) {
      # t = 1
      model <- model %>%
        add_constraint(
          S_jpt[j, p, 1] ==
            0 +
            sum_expr(X_ij[i, j, p, 1], i = I) -
            sum_expr(X_jk[j, k, p, 1], k = K)
        )
      # t >= 2
      for (t in Tset[Tset >= 2]) {
        model <- model %>%
          add_constraint(
            S_jpt[j, p, t] ==
              S_jpt[j, p, t - 1] +
              sum_expr(X_ij[i, j, p, t], i = I) -
              sum_expr(X_jk[j, k, p, t], k = K)
          )
      }
    }
  }
  
  # Storage capacity (C7)
  # sum_p S_{j,p,t} <= CAP_stor_j
  model <- model %>%
    add_constraint(
      sum_expr(S_jpt[j, p, t], p = Pset) <= instance$storages$CAP_stor[j],
      j = J, t = Tset
    )
  
  # # Processing capacity (C8)
  # # sum_{i,p} X_ijpt <= CAP_proc_j
   model <- model %>%
     add_constraint(
       sum_expr(X_ij[i, j, p, t], i = I, p = Pset) <= instance$storages$CAP_proc[j],
       j = J, t = Tset
     )
  # 
  # Demand satisfaction with cascade (C9)
  # For each (k,p,t), sum over products q that can satisfy demand p:
  #   sum_j sum_{q<=p} X_jkqt >= D_max_{k,p,t}
  # Implement cascade: 1 can satisfy 1,2,3; 2 -> 2,3; 3 -> 3.
  # C9: demand with cascade
  model <- model %>% add_constraint(
    sum_expr(X_jk[j, k, p_prime, t], j = J, p_prime = 1:p) <=
      {
        dem <- instance$demand %>%
          filter(consumer_id == k, product == p, period == t) %>%
          pull(D_max)
        if (length(dem) == 0) 0 else dem[1]
      },
    k = K, p = Pset, t = Tset
  ) 
  
  

  
  if (verbose) cat("Solving V8a model with Gurobi...\n")
  
  result <- solve_model(
    model,
    with_ROI(
      solver = "gurobi",
      control = list(
        TimeLimit = time_limit,
        MIPGap = mip_gap,
        OutputFlag = ifelse(verbose, 1, 0)
      )
    )
  )
  
  list(
    model = model,
    result = result,
    objective = objective_value(result),
    status = result$status
  )
}



################################################################################
# direct formulation
################################################################################


library(Matrix)
library(ROI)
library(ROI.plugin.gurobi)

build_agroforestry_lp_sparse <- function(instance) {
  ns  <- instance$n_sites
  nj  <- instance$n_storages
  nk  <- instance$n_consumers
  Tm  <- instance$n_periods
  P   <- instance$n_products
  I   <- 1:ns
  J   <- 1:nj
  K   <- 1:nk
  Tset <- 1:Tm
  Pset <- 1:P
  T_ext <- 0:(Tm + 1)
  
  # 1. Precompute yield matrix p x age (dense but once)
  # yields_by_age: columns product, age, yield_ha
  Ywide <- xtabs(yield_ha ~ product + age, data = instance$yields_by_age)
  yield_matrix <- matrix(0, nrow = P, ncol = Tm)
  yield_matrix[as.integer(rownames(Ywide)),
               as.integer(colnames(Ywide))] <- Ywide
  
  area_ha <- setNames(instance$sites$area_ha, seq_len(ns))
  
  # 2. Create global indices for each variable block
  # z[i,s,t] for s<t with s,t in T_ext
  z_tuples <- which(outer(T_ext, T_ext, "<"), arr.ind = TRUE)
  colnames(z_tuples) <- c("s_idx", "t_idx")
  z_df <- expand.grid(i = I, row = seq_len(nrow(z_tuples)))
  z_df$s <- T_ext[z_tuples[z_df$row, "s_idx"]]
  z_df$t <- T_ext[z_tuples[z_df$row, "t_idx"]]
  z_df$row <- NULL
  
  n_z <- nrow(z_df)
  z_df$col <- seq_len(n_z)
  
  # Y[i,p,t]
  Y_df <- expand.grid(i = I, p = Pset, t = Tset)
  n_Y <- nrow(Y_df)
  Y_df$col <- n_z + seq_len(n_Y)
  
  # X_ij[i,j,p,t]
  Xij_df <- expand.grid(i = I, j = J, p = Pset, t = Tset)
  n_Xij <- nrow(Xij_df)
  Xij_df$col <- max(Y_df$col) + seq_len(n_Xij)
  
  # X_jk[j,k,p,t]
  Xjk_df <- expand.grid(j = J, k = K, p = Pset, t = Tset)
  n_Xjk <- nrow(Xjk_df)
  Xjk_df$col <- max(Xij_df$col) + seq_len(n_Xjk)
  
  # S_jpt[j,p,t]
  S_df <- expand.grid(j = J, p = Pset, t = Tset)
  n_S <- nrow(S_df)
  S_df$col <- max(Xjk_df$col) + seq_len(n_S)
  
  n_vars <- max(S_df$col)
  
  # 3. Objective function coefficients
  
  c_vec <- numeric(n_vars)
  
  # Revenue: R_p[p] * X_jk[j,k,p,t]
  idx <- Xjk_df$col
  c_vec[idx] <- instance$R_p[Xjk_df$p]
  
  # Establishment: C_est[i] * sum_t z[i,0,t]
  # Only arcs with s==0, t in 1..Tm
  z_est <- subset(z_df, s == 0 & t %in% Tset)
  c_vec[z_est$col] <- c_vec[z_est$col] - instance$sites$C_est[z_est$i]
  
  # Harvest cost: C_harv[i] * area_ha[i] * sum_{s=1}^{t-1} z[i,s,t], t>=2
  z_harv <- subset(z_df, s >= 1 & t %in% Tset & t > s)
  harv_cost <- instance$sites$C_harv[z_harv$i] * area_ha[as.character(z_harv$i)]
  c_vec[z_harv$col] <- c_vec[z_harv$col] - harv_cost
  
  # Transport site->hub: c_tr_raw * d_ij[i,j] * X_ij
  idx <- Xij_df$col
  c_vec[idx] <- c_vec[idx] - instance$c_tr_raw *
    instance$d_ij[cbind(Xij_df$i, Xij_df$j)]
  
  # Transport hub->consumer: c_tr_pre * d_jk[j,k] * X_jk
  idx <- Xjk_df$col
  c_vec[idx] <- c_vec[idx] - instance$c_tr_pre *
    instance$d_jk[cbind(Xjk_df$j, Xjk_df$k)]
  
  # Storage cost: c_stor[j] * S_jpt
  idx <- S_df$col
  c_vec[idx] <- c_vec[idx] - instance$storages$c_stor[S_df$j]
  
  # 4. Constraints
  row_idx <- integer(0)
  col_idx <- integer(0)
  val     <- numeric(0)
  rhs     <- numeric(0)
  dir     <- character(0)

  # add row function
    
  add_row <- function(cols, vals, sense, rhs_val) {
    r <- length(rhs) + 1L
    row_idx <<- c(row_idx, rep.int(r, length(cols)))
    col_idx <<- c(col_idx, cols)
    val     <<- c(val, vals)
    dir     <<- c(dir, sense)
    rhs     <<- c(rhs, rhs_val)
  }
  
  # C1: path establishment
  
  for (i in I) {
    idx <- with(z_df, col[s == 0 & t %in% Tset & i == !!i])
    if (length(idx) > 0)
      add_row(idx, rep(1, length(idx)), "<=", 1)
  }
  
  # C2: path connectivity
  for (i in I) {
    for (t in Tset) {
      left  <- with(z_df, col[i == !!i & t == !!t & s < t])
      right <- with(z_df, col[i == !!i & s == !!t & t > !!t])
      if (length(left) + length(right) == 0) next
      cols <- c(left, right)
      vals <- c(rep(1, length(left)), rep(-1, length(right)))
      add_row(cols, vals, "==", 0)
    }
  }
  # C3a: biomass yield calculation
  for (i in I) {
    for (p in Pset) {
      for (t in Tset[Tset >= 2]) {
        s_vals <- 1:(t - 1)
        z_rows <- z_df[z_df$i == i & z_df$s %in% s_vals & z_df$t == t, ]
        if (nrow(z_rows) == 0) next
        cols <- c(Y_df$col[Y_df$i == i & Y_df$p == p & Y_df$t == t],
                  z_rows$col)
        vals <- c(1,
                  -yield_matrix[p, t - z_rows$s] * area_ha[as.character(i)])
        add_row(cols, vals, "<=", 0)
      }
    }
  }
  # C4 
  for (i in I) for (p in Pset) for (t in Tset) {
    x_rows <- Xij_df[Xij_df$i == i & Xij_df$p == p & Xij_df$t == t, ]
    if (nrow(x_rows) == 0) next
    y_col <- Y_df$col[Y_df$i == i & Y_df$p == p & Y_df$t == t]
    cols <- c(x_rows$col, y_col)
    vals <- c(rep(1, nrow(x_rows)), -1)
    add_row(cols, vals, "<=", 0)
  }
  
  # C6: inventory balance at storage hubs
  
  for (j in J) for (p in Pset) {
    # t = 1
    s_col <- S_df$col[S_df$j == j & S_df$p == p & S_df$t == 1]
    x_in  <- Xij_df$col[Xij_df$j == j & Xij_df$p == p & Xij_df$t == 1]
    x_out <- Xjk_df$col[Xjk_df$j == j & Xjk_df$p == p & Xjk_df$t == 1]
    cols <- c(s_col, x_in, x_out)
    vals <- c(1, rep(1, length(x_in)), rep(-1, length(x_out)))
    add_row(cols, vals, "==", 0)
    
    # t >= 2
    for (t in Tset[Tset >= 2]) {
      s_t  <- S_df$col[S_df$j == j & S_df$p == p & S_df$t == t]
      s_tm1<- S_df$col[S_df$j == j & S_df$p == p & S_df$t == t - 1]
      x_in <- Xij_df$col[Xij_df$j == j & Xij_df$p == p & Xij_df$t == t]
      x_out<- Xjk_df$col[Xjk_df$j == j & Xjk_df$p == p & Xjk_df$t == t]
      cols <- c(s_t, s_tm1, x_in, x_out)
      vals <- c(1, -1, rep(1, length(x_in)), rep(-1, length(x_out)))
      add_row(cols, vals, "==", 0)
    }
  }
  
  
  # C7: sum_p S_{j,p,t} <= CAP_stor_j
  for (j in J) for (t in Tset) {
    s_rows <- S_df[S_df$j == j & S_df$t == t, ]
    if (nrow(s_rows) == 0) next
    cols <- s_rows$col
    vals <- rep(1, nrow(s_rows))
    add_row(cols, vals, "<=", instance$storages$CAP_stor[j])
  }
  
  # C8: processing capacity
  for (j in J) for (t in Tset) {
    x_rows <- Xij_df[Xij_df$j == j & Xij_df$t == t, ]
    if (nrow(x_rows) == 0) next
    cols <- x_rows$col
    vals <- rep(1, nrow(x_rows))
    add_row(cols, vals, "<=", instance$storages$CAP_proc[j])
  }
  
  # C9: demand with cascade
  # For each consumer k, demanded product p, period t:
  # sum_j sum_{q <= p} X_jk[j,k,q,t] >= D_max(k,p,t)
  for (k in K) for (p in Pset) for (t in Tset) {
    # get demand value for (k,p,t)
    dem_row <- instance$demand[
      instance$demand$consumer_id == k &
        instance$demand$product     == p &
        instance$demand$period      == t,
    ]
    if (nrow(dem_row) == 0) next
    D_max_kpt <- dem_row$D_max[1]
    
    x_rows <- Xjk_df[
      Xjk_df$k == k &
        Xjk_df$t == t &
        Xjk_df$p <= p,
    ]
    if (nrow(x_rows) == 0) next
    
    cols <- x_rows$col
    vals <- rep(1, nrow(x_rows))
    add_row(cols, vals, ">=", D_max_kpt)
  }
  
  # 5. Build sparse constraint matrix
  
  A <- sparseMatrix(i = row_idx,
                    j = col_idx,
                    x = val,
                    dims = c(length(rhs), n_vars))
  
  # Split into <=, ==, >= for ROI
  leq <- dir == "<="
  eq  <- dir == "=="
  geq <- dir == ">="
  
  
  A <- sparseMatrix(i = row_idx, j = col_idx, x = val,
                    dims = c(length(rhs), n_vars))
  
  # solving via GUROBI directly
  
  model            <- list()
  model$A          <- A
  model$obj        <- c_vec
  model$modelsense <- 'max'
  model$rhs        <- rhs
  model$sense      <- dir
  

  result <- gurobi(model)
  
  # solving via ROI
  
  constr_list <- vector("list", length(rhs))
  for (r in seq_along(rhs)) {
    Lr  <- as.L_term(A[r, , drop = FALSE])  # convert row to L_term
    constr_list[[r]] <- L_constraint(L = Lr, dir = dir[r], rhs = rhs[r])
  }
  
  model <- OP(
    objective   = c_vec,
    constraints = L_constraint(
      L=A,
      dir=dir,
      rhs=rhs),
    max     = TRUE
  )
  
  
  
  L_model <- OP(
    objective = c_vec,
    constraints = rbind(
      L_constraint(L = A[leq, , drop = FALSE],
                   dir = rep("<=", sum(leq)),
                   rhs = rhs[leq]),
      L_constraint(L = A[eq, , drop = FALSE],
                   dir = rep("=", sum(eq)),
                   rhs = rhs[eq]),
      L_constraint(L = A[geq, , drop = FALSE],
                   dir = rep(">=", sum(geq)),
                   rhs = rhs[geq])
    ),
    maximum = TRUE
  )
  
  
  
  res <- ROI_solve(
    built$L_model,
    solver = "gurobi",
    control = list(TimeLimit = time_limit,
                   MIPGap = mip_gap,
                   OutputFlag = ifelse(verbose, 1, 0))
  )
  list(
    objective = res$objval,
    solution  = res$solution,
    maps      = built$var_maps
  )
  
  
}

  
    




