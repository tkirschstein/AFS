# build_lp_rcpp_wrapper.R
# Wrapper: compiles build_lp_rcpp.cpp and assembles the ROI OP object.
# Run Rcpp::sourceCpp("build_lp_rcpp.cpp") once per session before calling this.
#
# Instance must be produced by build_optimization_instance() in
# !helper_instance_builder_v8a.R (or generate_instance_v8_final()).
# Required fields:
#   $n_sites, $n_storages, $n_consumers, $n_periods, $n_products, $max_age, $min_age
#   $sites        — data.frame with: site_id, area_ha, C_est, C_harv, C_main, C_opp
#   $storages     — data.frame with: storage_id, CAP_stor, CAP_proc, c_stor
#   $consumer_prices — data.frame with: consumer_id, product, price
#   $demand       — data.frame with: consumer_id, product, period, D_max
#   $d_ij         — numeric matrix [n_sites x n_storages]
#   $d_jk         — numeric matrix [n_storages x n_consumers]
#   $yields_by_age — data.frame with: product (1..P), age (1..Tm), yield_ha
#   $c_tr_raw, $c_tr_pre

build_agroforestry_lp_rcpp_ROI <- function(instance) {
  cat("Calling Rcpp LP builder (v10-no-Y, matching !build_AFS_milp.R)...\n")

  # ── Validate required site columns ───────────────────────────────────────
  required_site_cols <- c("area_ha", "C_est", "C_harv", "C_main", "C_opp")
  missing_cols <- setdiff(required_site_cols, colnames(instance$sites))
  if (length(missing_cols) > 0)
    stop(paste("sites df missing columns required by Rcpp builder:",
               paste(missing_cols, collapse = ", ")))

  # ── Pre-compute yield_matrix [P x Tm] for the C++ side ───────────────────
  # yield_matrix[p, age] = eta_{p, age}  (1-based age → R matrix column = age)
  # C++ reads yield_matrix(p, age-1) with 0-based row/col.
  P  <- instance$n_products
  Tm <- instance$n_periods
  ym <- matrix(0.0, nrow = P, ncol = Tm)
  yba <- instance$yields_by_age
  for (r in seq_len(nrow(yba))) {
    p   <- yba$product[r]
    age <- yba$age[r]
    if (p >= 1L && p <= P && age >= 1L && age <= Tm)
      ym[p, age] <- yba$yield_ha[r]
  }
  instance$yield_matrix <- ym

  # ── Call C++ builder ──────────────────────────────────────────────────────
  cpp_out <- build_lp_rcpp(instance)

  # ── Assemble sparse constraint matrix ────────────────────────────────────
  A <- slam::simple_triplet_matrix(
    i    = cpp_out$row_idx,
    j    = cpp_out$col_idx,
    v    = cpp_out$val,
    nrow = cpp_out$n_constrs,
    ncol = cpp_out$n_vars
  )

  # ── Build ROI OP object ───────────────────────────────────────────────────
  model <- ROI::OP(
    objective   = ROI::L_objective(cpp_out$c_vec),
    constraints = ROI::L_constraint(
      L   = A,
      dir = cpp_out$sense,
      rhs = cpp_out$rhs
    ),
    maximum = TRUE,
    bounds  = ROI::V_bound(
      li   = seq_len(cpp_out$n_vars),
      ui   = seq_len(cpp_out$n_vars),
      lb   = rep(0, cpp_out$n_vars),
      ub   = cpp_out$ub,
      nobj = cpp_out$n_vars
    ),
    types = cpp_out$types
  )

  cat(sprintf("ROI model built: %d vars, %d constraints, %d nnz\n",
              cpp_out$n_vars, cpp_out$n_constrs, length(cpp_out$val)))

  list(
    model         = model,
    instance_info = cpp_out[c("n_vars", "n_constrs",
                               "n_z", "n_Xij", "n_S", "n_Xjk")]
  )
}



build_agroforestry_lp_rcpp <- function(instance) {
  
  # ── Validate required site columns ───────────────────────────────────────
  required_site_cols <- c("area_ha", "C_est", "C_harv", "C_main", "C_opp")
  missing_cols <- setdiff(required_site_cols, colnames(instance$sites))
  if (length(missing_cols) > 0)
    stop(paste("sites df missing columns required by Rcpp builder:",
               paste(missing_cols, collapse = ", ")))
  
  # ── Pre-compute yield_matrix [P x Tm] for the C++ side ───────────────────
  # yield_matrix[p, age] = eta_{p, age}  (1-based age → R matrix column = age)
  # C++ reads yield_matrix(p, age-1) with 0-based row/col.
  P  <- instance$n_products
  Tm <- instance$n_periods
  ym <- matrix(0.0, nrow = P, ncol = Tm)
  yba <- instance$yields_by_age
  for (r in seq_len(nrow(yba))) {
    p   <- yba$product[r]
    age <- yba$age[r]
    if (p >= 1L && p <= P && age >= 1L && age <= Tm)
      ym[p, age] <- yba$yield_ha[r]
  }
  instance$yield_matrix <- ym
  
  # ── Call C++ builder ──────────────────────────────────────────────────────
  cpp_out <- build_lp_rcpp(instance)
  
  
  return(cpp_out)
}
