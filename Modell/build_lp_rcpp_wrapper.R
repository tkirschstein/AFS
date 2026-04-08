# build_lp_rcpp_wrapper.R
# Wrapper: compiles build_lp_rcpp.cpp and assembles the ROI OP object.
# Run Rcpp::sourceCpp("build_lp_rcpp.cpp") once per session before calling this.
#
# Matches build_agroforestry_lp_sparse_v10_optimized() in !build_AFS_milp.R:
#   - No Y variable
#   - C3 combined yield+shipping bound (<=)
#   - C6 inventory balance over Tharv only
#   - C7/C8 over Tharv only
#   - C9 demand cascade: delivered product filters pp_pairs[pi].second == p

build_agroforestry_lp_rcpp <- function(instance) {
  cat("Calling Rcpp LP builder (v10-no-Y)...\n")

  # ── Pre-compute yield_matrix [P x Tm] for the C++ side ───────────────────
  # yield_matrix[p, age-1] = eta_{p, age}  (0-based column = age - 1)
  # instance$yields_by_age must have columns: product (1..P), age (1..Tm), yield_ha
  P  <- instance$n_products
  Tm <- instance$n_periods
  ym <- matrix(0.0, nrow = P, ncol = Tm)
  yba <- instance$yields_by_age
  for (r in seq_len(nrow(yba))) {
    p   <- yba$product[r]
    age <- yba$age[r]
    if (p >= 1 && p <= P && age >= 1 && age <= Tm)
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
