# build_lp_rcpp_wrapper.R
# Compiles and uses the Rcpp LP builder
# Rcpp::sourceCpp("build_lp_rcpp.cpp")  # run once per session

build_agroforestry_lp_rcpp <- function(instance) {
  cat("Calling Rcpp LP builder...\n")
  
  # Pre-compute yield_matrix once in R (clean interface boundary)
  instance$yield_matrix <- local({
    m <- matrix(0, nrow = instance$n_products, ncol = instance$n_periods)
    yba <- instance$yields_by_age
    for (p in seq_len(instance$n_products))
      for (age in seq_len(instance$n_periods)) {
        r <- yba[yba$product == p & yba$age == age, "yield_ha"]
        if (length(r) > 0) m[p, age] <- r[1]
      }
    m
  })
  
  # Call C++ builder
  cpp_out <- build_lp_rcpp(instance)
  
  # Assemble sparse constraint matrix
  A <- slam::simple_triplet_matrix(
    i    = cpp_out$row_idx,
    j    = cpp_out$col_idx,
    v    = cpp_out$val,
    nrow = cpp_out$n_constrs,
    ncol = cpp_out$n_vars
  )
  
  # Build ROI OP object
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
    model        = model,
    instance_info = cpp_out[c("n_vars","n_constrs","n_z","n_Y","n_Xij","n_S","n_Xjk")]
  )
}