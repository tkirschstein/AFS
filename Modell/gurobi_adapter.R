# ==============================================================================
# gurobi_adapter.R
#
# Native Gurobi R-API adapter for the AFS agroforestry MILP.
#
# Converts the LP object produced by build_lp_rcpp() / build_AFS_milp() into
# the list format expected by gurobi::gurobi(), optionally injects a MIP start,
# solves, and returns a structured result object.
#
# Functions exported:
#   lp_to_gurobi_model()      — pure conversion, no solve
#   solve_lp_with_gurobi()    — convert + solve with full parameter control
#   solve_afs_gurobi()        — one-call pipeline from instance → result
#
# See gurobi_warmstart.R  for make_afs_start()  (start vector construction)
# See gurobi_feascheck.R  for check_start_feasibility()  (pre-solve validation)
#
# Dependencies: Matrix, gurobi (installed from Gurobi distribution)
# Authors: tkirschstein / SmartAgroforst 2026
# ==============================================================================

# ── Required packages ─────────────────────────────────────────────────────────
.check_pkgs <- function() {
  if (!requireNamespace("Matrix",  quietly = TRUE)) stop("Package 'Matrix' required.")
  if (!requireNamespace("gurobi", quietly = TRUE))
    stop("Package 'gurobi' required. Install via: install.packages('gurobi', repos = <gurobi_R_path>)")
}


# ==============================================================================
# lp_to_gurobi_model()
#
# Converts a raw LP triplet object (output of build_lp_rcpp or build_AFS_milp)
# into a native Gurobi model list.  Does NOT call gurobi::gurobi().
#
# Args:
#   lp          — list with fields: row_idx, col_idx, val, rhs, sense, c_vec,
#                 ub, types, n_vars, n_constrs.  Optionally: lb, varnames.
#   start       — numeric vector length n_vars, or NULL.  NA entries = Gurobi
#                 decides those variables automatically (partial MIP start).
#   modelsense  — "max" (default) or "min"
#   add_names   — if TRUE, auto-generate varnames and constrnames for debugging
#
# Returns: named list suitable for gurobi::gurobi(model, params)
# ==============================================================================
lp_to_gurobi_model <- function(lp,
                                start      = NULL,
                                modelsense = c("max", "min"),
                                add_names  = FALSE) {
  .check_pkgs()
  modelsense <- match.arg(modelsense)

  # ── Field validation ────────────────────────────────────────────────────
  required <- c("row_idx", "col_idx", "val",
                "rhs", "sense", "c_vec",
                "ub", "types", "n_vars", "n_constrs")
  miss <- setdiff(required, names(lp))
  if (length(miss) > 0L)
    stop("[lp_to_gurobi_model] LP object missing fields: ", paste(miss, collapse = ", "))

  nv <- as.integer(lp$n_vars)
  nc <- as.integer(lp$n_constrs)

  if (length(lp$c_vec) != nv) stop("length(c_vec) != n_vars")
  if (length(lp$ub)    != nv) stop("length(ub) != n_vars")
  if (length(lp$types) != nv) stop("length(types) != n_vars")
  if (length(lp$rhs)   != nc) stop("length(rhs) != n_constrs")
  if (length(lp$sense) != nc) stop("length(sense) != n_constrs")
  nnz <- length(lp$val)
  if (length(lp$row_idx) != nnz || length(lp$col_idx) != nnz)
    stop("row_idx / col_idx / val must have equal length")

  # ── Sparse constraint matrix ─────────────────────────────────────────────
  A <- Matrix::sparseMatrix(
    i    = as.integer(lp$row_idx),
    j    = as.integer(lp$col_idx),
    x    = as.numeric(lp$val),
    dims = c(nc, nv)
  )

  # ── Constraint senses ────────────────────────────────────────────────────
  # Gurobi R API uses single characters: '<', '>', '='
  sense_map <- c("<="  = "<",
                 ">="  = ">",
                 "=="  = "=",
                 "<"   = "<",
                 ">"   = ">",
                 "="   = "=")
  sense <- sense_map[as.character(lp$sense)]
  bad   <- is.na(sense)
  if (any(bad))
    stop("Unsupported sense values: ",
         paste(unique(lp$sense[bad]), collapse = ", "))

  # ── Variable types ───────────────────────────────────────────────────────
  # Gurobi accepts: 'C' continuous, 'B' binary, 'I' integer
  vtype <- as.character(lp$types)
  bad_v <- setdiff(unique(vtype), c("C", "B", "I", "S", "N"))
  if (length(bad_v) > 0L)
    stop("Unsupported vtype: ", paste(bad_v, collapse = ", "))

  # ── Lower bounds ────────────────────────────────────────────────────────
  lb <- if (!is.null(lp$lb)) as.numeric(lp$lb) else rep(0.0, nv)

  # ── Assemble model list ──────────────────────────────────────────────────
  model <- list(
    A          = A,
    obj        = as.numeric(lp$c_vec),
    rhs        = as.numeric(lp$rhs),
    sense      = unname(sense),
    vtype      = vtype,
    lb         = lb,
    ub         = as.numeric(lp$ub),
    modelsense = modelsense
  )

  # ── Optional debugging names ─────────────────────────────────────────────
  if (isTRUE(add_names)) {
    model$varnames   <- sprintf("x_%d", seq_len(nv))
    model$constrnames <- sprintf("c_%d", seq_len(nc))
  }

  # ── MIP start ────────────────────────────────────────────────────────────
  # NA entries are passed as-is; Gurobi treats them as GRB_UNDEFINED and
  # completes the start automatically.  Full starts (no NA) are most effective.
  if (!is.null(start)) {
    if (length(start) != nv)
      stop("length(start) must equal n_vars (", nv, ")")
    model$start <- as.numeric(start)
  }

  model
}


# ==============================================================================
# solve_lp_with_gurobi()
#
# Converts an LP object and solves it with the native Gurobi R API.
#
# Args:
#   lp         — LP triplet object (see lp_to_gurobi_model)
#   start      — optional warm-start vector (see make_afs_start in
#                gurobi_warmstart.R)
#   params     — named list of Gurobi parameters (see Gurobi reference manual)
#                Useful keys: TimeLimit, MIPFocus, MIPGap, Threads,
#                             Heuristics, Cuts, Presolve, NumericFocus
#   modelsense — "max" or "min"
#   add_names  — forward to lp_to_gurobi_model (useful for IIS debugging)
#   verbose    — if FALSE, sets OutputFlag = 0
#
# Returns: list of class "afs_gurobi_result" with:
#   $model     — the Gurobi model list (for re-use / inspection)
#   $params    — the params list actually used
#   $result    — raw gurobi() output
#   $x         — solution vector (NA if infeasible / no incumbent)
#   $objval    — objective value
#   $status    — Gurobi status string
#   $runtime   — wall-clock seconds
# ==============================================================================
solve_lp_with_gurobi <- function(lp,
                                  start      = NULL,
                                  params     = list(),
                                  modelsense = "max",
                                  add_names  = FALSE,
                                  verbose    = F) {
  .check_pkgs()

  model <- lp_to_gurobi_model(
    lp, start = start, modelsense = modelsense, add_names = add_names
  )

  # ── Merge default params with user-supplied ──────────────────────────────
  defaults <- list(OutputFlag = if (verbose) 1L else 0L)
  params   <- utils::modifyList(defaults, params)

  cat(sprintf(
    "[solve_lp_with_gurobi] %d vars | %d constrs | %d nnz | sense=%s | start=%s\n",
    lp$n_vars, lp$n_constrs, length(lp$val),
    modelsense,
    if (is.null(start)) "none" else sprintf("partial(%d NA)", sum(is.na(start)))
  ))

  t0  <- proc.time()["elapsed"]
  res <- gurobi::gurobi(model, params)
  rt  <- as.numeric(proc.time()["elapsed"] - t0)

  cat(sprintf("[solve_lp_with_gurobi] Status: %s | ObjVal: %s | Runtime: %.1f s\n",
              res$status,
              if (!is.null(res$objval)) sprintf("%.4f", res$objval) else "N/A",
              rt))

  structure(
    list(
      model   = model,
      params  = params,
      result  = res,
      x       = res$x,
      objval  = res$objval,
      status  = res$status,
      runtime = rt
    ),
    class = "afs_gurobi_result"
  )
}


# ==============================================================================
# print.afs_gurobi_result()  —  S3 print method
# ==============================================================================
print.afs_gurobi_result <- function(x, ...) {
  cat("── AFS Gurobi Result ──────────────────────────────\n")
  cat(sprintf("  Status  : %s\n",  x$status))
  cat(sprintf("  ObjVal  : %s\n",  if (!is.null(x$objval)) sprintf("%.6f", x$objval) else "N/A"))
  cat(sprintf("  Runtime : %.1f s\n", x$runtime))
  if (!is.null(x$result$mipgap))
    cat(sprintf("  MIPGap  : %.4f%%\n", x$result$mipgap * 100))
  if (!is.null(x$result$nodecount))
    cat(sprintf("  Nodes   : %.0f\n",   x$result$nodecount))
  cat("────────────────────────────────────────────────────\n")
  invisible(x)
}


# ==============================================================================
# solve_afs_gurobi()
#
# Convenience one-liner: build LP from instance, optionally inject warm start,
# solve with native Gurobi.
#
# Args:
#   instance   — output of build_optimization_instance()
#   builder    — "rcpp" (default, fast) or "r" (pure-R reference)
#   start      — optional warm-start vector (use make_afs_start() to build)
#   params     — Gurobi parameter list
#   modelsense — "max" or "min"
#   verbose    — passed to solve_lp_with_gurobi()
#
# Returns: afs_gurobi_result
# ==============================================================================
solve_afs_gurobi <- function(instance,
                              builder    = c("rcpp", "r"),
                              start      = NULL,
                              params     = list(),
                              modelsense = "max",
                              verbose    = TRUE) {
  builder <- match.arg(builder)

  cat(sprintf("[solve_afs_gurobi] Building LP via '%s' builder ...\n", builder))

  lp_obj <- switch(
    builder,
    rcpp = {
      if (!exists("build_agroforestry_lp_rcpp", mode = "function"))
        stop("build_agroforestry_lp_rcpp() not found. Run: ",
             "Rcpp::sourceCpp('build_lp_rcpp.cpp') and ",
             "source('build_lp_rcpp_wrapper.R')")
      out <- build_agroforestry_lp_rcpp(instance)
      # Unwrap: the wrapper returns list(model=ROI_OP, instance_info=...)
      # We need the raw cpp_out; re-call the C++ builder directly.
      # Rebuild via the C++ function to get the raw triplet LP.
      cpp_out <- build_lp_rcpp(instance)
      cpp_out
    },
    r = {
      if (!exists("build_AFS_milp", mode = "function"))
        stop("build_AFS_milp() not found. source('!build_AFS_milp.R')")
      build_AFS_milp(instance)
    }
  )

  solve_lp_with_gurobi(
    lp         = lp_obj,
    start      = start,
    params     = params,
    modelsense = modelsense,
    verbose    = verbose
  )
}
