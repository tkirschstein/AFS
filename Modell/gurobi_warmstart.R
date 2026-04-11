# ==============================================================================
# gurobi_warmstart.R
#
# Helpers to construct and validate MIP start vectors for the AFS MILP.
#
# Variable layout in the LP vector (1-based):
#
#   [ z_1 ... z_{n_z} | Xij_1 ... Xij_{n_Xij} | S_1 ... S_{n_S} | Xjk_1 ... Xjk_{n_Xjk} ]
#   |<---- block 1 --->|<------- block 2 ------->|<-- block 3 --->|<------- block 4 ------->|
#
# Functions exported:
#   make_afs_start()         — build a start vector from named block vectors
#   start_from_solution()    — extract a start from a previous solve result
#   trivial_zero_start()     — all-zero start (usually feasible for MILPs
#                              with full Big-M relaxation)
#
# Authors: tkirschstein / SmartAgroforst 2026
# ==============================================================================


# ==============================================================================
# make_afs_start()
#
# Constructs a full-length start vector from separate block vectors.
# Any block left as NULL is filled with `fill` (default NA = Gurobi decides).
#
# Args:
#   lp   — raw LP triplet list; must have n_vars, n_z, n_Xij, n_S, n_Xjk
#   z    — numeric vector length n_z,   or NULL
#   Xij  — numeric vector length n_Xij, or NULL
#   S    — numeric vector length n_S,   or NULL
#   Xjk  — numeric vector length n_Xjk, or NULL
#   fill — scalar used for unspecified blocks (NA or 0)
#
# Returns: numeric vector length n_vars
# ==============================================================================
make_afs_start <- function(lp,
                            z    = NULL,
                            Xij  = NULL,
                            S    = NULL,
                            Xjk  = NULL,
                            fill = NA_real_) {

  # ── Validate layout fields ───────────────────────────────────────────────
  req <- c("n_vars", "n_z", "n_Xij", "n_S", "n_Xjk")
  miss <- setdiff(req, names(lp))
  if (length(miss) > 0L)
    stop("[make_afs_start] LP object missing block-size fields: ",
         paste(miss, collapse = ", "), "\n",
         "  These are set by build_lp_rcpp(); check your builder output.")

  nv    <- as.integer(lp$n_vars)
  n_z   <- as.integer(lp$n_z)
  n_Xij <- as.integer(lp$n_Xij)
  n_S   <- as.integer(lp$n_S)
  n_Xjk <- as.integer(lp$n_Xjk)

  declared_total <- n_z + n_Xij + n_S + n_Xjk
  if (declared_total != nv)
    warning(sprintf(
      "[make_afs_start] n_z+n_Xij+n_S+n_Xjk = %d != n_vars = %d.\n",
      declared_total, nv
    ))

  # ── Block offsets (1-based) ──────────────────────────────────────────────
  off_z   <- 0L
  off_Xij <- off_z   + n_z
  off_S   <- off_Xij + n_Xij
  off_Xjk <- off_S   + n_S

  x0 <- rep(fill, nv)

  .set_block <- function(x0, offset, n, vec, name) {
    if (is.null(vec)) return(x0)
    if (length(vec) != n)
      stop(sprintf("[make_afs_start] length(%s) = %d but n_%s = %d",
                   name, length(vec), name, n))
    idx <- seq.int(offset + 1L, offset + n)
    x0[idx] <- as.numeric(vec)
    x0
  }

  x0 <- .set_block(x0, off_z,   n_z,   z,   "z")
  x0 <- .set_block(x0, off_Xij, n_Xij, Xij, "Xij")
  x0 <- .set_block(x0, off_S,   n_S,   S,   "S")
  x0 <- .set_block(x0, off_Xjk, n_Xjk, Xjk, "Xjk")

  cat(sprintf(
    "[make_afs_start] start vector: %d vars | %d specified | %d NA (Gurobi auto)\n",
    nv,
    sum(!is.na(x0)),
    sum(is.na(x0))
  ))

  x0
}


# ==============================================================================
# trivial_zero_start()
#
# Returns an all-zero start vector: all paths inactive, no flows.
# This is always primal-feasible when all RHS are >= 0 and bounds include 0,
# which holds for the AFS MILP (no-plant = valid, zero-profit decision).
#
# Args:
#   lp — LP triplet object
#
# Returns: numeric vector of zeros, length n_vars
# ==============================================================================
trivial_zero_start <- function(lp) {
  make_afs_start(
    lp,
    z    = rep(0.0, lp$n_z),
    Xij  = rep(0.0, lp$n_Xij),
    S    = rep(0.0, lp$n_S),
    Xjk  = rep(0.0, lp$n_Xjk),
    fill = 0.0
  )
}


# ==============================================================================
# start_from_solution()
#
# Extracts a warm start from a previous afs_gurobi_result.
# Useful for iterative re-solves with modified parameters or constraints.
#
# Args:
#   res — object of class afs_gurobi_result (from solve_lp_with_gurobi)
#
# Returns: numeric vector length n_vars, or NULL if no primal solution exists
# ==============================================================================
start_from_solution <- function(res) {
  if (!inherits(res, "afs_gurobi_result"))
    stop("[start_from_solution] expected object of class 'afs_gurobi_result'")
  if (is.null(res$x)) {
    warning("[start_from_solution] result has no primal solution (status: ",
            res$status, "). Returning NULL.")
    return(NULL)
  }
  cat(sprintf("[start_from_solution] Reusing solution from status=%s, objval=%.4f\n",
              res$status, res$objval))
  res$x
}
