# ==============================================================================
# gurobi_feascheck.R
#
# Pre-solve feasibility and quality checks for AFS MILP warm-start vectors.
#
# Validates a candidate solution x0 against the LP object before submitting
# it to Gurobi.  Catches bugs in the start-vector construction (wrong offsets,
# violated bounds, binding equality constraints, C3 violations) without
# incurring any Gurobi API call.
#
# Functions exported:
#   check_start_feasibility()   — full constraint + bound check, returns report
#   check_c3_linkage()          — targeted check of the harvest-linkage
#                                 constraint block (the currently broken C3)
#   summarise_lp()              — compact structural summary of an LP object
#
# Authors: tkirschstein / SmartAgroforst 2026
# ==============================================================================


# ==============================================================================
# summarise_lp()
#
# Prints a compact structural summary of an LP triplet object.
# Useful for quick sanity-checks after building.
#
# Args:
#   lp      — LP triplet list
#   detail  — if TRUE, prints sense/type breakdown
# ==============================================================================
summarise_lp <- function(lp, detail = TRUE) {
  cat("── LP Object Summary ────────────────────────────────\n")
  cat(sprintf("  Variables    : %d\n", lp$n_vars))
  cat(sprintf("  Constraints  : %d\n", lp$n_constrs))
  cat(sprintf("  Non-zeros    : %d  (density %.4f%%)\n",
              length(lp$val),
              100 * length(lp$val) / (lp$n_vars * lp$n_constrs)))
  cat(sprintf("  Obj range    : [%.4g, %.4g]\n",
              min(lp$c_vec), max(lp$c_vec)))
  cat(sprintf("  UB range     : [%.4g, %.4g]\n",
              min(lp$ub), max(lp$ub)))
  cat(sprintf("  RHS range    : [%.4g, %.4g]\n",
              min(lp$rhs), max(lp$rhs)))

  if (isTRUE(detail)) {
    s_tbl <- table(lp$sense)
    cat("  Sense counts :")
    for (nm in names(s_tbl))
      cat(sprintf(" %s=%d", nm, s_tbl[[nm]]))
    cat("\n")

    v_tbl <- table(lp$types)
    cat("  Type counts  :")
    for (nm in names(v_tbl))
      cat(sprintf(" %s=%d", nm, v_tbl[[nm]]))
    cat("\n")

    if (!is.null(lp$n_z)) {
      cat(sprintf("  Block n_z    : %d\n",   lp$n_z))
      cat(sprintf("  Block n_Xij  : %d\n",   lp$n_Xij))
      cat(sprintf("  Block n_S    : %d\n",   lp$n_S))
      cat(sprintf("  Block n_Xjk  : %d\n",   lp$n_Xjk))
    }
  }
  cat("──────────────────────────────────────────────────────\n")
  invisible(lp)
}


# ==============================================================================
# check_start_feasibility()
#
# Checks whether a candidate start vector x0 satisfies:
#   1. Bound feasibility:    lb <= x0 <= ub
#   2. Integrality:          binary / integer vars are (near-)integer
#   3. Constraint residuals: A x0 vs rhs for each sense
#
# Returns a list (invisible) with:
#   $feasible      — TRUE / FALSE
#   $bound_viols   — data.frame of bound violations
#   $int_viols     — data.frame of integrality violations
#   $constr_viols  — data.frame of constraint violations
#   $summary       — named numeric: counts of each violation type
#
# Args:
#   lp        — LP triplet object (must have row_idx, col_idx, val, rhs,
#               sense, ub, lb (optional), types, n_vars, n_constrs)
#   x0        — numeric vector length n_vars
#   tol       — feasibility tolerance (default 1e-6)
#   max_print — max rows printed per violation table
# ==============================================================================
check_start_feasibility <- function(lp, x0,
                                     tol       = 1e-6,
                                     max_print = 20L) {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' required.")

  nv <- as.integer(lp$n_vars)
  nc <- as.integer(lp$n_constrs)

  if (length(x0) != nv)
    stop(sprintf("length(x0) = %d but n_vars = %d", length(x0), nv))

  lb <- if (!is.null(lp$lb)) as.numeric(lp$lb) else rep(0.0, nv)
  ub <- as.numeric(lp$ub)

  # ── 1. Bound feasibility ─────────────────────────────────────────────────
  lb_viol  <- which(x0 < lb - tol)
  ub_viol  <- which(x0 > ub + tol)
  bv_idx   <- union(lb_viol, ub_viol)

  bound_df <- if (length(bv_idx) > 0L) {
    data.frame(
      var_idx   = bv_idx,
      x0_val    = x0[bv_idx],
      lb        = lb[bv_idx],
      ub        = ub[bv_idx],
      viol_type = ifelse(bv_idx %in% lb_viol, "lb", "ub")
    )
  } else {
    data.frame(var_idx=integer(0), x0_val=numeric(0),
               lb=numeric(0), ub=numeric(0), viol_type=character(0))
  }

  # ── 2. Integrality ───────────────────────────────────────────────────────
  int_mask  <- lp$types %in% c("B", "I")
  frac_vals <- abs(x0[int_mask] - round(x0[int_mask]))
  int_viol_idx <- which(int_mask)[which(frac_vals > tol)]

  int_df <- if (length(int_viol_idx) > 0L) {
    data.frame(
      var_idx  = int_viol_idx,
      x0_val   = x0[int_viol_idx],
      frac_gap = frac_vals[match(int_viol_idx, which(int_mask))],
      vtype    = lp$types[int_viol_idx]
    )
  } else {
    data.frame(var_idx=integer(0), x0_val=numeric(0),
               frac_gap=numeric(0), vtype=character(0))
  }

  # ── 3. Constraint residuals ──────────────────────────────────────────────
  A <- Matrix::sparseMatrix(
    i    = as.integer(lp$row_idx),
    j    = as.integer(lp$col_idx),
    x    = as.numeric(lp$val),
    dims = c(nc, nv)
  )
  Ax  <- as.numeric(A %*% x0)
  rhs <- as.numeric(lp$rhs)

  sense_map <- as.character(lp$sense)
  sense_map[sense_map == "<="] <- "<"
  sense_map[sense_map == ">="] <- ">"
  sense_map[sense_map == "=="] <- "="

  resid <- numeric(nc)
  viol  <- logical(nc)
  for (r in seq_len(nc)) {
    resid[r] <- Ax[r] - rhs[r]
    viol[r]  <- switch(
      sense_map[r],
      "<" = resid[r] >  tol,
      ">" = resid[r] < -tol,
      "=" = abs(resid[r]) > tol,
      FALSE
    )
  }

  c_idx <- which(viol)
  constr_df <- if (length(c_idx) > 0L) {
    data.frame(
      constr_idx = c_idx,
      Ax_val     = Ax[c_idx],
      rhs_val    = rhs[c_idx],
      residual   = resid[c_idx],
      sense      = sense_map[c_idx]
    )
  } else {
    data.frame(constr_idx=integer(0), Ax_val=numeric(0),
               rhs_val=numeric(0), residual=numeric(0), sense=character(0))
  }

  # ── Report ───────────────────────────────────────────────────────────────
  n_bv  <- nrow(bound_df)
  n_iv  <- nrow(int_df)
  n_cv  <- nrow(constr_df)
  feasible <- (n_bv + n_iv + n_cv) == 0L

  cat("── Feasibility Check Report ─────────────────────────\n")
  cat(sprintf("  Variables    : %d\n", nv))
  cat(sprintf("  Constraints  : %d\n", nc))
  cat(sprintf("  Tolerance    : %.2e\n", tol))
  cat(sprintf("  Bound viols  : %d\n", n_bv))
  cat(sprintf("  Int viols    : %d\n", n_iv))
  cat(sprintf("  Constr viols : %d\n", n_cv))
  cat(sprintf("  FEASIBLE     : %s\n", if (feasible) "YES" else "NO"))
  cat("─────────────────────────────────────────────────────\n")

  if (n_bv > 0L) {
    cat("\n  Bound violations (first", min(max_print, n_bv), "):\n")
    print(utils::head(bound_df, max_print), row.names = FALSE)
  }
  if (n_iv > 0L) {
    cat("\n  Integrality violations (first", min(max_print, n_iv), "):\n")
    print(utils::head(int_df, max_print), row.names = FALSE)
  }
  if (n_cv > 0L) {
    cat("\n  Constraint violations (first", min(max_print, n_cv), "):\n")
    # Sort by absolute residual for diagnostic priority
    constr_df <- constr_df[order(-abs(constr_df$residual)), ]
    print(utils::head(constr_df, max_print), row.names = FALSE)
  }

  invisible(list(
    feasible     = feasible,
    bound_viols  = bound_df,
    int_viols    = int_df,
    constr_viols = constr_df,
    summary      = c(
      bound_violations      = n_bv,
      integrality_violations = n_iv,
      constraint_violations  = n_cv
    )
  ))
}


# ==============================================================================
# check_c3_linkage()
#
# Targeted diagnostic for the harvest-linkage constraint (C3):
#   sum_j X_{ijpt} <= sum_{(s,t) in S^harv} eta_{p,t-s} * AREA_i * z_{ist}
#
# Identifies (i, p, t) triples where X > 0 but the z-weighted RHS = 0,
# which is exactly the bug observed in the current implementation.
#
# Args:
#   lp         — LP triplet object
#   x0         — solution or start vector, length n_vars
#   c3_row_idx — integer vector of C3 constraint row indices in the LP
#                (set this from your builder's constraint-label output;
#                 if NULL the function attempts auto-detection via LHS structure)
#   tol        — violation tolerance
#
# Returns: data.frame with violating C3 rows (invisible)
# ==============================================================================
check_c3_linkage <- function(lp, x0,
                              c3_row_idx = NULL,
                              tol        = 1e-6) {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' required.")

  nv <- as.integer(lp$n_vars)
  nc <- as.integer(lp$n_constrs)

  if (length(x0) != nv)
    stop(sprintf("length(x0) = %d != n_vars = %d", length(x0), nv))

  # ── Auto-detect C3 rows if not supplied ─────────────────────────────────
  # C3 rows have structure:  +1*X_{...} - coeff*z_{...} <= 0
  # Heuristic: look for '<=' rows whose RHS = 0 and which have at least
  # one positive-coef variable in the Xij block and one negative-coef
  # variable in the z block.
  if (is.null(c3_row_idx)) {
    cat("[check_c3_linkage] c3_row_idx not supplied; attempting auto-detection ...\n")

    sense_norm <- as.character(lp$sense)
    sense_norm[sense_norm == "<="] <- "<"
    leq_rows_zero <- which(sense_norm == "<" & abs(lp$rhs) < tol)

    n_z   <- as.integer(lp$n_z)
    n_Xij <- as.integer(lp$n_Xij)
    Xij_range <- seq.int(n_z + 1L, n_z + n_Xij)
    z_range   <- seq.int(1L, n_z)

    A <- Matrix::sparseMatrix(
      i    = as.integer(lp$row_idx),
      j    = as.integer(lp$col_idx),
      x    = as.numeric(lp$val),
      dims = c(nc, nv)
    )

    is_c3 <- vapply(leq_rows_zero, function(r) {
      cols <- A@j[A@p[r]:(A@p[r + 1L] - 1L)] + 1L   # 1-based
      vals <- A@x[A@p[r]:(A@p[r + 1L] - 1L)]
      has_pos_Xij <- any(cols %in% Xij_range & vals > tol)
      has_neg_z   <- any(cols %in% z_range   & vals < -tol)
      has_pos_Xij && has_neg_z
    }, logical(1L))

    c3_row_idx <- leq_rows_zero[is_c3]
    cat(sprintf("  Auto-detected %d C3 rows.\n", length(c3_row_idx)))
  }

  if (length(c3_row_idx) == 0L) {
    cat("[check_c3_linkage] No C3 rows found. Check builder output.\n")
    return(invisible(data.frame()))
  }

  # ── Evaluate C3 rows ─────────────────────────────────────────────────────
  A <- Matrix::sparseMatrix(
    i    = as.integer(lp$row_idx),
    j    = as.integer(lp$col_idx),
    x    = as.numeric(lp$val),
    dims = c(nc, nv)
  )
  Ax  <- as.numeric(A[c3_row_idx, , drop = FALSE] %*% x0)
  rhs <- as.numeric(lp$rhs[c3_row_idx])

  viol_mask <- Ax > rhs + tol

  # Also flag rows where LHS > 0 but z-contribution = 0 (the known bug)
  n_z      <- as.integer(lp$n_z)
  z_range  <- seq.int(1L, n_z)
  z_part   <- as.numeric(A[c3_row_idx, z_range, drop = FALSE] %*% x0[z_range])
  Xij_part <- Ax - z_part   # the positive-coef (X) part of the LHS

  z_zero_X_pos <- Xij_part > tol & abs(z_part) < tol

  problem_mask <- viol_mask | z_zero_X_pos

  result_df <- data.frame(
    row_idx       = c3_row_idx[problem_mask],
    Ax            = Ax[problem_mask],
    rhs           = rhs[problem_mask],
    residual      = Ax[problem_mask] - rhs[problem_mask],
    Xij_part      = Xij_part[problem_mask],
    z_part        = z_part[problem_mask],
    violated_C3   = viol_mask[problem_mask],
    z_zero_X_pos  = z_zero_X_pos[problem_mask]
  )

  cat(sprintf("[check_c3_linkage] C3 rows checked: %d | Problems: %d\n",
              length(c3_row_idx), nrow(result_df)))
  if (nrow(result_df) > 0L) {
    cat("  Rows where X > 0 but z = 0 (C3 not enforced):\n")
    print(utils::head(result_df[order(-result_df$Xij_part), ], 20L),
          row.names = FALSE)
  } else {
    cat("  No C3 violations detected.\n")
  }

  invisible(result_df)
}
