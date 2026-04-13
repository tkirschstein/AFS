suppressPackageStartupMessages({
  library(data.table)
})


# =============================================================================
# LP equivalence diagnostics for:
#   build_AFS_milp(milp_instance)
#   build_agroforestry_lp_rcpp(milp_instance)
#
# Assumes both functions return a list with at least:
#   row_idx, col_idx, val, rhs, sense, c_vec, ub, types,
#   n_vars, n_constrs, n_z, n_Xij, n_S, n_Xjk
#
# Optional:
#   each builder may also attach a solution vector / objective externally
# =============================================================================

tol <- 1e-9

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
cat_rule <- function(x = "-", n = 78) cat(paste(rep(x, n), collapse = ""), "\n", sep = "")

assert_same_scalar <- function(a, b, nm) {
  if (!identical(a, b)) stop(sprintf("Mismatch in %s: %s vs %s", nm, a, b), call. = FALSE)
}

safe_max_abs_diff <- function(x, y) {
  if (length(x) == 0L && length(y) == 0L) return(0)
  max(abs(x - y), na.rm = TRUE)
}

same_num <- function(x, y, tol = 1e-9) {
  same_inf <- is.infinite(x) & is.infinite(y) & sign(x) == sign(y)
  same_fin <- (!is.infinite(x) & !is.infinite(y) & abs(x - y) <= tol)
  same_inf | same_fin
}

col_family <- function(lp, col) {
  if (col >= 1 && col <= lp$instance_info$n_z) return("z")
  if (col > lp$instance_info$n_z && col <= lp$instance_info$n_z + lp$instance_info$n_Xij) return("Xij")
  if (col > lp$instance_info$n_z + lp$instance_info$n_Xij && col <= lp$instance_info$n_z + lp$instance_info$n_Xij + lp$instance_info$n_S) return("S")
  if (col > lp$instance_info$n_z + lp$instance_info$n_Xij + lp$instance_info$n_S && col <= lp$instance_info$n_vars) return("Xjk")
  "UNKNOWN"
}

annotate_cols <- function(lp, cols) {
  data.table(
    col = cols,
    family = vapply(cols, function(cc) col_family(lp, cc), character(1))
  )
}

# Reconstruct row family labels using the same row ordering as build_lp_rcpp.cpp
# C1: ns
# C2: ns * Tm
# C3: generated only when harv_arcs exist; for current formulation this should
#     usually be ns * P * nTh, but we still compute from instance robustly.
# C4: nj * P * nTh
# C5: nj * nTh
# C6: nj * nTh
# C7: depends on existing demand rows in the instance and Tharv
#
# This function is conservative: it derives C7 from demand lookup and uses the
# standard sizes for C1..C6 under current logic.
row_family_map <- function(lp, milp_instance) {
  ns <- milp_instance$n_sites
  nj <- milp_instance$n_storages
  nk <- milp_instance$n_consumers
  Tm <- milp_instance$n_periods
  P  <- milp_instance$n_products
  Amin <- milp_instance$min_age
  nTh <- length(seq.int(max(1, Amin) + 1, Tm))
  
  demand <- as.data.table(milp_instance$demand)
  demand <- demand[period %in% seq.int(max(1, Amin) + 1, Tm)]
  nC7 <- nrow(unique(demand[, .(consumer_id, product, period)]))
  
  nC1 <- ns
  nC2 <- ns * Tm
  nC3 <- ns * P * nTh
  nC4 <- nj * P * nTh
  nC5 <- nj * nTh
  nC6 <- nj * nTh
  
  fam <- c(
    rep("C1_path_est", nC1),
    rep("C2_path_flow", nC2),
    rep("C3_yield",    nC3),
    rep("C4_inventory", nC4),
    rep("C5_storage_cap", nC5),
    rep("C6_proc_cap", nC6),
    rep("C7_demand", nC7)
  )
  
  if (length(fam) != lp$n_constrs) {
    warning(sprintf(
      paste0("Derived row family map has length %d, but lp$n_constrs = %d. ",
             "Proceeding with partial labels only."),
      length(fam), lp$n_constrs
    ))
    if (length(fam) < lp$n_constrs) {
      fam <- c(fam, rep("UNMAPPED", lp$n_constrs - length(fam)))
    } else {
      fam <- fam[seq_len(lp$n_constrs)]
    }
  }
  
  data.table(row = seq_len(lp$n_constrs), row_family = fam)
}

make_triplet_dt <- function(lp) {
  dt <- data.table(
    row = as.integer(lp$row_idx),
    col = as.integer(lp$col_idx),
    val = as.numeric(lp$val)
  )
  dt[, idx := .I]
  setkey(dt, row, col)
  dt
}

objective_diff <- function(lp_r, lp_cpp, outdir) {
  d <- data.table(
    col_r = lp_r$model$objective$L$j,
    col_cpp = lp_cpp$model$objective$L$j,
    c_r = as.numeric(lp_r$model$objective$L$v),
    c_cpp = as.numeric(lp_cpp$model$objective$L$v)
  )
  d[, diff := c_r - c_cpp]
  d[, diff_ind := col_r - col_cpp]
  
  d <- d[abs(diff_ind) > tol]
  
  if (nrow(d)) {
    d[, family_r := annotate_cols(lp_r, d$col_r)$family]
    d[, family_cpp := annotate_cols(lp_cpp, d$col_cpp)$family]
    
    fwrite(d, file.path(outdir, "objective_diffs.csv"))
  }
  d
}

bound_diff <- function(lp_r, lp_cpp, outdir) {
  d <- data.table(
    col = seq_len(lp_r$n_vars),
    ub_r = as.numeric(lp_r$ub),
    ub_cpp = as.numeric(lp_cpp$ub)
  )
  d <- d[!same_num(ub_r, ub_cpp, tol)]
  if (nrow(d)) {
    d <- merge(d, annotate_cols(lp_r, d$col), by = "col", all.x = TRUE)
    fwrite(d, file.path(outdir, "ub_diffs.csv"))
  }
  d
}

type_diff <- function(lp_r, lp_cpp, outdir) {
  d <- data.table(
    col = seq_len(lp_r$n_vars),
    type_r = as.character(lp_r$types),
    type_cpp = as.character(lp_cpp$types)
  )
  d <- d[type_r != type_cpp]
  if (nrow(d)) {
    d <- merge(d, annotate_cols(lp_r, d$col), by = "col", all.x = TRUE)
    fwrite(d, file.path(outdir, "type_diffs.csv"))
  }
  d
}

rhs_sense_diff <- function(lp_r, lp_cpp, milp_instance, outdir) {
  fam <- row_family_map(lp_r, milp_instance)
  
  rhs <- data.table(
    row = seq_len(lp_r$n_constrs),
    rhs_r = as.numeric(lp_r$rhs),
    rhs_cpp = as.numeric(lp_cpp$rhs)
  )
  rhs[, diff := rhs_r - rhs_cpp]
  rhs <- rhs[abs(diff) > tol]
  if (nrow(rhs)) {
    rhs <- merge(rhs, fam, by = "row", all.x = TRUE)
    fwrite(rhs, file.path(outdir, "rhs_diffs.csv"))
  }
  
  sns <- data.table(
    row = seq_len(lp_r$n_constrs),
    sense_r = as.character(lp_r$sense),
    sense_cpp = as.character(lp_cpp$sense)
  )
  sns <- sns[sense_r != sense_cpp]
  if (nrow(sns)) {
    sns <- merge(sns, fam, by = "row", all.x = TRUE)
    fwrite(sns, file.path(outdir, "sense_diffs.csv"))
  }
  
  list(rhs = rhs, sense = sns)
}

matrix_diff <- function(lp_r, lp_cpp, milp_instance, outdir) {
  A_r <- make_triplet_dt(lp_r)
  setnames(A_r, "val", "val_r")
  
  A_c <- make_triplet_dt(lp_cpp)
  setnames(A_c, "val", "val_cpp")
  
  A <- merge(A_r[, .(row, col, val_r)], A_c[, .(row, col, val_cpp)],
             by = c("row", "col"), all = TRUE)
  
  A[is.na(val_r), val_r := 0]
  A[is.na(val_cpp), val_cpp := 0]
  A[, diff := val_r - val_cpp]
  
  fam_r <- row_family_map(lp_r, milp_instance)
  A <- merge(A, fam_r, by = "row", all.x = TRUE)
  A <- merge(A, annotate_cols(lp_r, A$col), by = "col", all.x = TRUE)
  
  D <- A[abs(diff) > tol]
  if (nrow(D)) fwrite(D, file.path(outdir, "matrix_triplet_diffs.csv"))
  
  list(all = A, diff = D)
}

summarize_by_family <- function(dt_diff, outdir, filename) {
  if (!nrow(dt_diff)) return(invisible(NULL))
  s <- dt_diff[, .(
    n = .N,
    max_abs_diff = max(abs(diff), na.rm = TRUE)
  ), by = .(row_family, family)][order(-n, row_family, family)]
  fwrite(s, file.path(outdir, filename))
  s
}

# Evaluate a candidate x-vector against an LP object
eval_lp_solution <- function(lp, x) {
  stopifnot(length(x) == lp$n_vars)
  
  obj <- sum(lp$c_vec * x)
  
  A <- make_triplet_dt(lp)
  lhs <- A[, .(lhs = sum(val * x[col])), by = row][order(row)]
  
  sense <- as.character(lp$sense)
  rhs <- as.numeric(lp$rhs)
  
  viol <- numeric(length(rhs))
  viol[sense == "<="] <- pmax(lhs$lhs[sense == "<="] - rhs[sense == "<="], 0)
  viol[sense == "=="] <- abs(lhs$lhs[sense == "=="] - rhs[sense == "=="])
  viol[sense == ">="] <- pmax(rhs[sense == ">="] - lhs$lhs[sense == ">="], 0)
  
  list(
    objective = obj,
    max_violation = max(viol),
    n_violated = sum(viol > 1e-7),
    violation = viol
  )
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
diagnose_lp_mismatch <- function(milp_instance,
                                 build_r = build_AFS_milp,
                                 build_cpp = build_agroforestry_lp_rcpp,
                                 outdir = "lp_diagnostics",
                                 solve = FALSE,
                                 solver_fun = NULL) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  cat_rule("=")
  cat("Building LPs...\n")
  lp_r   <- build_r(milp_instance)
  lp_cpp <- build_cpp(milp_instance)
  
  saveRDS(lp_r, file.path(outdir, "lp_r.rds"))
  saveRDS(lp_cpp, file.path(outdir, "lp_cpp.rds"))
  
  cat_rule("=")
  cat("Basic dimensions\n")
  print(rbind(
    R = c(
      n_vars = lp_r$instance_info$n_vars, n_constrs = lp_r$instance_info$n_constrs,
      n_z = lp_r$instance_info$n_z, n_Xij = lp_r$instance_info$n_Xij, n_S = lp_r$instance_info$n_S, n_Xjk = lp_r$instance_info$n_Xjk
    ),
    CPP = c(
      n_vars = lp_cpp$instance_info$n_vars, n_constrs = lp_cpp$instance_info$n_constrs,
      n_z = lp_cpp$instance_info$n_z, n_Xij = lp_cpp$instance_info$n_Xij, n_S = lp_cpp$instance_info$n_S, n_Xjk = lp_cpp$instance_info$n_Xjk
    )
  ))
  
  assert_same_scalar(lp_r$instance_info$n_vars, lp_cpp$instance_info$n_vars, "n_vars")
  assert_same_scalar(lp_r$instance_info$n_constrs, lp_cpp$instance_info$n_constrs, "n_constrs")
  
  cat_rule("=")
  cat("Comparing objective, bounds, types, rhs, senses, matrix...\n")
  
  d_obj <- objective_diff(lp_r, lp_cpp, outdir)
  d_ub  <- bound_diff(lp_r, lp_cpp, outdir)
  d_ty  <- type_diff(lp_r, lp_cpp, outdir)
  d_rs  <- rhs_sense_diff(lp_r, lp_cpp, milp_instance, outdir)
  d_A   <- matrix_diff(lp_r, lp_cpp, milp_instance, outdir)
  
  cat("\nObjective diffs :", nrow(d_obj), "\n")
  cat("UB diffs        :", nrow(d_ub), "\n")
  cat("Type diffs      :", nrow(d_ty), "\n")
  cat("RHS diffs       :", nrow(d_rs$rhs), "\n")
  cat("Sense diffs     :", nrow(d_rs$sense), "\n")
  cat("Matrix diffs    :", nrow(d_A$diff), "\n")
  
  if (nrow(d_obj)) {
    cat("\nTop objective diffs:\n")
    print(d_obj[1:min(20, .N)])
    print(d_obj[, .N, by = family][order(-N)])
  }
  
  if (nrow(d_rs$rhs)) {
    cat("\nTop rhs diffs:\n")
    print(d_rs$rhs[1:min(20, .N)])
  }
  
  if (nrow(d_rs$sense)) {
    cat("\nTop sense diffs:\n")
    print(d_rs$sense[1:min(20, .N)])
  }
  
  if (nrow(d_A$diff)) {
    cat("\nTop matrix diffs:\n")
    print(d_A$diff[1:min(50, .N)])
    byfam <- summarize_by_family(d_A$diff, outdir, "matrix_diff_summary.csv")
    if (!is.null(byfam)) {
      cat("\nMatrix diffs by row/col family:\n")
      print(byfam[1:min(30, .N)])
    }
  }
  
  # Optional solving / cross-evaluation
  if (isTRUE(solve)) {
    if (is.null(solver_fun)) {
      warning("solve=TRUE but solver_fun is NULL; skipping solve step.")
    } else {
      cat_rule("=")
      cat("Solving both LPs via solver_fun...\n")
      sol_r <- ROI_solve(
        lp_r,
        solver = "highs",
        control = list(
          time_limit  = 60,
          mip_rel_gap  = .05,
          log_to_console  = TRUE
        )
      )
        
      sol_c <- ROI_solve(
        lp_cpp,
        solver = "highs",
        control = list(
          time_limit  = 60,
          mip_rel_gap  = .05,
          log_to_console  = TRUE
        )
      )
        
      
      saveRDS(sol_r, file.path(outdir, "sol_r.rds"))
      saveRDS(sol_c, file.path(outdir, "sol_cpp.rds"))
      
      if (!is.null(sol_r$x) && !is.null(sol_c$x)) {
        cat("\nCross-evaluating solutions...\n")
        ev_rr <- eval_lp_solution(lp_r,   sol_r$x)
        ev_rc <- eval_lp_solution(lp_cpp, sol_r$x)
        ev_cc <- eval_lp_solution(lp_cpp, sol_c$x)
        ev_cr <- eval_lp_solution(lp_r,   sol_c$x)
        
        evtab <- rbindlist(list(
          data.table(model = "R",   solution = "R",   objective = ev_rr$objective, max_violation = ev_rr$max_violation, n_violated = ev_rr$n_violated),
          data.table(model = "CPP", solution = "R",   objective = ev_rc$objective, max_violation = ev_rc$max_violation, n_violated = ev_rc$n_violated),
          data.table(model = "CPP", solution = "CPP", objective = ev_cc$objective, max_violation = ev_cc$max_violation, n_violated = ev_cc$n_violated),
          data.table(model = "R",   solution = "CPP", objective = ev_cr$objective, max_violation = ev_cr$max_violation, n_violated = ev_cr$n_violated)
        ))
        
        fwrite(evtab, file.path(outdir, "cross_evaluation.csv"))
        print(evtab)
      }
    }
  }
  
  cat_rule("=")
  cat("Artifacts written to:", normalizePath(outdir), "\n")
  
  invisible(list(
    lp_r = lp_r,
    lp_cpp = lp_cpp,
    objective_diffs = d_obj,
    ub_diffs = d_ub,
    type_diffs = d_ty,
    rhs_diffs = d_rs$rhs,
    sense_diffs = d_rs$sense,
    matrix_diffs = d_A$diff,
    matrix_all = d_A$all
  ))
}

# =============================================================================
# Example usage
# =============================================================================
result <- diagnose_lp_mismatch(milp_instance)
#
# If you already have a solver wrapper returning list(x = ..., obj = ...), e.g.
# my_solver <- function(lp) {
#   sol <- solve_lp_with_lpsolve(lp)   # adapt to your project
#   list(x = sol$solution, obj = sol$objval)
# }
# result <- diagnose_lp_mismatch(milp_instance, solve = TRUE, solver_fun = my_solver)