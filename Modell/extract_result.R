# ==============================================================================
# extract_result.R
#
# Decodes an afs_gurobi_result object into human-readable, analysis-ready
# data frames.  Mirrors the EXACT variable layout of build_lp_rcpp.cpp v10:
#
#   Block       Offset         Size
#   ─────────────────────────────────────────────────
#   z           off_z          ns * n_arcs   (binary)
#   Xij         off_Xij        ns*nj*P*nTh   (continuous)
#   S           off_S          nj*P*nTh      (continuous)
#   Xjk         off_Xjk        nj*nk*n_ppairs*nTh (continuous)
#
# Arc enumeration and Tharv construction MUST match build_lp_rcpp.cpp exactly
# (see inline comments).
#
# Functions exported:
#   build_index_tables(instance)        — re-derive all index structures
#   extract_z(res, idx)                 — planting / harvest decisions
#   extract_Xij(res, idx)               — site → hub biomass flows
#   extract_S(res, idx)                 — hub inventory levels
#   extract_Xjk(res, idx)               — hub → consumer deliveries
#   extract_result(res, instance)       — master extractor → named list
#   summarise_result(extracted)         — one-page KPI console print
#   validate_extracted(extracted, tol)  — mass-balance + C3 checks on output
#
# Dependencies: dplyr (tidyverse), tibble
# Authors:      tkirschstein / SmartAgroforst 2026
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})


# ==============================================================================
# build_index_tables()
#
# Reconstructs every index structure that build_lp_rcpp.cpp builds internally,
# so result extraction does not require re-running the C++ builder.
#
# Args:
#   instance — the same list passed to build_lp_rcpp()
#
# Returns: list with components
#   $ns, $nj, $nk, $Tm, $P, $Amin, $Amax
#   $arcs        — data.frame(arc_id, s, t)          — same order as cpp
#   $n_arcs
#   $Tharv       — integer vector of harvest periods
#   $nTh
#   $pp_pairs    — data.frame(pi, src_prod, del_prod) — 0-based products
#   $n_ppairs
#   $off_z, $off_Xij, $off_S, $off_Xjk
#   $n_z, $n_Xij, $n_S, $n_Xjk, $n_vars
#   $sites, $storages
# ==============================================================================
build_index_tables <- function(instance) {

  ns   <- as.integer(instance$n_sites)
  nj   <- as.integer(instance$n_storages)
  nk   <- as.integer(instance$n_consumers)
  Tm   <- as.integer(instance$n_periods)
  P    <- as.integer(instance$n_products)
  Amin <- as.integer(instance$min_age)
  Amax <- as.integer(instance$max_age)

  # ── Arc enumeration (mirrors build_lp_rcpp.cpp exactly) ──────────────────
  # s = 0 arcs: t in 1..(Tm+1), no age filter — includes "never-plant" (0, Tm+1)
  # s >= 1 arcs: kept iff Amin <= (t-s) <= Amax
  arcs <- vector("list", 0L)
  for (s in 0:(Tm + 1)) {
    for (t in (s + 1):(Tm + 1)) {
      if (t > Tm + 1) next
      if (s == 0) {
        arcs <- c(arcs, list(c(s = s, t = t)))
      } else {
        len <- t - s
        if (len >= Amin && len <= Amax)
          arcs <- c(arcs, list(c(s = s, t = t)))
      }
    }
  }
  arcs_df <- as.data.frame(do.call(rbind, arcs))
  arcs_df$arc_id <- seq_len(nrow(arcs_df)) - 1L   # 0-based arc index
  n_arcs <- nrow(arcs_df)

  # ── Harvest periods — mirrors cpp: t_min_harv = max(1, Amin) + 1 ─────────
  t_min_harv <- max(1L, Amin) + 1L
  Tharv <- t_min_harv:Tm
  nTh   <- length(Tharv)

  # ── Product pairs: (src, del) with del >= src (both 0-based) ─────────────
  pp_pairs <- data.frame(
    pi       = integer(0),
    src_prod = integer(0),
    del_prod = integer(0)
  )
  pi <- 0L
  for (p in 0:(P - 1)) {
    for (pp in p:(P - 1)) {
      pp_pairs <- rbind(pp_pairs,
                        data.frame(pi = pi, src_prod = p, del_prod = pp))
      pi <- pi + 1L
    }
  }
  n_ppairs <- nrow(pp_pairs)

  # ── Variable block sizes ─────────────────────────────────────────────────
  n_z   <- ns * n_arcs
  n_Xij <- ns * nj * P * nTh
  n_S   <- nj * P * nTh
  n_Xjk <- nj * nk * n_ppairs * nTh

  off_z   <- 0L
  off_Xij <- off_z   + n_z
  off_S   <- off_Xij + n_Xij
  off_Xjk <- off_S   + n_S
  n_vars  <- off_Xjk + n_Xjk

  list(
    ns = ns, nj = nj, nk = nk, Tm = Tm, P = P, Amin = Amin, Amax = Amax,
    arcs     = arcs_df,   n_arcs   = n_arcs,
    Tharv    = Tharv,     nTh      = nTh,
    pp_pairs = pp_pairs,  n_ppairs = n_ppairs,
    off_z    = off_z,     off_Xij  = off_Xij,
    off_S    = off_S,     off_Xjk  = off_Xjk,
    n_z      = n_z, n_Xij = n_Xij, n_S = n_S, n_Xjk = n_Xjk,
    n_vars   = n_vars,
    sites    = as.data.frame(instance$sites),
    storages = as.data.frame(instance$storages)
  )
}


# ==============================================================================
# .col_z / .col_Xij / .col_S / .col_Xjk  — vectorised 1-based column accessors
# (all inputs are 0-based indices, matching cpp col_* lambdas)
# ==============================================================================
.col_z <- function(idx, i, arc) {
  idx$off_z + i * idx$n_arcs + arc + 1L   # → 1-based R index
}
.col_Xij <- function(idx, i, j, p, th) {
  idx$off_Xij + i*(idx$nj*idx$P*idx$nTh) + j*(idx$P*idx$nTh) +
    p*idx$nTh + th + 1L
}
.col_S <- function(idx, j, p, th) {
  idx$off_S + j*(idx$P*idx$nTh) + p*idx$nTh + th + 1L
}
.col_Xjk <- function(idx, j, k, pi, th) {
  idx$off_Xjk + j*(idx$nk*idx$n_ppairs*idx$nTh) + k*(idx$n_ppairs*idx$nTh) +
    pi*idx$nTh + th + 1L
}


# ==============================================================================
# extract_z()
#
# Returns all active (value > tol) arc decisions as a tidy tibble.
#
# Output columns:
#   site_id (1-based), arc_id (0-based), s, t, arc_type, age_len,
#   value, area_ha
#   arc_type:
#     "establishment"  s=0,  t in 1..Tm
#     "never_plant"    s=0,  t = Tm+1
#     "harvest"        s>=1, t in 1..Tm   (active plantation → harvest)
#     "termination"    s>=1, t = Tm+1     (plantation ends at model horizon)
# ==============================================================================
extract_z <- function(res, idx, tol = 0.5) {
  x <- res$x
  rows <- list()
  for (i in 0:(idx$ns - 1)) {
    for (a in 0:(idx$n_arcs - 1)) {
      col <- .col_z(idx, i, a)
      val <- x[col]
      if (is.na(val) || val < tol) next
      arc <- idx$arcs[a + 1L, ]
      s   <- arc$s;  t <- arc$t
      rows[[length(rows) + 1L]] <- data.frame(
        site_id  = i + 1L,
        arc_id   = a,
        s        = s,
        t        = t,
        arc_type = dplyr::case_when(
          s == 0  & t == idx$Tm + 1 ~ "never_plant",
          s == 0                    ~ "establishment",
          t == idx$Tm + 1           ~ "termination",
          TRUE                      ~ "harvest"
        ),
        age_len  = t - s,
        value    = val
      )
    }
  }
  if (length(rows) == 0L) {
    return(tibble::tibble(
      site_id=integer(), arc_id=integer(), s=integer(), t=integer(),
      arc_type=character(), age_len=integer(), value=double(), area_ha=double()
    ))
  }
  df <- dplyr::bind_rows(rows)
  # Attach site area_ha for downstream KPI aggregation
  site_meta <- idx$sites %>%
    dplyr::mutate(site_id = dplyr::row_number()) %>%
    dplyr::select(site_id, area_ha)
  df <- dplyr::left_join(df, site_meta, by = "site_id")
  tibble::as_tibble(df)
}


# ==============================================================================
# extract_Xij()
#
# Returns all positive biomass flows from sites (i) to storage hubs (j).
#
# Output columns:
#   site_id (1-based), hub_id (1-based), product (1-based), period,
#   th (0-based index into Tharv), value
# ==============================================================================
extract_Xij <- function(res, idx, tol = 1e-6) {
  x <- res$x
  rows <- list()
  for (i in 0:(idx$ns - 1))
    for (j in 0:(idx$nj - 1))
      for (p in 0:(idx$P - 1))
        for (th in 0:(idx$nTh - 1)) {
          col <- .col_Xij(idx, i, j, p, th)
          val <- x[col]
          if (is.na(val) || val < tol) next
          rows[[length(rows) + 1L]] <- data.frame(
            site_id = i + 1L,
            hub_id  = j + 1L,
            product = p + 1L,
            period  = idx$Tharv[th + 1L],
            th      = th,
            value   = val
          )
        }
  if (length(rows) == 0L)
    return(tibble::tibble(
      site_id=integer(), hub_id=integer(), product=integer(),
      period=integer(), th=integer(), value=double()
    ))
  tibble::as_tibble(dplyr::bind_rows(rows))
}


# ==============================================================================
# extract_S()
#
# Returns hub inventory levels S(j, p, th) above tol.
#
# Output columns:
#   hub_id (1-based), product (1-based), period, th (0-based), value
# ==============================================================================
extract_S <- function(res, idx, tol = 1e-6) {
  x <- res$x
  rows <- list()
  for (j in 0:(idx$nj - 1))
    for (p in 0:(idx$P - 1))
      for (th in 0:(idx$nTh - 1)) {
        col <- .col_S(idx, j, p, th)
        val <- x[col]
        if (is.na(val) || val < tol) next
        rows[[length(rows) + 1L]] <- data.frame(
          hub_id  = j + 1L,
          product = p + 1L,
          period  = idx$Tharv[th + 1L],
          th      = th,
          value   = val
        )
      }
  if (length(rows) == 0L)
    return(tibble::tibble(
      hub_id=integer(), product=integer(),
      period=integer(), th=integer(), value=double()
    ))
  tibble::as_tibble(dplyr::bind_rows(rows))
}


# ==============================================================================
# extract_Xjk()
#
# Returns all positive deliveries from hubs (j) to consumers (k).
#
# Output columns:
#   hub_id, consumer_id (both 1-based), pi (0-based pair index in pp_pairs),
#   src_product, del_product (both 1-based), period, th (0-based), value
# ==============================================================================
extract_Xjk <- function(res, idx, tol = 1e-6) {
  x <- res$x
  rows <- list()
  for (j in 0:(idx$nj - 1))
    for (k in 0:(idx$nk - 1))
      for (pi in 0:(idx$n_ppairs - 1))
        for (th in 0:(idx$nTh - 1)) {
          col <- .col_Xjk(idx, j, k, pi, th)
          val <- x[col]
          if (is.na(val) || val < tol) next
          rows[[length(rows) + 1L]] <- data.frame(
            hub_id      = j + 1L,
            consumer_id = k + 1L,
            pi          = pi,
            src_product = idx$pp_pairs$src_prod[pi + 1L] + 1L,
            del_product = idx$pp_pairs$del_prod[pi + 1L] + 1L,
            period      = idx$Tharv[th + 1L],
            th          = th,
            value       = val
          )
        }
  if (length(rows) == 0L)
    return(tibble::tibble(
      hub_id=integer(), consumer_id=integer(), pi=integer(),
      src_product=integer(), del_product=integer(),
      period=integer(), th=integer(), value=double()
    ))
  tibble::as_tibble(dplyr::bind_rows(rows))
}


# ==============================================================================
# extract_result()
#
# Master extractor.  Decodes all four variable blocks and computes scalar KPIs.
#
# Args:
#   res      — afs_gurobi_result (from solve_lp_with_gurobi)
#   instance — the same instance list used to build the LP
#   tol      — threshold below which continuous values are treated as zero
#
# Returns: named list of class "afs_extracted_result" with:
#   $status, $objval, $mipgap, $runtime
#   $idx          — index tables from build_index_tables()
#   $z            — tibble: planting/harvest arc decisions
#   $Xij          — tibble: site → hub biomass flows
#   $S            — tibble: hub inventory levels
#   $Xjk          — tibble: hub → consumer deliveries
#   $kpi          — named list of scalar performance indicators
# ==============================================================================
extract_result <- function(res, instance, tol = 1e-6) {

  if (!inherits(res, "afs_gurobi_result"))
    stop("[extract_result] 'res' must be an afs_gurobi_result object.")

  feasible_status <- c("OPTIMAL", "SUBOPTIMAL", "TIME_LIMIT")
  if (!res$status %in% feasible_status || is.null(res$x)) {
    warning("[extract_result] No incumbent solution (status: ", res$status,
            "). Returning skeleton.")
    idx <- build_index_tables(instance)
    return(structure(
      list(status = res$status, objval = NA_real_,
           mipgap = res$result$mipgap, runtime = res$runtime,
           idx = idx, z = NULL, Xij = NULL, S = NULL, Xjk = NULL, kpi = list()),
      class = "afs_extracted_result"
    ))
  }

  idx <- build_index_tables(instance)

  z   <- extract_z(res,   idx, tol = 0.5)
  Xij <- extract_Xij(res, idx, tol)
  S   <- extract_S(res,   idx, tol)
  Xjk <- extract_Xjk(res, idx, tol)

  # ── Derived KPIs ───────────────────────────────────────────────────────────
  n_planted     <- z %>% dplyr::filter(arc_type == "establishment") %>%
                     dplyr::pull(site_id) %>% unique() %>% length()
  n_never_plant <- z %>% dplyr::filter(arc_type == "never_plant")  %>% nrow()
  n_harvests    <- z %>% dplyr::filter(arc_type == "harvest")      %>% nrow()

  total_area_ha <- z %>%
    dplyr::filter(arc_type == "establishment") %>%
    dplyr::group_by(site_id) %>%
    dplyr::summarise(area_ha = dplyr::first(area_ha), .groups = "drop") %>%
    dplyr::pull(area_ha) %>% sum(na.rm = TRUE)

  total_biomass_ij <- sum(Xij$value, na.rm = TRUE)
  total_delivery   <- sum(Xjk$value, na.rm = TRUE)
  total_inventory  <- sum(S$value,   na.rm = TRUE)

  # Revenue: price[k, del_product] * Xjk value
  total_revenue <- tryCatch({
    cp <- as.data.frame(instance$consumer_prices)
    if (ncol(cp) >= 3 && nrow(Xjk) > 0) {
      names(cp)[1:3] <- c("consumer_id", "del_product", "price")
      cp$consumer_id <- as.integer(cp$consumer_id)
      cp$del_product <- as.integer(cp$del_product)
      rev_df <- Xjk %>%
        dplyr::left_join(cp, by = c("consumer_id", "del_product")) %>%
        dplyr::mutate(rev = value * dplyr::coalesce(price, 0))
      sum(rev_df$rev, na.rm = TRUE)
    } else NA_real_
  }, error = function(e) NA_real_)

  kpi <- list(
    n_sites_planted   = n_planted,
    n_sites_never     = n_never_plant,
    n_harvest_events  = n_harvests,
    total_area_ha     = total_area_ha,
    total_biomass_ij  = total_biomass_ij,
    total_delivery    = total_delivery,
    total_inventory   = total_inventory,
    total_revenue     = total_revenue,
    objval            = res$objval,
    mipgap_pct        = if (!is.null(res$result$mipgap))
                          res$result$mipgap * 100 else NA_real_,
    runtime_s         = res$runtime,
    n_vars            = idx$n_vars,
    n_z               = idx$n_z,
    n_Xij             = idx$n_Xij,
    n_S               = idx$n_S,
    n_Xjk             = idx$n_Xjk
  )

  structure(
    list(
      status  = res$status,
      objval  = res$objval,
      mipgap  = res$result$mipgap,
      runtime = res$runtime,
      idx     = idx,
      z       = z,
      Xij     = Xij,
      S       = S,
      Xjk     = Xjk,
      kpi     = kpi
    ),
    class = "afs_extracted_result"
  )
}


# ==============================================================================
# print.afs_extracted_result()  — S3 method
# ==============================================================================
print.afs_extracted_result <- function(x, ...) {
  summarise_result(x)
  invisible(x)
}


# ==============================================================================
# summarise_result()
#
# Prints a structured one-page KPI summary to the console.
# ==============================================================================
summarise_result <- function(extracted) {
  k <- extracted$kpi
  cat("══════════════════════════════════════════════════════\n")
  cat("  AFS MILP — Solution Summary\n")
  cat("══════════════════════════════════════════════════════\n")
  cat(sprintf("  Status       : %s\n",     extracted$status))
  cat(sprintf("  Objective    : %.2f\n",   k$objval   %||% NA))
  cat(sprintf("  MIP gap      : %.3f%%\n", k$mipgap_pct %||% NA))
  cat(sprintf("  Runtime      : %.1f s\n", k$runtime_s %||% NA))
  cat("──────────────────────────────────────────────────────\n")
  cat("  Planting decisions\n")
  cat(sprintf("    Sites planted  : %d  (never-plant: %d)\n",
              k$n_sites_planted %||% 0L, k$n_sites_never %||% 0L))
  cat(sprintf("    Harvest events : %d\n", k$n_harvest_events %||% 0L))
  cat(sprintf("    Total area     : %.1f ha\n", k$total_area_ha %||% 0))
  cat("  Biomass flows\n")
  cat(sprintf("    Site→Hub  Xij  : %.1f t (total raw shipped)\n",
              k$total_biomass_ij %||% 0))
  cat(sprintf("    Hub→Consumer   : %.1f t (total delivered)\n",
              k$total_delivery   %||% 0))
  cat(sprintf("    Peak inventory : %.1f t (sum S over all j×p×t)\n",
              k$total_inventory  %||% 0))
  if (!is.na(k$total_revenue %||% NA_real_))
    cat(sprintf("    Total revenue  : %.2f\n", k$total_revenue))
  cat("  Model dimensions\n")
  cat(sprintf("    n_vars: %d  (z:%d  Xij:%d  S:%d  Xjk:%d)\n",
              k$n_vars %||% 0L, k$n_z %||% 0L, k$n_Xij %||% 0L,
              k$n_S %||% 0L, k$n_Xjk %||% 0L))
  cat("══════════════════════════════════════════════════════\n")
  invisible(extracted)
}

# Null-coalescing helper — avoids purrr dependency
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L) a else b


# ==============================================================================
# validate_extracted()
#
# Post-solve consistency checks on the decoded tibbles.
#
# Checks implemented:
#   V1  C3 linkage   — every (site, period) in Xij has a matching harvest arc z
#   V2  C4 balance   — S(t) == S(t-1) + inflow(Xij) - outflow(Xjk) per hub×product×t
#   V4  Cascade dir  — del_product >= src_product in all Xjk rows
#   V5  Cap storage  — sum_p S(j,p,t) <= CAP_stor[j] for all j, t
#
# Args:
#   extracted  — afs_extracted_result from extract_result()
#   tol        — numerical tolerance for equality / inequality checks
#
# Returns:
#   list(pass = logical,
#        violations = named list of data frames — empty if all pass)
# ==============================================================================
validate_extracted <- function(extracted, tol = 1e-4) {

  violations <- list()
  idx <- extracted$idx
  z   <- extracted$z
  Xij <- extracted$Xij
  S   <- extracted$S
  Xjk <- extracted$Xjk

  # ── V1: C3 linkage — every Xij (site, period) must have a harvest arc ────
  if (!is.null(Xij) && nrow(Xij) > 0 && !is.null(z) && nrow(z) > 0) {
    harvest_pairs <- z %>%
      dplyr::filter(arc_type == "harvest") %>%
      dplyr::select(site_id, period = t) %>%
      dplyr::distinct()
    Xij_sp <- Xij %>%
      dplyr::group_by(site_id, period) %>%
      dplyr::summarise(total_flow = sum(value), .groups = "drop")
    orphan <- dplyr::anti_join(Xij_sp, harvest_pairs,
                               by = c("site_id", "period"))
    if (nrow(orphan) > 0)
      violations$V1_orphan_Xij <- orphan
  }

  # ── V2: Hub inventory mass balance (C4) ──────────────────────────────────
  if (!is.null(S) && nrow(S) > 0) {
    bal_rows <- list()
    for (j in seq_len(idx$nj)) {
      for (p in seq_len(idx$P)) {
        for (th_idx in seq_along(idx$Tharv)) {
          t      <- idx$Tharv[th_idx]
          s_now  <- sum(S$value[S$hub_id == j & S$product == p & S$period == t],
                        na.rm = TRUE)
          s_prev <- if (th_idx > 1) {
            t_prev <- idx$Tharv[th_idx - 1L]
            sum(S$value[S$hub_id == j & S$product == p & S$period == t_prev],
                na.rm = TRUE)
          } else 0.0
          inflow <- if (!is.null(Xij) && nrow(Xij) > 0)
            sum(Xij$value[Xij$hub_id == j & Xij$product == p & Xij$period == t],
                na.rm = TRUE) else 0.0
          outflow <- if (!is.null(Xjk) && nrow(Xjk) > 0)
            sum(Xjk$value[Xjk$hub_id == j & Xjk$src_product == p & Xjk$period == t],
                na.rm = TRUE) else 0.0
          residual <- s_now - s_prev - inflow + outflow
          if (abs(residual) > tol)
            bal_rows[[length(bal_rows) + 1L]] <- data.frame(
              hub_id = j, product = p, period = t,
              S_now = s_now, S_prev = s_prev,
              inflow = inflow, outflow = outflow, residual = residual
            )
        }
      }
    }
    if (length(bal_rows) > 0)
      violations$V2_inventory_balance <- dplyr::bind_rows(bal_rows)
  }

  # ── V4: Cascade direction — delivered product index >= source product ─────
  if (!is.null(Xjk) && nrow(Xjk) > 0) {
    bad_dir <- Xjk %>% dplyr::filter(del_product < src_product)
    if (nrow(bad_dir) > 0)
      violations$V4_cascade_direction <- bad_dir
  }

  # ── V5: Storage capacity feasibility ─────────────────────────────────────
  if (!is.null(S) && nrow(S) > 0 && "CAP_stor" %in% names(idx$storages)) {
    cap_df <- idx$storages %>%
      dplyr::mutate(hub_id = dplyr::row_number()) %>%
      dplyr::select(hub_id, CAP_stor)
    cap_viol <- S %>%
      dplyr::group_by(hub_id, period) %>%
      dplyr::summarise(total_S = sum(value), .groups = "drop") %>%
      dplyr::left_join(cap_df, by = "hub_id") %>%
      dplyr::filter(total_S > CAP_stor + tol)
    if (nrow(cap_viol) > 0)
      violations$V5_storage_capacity <- cap_viol
  }

  pass <- length(violations) == 0L

  if (pass) {
    cat("[validate_extracted] \u2713 All checks passed.\n")
  } else {
    cat(sprintf("[validate_extracted] \u2717 %d check(s) failed:\n",
                length(violations)))
    for (nm in names(violations))
      cat(sprintf("  \u2022 %s  (%d row(s))\n", nm, nrow(violations[[nm]])))
  }

  list(pass = pass, violations = violations)
}
