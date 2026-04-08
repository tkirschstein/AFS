// [[file: build_lp_rcpp.cpp]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <algorithm>
using namespace Rcpp;

// ============================================================================
// build_lp_rcpp.cpp
// Rcpp LP builder — mirrors build_agroforestry_lp_sparse_v10_optimized() in
// !build_AFS_milp.R exactly.
//
// Variable layout (0-based column offsets):
//   [off_z  .. off_z  + n_z  - 1]  z[i, arc]       binary
//   [off_Xij.. off_Xij + n_Xij-1]  Xij[i,j,p,th]   continuous
//   [off_S  .. off_S  + n_S  - 1]  S[j,p,th]        continuous
//   [off_Xjk.. off_Xjk + n_Xjk-1] Xjk[j,k,pi,th]  continuous
//
// NOTE: Variable Y has been removed. C3 directly constrains
//       sum_j Xij(i,j,p,t) <= sum_{s>=1,s<t} eta_{p,t-s} * AREA_i * z(i,s,t)
//
// Constraints (labels match !build_AFS_milp.R):
//   C1  Path establishment
//   C2  Path connectivity
//   C3  Biomass yield + shipping bound  (merged, <=)
//   C6  Inventory balance               (== , only Tharv)
//   C7  Storage capacity                (<= , only Tharv)
//   C8  Processing capacity             (<= , only Tharv)
//   C9  Demand with cascade             (<= )
//
// Authors: SmartAgroforst 2026
// ============================================================================

// [[Rcpp::export]]
List build_lp_rcpp(List instance) {

  // ── Scalars ────────────────────────────────────────────────────────────────
  int    ns      = as<int>(instance["n_sites"]);
  int    nj      = as<int>(instance["n_storages"]);
  int    nk      = as<int>(instance["n_consumers"]);
  int    Tm      = as<int>(instance["n_periods"]);
  int    P       = as<int>(instance["n_products"]);
  int    Amax    = as<int>(instance["max_age"]);
  int    Amin    = as<int>(instance["min_age"]);
  double Copp     = as<double>(instance["c_opp"]);
  double c_tr_raw = as<double>(instance["c_tr_raw"]);
  double c_tr_pre = as<double>(instance["c_tr_pre"]);

  Rcpp::Rcout << "Building sparse LP (Rcpp v10-no-Y)...\n";

  // ── Site / storage vectors ─────────────────────────────────────────────────
  DataFrame sites_df  = as<DataFrame>(instance["sites"]);
  NumericVector area_ha = sites_df["area_ha"];
  NumericVector C_est   = sites_df["C_est"];
  NumericVector C_harv  = sites_df["C_harv"];

  DataFrame storages_df = as<DataFrame>(instance["storages"]);
  NumericVector CAP_stor = storages_df["CAP_stor"];
  NumericVector CAP_proc = storages_df["CAP_proc"];
  NumericVector c_stor   = storages_df["c_stor"];

  NumericMatrix d_ij = as<NumericMatrix>(instance["d_ij"]);  // [ns x nj]
  NumericMatrix d_jk = as<NumericMatrix>(instance["d_jk"]);  // [nj x nk]

  // ── Yield matrix [P x Tm] (1-based age → column age-1) ────────────────────
  // Passed from R wrapper as matrix(nrow=P, ncol=Tm)
  NumericMatrix yield_matrix = as<NumericMatrix>(instance["yield_matrix"]);

  // ── Consumer prices look-up (k-1)*P + (pp-1) → price ─────────────────────
  DataFrame cp       = as<DataFrame>(instance["consumer_prices"]);
  IntegerVector cp_k  = cp[0];
  IntegerVector cp_pp = cp[1];
  NumericVector cp_pr = cp[2];
  NumericVector price_lut(nk * P, 0.0);
  for (int r = 0; r < cp_k.size(); r++)
    price_lut[(cp_k[r]-1) * P + (cp_pp[r]-1)] = cp_pr[r];

  // ── Demand look-up (k-1)*P*Tm + (p-1)*Tm + (t-1) → D_max ────────────────
  DataFrame dem_df     = as<DataFrame>(instance["demand"]);
  IntegerVector dem_k  = dem_df["consumer_id"];
  IntegerVector dem_p  = dem_df["product"];
  IntegerVector dem_t  = dem_df["period"];
  NumericVector dem_Dmax = dem_df["D_max"];
  NumericVector dem_lut(nk * P * Tm, -1.0);
  for (int r = 0; r < dem_k.size(); r++)
    dem_lut[(dem_k[r]-1)*P*Tm + (dem_p[r]-1)*Tm + (dem_t[r]-1)] = dem_Dmax[r];

  // ==========================================================================
  // STEP 1: VARIABLE INDEXING
  // ==========================================================================

  // ── Valid arcs (s,t) ───────────────────────────────────────────────────────
  // Mirrors z_tuples construction in !build_AFS_milp.R:
  //   s == 0 && t in 1..Tm                         → establishment arcs
  //   s >= 1 && t == Tm+1                          → termination arcs
  //   s >= 1 && t in 1..Tm && (t-s) in [Amin,Amax]→ harvest arcs
  struct Arc { int s, t; };
  std::vector<Arc> arcs;
  arcs.reserve((Tm + 2) * (Amax - Amin + 2));
  for (int s = 0; s <= Tm + 1; s++) {
    for (int t = s + 1; t <= Tm + 1; t++) {
      int len = t - s;
      if (s == 0 && t >= 1 && t <= Tm) {
        arcs.push_back({s, t}); continue;
      }
      if (t == Tm + 1 && s >= 1) {
        arcs.push_back({s, t}); continue;
      }
      if (s >= 1 && t >= 1 && t <= Tm && len >= Amin && len <= Amax) {
        arcs.push_back({s, t});
      }
    }
  }
  int n_arcs = (int)arcs.size();

  // ── Harvest periods: Tset[Tset > max(1, Amin)] ────────────────────────────
  std::vector<int> Tharv;
  int t_min_harv = std::max(1, Amin) + 1;
  for (int t = t_min_harv; t <= Tm; t++) Tharv.push_back(t);
  int nTh = (int)Tharv.size();

  // ── Variable counts ───────────────────────────────────────────────────────
  // No Y variable — matches !build_AFS_milp.R
  int n_z   = ns * n_arcs;
  int n_Xij = ns * nj * P * nTh;
  int n_S   = nj * P * nTh;          // indexed over Tharv only (matches R)

  // Xjk: only pp >= p (0-based: second >= first)
  std::vector<std::pair<int,int>> pp_pairs;
  for (int p = 0; p < P; p++)
    for (int pp = p; pp < P; pp++)
      pp_pairs.push_back({p, pp});
  int n_ppairs = (int)pp_pairs.size();
  int n_Xjk = nj * nk * n_ppairs * nTh;

  int off_z   = 0;
  int off_Xij = off_z   + n_z;
  int off_S   = off_Xij + n_Xij;
  int off_Xjk = off_S   + n_S;
  int n_vars  = off_Xjk + n_Xjk;

  Rcpp::Rcout << "  n_vars=" << n_vars
              << " (z:" << n_z
              << " Xij:" << n_Xij
              << " S:" << n_S
              << " Xjk:" << n_Xjk << ")\n";

  // ── Column index lambdas (0-based) ────────────────────────────────────────
  auto col_z = [&](int i, int arc) {
    return off_z + i * n_arcs + arc;
  };
  auto col_Xij = [&](int i, int j, int p, int th) {
    return off_Xij + i*(nj*P*nTh) + j*(P*nTh) + p*nTh + th;
  };
  // S indexed over th (position in Tharv), not raw period
  auto col_S = [&](int j, int p, int th) {
    return off_S + j*(P*nTh) + p*nTh + th;
  };
  auto col_Xjk = [&](int j, int k, int pi, int th) {
    return off_Xjk + j*(nk*n_ppairs*nTh) + k*(n_ppairs*nTh) + pi*nTh + th;
  };

  // Helper: find th index for a given 1-based period; returns -1 if not in Tharv
  auto find_th = [&](int t1) -> int {
    for (int th = 0; th < nTh; th++)
      if (Tharv[th] == t1) return th;
    return -1;
  };

  // ==========================================================================
  // STEP 2: OBJECTIVE VECTOR
  // Mirrors Step 2 of !build_AFS_milp.R exactly
  // ==========================================================================
  std::vector<double> c_vec(n_vars, 0.0);

  // Revenue: price[k, pp] on Xjk (pp = pp_pairs[pi].second, 0-based)
  for (int j = 0; j < nj; j++)
    for (int k = 0; k < nk; k++)
      for (int pi = 0; pi < n_ppairs; pi++) {
        int pp = pp_pairs[pi].second;  // 0-based product index
        double pr = price_lut[k * P + pp];
        for (int th = 0; th < nTh; th++)
          c_vec[col_Xjk(j, k, pi, th)] += pr;
      }

  // Establishment cost + opportunity cost: z(i, 0, t), t in 1..Tm
  for (int a = 0; a < n_arcs; a++) {
    if (arcs[a].s != 0) continue;
    int t = arcs[a].t;
    if (t < 1 || t > Tm) continue;
    for (int i = 0; i < ns; i++) {
      int c = col_z(i, a);
      c_vec[c] -= C_est[i];                           // establishment cost
      c_vec[c] -= Copp * area_ha[i] * (Tm - t);      // opportunity cost
    }
  }

  // Harvest cost: z(i, s>=1, t in 1..Tm)
  for (int a = 0; a < n_arcs; a++) {
    if (arcs[a].s < 1 || arcs[a].t > Tm) continue;
    for (int i = 0; i < ns; i++)
      c_vec[col_z(i, a)] -= C_harv[i] * area_ha[i];
  }

  // Transport site → hub: -c_tr_raw * d_ij[i,j]
  for (int i = 0; i < ns; i++)
    for (int j = 0; j < nj; j++) {
      double cost = c_tr_raw * d_ij(i, j);
      for (int p = 0; p < P; p++)
        for (int th = 0; th < nTh; th++)
          c_vec[col_Xij(i, j, p, th)] -= cost;
    }

  // Transport hub → consumer: -c_tr_pre * d_jk[j,k]
  for (int j = 0; j < nj; j++)
    for (int k = 0; k < nk; k++) {
      double cost = c_tr_pre * d_jk(j, k);
      for (int pi = 0; pi < n_ppairs; pi++)
        for (int th = 0; th < nTh; th++)
          c_vec[col_Xjk(j, k, pi, th)] -= cost;
    }

  // Storage cost: -c_stor[j] per unit S(j,p,th)
  for (int j = 0; j < nj; j++)
    for (int p = 0; p < P; p++)
      for (int th = 0; th < nTh; th++)
        c_vec[col_S(j, p, th)] -= c_stor[j];

  // ==========================================================================
  // STEP 3: CONSTRAINTS
  // ==========================================================================
  int est_nnz = n_z * 4 + n_Xij * 3 + n_S * 4 + n_Xjk * 2;
  std::vector<int>    row_v; row_v.reserve(est_nnz);
  std::vector<int>    col_v; col_v.reserve(est_nnz);
  std::vector<double> val_v; val_v.reserve(est_nnz);
  std::vector<double> rhs_v;
  std::vector<std::string> sense_v;
  int nrow = 0;  // constraint counter (0-based)

  // push one nonzero (1-based indices for R slam)
  auto push = [&](int r, int c, double v) {
    row_v.push_back(r + 1);
    col_v.push_back(c + 1);
    val_v.push_back(v);
  };
  auto add_con = [&](const std::string& sense, double rhs) {
    sense_v.push_back(sense);
    rhs_v.push_back(rhs);
    nrow++;
  };

  // ── C1: Path establishment: sum_t z(i,0,t) <= 1 ──────────────────────────
  for (int i = 0; i < ns; i++) {
    bool any = false;
    for (int a = 0; a < n_arcs; a++) {
      if (arcs[a].s != 0 || arcs[a].t < 1 || arcs[a].t > Tm) continue;
      push(nrow, col_z(i, a), 1.0);
      any = true;
    }
    if (any) add_con("<=", 1.0);
  }
  Rcpp::Rcout << "  C1 done\n";

  // ── C2: Path connectivity: sum_{s<t} z(i,s,t) = sum_{u>t} z(i,t,u) ───────
  for (int i = 0; i < ns; i++) {
    for (int t = 1; t <= Tm; t++) {
      bool any = false;
      for (int a = 0; a < n_arcs; a++) {
        if (arcs[a].t == t && arcs[a].s < t) {
          push(nrow, col_z(i, a), +1.0); any = true;
        } else if (arcs[a].s == t && arcs[a].t > t) {
          push(nrow, col_z(i, a), -1.0); any = true;
        }
      }
      if (any) add_con("==", 0.0);
    }
  }
  Rcpp::Rcout << "  C2 done\n";

  // ── C3: Combined yield + shipping bound (NO Y variable) ───────────────────
  // Mirrors !build_AFS_milp.R C3 exactly:
  //   sum_j Xij(i,j,p,t) - sum_{s>=1, s<t} eta_{p, t-s} * AREA_i * z(i,s,t) <= 0
  // Age index: age = t - s  (1-based), yield_matrix is [P x Tm] (1-based cols)
  for (int i = 0; i < ns; i++) {
    for (int p = 0; p < P; p++) {
      for (int th = 0; th < nTh; th++) {
        int t = Tharv[th];
        // Collect harvest arcs s>=1, s<t arriving at t
        std::vector<int> harv_arcs;
        for (int a = 0; a < n_arcs; a++) {
          if (arcs[a].t == t && arcs[a].s >= 1 && arcs[a].s < t)
            harv_arcs.push_back(a);
        }
        if (harv_arcs.empty()) continue;

        // + sum_j Xij coefficients
        for (int j = 0; j < nj; j++)
          push(nrow, col_Xij(i, j, p, th), 1.0);

        // - eta_{p, age} * area_ha * z coefficients
        for (int a : harv_arcs) {
          int age = t - arcs[a].s;        // 1-based age (1..Amax)
          if (age < 1 || age > Tm) continue;
          double coef = yield_matrix(p, age - 1) * area_ha[i];  // age-1 = 0-based col
          if (coef > 0.0)
            push(nrow, col_z(i, a), -coef);
        }
        add_con("<=", 0.0);
      }
    }
  }
  Rcpp::Rcout << "  C3 done\n";

  // ── C6: Inventory balance ─────────────────────────────────────────────────
  // S(j,p,th) = S(j,p,th-1) + sum_i Xij(i,j,p,th) - sum_k sum_{pp>=p} Xjk
  // Loop only over Tharv (matches R); S boundary at th=0 is implicitly zero.
  for (int j = 0; j < nj; j++) {
    for (int p = 0; p < P; p++) {
      for (int th = 0; th < nTh; th++) {
        // + S(j,p,th)
        push(nrow, col_S(j, p, th), 1.0);
        // - S(j,p,th-1)  if th > 0
        if (th > 0)
          push(nrow, col_S(j, p, th - 1), -1.0);
        // - sum_i Xij(i,j,p,th)
        for (int i = 0; i < ns; i++)
          push(nrow, col_Xij(i, j, p, th), -1.0);
        // + sum_k sum_{pi: pp_pairs[pi].first == p} Xjk
        for (int k = 0; k < nk; k++)
          for (int pi = 0; pi < n_ppairs; pi++) {
            if (pp_pairs[pi].first != p) continue;
            push(nrow, col_Xjk(j, k, pi, th), +1.0);
          }
        add_con("==", 0.0);
      }
    }
  }
  Rcpp::Rcout << "  C6 done\n";

  // ── C7: Storage capacity: sum_p S(j,p,th) <= CAP_stor[j] ─────────────────
  // Loop over Tharv only (matches R — S_tuples uses Tharv)
  for (int j = 0; j < nj; j++) {
    for (int th = 0; th < nTh; th++) {
      for (int p = 0; p < P; p++)
        push(nrow, col_S(j, p, th), 1.0);
      add_con("<=", CAP_stor[j]);
    }
  }
  Rcpp::Rcout << "  C7 done\n";

  // ── C8: Processing capacity: sum_{i,p} Xij(i,j,p,th) <= CAP_proc[j] ──────
  for (int j = 0; j < nj; j++) {
    for (int th = 0; th < nTh; th++) {
      for (int i = 0; i < ns; i++)
        for (int p = 0; p < P; p++)
          push(nrow, col_Xij(i, j, p, th), 1.0);
      add_con("<=", CAP_proc[j]);
    }
  }
  Rcpp::Rcout << "  C8 done\n";

  // ── C9: Demand with cascade ───────────────────────────────────────────────
  // Paper C8: sum_{j, p' <= p} Xjk(j,k,p',p,t) <= D_max(k,p,t)
  // In pp_pairs: pair (p', p) means source=p', delivered=p.
  // Cascade: delivered product pp_pairs[pi].second must equal demanded product p;
  //          source product pp_pairs[pi].first <= p is implicitly guaranteed
  //          by the pp_pairs construction (first <= second).
  for (int k = 0; k < nk; k++) {
    for (int p = 0; p < P; p++) {
      for (int th = 0; th < nTh; th++) {
        int t1   = Tharv[th];
        double Dmax = dem_lut[k * P * Tm + p * Tm + (t1 - 1)];
        if (Dmax < 0.0) continue;  // no demand record → skip
        bool any = false;
        for (int j = 0; j < nj; j++)
          for (int pi = 0; pi < n_ppairs; pi++) {
            if (pp_pairs[pi].second != p) continue;  // delivered product must match
            push(nrow, col_Xjk(j, k, pi, th), 1.0);
            any = true;
          }
        if (any) add_con("<=", Dmax);
      }
    }
  }
  Rcpp::Rcout << "  C9 done\n";
  Rcpp::Rcout << "  Total constraints: " << nrow
              << ", nnz: " << (int)row_v.size() << "\n";

  // ==========================================================================
  // STEP 4: BOUNDS AND TYPES
  // ==========================================================================
  NumericVector ub(n_vars, R_PosInf);
  // z: binary → ub = 1
  for (int c = off_z; c < off_z + n_z; c++) ub[c] = 1.0;
  // S: per-product CAP_stor upper bound
  for (int j = 0; j < nj; j++)
    for (int p = 0; p < P; p++)
      for (int th = 0; th < nTh; th++)
        ub[col_S(j, p, th)] = CAP_stor[j];

  CharacterVector types(n_vars, "C");
  for (int c = off_z; c < off_z + n_z; c++) types[c] = "B";

  // ==========================================================================
  // STEP 5: RETURN
  // ==========================================================================
  return List::create(
    Named("row_idx")   = wrap(row_v),
    Named("col_idx")   = wrap(col_v),
    Named("val")       = wrap(val_v),
    Named("rhs")       = wrap(rhs_v),
    Named("sense")     = wrap(sense_v),
    Named("c_vec")     = wrap(c_vec),
    Named("ub")        = ub,
    Named("types")     = types,
    Named("n_vars")    = n_vars,
    Named("n_constrs") = nrow,
    Named("n_z")       = n_z,
    Named("n_Xij")     = n_Xij,
    Named("n_S")       = n_S,
    Named("n_Xjk")     = n_Xjk
  );
}
