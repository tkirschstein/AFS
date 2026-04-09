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
// Rcpp LP builder — mirrors build_AFS_milp() in !build_AFS_milp.R exactly.
//
// Variable layout (0-based column offsets):
//   [off_z  .. off_z  + n_z  - 1]  z[i, arc]       binary
//   [off_Xij.. off_Xij + n_Xij-1]  Xij[i,j,p,th]   continuous
//   [off_S  .. off_S  + n_S  - 1]  S[j,p,th]        continuous  (Tharv only)
//   [off_Xjk.. off_Xjk + n_Xjk-1] Xjk[j,k,pi,th]  continuous
//     where pp_pairs[pi] = (source_product, delivered_product), delivered >= source
//
// NOTE: Variable Y has been removed. C3 directly constrains
//       sum_j Xij(i,j,p,t) <= sum_{s>=1,s<t} eta_{p,t-s} * AREA_i * z(i,s,t)
//
// Constraints (labels match !build_AFS_milp.R):
//   C1  Path establishment       sum_t z(i,0,t) <= 1
//   C2  Path connectivity        flow conservation at each node t
//   C3  Biomass yield bound      sum_j Xij <= sum_s eta*area*z  (<=, Tharv)
//   C4  Inventory balance        S(t) = S(t-1) + Xij_in - Xjk_out  (==, Tharv)
//   C5  Storage capacity         sum_p S(j,p,t) <= CAP_stor[j]  (<=, Tharv)
//   C6  Processing capacity      sum_{i,p} Xij(i,j,p,t) <= CAP_proc[j]  (<=, Tharv)
//   C7  Demand with cascade      sum_{j,pprod<=p} Xjk(j,k,pprod,p,t) <= D_max  (<=, Tharv)
//
// Objective cost terms (matching !build_AFS_milp.R):
//   + price[k,pp]  * Xjk        revenue (delivered product pp)
//   - C_est[i]     * z(i,0,t)   establishment cost
//   - C_opp[i] * area_ha[i] * (Tm-t) * z(i,0,t)   opportunity cost  (PER SITE)
//   - C_main[i] * area_ha[i] * (t-s) * z(i,s,t)   maintenance cost  (PER SITE, s>=1)
//   - C_harv[i] * area_ha[i]   * z(i,s,t)          harvest cost      (s>=1)
//   - c_tr_raw * d_ij[i,j]     * Xij               raw transport cost
//   - c_tr_pre * d_jk[j,k]     * Xjk               pre-processed transport cost
//   - c_stor[j]                * S(j,p,t)          storage cost
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
  double c_tr_raw = as<double>(instance["c_tr_raw"]);
  double c_tr_pre = as<double>(instance["c_tr_pre"]);

  Rcpp::Rcout << "Building sparse LP (Rcpp v10-no-Y, matching !build_AFS_milp.R)...\n";

  // ── Site vectors ───────────────────────────────────────────────────────────
  DataFrame sites_df  = as<DataFrame>(instance["sites"]);
  NumericVector area_ha = sites_df["area_ha"];
  NumericVector C_est   = sites_df["C_est"];
  NumericVector C_harv  = sites_df["C_harv"];
  NumericVector C_main  = sites_df["C_main"];   // maintenance cost per ha per period
  NumericVector C_opp   = sites_df["C_opp"];    // per-site opportunity cost per ha per period

  // ── Storage vectors ────────────────────────────────────────────────────────
  DataFrame storages_df = as<DataFrame>(instance["storages"]);
  NumericVector CAP_stor = storages_df["CAP_stor"];
  NumericVector CAP_proc = storages_df["CAP_proc"];
  NumericVector c_stor   = storages_df["c_stor"];

  NumericMatrix d_ij = as<NumericMatrix>(instance["d_ij"]);  // [ns x nj]
  NumericMatrix d_jk = as<NumericMatrix>(instance["d_jk"]);  // [nj x nk]

  // ── Yield matrix [P x Tm] (1-based age → column age-1) ────────────────────
  // Passed from R wrapper as matrix(nrow=P, ncol=Tm)
  NumericMatrix yield_matrix = as<NumericMatrix>(instance["yield_matrix"]);

  // ── Consumer prices look-up (0-based): price_lut[k*P + pp] ───────────────
  DataFrame cp       = as<DataFrame>(instance["consumer_prices"]);
  IntegerVector cp_k  = cp[0];   // 1-based consumer id
  IntegerVector cp_pp = cp[1];   // 1-based product
  NumericVector cp_pr = cp[2];   // price
  NumericVector price_lut(nk * P, 0.0);
  for (int r = 0; r < cp_k.size(); r++)
    price_lut[(cp_k[r]-1) * P + (cp_pp[r]-1)] = cp_pr[r];

  // ── Demand look-up: dem_lut[(k-1)*P*Tm + (p-1)*Tm + (t-1)] = D_max ──────
  DataFrame dem_df       = as<DataFrame>(instance["demand"]);
  IntegerVector dem_k    = dem_df["consumer_id"];
  IntegerVector dem_p    = dem_df["product"];
  IntegerVector dem_t    = dem_df["period"];
  NumericVector dem_Dmax = dem_df["D_max"];
  NumericVector dem_lut(nk * P * Tm, -1.0);
  for (int r = 0; r < dem_k.size(); r++)
    dem_lut[(dem_k[r]-1)*P*Tm + (dem_p[r]-1)*Tm + (dem_t[r]-1)] = dem_Dmax[r];

  // ==========================================================================
  // STEP 1: VARIABLE INDEXING
  // ==========================================================================

  // ── Valid arcs (s,t) ───────────────────────────────────────────────────────
  // Mirrors z_tuples construction in !build_AFS_milp.R:
  //   s == 0 && t in 1..Tm                              → establishment arcs
  //   s >= 1 && t == Tm+1                               → termination arcs
  //   s >= 1 && t in 1..Tm && (t-s) in [Amin,Amax]     → harvest arcs
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

  // ── Variable counts (no Y variable) ──────────────────────────────────────
  int n_z   = ns * n_arcs;
  int n_Xij = ns * nj * P * nTh;
  int n_S   = nj * P * nTh;          // indexed over Tharv only

  // Xjk: pp_pairs where delivered >= source (both 0-based)
  std::vector<std::pair<int,int>> pp_pairs;  // (source, delivered)
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
  auto col_S = [&](int j, int p, int th) {
    return off_S + j*(P*nTh) + p*nTh + th;
  };
  auto col_Xjk = [&](int j, int k, int pi, int th) {
    return off_Xjk + j*(nk*n_ppairs*nTh) + k*(n_ppairs*nTh) + pi*nTh + th;
  };

  // ==========================================================================
  // STEP 2: OBJECTIVE VECTOR (mirrors !build_AFS_milp.R Step 2 exactly)
  // ==========================================================================
  std::vector<double> c_vec(n_vars, 0.0);

  // Revenue: price[k, pp] on Xjk — joined on DELIVERED product (pp_pairs[pi].second)
  for (int j = 0; j < nj; j++)
    for (int k = 0; k < nk; k++)
      for (int pi = 0; pi < n_ppairs; pi++) {
        int pp = pp_pairs[pi].second;  // 0-based delivered product
        double pr = price_lut[k * P + pp];
        for (int th = 0; th < nTh; th++)
          c_vec[col_Xjk(j, k, pi, th)] += pr;
      }

  // Establishment cost + opportunity cost: z(i, 0, t), t in 1..Tm
  // R: -C_est[ii]  and  -C_opp[ii] * area_ha[ii] * (Tm - t)   (PER SITE)
  for (int a = 0; a < n_arcs; a++) {
    if (arcs[a].s != 0) continue;
    int t = arcs[a].t;
    if (t < 1 || t > Tm) continue;
    for (int i = 0; i < ns; i++) {
      int c = col_z(i, a);
      c_vec[c] -= C_est[i];                              // establishment cost
      c_vec[c] -= C_opp[i] * area_ha[i] * (Tm - t);    // opportunity cost (per-site)
    }
  }

  // Maintenance cost + harvest cost: z(i, s>=1, t in 1..Tm), t > s
  // R: -C_main[ii] * area_ha[ii] * (t-s)   (maintenance, per arc length)
  //    -C_harv[ii] * area_ha[ii]            (harvest, per harvest event)
  for (int a = 0; a < n_arcs; a++) {
    if (arcs[a].s < 1 || arcs[a].t > Tm) continue;
    int arc_len = arcs[a].t - arcs[a].s;   // = t - s
    for (int i = 0; i < ns; i++) {
      int c = col_z(i, a);
      c_vec[c] -= C_main[i] * area_ha[i] * arc_len;    // maintenance cost
      c_vec[c] -= C_harv[i] * area_ha[i];              // harvest cost
    }
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
  int nrow = 0;

  // push one nonzero entry (1-based indices for R slam)
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

  // ── C1: Path establishment: sum_t z(i,0,t) <= 1  (one plantation per site)
  for (int i = 0; i < ns; i++) {
    bool any = false;
    for (int a = 0; a < n_arcs; a++) {
      if (arcs[a].s != 0 || arcs[a].t < 1 || arcs[a].t > Tm) continue;
      push(nrow, col_z(i, a), 1.0);
      any = true;
    }
    if (any) add_con("<=", 1.0);
  }
  Rcpp::Rcout << "  C1: Path establishment\n";

  // ── C2: Path connectivity: sum_{s<t} z(i,s,t) = sum_{u>t} z(i,t,u)
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
  Rcpp::Rcout << "  C2: Path connectivity\n";

  // ── C3: Biomass yield bound (NO Y variable) ───────────────────────────────
  // sum_j Xij(i,j,p,t) - sum_{s>=1,s<t} eta_{p,t-s} * area_i * z(i,s,t) <= 0
  // Age index: age = t - s  (1-based), yield_matrix[p, age-1] (0-based column)
  for (int i = 0; i < ns; i++) {
    for (int p = 0; p < P; p++) {
      for (int th = 0; th < nTh; th++) {
        int t = Tharv[th];
        std::vector<int> harv_arcs;
        for (int a = 0; a < n_arcs; a++) {
          if (arcs[a].t == t && arcs[a].s >= 1 && arcs[a].s < t)
            harv_arcs.push_back(a);
        }
        if (harv_arcs.empty()) continue;
        for (int j = 0; j < nj; j++)
          push(nrow, col_Xij(i, j, p, th), 1.0);
        for (int a : harv_arcs) {
          int age = t - arcs[a].s;          // 1-based
          if (age < 1 || age > Tm) continue;
          double coef = yield_matrix(p, age - 1) * area_ha[i];
          if (coef > 0.0)
            push(nrow, col_z(i, a), -coef);
        }
        add_con("<=", 0.0);
      }
    }
  }
  Rcpp::Rcout << "  C3: Biomass yield\n";

  // ── C4: Inventory balance ─────────────────────────────────────────────────
  // S(j,p,th) - S(j,p,th-1) - sum_i Xij(i,j,p,th) + sum_{k,pi:src=p} Xjk == 0
  // At th=0 there is no previous S (implicit zero boundary), matching R's
  // explicit first-period block for t = Amin+1.
  for (int j = 0; j < nj; j++) {
    for (int p = 0; p < P; p++) {
      for (int th = 0; th < nTh; th++) {
        push(nrow, col_S(j, p, th), +1.0);
        if (th > 0)
          push(nrow, col_S(j, p, th - 1), -1.0);
        for (int i = 0; i < ns; i++)
          push(nrow, col_Xij(i, j, p, th), -1.0);
        // outflow: Xjk where SOURCE product == p  (pp_pairs[pi].first == p)
        for (int k = 0; k < nk; k++)
          for (int pi = 0; pi < n_ppairs; pi++) {
            if (pp_pairs[pi].first != p) continue;
            push(nrow, col_Xjk(j, k, pi, th), +1.0);
          }
        add_con("==", 0.0);
      }
    }
  }
  Rcpp::Rcout << "  C4: Inventory balance\n";

  // ── C5: Storage capacity: sum_p S(j,p,th) <= CAP_stor[j]  (Tharv only)
  for (int j = 0; j < nj; j++) {
    for (int th = 0; th < nTh; th++) {
      for (int p = 0; p < P; p++)
        push(nrow, col_S(j, p, th), 1.0);
      add_con("<=", CAP_stor[j]);
    }
  }
  Rcpp::Rcout << "  C5: Storage capacity\n";

  // ── C6: Processing capacity: sum_{i,p} Xij(i,j,p,th) <= CAP_proc[j]
  for (int j = 0; j < nj; j++) {
    for (int th = 0; th < nTh; th++) {
      for (int i = 0; i < ns; i++)
        for (int p = 0; p < P; p++)
          push(nrow, col_Xij(i, j, p, th), 1.0);
      add_con("<=", CAP_proc[j]);
    }
  }
  Rcpp::Rcout << "  C6: Processing capacity\n";

  // ── C7: Demand with cascade ───────────────────────────────────────────────
  // sum_{j, pprod<=p} Xjk(j,k,pprod,p,t) <= D_max(k,p,t)
  // Filter: delivered product pp_pairs[pi].second == p  (0-based)
  //         source product pp_pairs[pi].first <= p is guaranteed by pp_pairs construction
  for (int k = 0; k < nk; k++) {
    for (int p = 0; p < P; p++) {
      for (int th = 0; th < nTh; th++) {
        int t1  = Tharv[th];
        double Dmax = dem_lut[k * P * Tm + p * Tm + (t1 - 1)];
        if (Dmax < 0.0) continue;
        bool any = false;
        for (int j = 0; j < nj; j++)
          for (int pi = 0; pi < n_ppairs; pi++) {
            if (pp_pairs[pi].second != p) continue;
            push(nrow, col_Xjk(j, k, pi, th), 1.0);
            any = true;
          }
        if (any) add_con("<=", Dmax);
      }
    }
  }
  Rcpp::Rcout << "  C7: Demand satisfaction\n";
  Rcpp::Rcout << "  Total constraints: " << nrow
              << ", nnz: " << (int)row_v.size() << "\n";

  // ==========================================================================
  // STEP 4: BOUNDS AND TYPES
  // ==========================================================================
  NumericVector ub(n_vars, R_PosInf);
  for (int c = off_z; c < off_z + n_z; c++) ub[c] = 1.0;
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
