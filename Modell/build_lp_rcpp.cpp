// [[file: build_lp_rcpp.cpp]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
#include <vector>
#include <string>
using namespace Rcpp;

// ============================================================================
// HELPER: flat index lookup into pre-built hash maps
//   Key encoding: (a,b,c,d) → a + A*(b + B*(c + C*d))
// ============================================================================
inline int idx3(int a, int b, int c, int B, int C) {
  return a + B * (b + C * c);
}
inline int idx4(int a, int b, int c, int d, int B, int C, int D) {
  return a + B * (b + C * (c + D * d));
}

// ============================================================================
// MAIN FUNCTION — mirrors build_agroforestry_lp_sparse_v10 in pure C++
//
// @param instance  Named R list (same structure as in the R version)
// @return          Named list: model components for ROI::OP()
// ============================================================================
// [[Rcpp::export]]
List build_lp_rcpp(List instance) {
  
  // ── Extract scalars ──────────────────────────────────────────────────────
  int ns   = as<int>(instance["n_sites"]);
  int nj   = as<int>(instance["n_storages"]);
  int nk   = as<int>(instance["n_consumers"]);
  int Tm   = as<int>(instance["n_periods"]);
  int P    = as<int>(instance["n_products"]);
  int Amax = as<int>(instance["max_age"]);
  int Amin = as<int>(instance["min_age"]);
  double Copp     = as<double>(instance["c_opp"]);
  double c_tr_raw = as<double>(instance["c_tr_raw"]);
  double c_tr_pre = as<double>(instance["c_tr_pre"]);
  
  Rcpp::Rcout << "Initializing scalars\n";
  
  // ── Extract vectors / matrices ───────────────────────────────────────────
  DataFrame sites_df   = as<DataFrame>(instance["sites"]);
  NumericVector area_ha = sites_df["area_ha"];
  NumericVector C_est   = sites_df["C_est"];
  NumericVector C_harv  = sites_df["C_harv"];
  
  DataFrame storages_df   = as<DataFrame>(instance["storages"]);
  NumericVector CAP_stor  = storages_df["CAP_stor"];
  NumericVector CAP_proc  = storages_df["CAP_proc"];
  NumericVector c_stor    = storages_df["c_stor"];
  
  NumericMatrix d_ij = as<NumericMatrix>(instance["d_ij"]);  // [ns x nj]
  NumericMatrix d_jk = as<NumericMatrix>(instance["d_jk"]);  // [nj x nk]
  
  // ── Yield matrix [P x Tm] ────────────────────────────────────────────────
  // Passed as a pre-computed matrix from R (avoids slow R loop in C++)
  NumericMatrix yield_matrix = as<NumericMatrix>(instance["yield_matrix"]);
  NumericVector yield_max(P);
  for (int p = 0; p < P; p++) {
    double mx = 0.0;
    for (int t = 0; t < Tm; t++) mx = std::max(mx, yield_matrix(p, t));
    yield_max[p] = mx;
  }
  
  // ── Consumer prices DataFrame [k, pp, price] ─────────────────────────────
  DataFrame cp        = as<DataFrame>(instance["consumer_prices"]);
  IntegerVector cp_k  = cp[0];
  IntegerVector cp_pp = cp[1];
  NumericVector cp_pr = cp[2];
  // Build price lookup: (k-1)*P + (pp-1) → price
  NumericVector price_lut(nk * P, 0.0);
  for (int r = 0; r < cp_k.size(); r++) {
    price_lut[(cp_k[r]-1) * P + (cp_pp[r]-1)] = cp_pr[r];
  }
  
  // ── Demand DataFrame [consumer_id, product, period, D_max] ───────────────
  DataFrame dem_df      = as<DataFrame>(instance["demand"]);
  IntegerVector dem_k   = dem_df["consumer_id"];
  IntegerVector dem_p   = dem_df["product"];
  IntegerVector dem_t   = dem_df["period"];
  NumericVector dem_Dmax= dem_df["D_max"];
  // Build demand lookup: (k-1)*P*Tm + (p-1)*Tm + (t-1) → D_max
  NumericVector dem_lut(nk * P * Tm, -1.0);
  for (int r = 0; r < dem_k.size(); r++) {
    dem_lut[(dem_k[r]-1)*P*Tm + (dem_p[r]-1)*Tm + (dem_t[r]-1)] = dem_Dmax[r];
  }
  
  //Rcpp::Rcout << "Building sparse LP (Rcpp)...\n";
  // Rcpp::Rcout << "  ns=" << ns << "  nj=" << nj << "  nk=" << nk
  //             << "  Tm=" << Tm << "  P=" << P << "\n";
  // 
  // ==========================================================================
  // STEP 1: VARIABLE INDEXING
  //   Layout (0-based column offsets):
  //     z_tuples : ns * n_arcs_per_site
  //     Y_tuples : ns * P * Tharv
  //     Xij      : ns * nj * P * Tharv
  //     S        : nj * P * Tm
  //     Xjk      : nj * nk * P * P * Tharv  (filtered: pp >= p)
  // ==========================================================================
  
  // ── Build valid arc list (s < t, arc length in [Amin,Amax]) ─────────────
  // T_ext = 0..Tm+1
  struct Arc { int s, t; };
  std::vector<Arc> arcs;
  arcs.reserve(Tm * Tm / 2);
  for (int s = 0; s <= Tm + 1; s++) {
    for (int t = s + 1; t <= Tm + 1; t++) {
      int len = t - s;
      // arc (0,t) with t in 1..Tm: establishment arcs
      if (s == 0 && t >= 1 && t <= Tm) { arcs.push_back({s, t}); continue; }
      // arc (s,Tm+1) with s >= 1: termination arcs
      if (t == Tm + 1 && s >= 1)      { arcs.push_back({s, t}); continue; }
      // harvest arcs: s >= 1, t in Tset, len in [Amin,Amax]
      if (s >= 1 && t >= 1 && t <= Tm && len >= Amin && len <= Amax) {
        arcs.push_back({s, t});
      }
    }
  }
  int n_arcs = arcs.size();
  
  // Harvest periods: t in 1..Tm where t > max(1, Amin)
  std::vector<int> Tharv;
  int t_min_harv = std::max(1, Amin) + 1;
  for (int t = t_min_harv; t <= Tm; t++) Tharv.push_back(t);
  int nTh = Tharv.size();
  
  // ── Contiguous column offsets ────────────────────────────────────────────
  int n_z   = ns * n_arcs;
  int n_Y   = ns * P * nTh;
  int n_Xij = ns * nj * P * nTh;
  int n_S   = nj * P * Tm;
  
  // Xjk: only pp >= p  →  count valid (p,pp) pairs
  std::vector<std::pair<int,int>> pp_pairs;  // (p, pp) 0-based
  for (int p = 0; p < P; p++)
    for (int pp = p; pp < P; pp++)
      pp_pairs.push_back({p, pp});
  int n_ppairs = pp_pairs.size();
  int n_Xjk = nj * nk * n_ppairs * nTh;
  
  int off_z   = 0;
  int off_Y   = off_z   + n_z;
  int off_Xij = off_Y   + n_Y;
  int off_S   = off_Xij + n_Xij;
  int off_Xjk = off_S   + n_S;
  int n_vars  = off_Xjk + n_Xjk;
  
  // Rcpp::Rcout << "  n_vars=" << n_vars
  //             << " (z:" << n_z << " Y:" << n_Y
  //             << " Xij:" << n_Xij << " S:" << n_S
  //             << " Xjk:" << n_Xjk << ")\n";
  // 
  // ── Inline column index helpers ──────────────────────────────────────────
  // All indices 0-based internally; +1 when stored for R's 1-based triplet
  auto col_z = [&](int i, int arc) { return off_z + i * n_arcs + arc; };
  auto col_Y = [&](int i, int p, int th) {
    return off_Y + i * (P * nTh) + p * nTh + th;
  };
  auto col_Xij = [&](int i, int j, int p, int th) {
    return off_Xij + i * (nj * P * nTh) + j * (P * nTh) + p * nTh + th;
  };
  auto col_S = [&](int j, int p, int t) {
    return off_S + j * (P * Tm) + p * Tm + t;
  };
  auto col_Xjk = [&](int j, int k, int pairIdx, int th) {
    return off_Xjk + j * (nk * n_ppairs * nTh)
    + k * (n_ppairs * nTh)
    + pairIdx * nTh + th;
  };
  
  // ==========================================================================
  // STEP 2: OBJECTIVE VECTOR  (pre-allocated)
  // ==========================================================================
  std::vector<double> c_vec(n_vars, 0.0);
  
  // Revenue: price[k,pp] for each Xjk
  for (int j = 0; j < nj; j++)
    for (int k = 0; k < nk; k++)
      for (int pi = 0; pi < n_ppairs; pi++) {
        int pp = pp_pairs[pi].second;
        double pr = price_lut[k * P + pp];
        for (int th = 0; th < nTh; th++)
          c_vec[col_Xjk(j, k, pi, th)] += pr;
      }
      
      // Establishment + opportunity cost on z(i, 0, t)
      for (int a = 0; a < n_arcs; a++) {
        if (arcs[a].s != 0) continue;
        int t = arcs[a].t;
        if (t < 1 || t > Tm) continue;
        for (int i = 0; i < ns; i++) {
          int col = col_z(i, a);
          c_vec[col] -= C_est[i];
          c_vec[col] -= Copp * area_ha[i] * (Tm - t);
        }
      }
      
      // Harvest cost on z(i, s>=1, t)
      for (int a = 0; a < n_arcs; a++) {
        if (arcs[a].s < 1 || arcs[a].t > Tm) continue;
        for (int i = 0; i < ns; i++)
          c_vec[col_z(i, a)] -= C_harv[i] * area_ha[i];
      }
      
      // Transport site→hub: -c_tr_raw * d_ij[i,j] for each Xij
      for (int i = 0; i < ns; i++)
        for (int j = 0; j < nj; j++) {
          double cost = c_tr_raw * d_ij(i, j);
          for (int p = 0; p < P; p++)
            for (int th = 0; th < nTh; th++)
              c_vec[col_Xij(i, j, p, th)] -= cost;
        }
        
        // Transport hub→consumer: -c_tr_pre * d_jk[j,k]
        for (int j = 0; j < nj; j++)
          for (int k = 0; k < nk; k++) {
            double cost = c_tr_pre * d_jk(j, k);
            for (int pi = 0; pi < n_ppairs; pi++)
              for (int th = 0; th < nTh; th++)
                c_vec[col_Xjk(j, k, pi, th)] -= cost;
          }
          
          // Storage cost
          for (int j = 0; j < nj; j++)
            for (int p = 0; p < P; p++)
              for (int t = 0; t < Tm; t++)
                c_vec[col_S(j, p, t)] -= c_stor[j];
  
  // ==========================================================================
  // STEP 3: CONSTRAINT ASSEMBLY  (pre-allocated STL vectors)
  // ==========================================================================
  // Generous pre-allocation: avoids any reallocation during fill
  int est_nnz = n_z * 4 + n_Y * 3 + n_Xij * 2 + n_S * 4 + n_Xjk * 2;
  std::vector<int>    row_v; row_v.reserve(est_nnz);
  std::vector<int>    col_v; col_v.reserve(est_nnz);
  std::vector<double> val_v; val_v.reserve(est_nnz);
  std::vector<double> rhs_v;
  std::vector<std::string> sense_v;
  int row = 0;  // 0-based; +1 when pushed to R
  
  auto push = [&](int r, int c, double v) {
    row_v.push_back(r + 1);   // 1-based for R
    col_v.push_back(c + 1);
    val_v.push_back(v);
  };
  
  // ── C1: Path establishment: sum_t z(i,0,t) <= 1 ─────────────────────────
  for (int i = 0; i < ns; i++) {
    bool any = false;
    for (int a = 0; a < n_arcs; a++) {
      if (arcs[a].s != 0 || arcs[a].t < 1 || arcs[a].t > Tm) continue;
      push(row, col_z(i, a), 1.0);
      any = true;
    }
    if (any) {
      rhs_v.push_back(1.0); sense_v.push_back("<="); row++;
    }
  }
  // Rcpp::Rcout << "  C1 done\n";
  
  // ── C2: Path connectivity: sum_{s<t} z(i,s,t) = sum_{u>t} z(i,t,u) ──────
  for (int i = 0; i < ns; i++) {
    for (int t = 1; t <= Tm; t++) {
      bool any = false;
      for (int a = 0; a < n_arcs; a++) {
        if (arcs[a].t == t && arcs[a].s < t) {
          push(row, col_z(i, a), +1.0); any = true;
        }
        if (arcs[a].s == t && arcs[a].t > t) {
          push(row, col_z(i, a), -1.0); any = true;
        }
      }
      if (any) {
        rhs_v.push_back(0.0); sense_v.push_back("=="); row++;
      }
    }
  }
  // Rcpp::Rcout << "  C2 done\n";
  
  // ── C3a: Biomass yield: Y(i,p,t) <= sum_s yield*area*z(i,s,t) ───────────
  for (int i = 0; i < ns; i++)
    for (int p = 0; p < P; p++)
      for (int th = 0; th < nTh; th++) {
        int t = Tharv[th];
        int Ycol = col_Y(i, p, th);
        bool any = false;
        push(row, Ycol, 1.0);
        for (int a = 0; a < n_arcs; a++) {
          if (arcs[a].t != t || arcs[a].s < 1 || arcs[a].s >= t) continue;
          int age = t - arcs[a].s - 1;  // 0-based index into yield_matrix
          if (age < 0 || age >= Tm) continue;
          double coef = yield_matrix(p, age) * area_ha[i];
          if (coef > 0.0) { push(row, col_z(i, a), -coef); any = true; }
        }
        if (any) {
          rhs_v.push_back(0.0); sense_v.push_back("<="); row++;
        } else {
          // Remove the Y entry we pushed if no arcs found
          row_v.pop_back(); col_v.pop_back(); val_v.pop_back();
        }
      }
      // Rcpp::Rcout << "  C3a done\n";
  
  // ── C4: Flow balance: sum_j Xij(i,j,p,t) <= Y(i,p,t) ───────────────────
  for (int i = 0; i < ns; i++)
    for (int p = 0; p < P; p++)
      for (int th = 0; th < nTh; th++) {
        push(row, col_Y(i, p, th), -1.0);
        for (int j = 0; j < nj; j++)
          push(row, col_Xij(i, j, p, th), 1.0);
        rhs_v.push_back(0.0); sense_v.push_back("<="); row++;
      }
      // Rcpp::Rcout << "  C4 done\n";
  
  // ── C6: Inventory balance: S(j,p,t) = S(j,p,t-1) + Xij_in - Xjk_out ────
  for (int j = 0; j < nj; j++)
    for (int p = 0; p < P; p++)
      for (int t = 0; t < Tm; t++) {   // t is 0-based (period t+1)
        int t1 = t + 1;                 // 1-based period
        // find th index for t1
        auto it = std::find(Tharv.begin(), Tharv.end(), t1);
        int th = (it != Tharv.end()) ? (int)(it - Tharv.begin()) : -1;
        
        // S(j,p,t)
        push(row, col_S(j, p, t), 1.0);
        // -S(j,p,t-1) if t > 0
        if (t > 0) push(row, col_S(j, p, t - 1), -1.0);
        // -sum_i Xij(i,j,p,th)
        if (th >= 0)
          for (int i = 0; i < ns; i++)
            push(row, col_Xij(i, j, p, th), -1.0);
        // +sum_k sum_{pp>=p} Xjk(j,k,p,pp,th)
        if (th >= 0)
          for (int k = 0; k < nk; k++)
            for (int pi = 0; pi < n_ppairs; pi++) {
              if (pp_pairs[pi].first != p) continue;  // filter: p == first of pair
              push(row, col_Xjk(j, k, pi, th), +1.0);
            }
            rhs_v.push_back(0.0); sense_v.push_back("=="); row++;
      }
      // Rcpp::Rcout << "  C6 done\n";
  
  // ── C7: Storage capacity: sum_p S(j,p,t) <= CAP_stor[j] ─────────────────
  for (int j = 0; j < nj; j++)
    for (int t = 0; t < Tm; t++) {
      for (int p = 0; p < P; p++) push(row, col_S(j, p, t), 1.0);
      rhs_v.push_back(CAP_stor[j]); sense_v.push_back("<="); row++;
    }
    // Rcpp::Rcout << "  C7 done\n";
  
  // ── C8: Processing capacity: sum_{i,p} Xij(i,j,p,t) <= CAP_proc[j] ──────
  for (int j = 0; j < nj; j++)
    for (int th = 0; th < nTh; th++) {
      for (int i = 0; i < ns; i++)
        for (int p = 0; p < P; p++)
          push(row, col_Xij(i, j, p, th), 1.0);
      rhs_v.push_back(CAP_proc[j]); sense_v.push_back("<="); row++;
    }
    // Rcpp::Rcout << "  C8 done\n";
  
  // ── C9: Demand with cascade: sum_j Xjk(j,k,q,p,t) <= D_max(k,p,t) ───────
  for (int k = 0; k < nk; k++)
    for (int p = 0; p < P; p++)
      for (int th = 0; th < nTh; th++) {
        int t1 = Tharv[th];
        double Dmax = dem_lut[k * P * Tm + p * Tm + (t1 - 1)];
        if (Dmax < 0.0) continue;
        for (int j = 0; j < nj; j++)
          for (int pi = 0; pi < n_ppairs; pi++) {
            // cascade: q (=pp_pairs[pi].first) <= p, delivered as product p
            if (pp_pairs[pi].second != p) continue;
            push(row, col_Xjk(j, k, pi, th), 1.0);
          }
          rhs_v.push_back(Dmax); sense_v.push_back("<="); row++;
      }
  //     Rcpp::Rcout << "  C9 done\n";
  // 
  // Rcpp::Rcout << "  Total constraints: " << row
  //             << ", nnz: " << row_v.size() << "\n";
  // 
  // ==========================================================================
  // STEP 4: RETURN LIST FOR ROI::OP() in R
  // ==========================================================================
  
  // Upper bounds vector
  NumericVector ub(n_vars, R_PosInf);
  // z: binary → ub = 1
  for (int c = off_z; c < off_z + n_z; c++) ub[c] = 1.0;
  // Y: area * yield_max per product
  for (int i = 0; i < ns; i++)
    for (int p = 0; p < P; p++)
      for (int th = 0; th < nTh; th++)
        ub[col_Y(i, p, th)] = area_ha[i] * yield_max[p];
  // S: CAP_stor
  for (int j = 0; j < nj; j++)
    for (int p = 0; p < P; p++)
      for (int t = 0; t < Tm; t++)
        ub[col_S(j, p, t)] = CAP_stor[j];
  
  // Variable types: "B" for z, "C" for rest
  CharacterVector types(n_vars, "C");
  for (int c = off_z; c < off_z + n_z; c++) types[c] = "B";
  
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
    Named("n_constrs") = row,
    Named("n_z")       = n_z,
    Named("n_Y")       = n_Y,
    Named("n_Xij")     = n_Xij,
    Named("n_S")       = n_S,
    Named("n_Xjk")     = n_Xjk
  );
}