// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
// ============================================================================
// build_and_solve_afs_lp_v12_highs.cpp  —  LP-Relaxation (MILP_v11.md)
//
// Ersetzt Gurobi durch den freien Hochleistungs-LP-Solver HiGHS
// (R-Paket: highs, https://cran.r-project.org/package=highs).
// HiGHS ist shinyapps.io-kompatibel (keine kommerzielle Lizenz nötig).
//
// Gegenüber v11 (Gurobi):
//   - Solver-Backend: gurobi::gurobi()  →  highs::highs_solve()
//   - Constraint-Matrix: arma::sp_mat (dgCMatrix) →
//     CSC-Format via R-Listen (highs::highs_solve erwartet sparse CSC)
//   - sense-Konvention: "<"/"="/">" → "<="/"=="/">=" (HiGHS-Notation)
//   - Objective-Sense: modelsense="max" → maximum=TRUE (highs)
//   - Output-Extraktion: result$primal_solution statt result$x
//   - Parameter-Interface: highs_params list (presolve, solver, time_limit …)
//
// Alle LP-Modellgleichungen (C1–C7), Variablen-Layout und
// Zielfunktion sind IDENTISCH zu v11.
//
// Variablen-Layout (0-basierte Spalten-Offsets):
//   [off_z  .. off_z   + n_z  - 1]  z[i,arc]       continuous [0, AREA_i]
//   [off_Xij.. off_Xij + n_Xij-1]  Xij[i,j,p,th]  continuous [0, +inf)
//   [off_S  .. off_S   + n_S  - 1]  S[j,p,th]      continuous [0, +inf)
//   [off_Xjk.. off_Xjk + n_Xjk-1]  Xjk[j,k,pi,th] continuous [0, +inf)
//
// Constraints:
//   C1  Sigma_{a in S^est} z[i,a] <= AREA_i                    for all i
//   C2  Sigma_{(s,t): t=tau} z[i,s,t] = Sigma_{(tau,u)} ...    for all i, tau
//   C3  Sigma_j Xij[i,j,p,t] <= Sigma_{harv-arcs} eta*z        for all i,p,t
//   C4  S[j,p,th] - S[j,p,th-1] - Sigma Xij + Sigma Xjk == 0  for all j,p,th
//   C5  Sigma_p S[j,p,th] <= CAP_stor[j]                       for all j,th
//   C6  Sigma_{i,p} Xij[i,j,p,th] <= CAP_proc[j]               for all j,th
//   C7  Sigma_{j,pi:p_low=p} Xjk[j,k,pi,th] <= D_max[k,p,t]   for all k,p,th
//
// SmartAgroforst 2026
// ============================================================================
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
List build_and_solve_afs_lp_v12_impl(List instance,
                                      List highs_params,
                                      bool verbose) {

  // ── 1. Skalare ─────────────────────────────────────────────────────────────
  const int    ns       = as<int>(instance["n_sites"]);
  const int    nj       = as<int>(instance["n_storages"]);
  const int    nk       = as<int>(instance["n_consumers"]);
  const int    Tm       = as<int>(instance["n_periods"]);
  const int    P        = as<int>(instance["n_products"]);
  const int    Amax     = as<int>(instance["max_age"]);
  const int    Amin     = as<int>(instance["min_age"]);
  const double c_tr_raw = as<double>(instance["c_tr_raw"]);
  const double c_tr_pre = as<double>(instance["c_tr_pre"]);

  if (verbose) {
    Rcpp::Rcout << "Building AFS-SCD LP (v12, HiGHS solver, shinyapps.io-compatible)...\n";
  }

  // ── 2. Site-Vektoren ───────────────────────────────────────────────────────
  DataFrame     sites_df = as<DataFrame>(instance["sites"]);
  NumericVector area_ha  = sites_df["area_ha"];
  NumericVector C_est_v  = sites_df["C_est"];
  NumericVector C_harv_v = sites_df["C_harv"];
  NumericVector C_main_v = sites_df["C_main"];
  NumericVector C_opp_v  = sites_df["C_opp"];

  // ── 3. Lager-Vektoren ──────────────────────────────────────────────────────
  DataFrame     stor_df  = as<DataFrame>(instance["storages"]);
  NumericVector CAP_stor = stor_df["CAP_stor"];
  NumericVector CAP_proc = stor_df["CAP_proc"];
  NumericVector c_stor_v = stor_df["c_stor"];

  NumericMatrix d_ij = as<NumericMatrix>(instance["d_ij"]);
  NumericMatrix d_jk = as<NumericMatrix>(instance["d_jk"]);

  // ── Ertragsmatrix: yield_mat(p, age-1) ────────────────────────────────────
  DataFrame     ydf_raw = as<DataFrame>(instance["yields_by_age"]);
  IntegerVector ya_v    = ydf_raw["age"];
  IntegerVector yp_v    = ydf_raw["product"];   // 1-basiert
  NumericVector yy_v    = ydf_raw["yield_ha"];
  NumericMatrix yield_mat(P, Tm);               // P x Tm, 0-initialisiert
  for (int r = 0; r < ya_v.size(); r++) {
    int col = ya_v[r] - 1;
    int row = yp_v[r] - 1;
    if (row >= 0 && row < P && col >= 0 && col < Tm) {
      yield_mat(row, col) = yy_v[r];
    }
  }

  // ── 4. Preis-LUT: price_lut[k*P + p] (beide 0-basiert) ───────────────────
  DataFrame     cp    = as<DataFrame>(instance["consumer_prices"]);
  IntegerVector cp_k  = cp[0];
  IntegerVector cp_p  = cp[1];
  NumericVector cp_pr = cp[2];
  NumericVector price_lut(nk * P, 0.0);
  for (int r = 0; r < cp_k.size(); r++) {
    price_lut[(cp_k[r] - 1) * P + (cp_p[r] - 1)] = cp_pr[r];
  }

  // ── 5. Nachfrage-LUT ──────────────────────────────────────────────────────
  DataFrame     dem_df   = as<DataFrame>(instance["demand"]);
  IntegerVector dem_k    = dem_df["consumer_id"];
  IntegerVector dem_p    = dem_df["product"];
  IntegerVector dem_t    = dem_df["period"];
  NumericVector dem_Dmax = dem_df["D_max"];
  NumericVector dem_lut(nk * P * Tm, -1.0);
  for (int r = 0; r < dem_k.size(); r++) {
    dem_lut[(dem_k[r]-1)*P*Tm + (dem_p[r]-1)*Tm + (dem_t[r]-1)] = dem_Dmax[r];
  }

  // ==========================================================================
  // SCHRITT 1: ARC-KONSTRUKTION  (identisch zu v11)
  // ==========================================================================
  struct Arc { int s, t; };
  std::vector<Arc> arcs;
  arcs.reserve(static_cast<size_t>((Tm + 2) * (Amax - Amin + 2)));

  // S^est: (0, t),  t = 1 .. Tm - Amin
  for (int t = 1; t <= Tm - Amin; ++t) {
    arcs.push_back({0, t});
  }
  // S^harv: (s, t),  t-s in [Amin, Amax], t <= Tm
  for (int s = 1; s <= Tm; ++s) {
    for (int len = Amin; len <= Amax; ++len) {
      int t = s + len;
      if (t <= Tm) {
        arcs.push_back({s, t});
      }
    }
  }
  // S^term: (s, Tm+1),  s = Amin+1 .. Tm
  for (int s = Amin + 1; s <= Tm; ++s) {
    arcs.push_back({s, Tm + 1});
  }

  const int n_arcs = static_cast<int>(arcs.size());

  // T_harv-Vektor
  std::vector<int> Tharv;
  for (int t = std::max(1, Amin) + 1; t <= Tm; t++) {
    Tharv.push_back(t);
  }
  const int nTh = static_cast<int>(Tharv.size());

  // pp_pairs: Q = {(p_high, p_low) | p_low >= p_high}  (0-basiert)
  std::vector<std::pair<int,int>> pp_pairs;
  for (int ph = 0; ph < P; ph++) {
    for (int pl = ph; pl < P; pl++) {
      pp_pairs.push_back({ph, pl});
    }
  }
  const int n_ppairs = static_cast<int>(pp_pairs.size());

  // Variablenanzahlen und Offsets
  const int n_z   = ns * n_arcs;
  const int n_Xij = ns * nj * P * nTh;
  const int n_S   = nj * P * nTh;
  const int n_Xjk = nj * nk * n_ppairs * nTh;

  const int off_z   = 0;
  const int off_Xij = off_z   + n_z;
  const int off_S   = off_Xij + n_Xij;
  const int off_Xjk = off_S   + n_S;
  const int n_vars  = off_Xjk + n_Xjk;

  if (verbose) {
    Rcpp::Rcout << "  n_arcs="   << n_arcs
                << "  n_ppairs=" << n_ppairs
                << "  nTh="      << nTh      << "\n";
    Rcpp::Rcout << "  n_vars="   << n_vars
                << " (z:"  << n_z
                << " Xij:" << n_Xij
                << " S:"   << n_S
                << " Xjk:" << n_Xjk << ")\n";
    Rcpp::Rcout << "  Solver: HiGHS (free, shinyapps.io-compatible)\n";
  }

  // Spalten-Index-Lambdas (0-basiert)
  auto col_z = [&](int i, int a) -> int {
    return off_z + i * n_arcs + a;
  };
  auto col_Xij = [&](int i, int j, int p, int th) -> int {
    return off_Xij + i*(nj*P*nTh) + j*(P*nTh) + p*nTh + th;
  };
  auto col_S = [&](int j, int p, int th) -> int {
    return off_S + j*(P*nTh) + p*nTh + th;
  };
  auto col_Xjk = [&](int j, int k, int pi, int th) -> int {
    return off_Xjk + j*(nk*n_ppairs*nTh) + k*(n_ppairs*nTh) + pi*nTh + th;
  };

  // ==========================================================================
  // SCHRITT 2: ZIELFUNKTIONSVEKTOR  (identisch zu v11)
  // HiGHS unterstützt maximum = TRUE direkt.
  // ==========================================================================
  std::vector<double> c_vec(n_vars, 0.0);

  // +price[k, p_low] * Xjk
  for (int j = 0; j < nj; j++) {
    for (int k = 0; k < nk; k++) {
      for (int pi = 0; pi < n_ppairs; pi++) {
        double pr = price_lut[k * P + pp_pairs[pi].second];
        for (int th = 0; th < nTh; th++) {
          c_vec[col_Xjk(j, k, pi, th)] += pr;
        }
      }
    }
  }

  // -C_est[i] * z(est-arcs)
  for (int a = 0; a < n_arcs; a++) {
    if (arcs[a].s == 0) {
      for (int i = 0; i < ns; i++) {
        c_vec[col_z(i, a)] -= C_est_v[i];
      }
    }
  }

  // -(C_main+C_opp)*(t-s)*z  und  -C_harv*z  (Harvest-Arcs)
  for (int a = 0; a < n_arcs; a++) {
    if (arcs[a].s >= 1 && arcs[a].t <= Tm) {
      int len = arcs[a].t - arcs[a].s;
      for (int i = 0; i < ns; i++) {
        c_vec[col_z(i, a)] -= (C_main_v[i] + C_opp_v[i]) * len;
        c_vec[col_z(i, a)] -= C_harv_v[i];
      }
    }
  }

  // -c_tr_raw * d_ij[i,j] * Xij
  for (int i = 0; i < ns; i++) {
    for (int j = 0; j < nj; j++) {
      double cost = c_tr_raw * d_ij(i, j);
      for (int p = 0; p < P; p++) {
        for (int th = 0; th < nTh; th++) {
          c_vec[col_Xij(i, j, p, th)] -= cost;
        }
      }
    }
  }

  // -c_tr_pre * d_jk[j,k] * Xjk
  for (int j = 0; j < nj; j++) {
    for (int k = 0; k < nk; k++) {
      double cost = c_tr_pre * d_jk(j, k);
      for (int pi = 0; pi < n_ppairs; pi++) {
        for (int th = 0; th < nTh; th++) {
          c_vec[col_Xjk(j, k, pi, th)] -= cost;
        }
      }
    }
  }

  // -c_stor[j] * S[j,p,th]
  for (int j = 0; j < nj; j++) {
    for (int p = 0; p < P; p++) {
      for (int th = 0; th < nTh; th++) {
        c_vec[col_S(j, p, th)] -= c_stor_v[j];
      }
    }
  }

  // ==========================================================================
  // SCHRITT 3: CONSTRAINTS (Triplet-Format: 0-basiert)
  // HiGHS erwartet lhs/rhs statt sense+rhs:
  //   "<=" : lhs = -Inf,      rhs = rhs_val
  //   "==" : lhs = rhs_val,   rhs = rhs_val
  //   ">=" : lhs = rhs_val,   rhs = +Inf
  // ==========================================================================
  std::vector<int>    row_v, col_v;
  std::vector<double> val_v, rhs_v;
  std::vector<std::string> sense_v;
  int nrow = 0;

  auto push = [&](int r, int c, double v) {
    row_v.push_back(r);     // 0-basiert
    col_v.push_back(c);     // 0-basiert
    val_v.push_back(v);
  };
  auto add_con = [&](const std::string& sense, double rhs) {
    sense_v.push_back(sense);
    rhs_v.push_back(rhs);
    nrow++;
  };

  // ── C1 (LP): Sigma_{a in S^est} z[i,a] <= AREA_i   for all i ─────────────
  for (int i = 0; i < ns; i++) {
    bool any = false;
    for (int a = 0; a < n_arcs; a++) {
      if (arcs[a].s == 0) {
        push(nrow, col_z(i, a), 1.0);
        any = true;
      }
    }
    if (any) {
      add_con("<=", area_ha[i]);
    }
  }
  if (verbose) { Rcpp::Rcout << "  C1 (LP, RHS=AREA_i): " << nrow << " rows\n"; }

  // ── C2: Pfadkonnektivitat   for all i, t in 1..Tm  ───────────────────────
  {
    int c2_start = nrow;
    for (int i = 0; i < ns; i++) {
      for (int t = 1; t <= Tm; t++) {
        bool any = false;
        for (int a = 0; a < n_arcs; a++) {
          if (arcs[a].t == t && arcs[a].s < t) {
            push(nrow, col_z(i, a), +1.0);
            any = true;
          } else if (arcs[a].s == t && arcs[a].t > t) {
            push(nrow, col_z(i, a), -1.0);
            any = true;
          }
        }
        if (any) {
          add_con("==", 0.0);
        }
      }
    }
    if (verbose) { Rcpp::Rcout << "  C2: " << nrow - c2_start << " rows\n"; }
  }

  // ── C3 (LP): Sigma_j Xij <= Sigma_{harv-arcs} eta*z   for all i,p,t ──────
  {
    int c3_start = nrow;
    for (int i = 0; i < ns; i++) {
      for (int p = 0; p < P; p++) {
        for (int th = 0; th < nTh; th++) {
          int t = Tharv[th];
          std::vector<int> harv_a;
          for (int a = 0; a < n_arcs; a++) {
            if (arcs[a].t == t && arcs[a].s >= 1 && arcs[a].s < t) {
              harv_a.push_back(a);
            }
          }
          if (harv_a.empty()) { continue; }
          for (int j = 0; j < nj; j++) {
            push(nrow, col_Xij(i, j, p, th), 1.0);
          }
          for (int a : harv_a) {
            int age = t - arcs[a].s;
            if (age < 1 || age > Tm) { continue; }
            double coef = yield_mat(p, age - 1);
            if (coef > 0.0) {
              push(nrow, col_z(i, a), -coef);
            }
          }
          add_con("<=", 0.0);
        }
      }
    }
    if (verbose) { Rcpp::Rcout << "  C3 (LP, no AREA_i): " << nrow - c3_start << " rows\n"; }
  }

  // ── C4: Inventarbilanz   for all j,p,th  ─────────────────────────────────
  {
    int c4_start = nrow;
    for (int j = 0; j < nj; j++) {
      for (int p = 0; p < P; p++) {
        for (int th = 0; th < nTh; th++) {
          push(nrow, col_S(j, p, th), +1.0);
          if (th > 0) {
            push(nrow, col_S(j, p, th - 1), -1.0);
          }
          for (int i = 0; i < ns; i++) {
            push(nrow, col_Xij(i, j, p, th), -1.0);
          }
          for (int k = 0; k < nk; k++) {
            for (int pi = 0; pi < n_ppairs; pi++) {
              if (pp_pairs[pi].first == p) {
                push(nrow, col_Xjk(j, k, pi, th), +1.0);
              }
            }
          }
          add_con("==", 0.0);
        }
      }
    }
    if (verbose) { Rcpp::Rcout << "  C4: " << nrow - c4_start << " rows\n"; }
  }

  // ── C5: Lagerkapazitat   for all j,th  ────────────────────────────────────
  {
    int c5_start = nrow;
    for (int j = 0; j < nj; j++) {
      for (int th = 0; th < nTh; th++) {
        for (int p = 0; p < P; p++) {
          push(nrow, col_S(j, p, th), 1.0);
        }
        add_con("<=", CAP_stor[j]);
      }
    }
    if (verbose) { Rcpp::Rcout << "  C5: " << nrow - c5_start << " rows\n"; }
  }

  // ── C6: Verarbeitungskapazitat   for all j,th  ───────────────────────────
  {
    int c6_start = nrow;
    for (int j = 0; j < nj; j++) {
      for (int th = 0; th < nTh; th++) {
        for (int i = 0; i < ns; i++) {
          for (int p = 0; p < P; p++) {
            push(nrow, col_Xij(i, j, p, th), 1.0);
          }
        }
        add_con("<=", CAP_proc[j]);
      }
    }
    if (verbose) { Rcpp::Rcout << "  C6: " << nrow - c6_start << " rows\n"; }
  }

  // ── C7: Nachfrage mit Quality Cascade   for all k,p,th  ──────────────────
  {
    int c7_start = nrow;
    for (int k = 0; k < nk; k++) {
      for (int p = 0; p < P; p++) {
        for (int th = 0; th < nTh; th++) {
          double Dmax = dem_lut[k * P * Tm + p * Tm + (Tharv[th] - 1)];
          if (Dmax < 0.0) { continue; }
          bool any = false;
          for (int j = 0; j < nj; j++) {
            for (int pi = 0; pi < n_ppairs; pi++) {
              if (pp_pairs[pi].second == p) {
                push(nrow, col_Xjk(j, k, pi, th), 1.0);
                any = true;
              }
            }
          }
          if (any) {
            add_con("<=", Dmax);
          }
        }
      }
    }
    if (verbose) { Rcpp::Rcout << "  C7: " << nrow - c7_start << " rows\n"; }
  }

  if (verbose) {
    Rcpp::Rcout << "  Total constraints: " << nrow
                << "  nnz: " << static_cast<int>(row_v.size()) << "\n";
  }

  // ==========================================================================
  // SCHRITT 4: BOUNDS (identisch zu v11)
  // ==========================================================================
  NumericVector lb(n_vars, 0.0);
  NumericVector ub(n_vars, R_PosInf);

  // z[i,a] in [0, AREA_i]
  for (int i = 0; i < ns; i++) {
    for (int a = 0; a < n_arcs; a++) {
      ub[col_z(i, a)] = area_ha[i];
    }
  }

  // ==========================================================================
  // SCHRITT 5: HIGHS-AUFRUF (ersetzt Gurobi vollständig)
  //
  // highs::highs_solve() API (R-Paket highs >= 0.1.0):
  //   highs_solve(
  //     Q      = NULL,       # keine quadratische Zielfunktion (LP)
  //     L      = c_vec,      # lineare Zielfunktionskoeffizienten
  //     lower  = lb,         # untere Variablen-Schranken
  //     upper  = ub,         # obere Variablen-Schranken
  //     A      = A_sparse,   # Constraint-Matrix (dgCMatrix, CSC)
  //     lhs    = lhs_vec,    # linke Seite: -Inf fuer <=, rhs fuer ==/>= 
  //     rhs    = rhs_vec,    # rechte Seite: rhs fuer <=, +Inf fuer >=
  //     types  = NULL,       # NULL = LP (alle Variablen kontinuierlich)
  //     maximum = TRUE,      # Maximierung
  //     control = list(…)    # Solver-Parameter
  //   )
  // ==========================================================================

  // COO → dgCMatrix via Armadillo (wird von Rcpp automatisch konvertiert)
  arma::umat locs(2, row_v.size());
  arma::vec  vals(val_v.size());
  for (size_t ii = 0; ii < row_v.size(); ii++) {
    locs(0, ii) = static_cast<arma::uword>(row_v[ii]);
    locs(1, ii) = static_cast<arma::uword>(col_v[ii]);
    vals[ii]    = val_v[ii];
  }
  arma::sp_mat A_arma(locs, vals, nrow, n_vars);
  RObject A_sparse = wrap(A_arma);  // sp_mat → dgCMatrix (CSC)

  // lhs/rhs-Vektoren aus sense+rhs aufbauen
  NumericVector lhs_vec(nrow);
  NumericVector rhs_vec(nrow);
  for (int r = 0; r < nrow; r++) {
    double rv = rhs_v[r];
    if        (sense_v[r] == "<=") {
      lhs_vec[r] = R_NegInf;
      rhs_vec[r] = rv;
    } else if (sense_v[r] == "==") {
      lhs_vec[r] = rv;
      rhs_vec[r] = rv;
    } else {                         // ">="
      lhs_vec[r] = rv;
      rhs_vec[r] = R_PosInf;
    }
  }

  // HiGHS-Kontrollparameter (Defaults + User-Overrides)
  List h_ctrl = List::create(
    Named("presolve")    = "on",
    Named("solver")      = "ipm",        // Interior-Point: schnellster fuer grosse LP
    Named("time_limit")  = 300.0,
    Named("output_flag") = (verbose ? true : false)
  );
  if (highs_params.size() > 0) {
    CharacterVector pnames = highs_params.names();
    for (int pi = 0; pi < pnames.size(); pi++) {
      h_ctrl[as<std::string>(pnames[pi])] = highs_params[pi];
    }
  }

  if (verbose) {
    Rcpp::Rcout << "  Calling highs::highs_solve() [LP mode]"
                << "  vars=" << n_vars
                << "  constrs=" << nrow << "\n";
  }

  // highs::highs_solve aufrufen
  Environment highs_env = Environment::namespace_env("highs");
  Function    highs_fn  = highs_env["highs_solve"];

  List highs_result = highs_fn(
    Named("Q")       = R_NilValue,
    Named("L")       = wrap(c_vec),
    Named("lower")   = lb,
    Named("upper")   = ub,
    Named("A")       = A_sparse,
    Named("lhs")     = lhs_vec,
    Named("rhs")     = rhs_vec,
    Named("types")   = R_NilValue,   // NULL = LP
    Named("maximum") = true,          // Maximierung
    Named("control") = h_ctrl
  );

  // ==========================================================================
  // SCHRITT 6: ERGEBNIS EXTRAHIEREN (HiGHS-Notation)
  //
  // highs_solve() gibt zurück:
  //   $primal_solution  — Lösungsvektor (= x in Gurobi)
  //   $objective_value  — Zielfunktionswert
  //   $status           — "Optimal", "Infeasible", "Time limit reached", ...
  //   $info             — Solver-Statistiken
  // ==========================================================================
  std::string status = "Unknown";
  if (highs_result.containsElementNamed("status")) {
    status = as<std::string>(highs_result["status"]);
  }

  double objval = 0.0;
  if (highs_result.containsElementNamed("objective_value")) {
    objval = as<double>(highs_result["objective_value"]);
  }

  if (verbose) {
    Rcpp::Rcout << "  Status: " << status
                << "  ObjVal: " << objval << "\n";
  }

  // Lösungsvektor extrahieren ($primal_solution statt $x)
  NumericVector x_sol(n_vars, 0.0);
  if (highs_result.containsElementNamed("primal_solution")) {
    SEXP ps = highs_result["primal_solution"];
    if (!Rf_isNull(ps)) {
      x_sol = as<NumericVector>(ps);
    }
  }

  // ── z-Lösung ──────────────────────────────────────────────────────────────
  std::vector<int>         z_i_v, z_a_v, z_s_v, z_t_v;
  std::vector<double>      z_val_v;
  std::vector<std::string> z_type_v;
  for (int i = 0; i < ns; i++) {
    for (int a = 0; a < n_arcs; a++) {
      double zv = x_sol[col_z(i, a)];
      if (zv > 1e-6) {
        z_i_v.push_back(i + 1);
        z_a_v.push_back(a + 1);
        z_s_v.push_back(arcs[a].s);
        z_t_v.push_back(arcs[a].t);
        z_val_v.push_back(zv);
        if      (arcs[a].s == 0)      { z_type_v.push_back("establishment"); }
        else if (arcs[a].t == Tm + 1) { z_type_v.push_back("termination");   }
        else                           { z_type_v.push_back("harvest");       }
      }
    }
  }

  // ── Xij-Lösung ────────────────────────────────────────────────────────────
  std::vector<int>    xij_i_v, xij_j_v, xij_p_v, xij_t_v;
  std::vector<double> xij_val_v;
  for (int i = 0; i < ns; i++) {
    for (int j = 0; j < nj; j++) {
      for (int p = 0; p < P; p++) {
        for (int th = 0; th < nTh; th++) {
          double v = x_sol[col_Xij(i, j, p, th)];
          if (v > 1e-6) {
            xij_i_v.push_back(i + 1);
            xij_j_v.push_back(j + 1);
            xij_p_v.push_back(p + 1);
            xij_t_v.push_back(Tharv[th]);
            xij_val_v.push_back(v);
          }
        }
      }
    }
  }

  // ── Xjk-Lösung ────────────────────────────────────────────────────────────
  std::vector<int>    xjk_j_v, xjk_k_v, xjk_ph_v, xjk_pl_v, xjk_t_v;
  std::vector<double> xjk_val_v;
  for (int j = 0; j < nj; j++) {
    for (int k = 0; k < nk; k++) {
      for (int pi = 0; pi < n_ppairs; pi++) {
        for (int th = 0; th < nTh; th++) {
          double v = x_sol[col_Xjk(j, k, pi, th)];
          if (v > 1e-6) {
            xjk_j_v.push_back(j + 1);
            xjk_k_v.push_back(k + 1);
            xjk_ph_v.push_back(pp_pairs[pi].first  + 1);
            xjk_pl_v.push_back(pp_pairs[pi].second + 1);
            xjk_t_v.push_back(Tharv[th]);
            xjk_val_v.push_back(v);
          }
        }
      }
    }
  }

  // ── S-Lösung ──────────────────────────────────────────────────────────────
  std::vector<int>    s_j_v, s_p_v, s_t_v;
  std::vector<double> s_val_v;
  for (int j = 0; j < nj; j++) {
    for (int p = 0; p < P; p++) {
      for (int th = 0; th < nTh; th++) {
        double v = x_sol[col_S(j, p, th)];
        if (v > 1e-6) {
          s_j_v.push_back(j + 1);
          s_p_v.push_back(p + 1);
          s_t_v.push_back(Tharv[th]);
          s_val_v.push_back(v);
        }
      }
    }
  }

  // DataFrames zusammenstellen
  DataFrame z_df = DataFrame::create(
    Named("i")    = wrap(z_i_v),
    Named("a")    = wrap(z_a_v),
    Named("s")    = wrap(z_s_v),
    Named("t")    = wrap(z_t_v),
    Named("type") = wrap(z_type_v),
    Named("value")= wrap(z_val_v)
  );
  DataFrame xij_df = DataFrame::create(
    Named("i")    = wrap(xij_i_v),
    Named("j")    = wrap(xij_j_v),
    Named("p")    = wrap(xij_p_v),
    Named("t")    = wrap(xij_t_v),
    Named("value")= wrap(xij_val_v)
  );
  DataFrame xjk_df = DataFrame::create(
    Named("j")    = wrap(xjk_j_v),
    Named("k")    = wrap(xjk_k_v),
    Named("p")    = wrap(xjk_ph_v),
    Named("q")    = wrap(xjk_pl_v),
    Named("t")    = wrap(xjk_t_v),
    Named("value")= wrap(xjk_val_v)
  );
  DataFrame s_df = DataFrame::create(
    Named("j")    = wrap(s_j_v),
    Named("p")    = wrap(s_p_v),
    Named("t")    = wrap(s_t_v),
    Named("value")= wrap(s_val_v)
  );

  return List::create(
    Named("status")          = status,
    Named("objective_value") = objval,
    Named("model_type")      = "LP",
    Named("solver")          = "HiGHS",
    Named("highs_result")    = highs_result,
    Named("solution_vector") = List::create(
      Named("z")   = z_df,
      Named("Xij") = xij_df,
      Named("Xjk") = xjk_df,
      Named("S")   = s_df
    )
  );
}

// [[Rcpp::export]]
List build_and_solve_afs_lp_v12(List instance,
                                  List highs_params,
                                  bool verbose = true) {
  return build_and_solve_afs_lp_v12_impl(instance, highs_params, verbose);
}

// [[Rcpp::export]]
List build_and_solve_afs_lp_v12_default(List instance,
                                         bool verbose = true) {
  List empty_params = List::create();
  return build_and_solve_afs_lp_v12_impl(instance, empty_params, verbose);
}
