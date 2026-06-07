// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
List build_and_solve_afs_milp(List instance,
                              List gurobi_params = List::create(),
                              bool verbose = true) {
  
  // ── 1. Skalare Parameter ──────────────────────────────────────────────────
  int ns   = as<int>(instance["n_sites"]);
  int nj   = as<int>(instance["n_storages"]);
  int nk   = as<int>(instance["n_consumers"]);
  int Tm   = as<int>(instance["n_periods"]);
  int P    = as<int>(instance["n_products"]);
  int Amax = as<int>(instance["max_age"]);
  int Amin = as<int>(instance["min_age"]);
  double c_tr_raw = as<double>(instance["c_tr_raw"]);
  double c_tr_pre = as<double>(instance["c_tr_pre"]);
  
  if (verbose)
    Rcpp::Rcout << "Building sparse LP (Rcpp v3 — OMPR-aligned, Gurobi direct)...\n";
  
  // ── 2. Site-Vektoren ──────────────────────────────────────────────────────
  DataFrame sites_df = as<DataFrame>(instance["sites"]);
  NumericVector area_ha = sites_df["area_ha"];
  NumericVector C_est   = sites_df["C_est"];
  NumericVector C_harv  = sites_df["C_harv"];
  NumericVector C_main  = sites_df["C_main"];
  NumericVector C_opp   = sites_df["C_opp"];
  
  // ── 3. Lager-Vektoren ─────────────────────────────────────────────────────
  DataFrame storages_df = as<DataFrame>(instance["storages"]);
  NumericVector CAP_stor = storages_df["CAP_stor"];
  NumericVector CAP_proc = storages_df["CAP_proc"];
  NumericVector c_stor   = storages_df["c_stor"];
  
  NumericMatrix d_ij = as<NumericMatrix>(instance["d_ij"]); // [ns x nj]
  NumericMatrix d_jk = as<NumericMatrix>(instance["d_jk"]); // [nj x nk]
  
  // ── 4. Ertragsmatrix [P x Tm] (age=1 → Spalte 0) ─────────────────────────
  NumericMatrix yield_matrix = as<NumericMatrix>(instance["yield_matrix"]);
  
  // ── 5. Preis-Lookup (0-basiert): price_lut[k*P + p] ─────────────────────
  DataFrame cp     = as<DataFrame>(instance["consumer_prices"]);
  IntegerVector cp_k  = cp[0]; // consumer_id (1-basiert)
  IntegerVector cp_pp = cp[1]; // product     (1-basiert)
  NumericVector cp_pr = cp[2]; // price
  NumericVector price_lut(nk * P, 0.0);
  for (int r = 0; r < cp_k.size(); r++)
    price_lut[(cp_k[r]-1) * P + (cp_pp[r]-1)] = cp_pr[r];
  
  // ── 6. Nachfrage-Lookup: dem_lut[(k-1)*P*Tm + (p-1)*Tm + (t-1)] ─────────
  DataFrame dem_df       = as<DataFrame>(instance["demand"]);
  IntegerVector dem_k    = dem_df["consumer_id"];
  IntegerVector dem_p    = dem_df["product"];
  IntegerVector dem_t    = dem_df["period"];
  NumericVector dem_Dmax = dem_df["D_max"];
  NumericVector dem_lut(nk * P * Tm, -1.0); // -1 = kein Eintrag
  for (int r = 0; r < dem_k.size(); r++)
    dem_lut[(dem_k[r]-1)*P*Tm + (dem_p[r]-1)*Tm + (dem_t[r]-1)] = dem_Dmax[r];
  
  // ==========================================================================
  // SCHRITT 1: VARIABLEN-INDIZIERUNG
  // ==========================================================================
  
  // ── Arc-Konstruktion (drei Kategorien wie in OMPR) ────────────────────────
  struct Arc { int s, t; };
  std::vector<Arc> arcs;
  arcs.reserve((Tm + 2) * (Amax - Amin + 2) * 3);
  
  // S^est: (0,t), t = 1 .. Tm - Amin
  for (int t = 1; t <= Tm - Amin; ++t)
    arcs.push_back({0, t});
  
  // S^harv: (s,t), s,t in 1..Tm, t-s in [Amin..Amax]
  for (int s = 1; s <= Tm; ++s)
    for (int len = Amin; len <= Amax; ++len) {
      int t = s + len;
      if (t >= 1 && t <= Tm)
        arcs.push_back({s, t});
    }
    
    // S^term: (s, Tm+1), s = Amin+1 .. Tm
    for (int s = Amin + 1; s <= Tm; ++s)
      arcs.push_back({s, Tm + 1});
  
  int n_arcs = (int)arcs.size();
  
  // ── T_harv: {max(1,Amin)+1 .. Tm} ────────────────────────────────────────
  std::vector<int> Tharv;
  int t_min_harv = std::max(1, Amin) + 1;
  for (int t = t_min_harv; t <= Tm; t++) Tharv.push_back(t);
  int nTh = (int)Tharv.size();
  
  // ── pp_pairs: (p_high, p_low) 0-basiert, p_low >= p_high ─────────────────
  // Entspricht OMPR Q-Tabelle: p_high kann p_low decken
  std::vector<std::pair<int,int>> pp_pairs;
  for (int ph = 0; ph < P; ph++)
    for (int pl = ph; pl < P; pl++)
      pp_pairs.push_back({ph, pl});
  int n_ppairs = (int)pp_pairs.size();
  
  // ── Variablenanzahlen und Offsets ─────────────────────────────────────────
  int n_z   = ns * n_arcs;
  int n_Xij = ns * nj * P * nTh;
  int n_S   = nj * P * nTh;
  int n_Xjk = nj * nk * n_ppairs * nTh;
  
  int off_z   = 0;
  int off_Xij = off_z   + n_z;
  int off_S   = off_Xij + n_Xij;
  int off_Xjk = off_S   + n_S;
  int n_vars  = off_Xjk + n_Xjk;
  
  if (verbose)
    Rcpp::Rcout << "  n_vars=" << n_vars
                << " (z:"   << n_z
                << " Xij:" << n_Xij
                << " S:"   << n_S
                << " Xjk:" << n_Xjk << ")\n";
    
    // ── Spalten-Index-Lambdas (0-basiert) ─────────────────────────────────────
    auto col_z = [&](int i, int arc) -> int {
      return off_z + i * n_arcs + arc;
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
    // SCHRITT 2: ZIELFUNKTIONSVEKTOR
    // ==========================================================================
    std::vector<double> c_vec(n_vars, 0.0);
    
    // Revenue: +price[k, p_low] * Xjk[j,k,pi,th]
    for (int j = 0; j < nj; j++)
      for (int k = 0; k < nk; k++)
        for (int pi = 0; pi < n_ppairs; pi++) {
          int pl = pp_pairs[pi].second; // p_low (0-basiert) = geliefertes Produkt
          double pr = price_lut[k * P + pl];
          for (int th = 0; th < nTh; th++)
            c_vec[col_Xjk(j, k, pi, th)] += pr;
        }
        
        // Etablierungskosten: -C_est[i] * area_ha[i] * z(est-arcs)
        for (int a = 0; a < n_arcs; a++) {
          if (arcs[a].s != 0) continue;
          for (int i = 0; i < ns; i++)
            c_vec[col_z(i, a)] -= C_est[i] * area_ha[i];
        }
        
        // Maintenance + Opp auf Ernte-Arcs: -(C_main[i]+C_opp[i])*area_ha[i]*(t-s)*z(harv)
        for (int a = 0; a < n_arcs; a++) {
          if (arcs[a].s < 1 || arcs[a].t > Tm) continue; // nur Harvest-Arcs
          int arc_len = arcs[a].t - arcs[a].s;
          for (int i = 0; i < ns; i++)
            c_vec[col_z(i, a)] -= (C_main[i] + C_opp[i]) * area_ha[i] * arc_len;
        }
        
        // Erntekosten: -C_harv[i] * area_ha[i] * z(harv-arcs)
        for (int a = 0; a < n_arcs; a++) {
          if (arcs[a].s < 1 || arcs[a].t > Tm) continue;
          for (int i = 0; i < ns; i++)
            c_vec[col_z(i, a)] -= C_harv[i] * area_ha[i];
        }
        
        // Rohtransportkosten: -c_tr_raw * d_ij[i,j] * Xij
        for (int i = 0; i < ns; i++)
          for (int j = 0; j < nj; j++) {
            double cost = c_tr_raw * d_ij(i, j);
            for (int p = 0; p < P; p++)
              for (int th = 0; th < nTh; th++)
                c_vec[col_Xij(i, j, p, th)] -= cost;
          }
          
          // Verarbeitete Transportkosten: -c_tr_pre * d_jk[j,k] * Xjk
          for (int j = 0; j < nj; j++)
            for (int k = 0; k < nk; k++) {
              double cost = c_tr_pre * d_jk(j, k);
              for (int pi = 0; pi < n_ppairs; pi++)
                for (int th = 0; th < nTh; th++)
                  c_vec[col_Xjk(j, k, pi, th)] -= cost;
            }
            
            // Lagerkosten: -c_stor[j] * S[j,p,th]
            for (int j = 0; j < nj; j++)
              for (int p = 0; p < P; p++)
                for (int th = 0; th < nTh; th++)
                  c_vec[col_S(j, p, th)] -= c_stor[j];
    
    // ==========================================================================
    // SCHRITT 3: CONSTRAINTS (sparse tripel-Format, 1-basiert für R/Gurobi)
    // ==========================================================================
    std::vector<int>    row_v, col_v;
    std::vector<double> val_v;
    std::vector<double> rhs_v;
    std::vector<std::string> sense_v;
    int nrow = 0;
    
    auto push = [&](int r, int c, double v) {
      row_v.push_back(r + 1); // 1-basiert
      col_v.push_back(c + 1); // 1-basiert
      val_v.push_back(v);
    };
    auto add_con = [&](const std::string& sense, double rhs) {
      sense_v.push_back(sense);
      rhs_v.push_back(rhs);
      nrow++;
    };
    
    // ── C1: sum_{a in S^est} z[i,a] <= 1  ∀ i ───────────────────────────────
    for (int i = 0; i < ns; i++) {
      bool any = false;
      for (int a = 0; a < n_arcs; a++) {
        if (arcs[a].s != 0) continue;
        push(nrow, col_z(i, a), 1.0);
        any = true;
      }
      if (any) add_con("<=", 1.0);
    }
    if (verbose) Rcpp::Rcout << "  C1: Path establishment\n";
    
    // ── C2: sum_{(s,t)} z[i,s,t] = sum_{(t,u)} z[i,t,u]  ∀ i, t ────────────
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
    if (verbose) Rcpp::Rcout << "  C2: Path connectivity\n";
    
    // ── C3: sum_j Xij[i,j,p,t] <= sum_{harv-arcs: t=t} eta*area*z  ∀ i,p,t ─
    for (int i = 0; i < ns; i++) {
      for (int p = 0; p < P; p++) {
        for (int th = 0; th < nTh; th++) {
          int t = Tharv[th];
          std::vector<int> harv_arcs;
          for (int a = 0; a < n_arcs; a++)
            if (arcs[a].t == t && arcs[a].s >= 1 && arcs[a].s < t)
              harv_arcs.push_back(a);
            if (harv_arcs.empty()) continue;
            
            for (int j = 0; j < nj; j++)
              push(nrow, col_Xij(i, j, p, th), 1.0);
            
            for (int a : harv_arcs) {
              int age = t - arcs[a].s;
              if (age < 1 || age > Tm) continue;
              double coef = yield_matrix(p, age - 1) * area_ha[i];
              if (coef > 0.0)
                push(nrow, col_z(i, a), -coef);
            }
            add_con("<=", 0.0);
        }
      }
    }
    if (verbose) Rcpp::Rcout << "  C3: Biomass yield\n";
    
    // ── C4: Inventarbalance  ─────────────────────────────────────────────────
    // S[j,p,t] = S[j,p,t-1] + sum_i Xij[i,j,p,t] - sum_{k,pi:p_high=p} Xjk[j,k,pi,t]
    // Umgestellt: S[j,p,th] - S[j,p,th-1] - sum_i Xij + sum_{k,pi} Xjk == 0
    // th=0: kein S[th-1] Term => implizit S[j,p,Amin]=0 (OMPR C5)
    for (int j = 0; j < nj; j++) {
      for (int p = 0; p < P; p++) {
        for (int th = 0; th < nTh; th++) {
          push(nrow, col_S(j, p, th), +1.0);
          if (th > 0)
            push(nrow, col_S(j, p, th-1), -1.0);
          for (int i = 0; i < ns; i++)
            push(nrow, col_Xij(i, j, p, th), -1.0);
          // sum_{k, pi: p_high == p} Xjk[j,k,pi,th]
          for (int k = 0; k < nk; k++)
            for (int pi = 0; pi < n_ppairs; pi++) {
              if (pp_pairs[pi].first != p) continue; // p_high == p
              push(nrow, col_Xjk(j, k, pi, th), +1.0);
            }
            add_con("==", 0.0);
        }
      }
    }
    if (verbose) Rcpp::Rcout << "  C4: Inventory balance\n";
    
    // ── C5: Lagerkapazität: sum_p S[j,p,th] <= CAP_stor[j]  ∀ j, th ─────────
    for (int j = 0; j < nj; j++) {
      for (int th = 0; th < nTh; th++) {
        for (int p = 0; p < P; p++)
          push(nrow, col_S(j, p, th), 1.0);
        add_con("<=", CAP_stor[j]);
      }
    }
    if (verbose) Rcpp::Rcout << "  C5: Storage capacity\n";
    
    // ── C6: Verarbeitungskapazität: sum_{i,p} Xij[i,j,p,th] <= CAP_proc[j] ──
    for (int j = 0; j < nj; j++) {
      for (int th = 0; th < nTh; th++) {
        for (int i = 0; i < ns; i++)
          for (int p = 0; p < P; p++)
            push(nrow, col_Xij(i, j, p, th), 1.0);
        add_con("<=", CAP_proc[j]);
      }
    }
    if (verbose) Rcpp::Rcout << "  C6: Processing capacity\n";
    
    // ── C7: Nachfrage mit Kaskade: sum_{j, pi: p_low==p} Xjk <= D_max[k,p,t] ─
    for (int k = 0; k < nk; k++) {
      for (int p = 0; p < P; p++) {
        for (int th = 0; th < nTh; th++) {
          int t1 = Tharv[th];
          double Dmax = dem_lut[k * P * Tm + p * Tm + (t1 - 1)];
          if (Dmax < 0.0) continue; // kein Eintrag in demand-Tabelle
          bool any = false;
          for (int j = 0; j < nj; j++)
            for (int pi = 0; pi < n_ppairs; pi++) {
              if (pp_pairs[pi].second != p) continue; // p_low == p
              push(nrow, col_Xjk(j, k, pi, th), 1.0);
              any = true;
            }
            if (any) add_con("<=", Dmax);
        }
      }
    }
    if (verbose)
      Rcpp::Rcout << "  C7: Demand satisfaction\n"
                  << "  Total constraints: " << nrow
                  << ", nnz: " << (int)row_v.size() << "\n";
      
      // ==========================================================================
      // SCHRITT 4: BOUNDS UND TYPEN
      // ==========================================================================
      NumericVector lb(n_vars, 0.0);
      NumericVector ub(n_vars, R_PosInf);
      
      // z-Variablen: [0,1] binär
      for (int c = off_z; c < off_z + n_z; c++) ub[c] = 1.0;
      
      CharacterVector types(n_vars, "C");
      for (int c = off_z; c < off_z + n_z; c++) types[c] = "B";
      
      // ==========================================================================
      // SCHRITT 5: GUROBI-MODELL AUFBAUEN UND LÖSEN
      // ==========================================================================
      
      // Sparse Matrix für Gurobi (Matrix-Format über sparseMatrix in R)
      // Rückgabe der LP-Daten für solve_lp_with_gurobi()-Wrapper in R
      // ODER: direkter Gurobi-Aufruf via R's gurobi-Paket (über Rcpp -> R-Aufruf)
      
      // Gurobi-Modell als R-Liste zusammenbauen (kompatibel mit gurobi::gurobi())
      // A-Matrix: sparse (nrow x n_vars) im R sparseMatrix-Format
      // Wird als tripel-Liste übergeben und in R via Matrix::sparseMatrix() assembliert
      
      // LP-Objekt zurückgeben (für Wrapper build_agroforestry_lp_rcpp / solve_lp_with_gurobi)
      List lp_obj = List::create(
        Named("row_idx")   = wrap(row_v),
        Named("col_idx")   = wrap(col_v),
        Named("val")       = wrap(val_v),
        Named("rhs")       = wrap(rhs_v),
        Named("sense")     = wrap(sense_v),
        Named("c_vec")     = wrap(c_vec),
        Named("lb")        = lb,
        Named("ub")        = ub,
        Named("types")     = types,
        Named("n_vars")    = n_vars,
        Named("n_constrs") = nrow,
        Named("n_z")       = n_z,
        Named("n_Xij")     = n_Xij,
        Named("n_S")       = n_S,
        Named("n_Xjk")     = n_Xjk,
        Named("n_arcs")    = n_arcs,
        Named("n_ppairs")  = n_ppairs,
        Named("Tharv")     = wrap(Tharv),
        Named("sense_model") = "max"
      );
      
      // Gurobi direkt aus Rcpp via R-Funktionsaufruf
      Environment gurobi_env = Environment::namespace_env("gurobi");
      Function gurobi_fn     = gurobi_env["gurobi"];
      
      // A-Matrix assemblieren (via Matrix-Paket in R)
      Environment matrix_env = Environment::namespace_env("Matrix");
      Function spMatrix_fn   = matrix_env["sparseMatrix"];
      
      List A_mat = spMatrix_fn(
        Named("i") = wrap(row_v),
        Named("j") = wrap(col_v),
        Named("x") = wrap(val_v),
        Named("dims") = IntegerVector::create(nrow, n_vars)
      );
      
      // Gurobi-Modell-Liste
      // sense-Vektor: "<=" → "<", "==" → "=", ">=" → ">"
      CharacterVector g_sense(nrow);
      for (int r = 0; r < nrow; r++) {
        if      (sense_v[r] == "<=") g_sense[r] = "<";
        else if (sense_v[r] == "==") g_sense[r] = "=";
        else                          g_sense[r] = ">";
      }
      
      // Gurobi-Params zusammenführen
      List g_params = List::create(
        Named("OutputFlag") = (verbose ? 1 : 0)
      );
      // User-Params übernehmen
      CharacterVector param_names = gurobi_params.names();
      for (int pi = 0; pi < param_names.size(); pi++)
        g_params[as<std::string>(param_names[pi])] = gurobi_params[pi];
      
      List gurobi_model = List::create(
        Named("A")      = A_mat,
        Named("obj")    = wrap(c_vec),
        Named("sense")  = g_sense,
        Named("rhs")    = wrap(rhs_v),
        Named("lb")     = lb,
        Named("ub")     = ub,
        Named("vtype")  = types,
        Named("modelsense") = "max"
      );
      
      if (verbose)
        Rcpp::Rcout << "  Calling gurobi() with " << n_vars
                    << " vars | " << nrow << " constrs\n";
        
        List gurobi_result = gurobi_fn(gurobi_model, g_params);
        
        // ==========================================================================
        // SCHRITT 6: ERGEBNIS EXTRAHIEREN
        // ==========================================================================
        std::string status = as<std::string>(gurobi_result["status"]);
        double objval = 0.0;
        if (gurobi_result.containsElementNamed("objval"))
          objval = as<double>(gurobi_result["objval"]);
        
        if (verbose)
          Rcpp::Rcout << "  Status: " << status
                      << " | ObjVal: " << objval << "\n";
          
          // Lösungsvektor
          NumericVector x_sol(n_vars, 0.0);
          if (gurobi_result.containsElementNamed("x"))
            x_sol = as<NumericVector>(gurobi_result["x"]);
          
          // z-Lösung extrahieren: aktive Arcs (z > 0.5)
          std::vector<int>    z_i_vec, z_a_vec, z_s_vec, z_t_vec;
          std::vector<double> z_val_vec;
          for (int i = 0; i < ns; i++)
            for (int a = 0; a < n_arcs; a++) {
              double zv = x_sol[col_z(i, a)];
              if (zv > 0.5) {
                z_i_vec.push_back(i + 1);   // 1-basiert
                z_a_vec.push_back(a + 1);
                z_s_vec.push_back(arcs[a].s);
                z_t_vec.push_back(arcs[a].t);
                z_val_vec.push_back(zv);
              }
            }
            
            // Xij-Lösung extrahieren
            std::vector<int>    xij_i, xij_j, xij_p, xij_t;
          std::vector<double> xij_val;
          for (int i = 0; i < ns; i++)
            for (int j = 0; j < nj; j++)
              for (int p = 0; p < P; p++)
                for (int th = 0; th < nTh; th++) {
                  double v = x_sol[col_Xij(i, j, p, th)];
                  if (v > 1e-6) {
                    xij_i.push_back(i + 1);
                    xij_j.push_back(j + 1);
                    xij_p.push_back(p + 1);
                    xij_t.push_back(Tharv[th]);
                    xij_val.push_back(v);
                  }
                }
                
                // Xjk-Lösung extrahieren
                std::vector<int>    xjk_j, xjk_k, xjk_ph, xjk_pl, xjk_t;
          std::vector<double> xjk_val;
          for (int j = 0; j < nj; j++)
            for (int k = 0; k < nk; k++)
              for (int pi = 0; pi < n_ppairs; pi++)
                for (int th = 0; th < nTh; th++) {
                  double v = x_sol[col_Xjk(j, k, pi, th)];
                  if (v > 1e-6) {
                    xjk_j.push_back(j + 1);
                    xjk_k.push_back(k + 1);
                    xjk_ph.push_back(pp_pairs[pi].first  + 1); // p_high 1-basiert
                    xjk_pl.push_back(pp_pairs[pi].second + 1); // p_low  1-basiert
                    xjk_t.push_back(Tharv[th]);
                    xjk_val.push_back(v);
                  }
                }
                
                // S-Lösung extrahieren
                std::vector<int>    s_j, s_p, s_t;
          std::vector<double> s_val;
          for (int j = 0; j < nj; j++)
            for (int p = 0; p < P; p++)
              for (int th = 0; th < nTh; th++) {
                double v = x_sol[col_S(j, p, th)];
                if (v > 1e-6) {
                  s_j.push_back(j + 1);
                  s_p.push_back(p + 1);
                  s_t.push_back(Tharv[th]);
                  s_val.push_back(v);
                }
              }
              
              // arc-Typ-Vektor für z-Lösung
              CharacterVector z_type(z_a_vec.size());
          for (int r = 0; r < (int)z_a_vec.size(); r++) {
            int a = z_a_vec[r] - 1;
            if      (arcs[a].s == 0)      z_type[r] = "establishment";
            else if (arcs[a].t == Tm + 1) z_type[r] = "termination";
            else                           z_type[r] = "harvest";
          }
          
          DataFrame z_df = DataFrame::create(
            Named("i")    = wrap(z_i_vec),
            Named("a")    = wrap(z_a_vec),
            Named("s")    = wrap(z_s_vec),
            Named("t")    = wrap(z_t_vec),
            Named("type") = z_type,
            Named("value")= wrap(z_val_vec)
          );
          DataFrame xij_df = DataFrame::create(
            Named("i")    = wrap(xij_i),
            Named("j")    = wrap(xij_j),
            Named("p")    = wrap(xij_p),
            Named("t")    = wrap(xij_t),
            Named("value")= wrap(xij_val)
          );
          DataFrame xjk_df = DataFrame::create(
            Named("j")     = wrap(xjk_j),
            Named("k")     = wrap(xjk_k),
            Named("p")     = wrap(xjk_ph),
            Named("q")     = wrap(xjk_pl),
            Named("t")     = wrap(xjk_t),
            Named("value") = wrap(xjk_val)
          );
          DataFrame s_df = DataFrame::create(
            Named("j")    = wrap(s_j),
            Named("p")    = wrap(s_p),
            Named("t")    = wrap(s_t),
            Named("value")= wrap(s_val)
          );
          
          return List::create(
            Named("status")          = status,
            Named("objective_value") = objval,
            Named("gurobi_result")   = gurobi_result,
            Named("lp_obj")          = lp_obj,
            Named("solution_vector") = List::create(
              Named("z")   = z_df,
              Named("Xij") = xij_df,
              Named("Xjk") = xjk_df,
              Named("S")   = s_df
            )
          );
}