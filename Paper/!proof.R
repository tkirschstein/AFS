library(tidyverse)
library(tibble)
library(stringr)
library(leaflet)
library(osrm)
library(ROI)
library(ROI.plugin.glpk)
library(gurobi)
library(ROI.plugin.gurobi)
library(ompr)
library(magrittr)
library(ompr.roi)
library(sf)
library(dplyr)

source("../Modell/!afs_biomass_setup.r")
source("../Modell/!helper_func.r")
source("../Modell/!helper_instance_builder_v8a.R")
source("../Modell/!build_AFS_milp.r")
source("../Modell/build_lp_rcpp_wrapper.R")
source("../Modell/gurobi_adapter.R")
source("../Modell/gurobi_warmstart.R")
source("../Modell/gurobi_feascheck.R")
source("../Modell/extract_result.R")
#Rcpp::sourceCpp("../Modell/build_lp_rcpp.cpp")
#Rcpp::sourceCpp("../Modell/build_lp_rcpp_v2.cpp")
#Rcpp::sourceCpp("../Modell/build_lp_rcpp_v3.cpp")
Rcpp::sourceCpp("../Modell/build_and_solve_afs_milp.cpp")
Rcpp::sourceCpp("../Modell/build_and_solve_afs_lp_v11.cpp")

build_ompr_model <- function(milp_instance) {
 
  
  # 1) Sets und Basisgrößen ----------------------------------------------
  
  I <- seq_len(milp_instance$n_sites)
  J <- seq_len(milp_instance$n_storages)
  K <- seq_len(milp_instance$n_consumers)
  P <- seq_len(milp_instance$n_products)
  
  nbI <- milp_instance$n_sites
  nbJ <- milp_instance$n_storages
  nbK <- milp_instance$n_consumers
  nbP <- milp_instance$n_products
  nbT <- milp_instance$n_periods
  
  Amax <- milp_instance$max_age
  Amin <- milp_instance$min_age
  
  T <- seq_len(nbT+1)
  T_harv <- seq.int(from = max(1, Amin) + 1, to = nbT)
  
  A_i <- milp_instance$sites$area_ha
  
  # Kosten-Parameter pro Site
  C_est  <- milp_instance$sites$C_est
  C_harv <- milp_instance$sites$C_harv
  C_main <- milp_instance$sites$C_main
  C_opp  <- milp_instance$sites$C_opp
  
  # Lager-/Prozesskapazitäten und Lagerkosten
  CAP_st <- milp_instance$storages$CAP_stor
  CAP_pr <- milp_instance$storages$CAP_proc
  c_stor <- milp_instance$storages$c_stor
  
  # Transportkostenraten
  c_tr_raw <- milp_instance$c_tr_raw
  c_tr_pre <- milp_instance$c_tr_pre
  
  # Distanzen
  d_ij <- milp_instance$d_ij
  d_jk <- milp_instance$d_jk
  
  # Yields: yields_by_age(age, product, yield_ha)
  yields_by_age <- milp_instance$yields_by_age
  
  # Nachfrage D_kpt aus milp_instance$demand
  # demand: consumer_id, product, period, D_max
  demand_tbl <- milp_instance$demand
  
  # Preise R_kp aus consumer_prices
  price_tbl <- milp_instance$consumer_prices
  
  # Produkt-Kaskade Q: hier einfach P1 -> P2 -> P3 (hochwertig -> niedrig)
  # Wenn du eine andere Kaskade brauchst, passe das Mapping an.
  Q <- tibble::tribble(
    ~p_high, ~p_low,
    1L,      1L,   # P1 deckt P1
    1L,      2L,   # P1 kann P2 decken
    1L,      3L,   # P1 kann P3 decken
    2L,      2L,   # P2 deckt P2
    2L,      3L,   # P2 kann P3 decken
    3L,      3L    # P3 deckt P3
  ) %>% tibble::rowid_to_column("pi")
    
  
  # 1) Arcs (s,t) wie in build_lp_rcpp
  arcs_harvest <- expand.grid(s = 1:(nbT), t = 1:(nbT)) |>
    dplyr::filter(t > s) |>
    dplyr::filter( (t - s) >= Amin & (t - s) <= Amax) 

  arcs_est <- data.frame(s = 0, t = 1:(nbT-Amin))

  arcs_term <- data.frame(s = (Amin+1):nbT , t = nbT+1)
  
  arcs <- rbind(arcs_harvest, arcs_est, arcs_term) %>% 
    tibble::rowid_to_column("arc_id") %>% 
    mutate(type = ifelse(s==0, "establishment", ifelse(t == nbT+1, "termination", "harvest") ) )
  
  # 2) Indizes als Hilfsfunktionen
  #   z[i, arc_id]; Xij[i,j,p,th]; S[j,p,th]; Xjk[j,k,pi,th]
  idx_z   <- expand.grid(i = 1:nbI, arc_id = arcs$arc_id)
  idx_Xij <- expand.grid(i = 1:nbI, j = 1:nbJ, p = 1:nbP, th = 1:length(T_harv))
  idx_S   <- expand.grid(j = 1:nbJ, p = 1:nbP, th = 1:length(T_harv))
  idx_Xjk <- expand.grid(j = 1:nbJ, k = 1:nbK,
                         pi = Q$pi, th = 1:length(T_harv))
  
  # Hilfsfunktionen für R_kp und D_kpt ----------------------------------
  
  get_price <- function(k, p) {
    val <- price_tbl %>%
      filter(consumer_id == k, product == p) %>%
      pull(price)
    if (length(val) == 0) 0 else val[1]
  }
  
  get_demand <- function(k, p, t) {
    val <- demand_tbl %>%
      filter(consumer_id == k, product == p, period == t) %>%
      pull(D_max)
    if (length(val) == 0) 0 else val[1]
  }
  
  get_yield <- function(p, len) {
    val <- yields_by_age %>%
      filter(product == p, age == len) %>%
      pull(yield_ha)
    if (length(val) == 0) 0 else val[1]
  }
  
  # 2) Variablen ---------------------------------------------------------
  
  model <- MIPModel() %>%
    # z_{i,t} ∈ {0,1}: "AFS aktiv" in Periode t an Site i
    add_variable(z[i, a], i = 1:nbI, a = arcs$arc_id, type = "binary") |>
    # Xij, S, Xjk continuous >= 0
    add_variable(Xij[i, j, p, th],
                 i = 1:nbI, j = 1:nbJ, p = 1:nbP, th = T_harv, lb=0) |>
    add_variable(S[j, p, th],
                 j = 1:nbJ, p = 1:nbP, th = Amin:nbT, lb=0) |>
    add_variable(Xjk[j, k, pi, th] ,
                 j = 1:nbJ, k = 1:nbK, pi = Q$pi, th = T_harv, lb=0) %>% 
    add_variable(rev, lb=0) %>% 
    add_variable(cost_est, lb=0) %>% 
    add_variable(cost_stor, lb=0) %>% 
    add_variable(cost_main, lb=0) %>% 
    add_variable(cost_harv, lb=0) %>% 
    add_variable(cost_trans_raw, lb=0) %>% 
    add_variable(cost_trans_pre, lb=0) 
  
  # 3) Zielfunktion ------------------------------------------------------
  
  model <- model %>%
    set_objective(
      rev - cost_est - cost_stor - cost_main - cost_harv - cost_trans_raw - cost_trans_pre
      , sense = "max" )
  
  # 4) C1: Establishment höchstens einmal je Site ------------------------
  model <- model %>%
    add_constraint(sum_expr(z[i, a], a = arcs$arc_id[arcs$type == "establishment" ] ) <= 1, i = I)
  
  # C2 harvest consecutiation
  model <- model %>%
    add_constraint(sum_expr(z[i, a], a = arcs$arc_id[arcs$t == t] ) == sum_expr(z[i, a], a = arcs$arc_id[arcs$s == t] ), i = I, t = 1:nbT)
  
  # 5) C3 + C4 (kombiniert): Yield und Shipping-Link --------------------
  model <- model %>%
    add_constraint(
      sum_expr(Xij[i, j, p ,t], j = J ) <= 
        sum_expr({
          ss <- arcs[arcs$arc_id == a, "s"]
          tt <- arcs[arcs$arc_id == a, "t"]
          age <- tt-ss
          z[i, a] * get_yield(p, age) * A_i[i]
          }, a = arcs$arc_id[arcs$t == t & arcs$type == "harvest"] ), 
      i = I, p = P, t = T_harv)
  
  # C5 Lagerinitialisierung ----------------------------------------------------
  model <- model %>%
    add_constraint(
      S[j,p,Amin] == 0 , j = J, p = P)
  
  model <- model %>%
    add_constraint(
      S[j,p,t] == S[j,p,t-1] + 
        sum_expr( Xij[i,j,p,t], i = I ) - 
        sum_expr( Xjk[j,k,q,t] , q = Q$pi[Q$p_high == p], k = K), 
      j = J, p = P, t = T_harv)
  
  # 7) C6 Lagerkapazität -------------------------------------------------
  
  model <- model %>%
    add_constraint(
      sum_expr(S[j, p, t], p = P) <= CAP_st[j],
      j = J, t = T_harv
    )
  
  # 8) C7 Verarbeitungskapazität ----------------------------------------
  
  model <- model %>%
    add_constraint(
      sum_expr(Xij[i, j, p, t], i = I, p = P) <= CAP_pr[j],
      j = J, t = T_harv
    )
  
  # 9) C8 Nachfrage + Qualitätskaskade ----------------------------------
  
  model <- model %>%
    add_constraint(
      get_demand(k,p,t) >= sum_expr( Xjk[j,k,q,t] , q = Q$pi[Q$p_low == p], j = J) , 
      k = K, p = P, t = T_harv)
  
  # revenue
  model <- model %>% 
    add_constraint(
      rev == sum_expr({
      src <- Q$p_high[Q$pi == pi]
      pp  <- Q$p_low[Q$pi == pi]
      get_price(k, pp) * Xjk[j, k, pi, t]
    }, j = J, k = K, pi = Q$pi , t = T_harv))
  # establishment cost
  model <- model %>% 
    add_constraint(
      cost_est == sum_expr({
      C_est[i] * A_i[i] * z[i, a]
    }, i = I, a = arcs$arc_id[arcs$type == "establishment" ] ))
  # maintenance cost                 
  model <- model %>% 
    add_constraint(
      cost_main ==  sum_expr({
        s <- arcs[arcs$arc_id == a, "s"]
        t <- arcs[arcs$arc_id == a, "t"]
        (C_main[i] + C_opp[i]) * A_i[i] * (t-s) * z[i, a]
      }, i = I, a = arcs$arc_id[arcs$type =="harvest"])
    )
                   
  model <- model %>% 
    add_constraint(
      cost_harv == sum_expr({
        C_harv[i] * A_i[i] * z[i,a]
      }, i = I, a = arcs$arc_id[arcs$type == "harvest" ])
    )
  
  model <- model %>% 
    add_constraint(
      cost_trans_raw ==
        sum_expr({
          c_tr_raw * d_ij[i, j] * Xij[i, j, p, t]
        }, i = I, j = J, p = P, t = T_harv)  
    )
  
  model <- model %>% 
    add_constraint(
      cost_trans_pre ==
        sum_expr({
          c_tr_pre * d_jk[j, k] * Xjk[j, k, pi, t]
        }, j = J, k = K, pi = Q$pi, t = T_harv)
    )
  
  
  model <- model %>% 
    add_constraint(
      cost_stor ==
        sum_expr({
      c_stor[j] * S[j, p, t]
    }, j = J, p = P, t = T_harv)
    ) 
  
  
  # solve model
  result <- solve_model(model, with_ROI(solver = "gurobi"))
  
  # Zielfunktionswert
  obj_val <- objective_value(result)
  
  # Variablenwerte als Named‑Vektor
  sol_vec_z <- get_solution(result, z[i,a]) %>% 
    filter(value >.5)
  
  sol_vec_z <- sol_vec_z %>% 
    left_join(arcs %>% rename(a = arc_id), by = "a")
  
  sol_vec_Xij <- get_solution(result, Xij[i,j,p,t]) %>% 
    filter(value >.005)
  
  sol_vec_Xjk <- get_solution(result, Xjk[j,k,pi,t]) %>% 
    filter(value >.005)
  
  sol_vec_Xjk <- sol_vec_Xjk %>% 
    left_join(Q , by = "pi") %>% 
    rename(p = p_high, q = p_low)
  
  sol_vec_S <- get_solution(result, S[j,p,t]) %>% 
    filter(value >.005)
  
  sol_vec_rev <- get_solution(result, rev)
  
  sol_vec_cost_est <- get_solution(result, cost_est)
  
  sol_vec_cost_main <- get_solution(result, cost_main)
  
  sol_vec_cost_harv <- get_solution(result, cost_harv)
  
  sol_vec_cost_stor <- get_solution(result, cost_stor)
  
  sol_vec_cost_trans_pre <- get_solution(result, cost_trans_pre)
  
  sol_vec_cost_trans_raw <- get_solution(result, cost_trans_raw)
  
  cost.tab <- c(rev = sol_vec_rev, cost_est = sol_vec_cost_est , cost_main = sol_vec_cost_main, 
                cost_harv = sol_vec_cost_harv, cost_stor = sol_vec_cost_stor, cost_trans_pre = sol_vec_cost_trans_pre, cost_trans_raw = sol_vec_cost_trans_raw)
  
  list(
    model    = model,
    result   = result,
    objective_value = obj_val,
    solution_vector = list(z = sol_vec_z, Xij = sol_vec_Xij, Xjk = sol_vec_Xjk, S = sol_vec_S, obj.vec = cost.tab),
    helper     = list(arcs= arcs, Q = Q)
  )
  
  
  
}

# ==============================================================================
# solve_lp_given_z()
#
# Fixiert die z-Entscheidung (Anbau-Schedule) auf die vorgegebene z_fixed-Tabelle
# und löst dann das verbleibende LP (Xij, S, Xjk) als reines LP.
#
# Argumente:
#   milp_instance  – Instanzobjekt wie in build_ompr_model()
#   z_fixed        – data.frame/tibble mit Spalten (i, a), d.h. die aktiven Arcs
#                    Typischerweise: result$solution_vector$z %>% select(i, a)
#
# Rückgabe: analog zu build_ompr_model()
# ==============================================================================

solve_lp_given_z <- function(milp_instance, z_fixed) {
  
  # 1) Sets und Basisgrößen ----------------------------------------------
  
  I <- seq_len(milp_instance$n_sites)
  J <- seq_len(milp_instance$n_storages)
  K <- seq_len(milp_instance$n_consumers)
  P <- seq_len(milp_instance$n_products)
  
  nbI <- milp_instance$n_sites
  nbJ <- milp_instance$n_storages
  nbK <- milp_instance$n_consumers
  nbP <- milp_instance$n_products
  nbT <- milp_instance$n_periods
  
  Amax <- milp_instance$max_age
  Amin <- milp_instance$min_age
  
  T <- seq_len(nbT+1)
  T_harv <- seq.int(from = max(1, Amin) + 1, to = nbT)
  
  A_i <- milp_instance$sites$area_ha
  
  # Kosten-Parameter pro Site
  C_est  <- milp_instance$sites$C_est
  C_harv <- milp_instance$sites$C_harv
  C_main <- milp_instance$sites$C_main
  C_opp  <- milp_instance$sites$C_opp
  
  # Lager-/Prozesskapazitäten und Lagerkosten
  CAP_st <- milp_instance$storages$CAP_stor
  CAP_pr <- milp_instance$storages$CAP_proc
  c_stor <- milp_instance$storages$c_stor
  
  # Transportkostenraten
  c_tr_raw <- milp_instance$c_tr_raw
  c_tr_pre <- milp_instance$c_tr_pre
  
  # Distanzen
  d_ij <- milp_instance$d_ij
  d_jk <- milp_instance$d_jk
  
  # Yields: yields_by_age(age, product, yield_ha)
  yields_by_age <- milp_instance$yields_by_age
  
  # Nachfrage D_kpt aus milp_instance$demand
  # demand: consumer_id, product, period, D_max
  demand_tbl <- milp_instance$demand
  
  # Preise R_kp aus consumer_prices
  price_tbl <- milp_instance$consumer_prices
  
  # Produkt-Kaskade Q: hier einfach P1 -> P2 -> P3 (hochwertig -> niedrig)
  # Wenn du eine andere Kaskade brauchst, passe das Mapping an.
  Q <- tibble::tribble(
    ~p_high, ~p_low,
    1L,      1L,   # P1 deckt P1
    1L,      2L,   # P1 kann P2 decken
    1L,      3L,   # P1 kann P3 decken
    2L,      2L,   # P2 deckt P2
    2L,      3L,   # P2 kann P3 decken
    3L,      3L    # P3 deckt P3
  ) %>% tibble::rowid_to_column("pi")
  
  
  # 1) Arcs (s,t) wie in build_lp_rcpp
  arcs_harvest <- expand.grid(s = 1:(nbT), t = 1:(nbT)) |>
    dplyr::filter(t > s) |>
    dplyr::filter( (t - s) >= Amin & (t - s) <= Amax) 
  
  arcs_est <- data.frame(s = 0, t = 1:(nbT-Amin))
  
  arcs_term <- data.frame(s = (Amin+1):nbT , t = nbT+1)
  
  arcs <- rbind(arcs_harvest, arcs_est, arcs_term) %>% 
    tibble::rowid_to_column("arc_id") %>% 
    mutate(type = ifelse(s==0, "establishment", ifelse(t == nbT+1, "termination", "harvest") ) )
  
  # 2) Indizes als Hilfsfunktionen
  #   z[i, arc_id]; Xij[i,j,p,th]; S[j,p,th]; Xjk[j,k,pi,th]
  idx_z   <- expand.grid(i = 1:nbI, arc_id = arcs$arc_id)
  idx_Xij <- expand.grid(i = 1:nbI, j = 1:nbJ, p = 1:nbP, th = 1:length(T_harv))
  idx_S   <- expand.grid(j = 1:nbJ, p = 1:nbP, th = 1:length(T_harv))
  idx_Xjk <- expand.grid(j = 1:nbJ, k = 1:nbK,
                         pi = Q$pi, th = 1:length(T_harv))
  
  # Identify active and inactive sites
  id.active.sites <- unique(z_fixed$i)
  id.inactive.sites <- setdiff(I, id.active.sites)
  
  # Identify active arcs 
  
  # Hilfsfunktionen für R_kp und D_kpt ----------------------------------
  
  get_price <- function(k, p) {
    val <- price_tbl %>%
      filter(consumer_id == k, product == p) %>%
      pull(price)
    if (length(val) == 0) 0 else val[1]
  }
  
  get_demand <- function(k, p, t) {
    val <- demand_tbl %>%
      filter(consumer_id == k, product == p, period == t) %>%
      pull(D_max)
    if (length(val) == 0) 0 else val[1]
  }
  
  get_yield <- function(p, len) {
    val <- yields_by_age %>%
      filter(product == p, age == len) %>%
      pull(yield_ha)
    if (length(val) == 0) 0 else val[1]
  }
  
  # 2) Variablen ---------------------------------------------------------
  
  model <- MIPModel() %>%
    # z_{i,t} ∈ {0,1}: "AFS aktiv" in Periode t an Site i
    add_variable(z[i, a], i = 1:nbI, a = arcs$arc_id, type = "binary") |>
    # Xij, S, Xjk continuous >= 0
    add_variable(Xij[i, j, p, th],
                 i = 1:nbI, j = 1:nbJ, p = 1:nbP, th = T_harv, lb=0) |>
    add_variable(S[j, p, th],
                 j = 1:nbJ, p = 1:nbP, th = Amin:nbT, lb=0) |>
    add_variable(Xjk[j, k, pi, th] ,
                 j = 1:nbJ, k = 1:nbK, pi = Q$pi, th = T_harv, lb=0) %>% 
    add_variable(rev, lb=0) %>% 
    add_variable(cost_est, lb=0) %>% 
    add_variable(cost_stor, lb=0) %>% 
    add_variable(cost_main, lb=0) %>% 
    add_variable(cost_harv, lb=0) %>% 
    add_variable(cost_trans_raw, lb=0) %>% 
    add_variable(cost_trans_pre, lb=0) 
  
  # 3) Zielfunktion ------------------------------------------------------
  
  model <- model %>%
    set_objective(
      rev - cost_est - cost_stor - cost_main - cost_harv - cost_trans_raw - cost_trans_pre
      , sense = "max" )
  
  # 4) C1: Establishment höchstens einmal je Site ------------------------
  model <- model %>%
    add_constraint(sum_expr(z[i, a], a = arcs$arc_id[arcs$type == "establishment" ] ) <= 1, i = I)

  # C2 harvest consecutiation
  model <- model %>%
    add_constraint(sum_expr(z[i, a], a = arcs$arc_id[arcs$t == t] ) == sum_expr(z[i, a], a = arcs$arc_id[arcs$s == t] ), i = I, t = 1:nbT)

  # 5) C3 + C4 (kombiniert): Yield und Shipping-Link --------------------
  model <- model %>%
    add_constraint(
      sum_expr(Xij[i, j, p ,t], j = J ) <= 
        sum_expr({
          ss <- arcs[arcs$arc_id == a, "s"]
          tt <- arcs[arcs$arc_id == a, "t"]
          age <- tt-ss
          z[i, a] * get_yield(p, age) * A_i[i]
        }, a = arcs$arc_id[arcs$t == t & arcs$type == "harvest"] ), 
      i = I, p = P, t = T_harv)
  
  # C5 Lagerinitialisierung ----------------------------------------------------
  model <- model %>%
    add_constraint(
      S[j,p,Amin] == 0 , j = J, p = P)
  
  model <- model %>%
    add_constraint(
      S[j,p,t] == S[j,p,t-1] + 
        sum_expr( Xij[i,j,p,t], i = I ) - 
        sum_expr( Xjk[j,k,q,t] , q = Q$pi[Q$p_high == p], k = K), 
      j = J, p = P, t = T_harv)
  
  # 7) C6 Lagerkapazität -------------------------------------------------
  
  model <- model %>%
    add_constraint(
      sum_expr(S[j, p, t], p = P) <= CAP_st[j],
      j = J, t = T_harv
    )
  
  # 8) C7 Verarbeitungskapazität ----------------------------------------
  
  model <- model %>%
    add_constraint(
      sum_expr(Xij[i, j, p, t], i = I, p = P) <= CAP_pr[j],
      j = J, t = T_harv
    )
  
  # 9) C8 Nachfrage + Qualitätskaskade ----------------------------------
  
  model <- model %>%
    add_constraint(
      get_demand(k,p,t) >= sum_expr( Xjk[j,k,q,t] , q = Q$pi[Q$p_low == p], j = J) , 
      k = K, p = P, t = T_harv)
  
  # revenue
  model <- model %>% 
    add_constraint(
      rev == sum_expr({
        src <- Q$p_high[Q$pi == pi]
        pp  <- Q$p_low[Q$pi == pi]
        get_price(k, pp) * Xjk[j, k, pi, t]
      }, j = J, k = K, pi = Q$pi , t = T_harv))
  # establishment cost
  model <- model %>% 
    add_constraint(
      cost_est == sum_expr({
        C_est[i] * A_i[i] * z[i, a]
      }, i = I, a = arcs$arc_id[arcs$type == "establishment" ] ))
  # maintenance cost                 
  model <- model %>% 
    add_constraint(
      cost_main ==  sum_expr({
        s <- arcs[arcs$arc_id == a, "s"]
        t <- arcs[arcs$arc_id == a, "t"]
        (C_main[i] + C_opp[i]) * A_i[i] * (t-s) * z[i, a]
      }, i = I, a = arcs$arc_id[arcs$type =="harvest"])
    )
  
  model <- model %>% 
    add_constraint(
      cost_harv == sum_expr({
        C_harv[i] * A_i[i] * z[i,a]
      }, i = I, a = arcs$arc_id[arcs$type == "harvest" ])
    )
  
  model <- model %>% 
    add_constraint(
      cost_trans_raw ==
        sum_expr({
          c_tr_raw * d_ij[i, j] * Xij[i, j, p, t]
        }, i = I, j = J, p = P, t = T_harv)  
    )
  
  model <- model %>% 
    add_constraint(
      cost_trans_pre ==
        sum_expr({
          c_tr_pre * d_jk[j, k] * Xjk[j, k, pi, t]
        }, j = J, k = K, pi = Q$pi, t = T_harv)
    )
  
  
  model <- model %>% 
    add_constraint(
      cost_stor ==
        sum_expr({
          c_stor[j] * S[j, p, t]
        }, j = J, p = P, t = T_harv)
    ) 
  
  # FIX SITES ################################################################
  # ## inactivate
  # model <- model %>%
  #   add_constraint(
  #     sum_expr( z[i,a], a = arcs$arc_id ) == 0, i = as.character(id.inactive.sites)
  #   )
  # ## activate
  # model <- model %>% 
  #   add_constraint(
  #     sum_expr( z[i,a], a = arcs$arc_id ) >= 1, i = id.active.sites  
  #   )
  # 
  # 
  # FIX arcs ################################################################
  for(ii in id.active.sites){
    tmp.z.fix <- z_fixed[z_fixed$i == ii, ]
    
    tmp.z.fix <- tmp.z.fix %>% 
      select(-i) %>% 
      filter(s > 0, t < nbT+1)
    
    tmp.id <- arcs %>%
      filter(s > 0, t < nbT+1) %>% 
      dplyr::anti_join(tmp.z.fix, by = c("s", "t")) %>% 
      select(arc_id)
    
    tmp.id <- as.character(tmp.id$arc_id)
    
    # find arc id not belonging to the arcs in tmp.z.fix
    
    model <- model %>%
      add_constraint(
        sum_expr( z[ii,a], a = tmp.id) == 0 )

    # find arc id  belonging to the arcs in tmp.z.fix
    
    arcs_included <- arcs %>%
      dplyr::semi_join(tmp.z.fix, by = c("s", "t"))
     
    tmp.id2 <- as.character(arcs_included$arc_id)
    # find arc id not belonging to the arcs in tmp.z.fix
    
    model <- model %>% 
       add_constraint( 
         sum_expr( z[ ii, a ], a = tmp.id2) >= 1  )
  }
  
  # solve model
  result <- solve_model(model, with_ROI(solver = "gurobi"))
  
  # Zielfunktionswert
  obj_val <- objective_value(result)
  
  # Variablenwerte als Named‑Vektor
  sol_vec_z <- get_solution(result, z[i,a]) %>% 
    filter(value >.5)
  
  sol_vec_z <- sol_vec_z %>% 
    left_join(arcs %>% rename(a = arc_id), by = "a")
  
  sol_vec_Xij <- get_solution(result, Xij[i,j,p,t]) %>% 
    filter(value >.005)
  
  sol_vec_Xjk <- get_solution(result, Xjk[j,k,pi,t]) %>% 
    filter(value >.005)
  
  sol_vec_Xjk <- sol_vec_Xjk %>% 
    left_join(Q , by = "pi") %>% 
    rename(p = p_high, q = p_low)
  
  sol_vec_S <- get_solution(result, S[j,p,t]) %>% 
    filter(value >.005)
  
  sol_vec_rev <- get_solution(result, rev)
  
  sol_vec_cost_est <- get_solution(result, cost_est)
  
  sol_vec_cost_main <- get_solution(result, cost_main)
  
  sol_vec_cost_harv <- get_solution(result, cost_harv)
  
  sol_vec_cost_stor <- get_solution(result, cost_stor)
  
  sol_vec_cost_trans_pre <- get_solution(result, cost_trans_pre)
  
  sol_vec_cost_trans_raw <- get_solution(result, cost_trans_raw)
  
  cost.tab <- c(rev = sol_vec_rev, cost_est = sol_vec_cost_est , cost_main = sol_vec_cost_main, 
                cost_harv = sol_vec_cost_harv, cost_stor = sol_vec_cost_stor, cost_trans_pre = sol_vec_cost_trans_pre, cost_trans_raw = sol_vec_cost_trans_raw)
  
  list(
    model    = model,
    result   = result,
    objective_value = obj_val,
    solution_vector = list(z = sol_vec_z, Xij = sol_vec_Xij, Xjk = sol_vec_Xjk, S = sol_vec_S, obj.vec = cost.tab),
    helper     = list(arcs= arcs, Q = Q)
  )
  
  
}


#load("../Modell/afs_workspace.RData")
load("afs_workspace_red.RData")

sites_sf   <- afs_workspace$sites_sf
sites      <- afs_workspace$sites
sites_clustered      <- afs_workspace$sites_clustered
storages   <- afs_workspace$storages
consumers  <- afs_workspace$consumers

# scale down industrial demand quantities by dividing by 5
consumers <- consumers %>%
  mutate(demand_P1 = demand_P1 / 5)

# drop small Biogas plants --> compute total demand and drop all consumer with D<2
consumers <- consumers %>%
  mutate(total_demand = demand_P1 + demand_P2 + demand_P3) %>%
  filter(total_demand >= 95) %>%
  select(-total_demand) 


yields_by_age <- build_scenario_ts(ages = 1:20 , N = 250, C_site = 6500, k = 0.194, t0 = 9.7, label = "Conservative (250 trees/ha)")
colnames(yields_by_age)

yields_by_age <- yields_by_age %>% 
  # drop column scenario from yields_by_age
  select(-scenario) %>% 
  # make long --> bring columns Merch. stem Merch. branch Residue to column "product" and values to "yield_ha"
  pivot_longer(-1 , names_to = "product", values_to = "yield_ha") %>% 
  # rename entries in "product" to 1:3
  mutate(product = as.integer(as.character(factor(product, labels = c(2,1,3) ))))

# scale yields to fresh biomass
yields_by_age$yield_ha <- yields_by_age$yield_ha / (1-0.5)


sites_filt <- afs_workspace$sites %>%
  filter(area > 400) %>%
  rename(area_ha = area)

set.seed(123578)
# generate random vector of opportunity cost for clustered sites in range [0,600]
c.opp.vec <- rnorm(nrow(sites_filt), 1, .75)


sites_filt <- sites_filt %>% 
  mutate(
  C_est  = 2500,   # €/ha  establishment (Section 2.1 mid-range)
  C_harv =  150,    # €/ha  harvest cost  (Section 2.1 mid-range)
  C_main = 10 + 150*c.opp.vec, #  €/ha  maintenance cost  (Section 2.1 mid-range)
  C_opp = 0, #  €/ha  opportunity cost  (Section 2.1 mid-range)
) 


# # change storage costs
# storages$c_stor <- storages$c_stor
storages <- storages[1:3,] 
  
# # change capacities
storages$CAP_stor <- storages$CAP_stor * 100000
storages$CAP_proc <- storages$CAP_proc * 100000
# 


# calculate distances anew
dist_ij_filt <- calculate_distance_matrix_osm(
  starts       = sites_filt %>% select(lat, lng),
  destinations = storages %>% select(lat, lng),
  max_entries  = 100
)


dist_jk <- afs_workspace$dist_jk[1:3,consumers$consumer_id]


data_obj <- list(
  sites    = sites_filt,
  storages = storages,
  consumers = consumers,
  dist_ij  = dist_ij_filt$distance_matrix_km,
  dist_jk  = dist_jk,
  yields_by_age = yields_by_age
)
# rescale demands in consumers to tons
data_obj$consumers[,c("demand_P1","demand_P2","demand_P3")]  <- data_obj$consumers[,c("demand_P1","demand_P2","demand_P3")] * 1000

# 
# # ── 3. PARAMS OBJECT — all scalar parameters ──────────────────────────────────
params_obj <- list(
  n_periods  = 15,
  max_age    = 6L,    # A_max: LRAFS upper rotation bound
  min_age    =  3L,    # A_min: SRAFS lower rotation bound
  c_tr_raw   = 0.08,   # €/(t·km)  raw biomass  site → storage
  c_tr_pre   = 0.06   # €/(t·km)  chips        storage → consumer
)
# 
# ── 4. generate instance ─────────────────────────────────────
milp_instance <- build_optimization_instance(
  data         = data_obj,
  params       = params_obj
) 
# 

result.cpp <- build_and_solve_afs_milp(milp_instance,
                                   gurobi_params = list(TimeLimit = 300, MIPGap = 0.01),
                                   verbose = TRUE)

result.lp <- build_and_solve_afs_lp_v11(milp_instance,
                                       gurobi_params = list(TimeLimit = 300),
                                       verbose = TRUE)

result.cpp <- result.lp

################################################################################
# recalculate objective value of CPP solution
################################################################################

# -- 1. Harvest- und Maintenance-Kosten ----------------------------------------
oper.cost.cpp <- result.cpp$solution_vector$z |>
  filter(type == "harvest") |>
  mutate(len = t - s) |>
  select(i, len)

oper.cost.cpp <- oper.cost.cpp |>
  left_join(milp_instance$sites |> rename(i = site_id), by = "i") |>
  mutate(
    cost_harv = C_harv * area_ha,
    cost_main = len * (C_main + C_opp) * area_ha
  ) |>
  select(cost_harv, cost_main)

harv.cost.cpp <- sum(oper.cost.cpp$cost_harv)
main.cost.cpp <- sum(oper.cost.cpp$cost_main)

# -- 2. Establishment-Kosten ---------------------------------------------------
id.active.sites.cpp <- result.cpp$solution_vector$z |>
  filter(type == "establishment") |>
  pull(i) |> sort() |> unique()

est.cost.cpp <- milp_instance$sites[id.active.sites.cpp, ] |>
  select(C_est, area_ha) |>
  mutate(cost = C_est * area_ha) |>
  pull(cost) |> sum()

# -- 3. Lagerkosten ------------------------------------------------------------
if (nrow(result.cpp$solution_vector$S) == 0) {
  stor.cost.cpp <- 0
} else {
  stor.cost.cpp <- result.cpp$solution_vector$S |>
    left_join(
      milp_instance$storages |> mutate(j = row_number()) |> select(j, c_stor),
      by = "j"
    ) |>
    mutate(cost = c_stor * value) |>
    pull(cost) |> sum()
}

# -- 4. Transportkosten roh (Site → Lager) ------------------------------------
tmp.cost.ij <- milp_instance$c_tr_raw * milp_instance$d_ij   # Matrix [ns x nj]

tmp.flows.ij <- result.cpp$solution_vector$Xij |>
  mutate(
    i = factor(i, levels = seq_len(nrow(tmp.cost.ij))),
    j = factor(j, levels = seq_len(ncol(tmp.cost.ij)))
  ) |>
  group_by(i, j, .drop = FALSE) |>
  summarise(tot_val = sum(value), .groups = "drop") |>
  pivot_wider(names_from = j, values_from = tot_val,
              values_fill = 0, names_expand = TRUE) |>
  column_to_rownames("i")

trans.cost.raw.cpp <- sum(tmp.flows.ij * tmp.cost.ij)

# -- 5. Transportkosten verarbeitet (Lager → Konsument) -----------------------
tmp.cost.jk <- milp_instance$c_tr_pre * milp_instance$d_jk   # Matrix [nj x nk]

tmp.flows.jk <- result.cpp$solution_vector$Xjk |>
  mutate(
    j = factor(j, levels = seq_len(nrow(tmp.cost.jk))),
    k = factor(k, levels = seq_len(ncol(tmp.cost.jk)))
  ) |>
  group_by(j, k, .drop = FALSE) |>
  summarise(tot_val = sum(value), .groups = "drop") |>
  pivot_wider(names_from = k, values_from = tot_val,
              values_fill = 0, names_expand = TRUE) |>
  column_to_rownames("j")

trans.cost.pre.cpp   <- sum(tmp.flows.jk * tmp.cost.jk)
trans.cost.total.cpp <- trans.cost.pre.cpp + trans.cost.raw.cpp

# -- 6. Erlöse ----------------------------------------------------------------
# Xjk enthält q = plow = geliefertes Produkt → Preis über (k, q)
tmp.flows.rev <- result.cpp$solution_vector$Xjk |>
  mutate(
    k = factor(k, levels = seq_len(ncol(tmp.cost.jk))),
    q = factor(q, levels = 1:3)
  ) |>
  group_by(k, q, .drop = FALSE) |>
  summarise(tot_val = sum(value), .groups = "drop") |>
  pivot_wider(names_from = q, values_from = tot_val,
              values_fill = 0, names_expand = TRUE) |>
  column_to_rownames("k")

tmp.rev.mat <- milp_instance$consumer_prices |>
  rename(k = consumer_id, q = product) |>
  mutate(
    k = factor(k, levels = seq_len(ncol(tmp.cost.jk))),
    q = factor(q, levels = 1:3)
  ) |>
  pivot_wider(names_from = q, values_from = price,
              values_fill = 0, names_expand = TRUE) |>
  column_to_rownames("k")

tot.rev.cpp <- sum(tmp.rev.mat * tmp.flows.rev)

# -- 7. Gesamtkosten & Zielfunktionswert --------------------------------------
tot.cost.cpp    <- harv.cost.cpp + main.cost.cpp + est.cost.cpp +
  stor.cost.cpp + trans.cost.total.cpp
obj.recalc.cpp  <- tot.rev.cpp - tot.cost.cpp

# -- 8. Ergebnistabelle -------------------------------------------------------
cost.tab.cpp <- c(
  rev            = tot.rev.cpp,
  cost_est       = est.cost.cpp,
  cost_main      = main.cost.cpp,
  cost_harv      = harv.cost.cpp,
  cost_stor      = stor.cost.cpp,
  cost_trans_pre = trans.cost.pre.cpp,
  cost_trans_raw = trans.cost.raw.cpp
)


obj.recalc.cpp - result.cpp$objective_value



############################################################################################
## OLD STUFF
#############################################################################################

{
  
  # Solve with OMPR Function
  
  result.ompr <- build_ompr_model(milp_instance)
  
  # Compare objective values ###################################################
  
  result.cpp$objective_value/10^3
  result.ompr$objective_value/10^3
  
  # total production 
  
  tmp.ompr <- result.ompr$solution_vector$Xij %>% 
    group_by(i,t) %>% 
    summarise(tot_val = sum(value))
  
  sum(tmp.ompr$tot_val) - sum(result.ompr$solution_vector$Xij$value)
  
  
  result.ompr$solution_vector$z %>% 
    left_join(tmp.ompr, by= c("i" , "t")) %>% 
    select(i,s,t,tot_val, type)
  
  # total shipments
  
  tmp.ompr <- result.ompr$solution_vector$Xij %>% 
    group_by(j,t) %>% 
    summarise(input = sum(value))
  
  tot.ship <- result.ompr$solution_vector$Xjk %>% 
    group_by(j,t) %>% 
    summarise(output = sum(value)) %>% 
    left_join(tmp.ompr, by =c("j","t"))
  
  sum(tot.ship$input)
  
  # recalculate objective value
  
  ## transport costs
  
  ### raw
  tmp.flows <- tmp.cost <- milp_instance$c_tr_raw * milp_instance$d_ij
  
  tmp.flows <- tmp.flows*0
  
  tmp.flows <- result.ompr$solution_vector$Xij %>% 
    mutate(i = factor(i, levels = 1:nrow(tmp.flows) ), j = factor(j, levels = 1:ncol(tmp.flows)) ) %>% 
    group_by(i,j, .drop = FALSE) %>% 
    summarise(tot_val = sum(value)) %>% 
    pivot_wider(names_from = c("j"), values_from = tot_val, values_fill = 0, names_expand =TRUE) %>% 
    column_to_rownames("i")
  
  trans.cost.raw <- sum(tmp.flows*tmp.cost)
  
  result.ompr$solution_vector$obj.vec["cost_trans_raw.cost_trans_raw"] - trans.cost.raw
  
  ### pre
  tmp.cost <- milp_instance$c_tr_pre * milp_instance$d_jk
  
  tmp.flows <- result.ompr$solution_vector$Xjk %>% 
    mutate(j = factor(j, levels = 1:nrow(tmp.cost) ), k = factor(k, levels = 1:ncol(tmp.cost)) ) %>% 
    group_by(j,k, .drop = FALSE) %>% 
    summarise(tot_val = sum(value)) %>% 
    pivot_wider(names_from = c("k"), values_from = tot_val, values_fill = 0, names_expand =TRUE) %>% 
    column_to_rownames("j")
  
  trans.cost.pre <- sum(tmp.flows*tmp.cost)
  
  result.ompr$solution_vector$obj.vec["cost_trans_pre.cost_trans_pre"] - trans.cost.pre
  
  ### total 
  trans.cost.total <- trans.cost.pre+ trans.cost.raw
  
  ## storage cost 
  if(nrow(result.ompr$solution_vector$S) == 0) stor.cost <- 0 else{
    
  }
  
  ## establishment cost
  
  id.active.sites <- result.ompr$solution_vector$z %>% 
    filter(type== "establishment") %>% 
    select(i) %>% 
    unlist() %>% 
    sort()
  
  est.cost <- milp_instance$sites[id.active.sites,] %>% 
    select(C_est, area_ha) %>% 
    mutate(cost = C_est*area_ha) %>% 
    select(cost)
  
  est.cost <- sum(est.cost)
  
  result.ompr$solution_vector$obj.vec["cost_est.cost_est"] - est.cost
  
  ## harvest & main + opp costs 
  
  oper.cost <- result.ompr$solution_vector$z %>% 
    filter(type== "harvest") %>% 
    mutate(len = t-s) %>% 
    select(i, len)
  
  oper.cost <- oper.cost %>% 
    left_join(milp_instance$sites %>% rename(i = site_id), by="i") %>% 
    mutate(cost_harv = (C_harv)*area_ha, cost_main = len * (C_main + C_opp) * area_ha ) %>% 
    select(cost_harv, cost_main)
  
  harv.cost <- sum(oper.cost$cost_harv)
  main.cost <- sum(oper.cost$cost_main)
  
  result.ompr$solution_vector$obj.vec["cost_main.cost_main"] - main.cost
  
  result.ompr$solution_vector$obj.vec["cost_harv.cost_harv"] - harv.cost
  
  ## total cost
  
  tot.cost <- (harv.cost + main.cost + est.cost + stor.cost + trans.cost.total)
  
  
  ## revenue
  tmp.flows <- result.ompr$solution_vector$Xjk %>% 
    mutate(k = factor(k, levels = 1:ncol(tmp.cost)), q = factor(q, levels = 1:3) ) %>% 
    group_by(k, q, .drop = FALSE) %>% 
    summarise(tot_val = sum(value)) %>% 
    pivot_wider(names_from = c("q"), values_from = tot_val, values_fill = 0, names_expand =TRUE) %>% 
    column_to_rownames("k")
  
  tmp.rev <- milp_instance$consumer_prices %>% 
    rename(k = consumer_id, q = product) %>% 
    mutate(k = factor(k, levels = 1:ncol(tmp.cost)), q = factor(q, levels = 1:3) ) %>% 
    pivot_wider(names_from = c("q"), values_from = price, values_fill = 0, names_expand =TRUE) %>% 
    column_to_rownames("k")
  
  
  tot.rev <- sum(tmp.rev * tmp.flows)
  
  result.ompr$solution_vector$obj.vec["rev.rev"] - tot.rev
  
  # objective value
  (tot.rev -  tot.cost) - result.ompr$objective_value
  
  
  # -- 9. Vergleich OMPR vs. C++ -----------------------------------------------
  cat("\n========== Nachkalkulation: OMPR vs. C++ ==========\n")
  comparison <- data.frame(
    Kostenart  = names(cost.tab.cpp),
    OMPR       = as.numeric(result.ompr$solution_vector$obj.vec[names(cost.tab.cpp)]),
    CPP_recalc = as.numeric(cost.tab.cpp),
    Differenz  = as.numeric(cost.tab.cpp) -
      as.numeric(result.ompr$solution_vector$obj.vec[names(cost.tab.cpp)])
  )
  print(comparison, digits = 6)
  
  cat("\n--- Zielfunktionswerte ---\n")
  cat(sprintf("OMPR  reported  : %12.2f\n", result.ompr$objective_value))
  cat(sprintf("OMPR  recalc    : %12.2f\n",
              result.ompr$solution_vector$obj.vec[["rev.rev"]] -
                sum(as.numeric(result.ompr$solution_vector$obj.vec[
                  c("cost_est.cost_est","cost_main.cost_main","cost_harv.cost_harv",
                    "cost_stor.cost_stor","cost_trans_pre.cost_trans_pre",
                    "cost_trans_raw.cost_trans_raw")]))))
  cat(sprintf("C++   reported  : %12.2f\n", result.cpp$objective_value))
  cat(sprintf("C++   recalc    : %12.2f\n", obj.recalc.cpp))
  cat(sprintf("Diff  (cpp-ompr): %12.2f\n",
              result.cpp$objective_value - result.ompr$objective_value))
  

  
  # z-Lösung aus dem MILP-Ergebnis extrahieren
  z_fix <- result.cpp$solution_vector$z %>% select(i, s, t)
  
  # LP mit fixiertem Anbaupfad lösen
  result.lp <- solve_lp_given_z(milp_instance, z_fix)
  
  # Vergleich
  result.cpp$objective_value / 1e3
  result.lp$objective_value  / 1e3   # sollte ≈ MILP-ObjFkt sein
  
  
    
}
