
# ── 1. Hilfsfunktion: Site-Level KPIs für ein Szenario extrahieren ──────────
extract_site_profit <- function(res, scenario_name = NA_character_) {
  
  milp_instance <- res$milp_instance
  ext           <- res$ext
  setting       <- res$setting
  
  # Szenario-Parameter
  scen_opp      <- if (!any(is.na(setting))) setting$opp      else 150
  scen_cost_log <- if (!any(is.na(setting))) setting$cost.log else 1.0   
  scen_cost_est <- if (!any(is.na(setting))) setting$cost.est else 1.0   
  scen_revenue  <- if (!any(is.na(setting))) setting$revenue  else 1.0
  
  # no open sites
  if(nrow(ext$z) == 0){
    
    result_df <- data.frame( 
      scenario = scenario_name, 
      opp_cost = scen_opp,
      cost_log_level       = scen_cost_log,   
      cost_est_level       = scen_cost_est,   
      revenue_level        = scen_revenue,
      profit_ha_yr = NA,
      opp_cost = NA, 
      opp_cost_site = NA, 
      avg_dist_hub_km = NA, 
      avg_dist_consumer_km = NA,
      avg_rotation_yr = NA, 
      n_harvests = NA, 
      share_p1 = NA,
      area_ha =NA, 
      area_afs = NA, 
      active_years = NA, 
      n_sites = NA)
    
  }
  else{
    # Aktive Standorte
    active_site_ids <- ext$z %>%
      filter(arc_type == "harvest") %>%
      pull(site_id) %>% unique()
    
    site_meta <- milp_instance$sites %>%
      mutate(site_id = row_number()) %>%
      filter(site_id %in% active_site_ids)
    
    # Preistabelle
    price_df <- as.data.frame(milp_instance$consumer_prices) %>%
      setNames(c("consumer_id", "del_product", "price")) %>%
      mutate(consumer_id = as.integer(consumer_id),
             del_product  = as.integer(del_product))
    
    # ── Distanzmatrizen ────────────────────────────────────────────────────────
    dist_ij_m <- as.matrix(milp_instance$d_ij)
    dist_jk_m <- as.matrix(milp_instance$d_jk)
    
    # Hub-Anteile pro Standort und Periode (für Kostenzurechnung)
    hub_site_share <- ext$Xij %>%
      group_by(hub_id, period) %>%
      mutate(hub_total = sum(value)) %>%
      ungroup() %>%
      mutate(share = ifelse(hub_total > 0, value / hub_total, 0)) %>%
      select(site_id, hub_id, period, share)
    
    # ── Kosten je Standort ─────────────────────────────────────────────────────
    sc_est <- ext$z %>% filter(arc_type == "establishment") %>%
      left_join(site_meta, by = "site_id") %>%
      mutate(cost = C_est * value) %>%
      group_by(site_id) %>% summarise(cost = sum(cost, na.rm = TRUE), .groups = "drop") %>%
      mutate(component = "Establishment")
    
    sc_harv <- ext$z %>% filter(arc_type == "harvest") %>%
      left_join(site_meta, by = "site_id") %>%
      mutate(cost = C_harv * value) %>%
      group_by(site_id) %>% summarise(cost = sum(cost, na.rm = TRUE), .groups = "drop") %>%
      mutate(component = "Harvesting")
    
    sc_main <- ext$z %>% filter(arc_type == "harvest") %>%
      mutate(age_len = t - s) %>%
      left_join(site_meta, by = "site_id") %>%
      mutate(cost = C_main * value * age_len) %>%
      group_by(site_id) %>% summarise(cost = sum(cost, na.rm = TRUE), .groups = "drop") %>%
      mutate(component = "Maintenance")
    
    sc_opp <- ext$z %>% filter(arc_type == "harvest") %>%
      mutate(age_len = t - s) %>%
      left_join(site_meta, by = "site_id") %>%
      mutate(cost = C_opp * value * age_len) %>%
      group_by(site_id) %>% summarise(cost = sum(cost, na.rm = TRUE), .groups = "drop") %>%
      mutate(component = "Opportunity")
    
    sc_tr_raw <- ext$Xij %>%
      mutate(dist_km = mapply(function(i, j) dist_ij_m[i, j], site_id, hub_id),
             cost    = milp_instance$c_tr_raw * dist_km * value) %>%
      group_by(site_id) %>% summarise(cost = sum(cost, na.rm = TRUE), .groups = "drop") %>%
      mutate(component = "Transport_raw")
    
    sc_tr_pre <- ext$Xjk %>%
      mutate(dist_km  = mapply(function(j, k) dist_jk_m[j, k], hub_id, consumer_id),
             cost_jkt = milp_instance$c_tr_pre * dist_km * value) %>%
      left_join(hub_site_share, by = c("hub_id", "period"), relationship = "many-to-many") %>%
      mutate(site_cost = cost_jkt * share) %>%
      group_by(site_id) %>% summarise(cost = sum(site_cost, na.rm = TRUE), .groups = "drop") %>%
      filter(!is.na(site_id)) %>%
      mutate(component = "Transport_pre")
    
    sc_stor <- ext$S %>%
      left_join(milp_instance$storages %>% mutate(hub_id = row_number()) %>% select(hub_id, c_stor),
                by = "hub_id") %>%
      mutate(stor_cost = c_stor * value) %>%
      select(hub_id, period, stor_cost) %>%
      left_join(hub_site_share, by = c("hub_id", "period"), relationship = "many-to-many") %>%
      mutate(site_cost = stor_cost * share) %>%
      group_by(site_id) %>% summarise(cost = sum(site_cost, na.rm = TRUE), .groups = "drop") %>%
      filter(!is.na(site_id)) %>%
      mutate(component = "Storage")
    
    # ── Erlöse je Standort ─────────────────────────────────────────────────────
    site_revenue <- ext$Xjk %>%
      left_join(price_df, by = c("consumer_id", "del_product")) %>%
      mutate(rev_jkt = value * coalesce(price, 0)) %>%
      left_join(hub_site_share, by = c("hub_id", "period"), relationship = "many-to-many") %>%
      mutate(site_rev = rev_jkt * share) %>%
      group_by(site_id) %>%
      summarise(revenue = sum(site_rev, na.rm = TRUE), .groups = "drop") %>%
      filter(!is.na(site_id))
    
    # ── Normierung: €/ha/a (nur aktive Produktionsjahre) ──────────────────────
    #est_p  <- ext$z %>% filter(arc_type == "establishment") %>% select(site_id, t_est = t)
    #term_p <- ext$z %>% filter(arc_type == "termination")  %>% select(site_id, t_term = s)
    active_p <- ext$z %>% filter(arc_type == "establishment") %>% 
      group_by(site_id) %>% 
      summarise(tot_value = sum(value))
    # 
    est_p  <- ext$z %>% filter(arc_type == "establishment") %>% group_by(site_id) %>% summarise(t_est = mean(sum(t*value)/sum(value)))
    term_p <- ext$z %>% filter(arc_type == "termination")  %>% group_by(site_id) %>% summarise(t_term = mean(sum(s*value)/sum(value)))
    
    site_norm <- site_meta %>%
      #select(site_id, area_ha) %>%
      left_join(active_p,  by = "site_id") %>%
      left_join(est_p,  by = "site_id") %>%
      left_join(term_p, by = "site_id") %>%
      mutate(t_term       = ifelse(is.na(t_term), milp_instance$n_periods, t_term),
             active_years = pmax(t_term - t_est, 1L),
             denom        = tot_value * active_years)
    
    all_costs <- bind_rows(sc_est, sc_harv, sc_main, sc_opp,
                           sc_tr_raw, sc_tr_pre, sc_stor) %>%
      group_by(site_id) %>%
      summarise(total_cost = sum(cost, na.rm = TRUE), .groups = "drop")
    
    # ── Durchschnittliche Distanzen je Standort ────────────────────────────────
    # (a) Ø Distanz Standort → Hub (gewichtet nach Liefervolumen)
    avg_dist_to_hub <- ext$Xij %>%
      mutate(dist_km = mapply(function(i, j) dist_ij_m[i, j], site_id, hub_id)) %>%
      group_by(site_id) %>%
      summarise(avg_dist_hub_km = weighted.mean(dist_km, w = pmax(value, 0), na.rm = TRUE),
                .groups = "drop")
    
    # (b) Ø Distanz Hub → Kunden (anteilig auf Standort zurückgerechnet)
    avg_dist_to_consumer <- ext$Xjk %>%
      mutate(dist_km = mapply(function(j, k) dist_jk_m[j, k], hub_id, consumer_id)) %>%
      left_join(hub_site_share, by = c("hub_id", "period"), relationship = "many-to-many") %>%
      filter(!is.na(site_id), share > 0) %>%
      group_by(site_id) %>%
      summarise(avg_dist_consumer_km = weighted.mean(dist_km, w = pmax(share * value, 0), na.rm = TRUE),
                .groups = "drop")
    
    # ── Weitere strukturelle Variablen ─────────────────────────────────────────
    # (c) Ø Erntezyklus-Länge (Rotationsalter) je Standort
    avg_rotation <- ext$z %>%
      filter(arc_type == "harvest") %>%
      mutate(cycle_len = t - s) %>%
      group_by(site_id) %>%
      summarise(avg_rotation_yr = sum(cycle_len*value/sum(value), na.rm = TRUE),
                n_harvests       = n(),
                .groups = "drop")
    
    # (d) Produktmix-Index: Anteil hochwertiger Biomasse (P1 Stämme) an Gesamtlieferung
    p1_share <- ext$Xij %>%
      group_by(site_id) %>%
      summarise(
        vol_total = sum(value, na.rm = TRUE),
        vol_p1    = sum(value[p == 1], na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      mutate(share_p1 = ifelse(vol_total > 0, vol_p1 / vol_total, NA_real_)) %>%
      select(site_id, share_p1, vol_total)
    
    # ── Zusammenführen ─────────────────────────────────────────────────────────
    result_df <- site_norm %>%
      left_join(site_revenue,          by = "site_id") %>%
      left_join(all_costs,             by = "site_id") %>%
      left_join(avg_dist_to_hub,       by = "site_id") %>%
      left_join(avg_dist_to_consumer,  by = "site_id") %>%
      left_join(avg_rotation,          by = "site_id") %>%
      left_join(p1_share,              by = "site_id") %>%
      mutate(
        profit_ha_yr         = (coalesce(revenue, 0) - coalesce(total_cost, 0)) / denom,
        scenario             = scenario_name,
        opp_cost             = scen_opp,
        cost_log_level       = scen_cost_log,   # renamed: was cost_level
        cost_est_level       = scen_cost_est,   # new
        revenue_level        = scen_revenue,
        opp_cost_site        = C_opp,
        area_afs             = tot_value
      ) %>%
      select(site_id, scenario, profit_ha_yr,
             opp_cost, opp_cost_site, cost_log_level, cost_est_level, revenue_level,
             avg_dist_hub_km, avg_dist_consumer_km,
             avg_rotation_yr, n_harvests, share_p1,
             area_ha, area_afs, active_years, n_sites)
    
  }
  
  
  
  return(result_df)
}
