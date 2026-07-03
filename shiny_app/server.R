library(shiny)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(plotly)
library(leaflet)
library(sf)
library(DT)
library(networkD3)
library(scales)
library(ggalluvial)
library(Rcpp)

# Define server logic required to draw a histogram
function(input, output, session) {
  # --------------------------------------------------------------------------
  # 1) INITIALISIERUNG
  # --------------------------------------------------------------------------
  rv <- reactiveValues(
    afs_workspace = NULL,
    sites_sf = NULL,
    sites = NULL,
    storages = NULL,
    consumers = NULL,
    yields_by_age = NULL,
    data_obj = NULL,
    params_obj = NULL,
    milp_instance = NULL,
    solve_result = NULL,
    ext = NULL,
    site_profit = NULL,
    sensitivity = NULL,
    status = "Idle"
  )
  
  
  observe({
    # Daten & Funktionen laden (beim Start)
    if (is.null(rv$afs_workspace)) {
      showNotification("Loading workspace and source files ...", type = "message")
      
      source("!afs_biomass_setup.r")
      source("!helper_func.r")
      source("!helper_instance_builder_v8a.R")
      source("build_lp_rcpp_wrapper.R")
      
      # Falls C++-Solver dynamisch genutzt werden soll:
      Rcpp::sourceCpp("build_and_solve_afs_lp_v11.cpp")
      
      load("afs_workspace_red.RData")
      rv$afs_workspace <- afs_workspace
      rv$sites_sf      <- afs_workspace$sites_sf
      rv$sites         <- afs_workspace$sites_clustered
      rv$storages      <- afs_workspace$storages
      rv$consumers     <- afs_workspace$consumers
      
      # Vorverarbeitung analog paper_lp.rmd
      rv$consumers <- rv$consumers %>%
        mutate(
          demand_P1 = demand_P1 / 5,
          P1 = P1 / 2,
          P2 = P2 / 2,
          P3 = P3 * 0.6
        ) %>%
        mutate(total_demand = demand_P1 + demand_P2 + demand_P3) %>%
        filter(total_demand >= 1.5) %>%
        select(-total_demand)
      
      rv$storages$CAP_stor <- rv$storages$CAP_stor * 100000
      rv$storages$CAP_proc <- rv$storages$CAP_proc * 100000
      
      rv$status <- "Base data loaded"
    }
  })
  
  
  # --------------------------------------------------------------------------
  # 2) HILFSFUNKTIONEN
  # --------------------------------------------------------------------------
  
  build_growth_df <- reactive({
    req(rv$afs_workspace)
    
    age.vec <- seq(0, max(20, input$max_age), by = 0.1)
    
    tmp <- build_scenario_ts(
      ages   = age.vec,
      N      = input$N_trees,
      C_site = input$C_site,
      k      = input$k_gomp,
      t0     = input$t0_gomp,
      label  = "Interactive scenario"
    )
    
    # Vereinheitlichung wie in paper_lp.rmd
    tmp_long <- tmp %>%
      rename(
        stem   = "Merch. stem",
        branch = "Merch. branch",
        residue = Residue
      ) %>%
      pivot_longer(c(stem, branch, residue), names_to = "component", values_to = "biomass_dm") %>%
      group_by(age) %>%
      mutate(total_dm = sum(biomass_dm, na.rm = TRUE),
             share = ifelse(total_dm > 0, biomass_dm / total_dm, 0)) %>%
      ungroup()
    
    tmp_long
  })
  
  build_yields_input <- reactive({
    req(rv$afs_workspace)
    
    y <- build_scenario_ts(
      ages   = 1:input$max_age,
      N      = input$N_trees,
      C_site = input$C_site,
      k      = input$k_gomp,
      t0     = input$t0_gomp,
      label  = "Interactive scenario"
    )
    
    y %>%
      select(-scenario) %>%
      pivot_longer(-1, names_to = "product", values_to = "yield_ha") %>%
      mutate(product = as.integer(as.character(factor(product, labels = c(2, 1, 3))))) %>%
      mutate(yield_ha = yield_ha / (1 - input$moisture))
  })
  
  build_model_data <- reactive({
    req(rv$sites, rv$storages, rv$consumers, rv$afs_workspace)
    
    set.seed(123578)
    c.opp.vec <- rnorm(nrow(rv$sites), 1, .3)
    
    sites_model <- rv$sites %>%
      mutate(
        C_est  = input$C_est * input$est_mult,
        C_harv = input$C_harv * input$log_mult,
        C_main = input$C_main,
        C_opp  = input$opp_mean * c.opp.vec
      )
    
    consumers_model <- rv$consumers %>%
      mutate(
        P1 = P1 * input$rev_mult,
        P2 = P2 * input$rev_mult,
        P3 = P3 * input$rev_mult
      )
    
    data_obj <- list(
      sites = sites_model,
      storages = rv$storages,
      consumers = consumers_model,
      dist_ij = rv$afs_workspace$dist_ij_clust,
      dist_jk = rv$afs_workspace$dist_jk[, consumers_model$consumer_id, drop = FALSE],
      yields_by_age = build_yields_input()
    )
    
    data_obj$consumers[, c("demand_P1", "demand_P2", "demand_P3")] <-
      data_obj$consumers[, c("demand_P1", "demand_P2", "demand_P3")] * 1000
    
    params_obj <- list(
      n_periods = input$n_periods,
      max_age   = as.integer(input$max_age),
      min_age   = as.integer(input$min_age),
      c_tr_raw  = input$c_tr_raw * input$log_mult,
      c_tr_pre  = input$c_tr_pre * input$log_mult
    )
    
    list(data_obj = data_obj, params_obj = params_obj)
  })
  
  prepare_solution_objects <- function(result) {
    ext <- result$solution_vector
    
    ext$Xjk <- ext$Xjk %>%
      rename(src_product = p, del_product = q, period = t, hub_id = j, consumer_id = k)
    
    ext$Xij <- ext$Xij %>%
      rename(period = t, site_id = i, hub_id = j)
    
    ext$z <- ext$z %>%
      rename(arc_type = type, site_id = i)
    
    ext$S <- ext$S %>%
      rename(hub_id = j, period = t)
    
    ext
  }
  
  build_site_profit <- reactive({
    req(rv$milp_instance, rv$ext)
    
    # vereinfachte Delegation an Funktion aus paper_lp.rmd
    extract_site_profit(list(
      setting = list(
        opp      = input$opp_mean,
        cost.log = input$log_mult,
        cost.est = input$est_mult,
        revenue  = input$rev_mult
      ),
      milp_instance = rv$milp_instance,
      ext = rv$ext
    ), scenario_name = "Interactive")
  })
  
  # --------------------------------------------------------------------------
  # 3) OPTIMIERUNG AUSLÖSEN
  # --------------------------------------------------------------------------
  
  observeEvent(input$run_opt, {
    req(build_model_data())
    rv$status <- "Building optimization instance"
    showNotification("Building optimization instance ...", type = "message")
    
    obj <- build_model_data()
    rv$data_obj   <- obj$data_obj
    rv$params_obj <- obj$params_obj
    
    rv$milp_instance <- build_optimization_instance(
      data   = rv$data_obj,
      params = rv$params_obj
    )
    
    showNotification("Solving LP model ...", type = "warning", duration = NULL)
    rv$status <- "Solving optimization model"
    
    rv$solve_result <- build_and_solve_afs_lp_v11(
       rv$milp_instance,
       gurobi_params = list(TimeLimit = 600),
       verbose = TRUE
     )
    
    rv$ext         <- prepare_solution_objects(rv$solve_result)
    rv$site_profit <- build_site_profit()
    rv$status      <- "Optimization complete"
    showNotification("Optimization finished.", type = "message")
  })
  
  # --------------------------------------------------------------------------
  # 5) KPI-OUTPUTS
  # --------------------------------------------------------------------------
  
  output$kpi_max_yield <- renderInfoBox({
    df <- build_growth_df()
    total_max <- max(df$total_dm, na.rm = TRUE)
    infoBox("Max biomass", paste0(round(total_max, 1), " t DM/ha"), icon = icon("chart-line"), color = "green")
  })
  
  output$kpi_stem_share_10 <- renderInfoBox({
    df <- build_growth_df() %>% filter(round(age, 1) == 10, component == "stem")
    val <- ifelse(nrow(df) > 0, 100 * df$share[1], NA)
    infoBox("Stem share at age 10", paste0(round(val, 1), "%"), icon = icon("tree"), color = "olive")
  })
  
  output$kpi_optimal_age <- renderInfoBox({
    df <- build_growth_df() %>% distinct(age, total_dm)
    best_age <- df$age[which.max(df$total_dm)]
    infoBox("Peak AGB age", paste0(round(best_age, 1), " yr"), icon = icon("bullseye"), color = "light-blue")
  })
  
  output$kpi_fresh_weight <- renderInfoBox({
    df <- build_growth_df() %>% distinct(age, total_dm)
    fw <- max(df$total_dm, na.rm = TRUE) / (1 - input$moisture)
    infoBox("Peak fresh weight", paste0(round(fw, 1), " t FW/ha"), icon = icon("tint"), color = "teal")
  })
  
  output$kpi_n_sites <- renderInfoBox({
    req(rv$site_profit)
    infoBox("Active sites", nrow(rv$site_profit), icon = icon("map-pin"), color = "green")
  })
  
  output$kpi_total_area <- renderInfoBox({
    req(rv$site_profit)
    infoBox("AFS area", paste0(round(sum(rv$site_profit$area_afs, na.rm = TRUE), 0), " ha"), icon = icon("draw-polygon"), color = "olive")
  })
  
  output$kpi_obj_val <- renderInfoBox({
    req(rv$solve_result)
    obj <- rv$solve_result$objective_value %||% rv$solve_result$objval %||% NA
    infoBox("Objective", comma(round(obj, 0)), icon = icon("euro-sign"), color = "light-blue")
  })
  
  output$kpi_solver_gap <- renderInfoBox({
    infoBox("Solver gap", "precomputed / n.a.", icon = icon("calculator"), color = "teal")
  })
  
  output$txt_solver_status <- renderText({ rv$status })
  output$txt_obj           <- renderText({ if (is.null(rv$solve_result)) "—" else format(rv$solve_result$objective_value, big.mark = ",") })
  output$txt_active_sites  <- renderText({ if (is.null(rv$site_profit)) "—" else nrow(rv$site_profit) })
  output$txt_total_area    <- renderText({ if (is.null(rv$site_profit)) "—" else round(sum(rv$site_profit$area_afs, na.rm = TRUE), 1) })
  output$txt_mean_profit   <- renderText({ if (is.null(rv$site_profit)) "—" else round(mean(rv$site_profit$profit_ha_yr, na.rm = TRUE), 1) })
  
  # --------------------------------------------------------------------------
  # 6) PLOTS — BIOMASS GROWTH
  # --------------------------------------------------------------------------
  
  output$plot_growth_stacked <- renderPlotly({
    df <- build_growth_df()
    
    browser()
    #text = paste0("Age: ", round(age,1), "<br>Component: ", component, "<br>Biomass: ", round(biomass_dm,2), " t DM/ha"))
    
    p <- ggplot(df, aes(x = age, y = df, fill = component)) +
      annotate("rect", xmin = input$min_age, xmax = input$max_age, ymin = 0, ymax = Inf, fill = "grey80", alpha = 0.25) +
      geom_area(color = "white", linewidth = 0.2, alpha = 0.75) +
      scale_fill_manual(values = c(stem = "#1b9e77", branch = "#d95f02", residue = "#7570b3")) +
      labs(x = "Stand age (years)", y = "Biomass (t DM/ha)", fill = "Fraction") +
      theme_minimal(base_size = 12)
    
    ggplotly(p, tooltip = "text")
  })
  
  output$plot_fraction_shares <- renderPlotly({
    df <- build_growth_df()
    p <- ggplot(df, aes(age, share, color = component,
                        text = paste0("Age: ", round(age,1),
                                      "<br>", component,
                                      "<br>Share: ", percent(share, accuracy = 0.1)))) +
      geom_line(linewidth = 1.1) +
      scale_y_continuous(labels = percent) +
      scale_color_manual(values = c(stem = "#1b9e77", branch = "#d95f02", residue = "#7570b3")) +
      labs(x = "Stand age (years)", y = "Share of total biomass", color = "Fraction") +
      theme_minimal(base_size = 12)
    ggplotly(p, tooltip = "text")
  })
  
  output$plot_growth_sensitivity <- renderPlotly({
    pars <- tidyr::crossing(k = c(input$k_gomp*0.8, input$k_gomp, input$k_gomp*1.2),
                            t0 = c(input$t0_gomp-1, input$t0_gomp, input$t0_gomp+1))
    
    df <- purrr::pmap_dfr(pars, function(k, t0) {
      build_scenario_ts(ages = seq(0, 20, 0.25), N = input$N_trees, C_site = input$C_site,
                        k = k, t0 = t0, label = paste0("k=", round(k,3), ", t0=", round(t0,1))) %>%
        mutate(total_dm = `Merch. stem` + `Merch. branch` + Residue,
               scenario = label)
    })
    
    p <- ggplot(df, aes(age, total_dm, color = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Age", y = "Total biomass (t DM/ha)", color = "Scenario") +
      theme_minimal(base_size = 11)
    
    ggplotly(p)
  })
  
  output$table_yields <- DT::renderDataTable({
    build_yields_input() %>%
      tidyr::pivot_wider(names_from = age, values_from = yield_ha) %>%
      arrange(product)
  })
  

}
