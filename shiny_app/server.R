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
library(Matrix)
library(slam) 
library(ROI)

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
      source("!helper_extract_site_profit.R")
      source("build_lp_rcpp_wrapper.R")
      
      # Falls C++-Solver dynamisch genutzt werden soll:
      Rcpp::sourceCpp("build_and_solve_afs_lp_v12_highs.cpp")
      
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
    
    tmp_long <- tmp %>%
      pivot_longer(
        cols      = c(`Merch. stem`, `Merch. branch`, Residue),
        names_to  = "component",
        values_to = "biomass_dm"
      ) %>%
      # Umbenennung NACH pivot_longer — jetzt sicher auf character-Werte
      mutate(component = recode(component,
                                "Merch. stem"   = "stem",
                                "Merch. branch" = "branch",
                                "Residue"       = "residue"
      )) %>%
      group_by(age) %>%
      mutate(
        total_dm = sum(biomass_dm, na.rm = TRUE),
        share    = ifelse(total_dm > 0, biomass_dm / total_dm, 0)
      ) %>%
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
        C_est  = input$C_est ,
        C_harv = input$C_harv ,
        C_main = input$C_main,
        C_opp  = input$opp_mean * c.opp.vec
      )
    
    consumers_model <- rv$consumers %>%
      mutate(
        P1 = P1 * input$rev_mult_P1,
        P2 = P2 * input$rev_mult_P2,
        P3 = P3 * input$rev_mult_P3
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
      c_tr_raw  = input$c_tr_raw ,
      c_tr_pre  = input$c_tr_pre 
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
        cost.log = input$c_tr_raw,
        cost.est = input$C_est,
        revenue  = mean(input$rev_mult_P1,input$rev_mult_P2,input$rev_mult_P3)
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
    
    #browser()
    
    rv$solve_result <- build_and_solve_afs_lp_v12(
       rv$milp_instance,
       highs_params = list(time_limit = 600, log_to_console = TRUE, presolve = "on"),
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

    p <- ggplot(df, aes(x = age, y = biomass_dm, fill = component)) +
      geom_area(color = "white", linewidth = 0.2, alpha = 0.75) +
      scale_fill_manual(values = c(stem = "#1b9e77", branch = "#d95f02", residue = "#7570b3")) +
      labs(x = "Stand age (years)", y = "Biomass (t DM/ha)", fill = "Fraction") +
      theme_minimal(base_size = 12)
    
    ggplotly(p)   # ← Fehler 3 behoben (text jetzt in aes)
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
  
  output$plot_growth_asymptotes <- renderPlotly({
    
    #browser()
    
    N_seq <- seq(20, 400, by = 1)
    df_a  <- bind_rows(
      tibble(N = N_seq,
             A = A_stand(N_seq, C_site = 4.000, beta=-.3),
             scenario = "advanteguous"),
      tibble(N = N_seq,
             A = A_stand(N_seq, C_site = 6.500, beta = -.5),
             scenario = "disadvantageous")
    )
    
    df_tree  <- bind_rows(
      tibble(N = N_seq,
             A = A_tree(N_seq, C_site = 4.000, beta=-.3),
             scenario = "advanteguous"),
      tibble(N = N_seq,
             A = A_tree(N_seq, C_site = 6.500, beta = -.5),
             scenario = "disadvantageous")
    )
    
    
    p_a <- ggplot(df_a, aes(x = N, y = A, linetype = scenario)) +
      geom_line(colour = "grey20", linewidth = 0.7) +
      scale_linetype_manual(
        name   = "Scenario",
        values = c("advanteguous" = "solid", "disadvantageous" = "dashed")) +
      labs(x = "Density N (trees/ha)",
           y = "A(N) (t DM/ha)") +
      theme_minimal(base_size = 8) +
      theme(legend.position  = c(0.72, 0.82),
            legend.key.size  = unit(0.4, "cm"),
            legend.text      = element_text(size = 7),
            legend.title     = element_text(size = 7, face = "bold"),
            legend.background = element_rect(fill = "white", colour = NA),
            panel.grid.minor = element_blank())
    
    p_b <- ggplot(df_tree, aes(x = N, y = A, linetype = scenario)) +
      geom_line(colour = "grey20", linewidth = 0.7) +
      scale_linetype_manual(
        name   = "Scenario",
        values = c("advanteguous" = "solid", "disadvantageous" = "dashed")) +
      labs(x = "Density N (trees/ha)",
           y = "A_tree(N) (t DM per tree)") +
       theme_minimal(base_size = 8) +
      theme(legend.position  = c(0.72, 0.82),
            legend.key.size  = unit(0.4, "cm"),
            legend.text      = element_text(size = 7),
            legend.title     = element_text(size = 7, face = "bold"),
            legend.background = element_rect(fill = "white", colour = NA),
            panel.grid.minor = element_blank())
    # ── Assemble with patchwork ─────────────────────────────────────────────────
    # Beide zu Plotly konvertieren, dann untereinander kombinieren
    subplot(
      ggplotly(p_a, tooltip = "none"),
      ggplotly(p_b, tooltip = "none"),
      nrows       = 2,          # ← untereinander
      shareX      = TRUE,       # ← gemeinsame X-Achse (beide haben N)
      titleY      = TRUE,       # ← Y-Achsentitel behalten
      margin      = 0.08        # ← Abstand zwischen den Plots
    )%>%
      layout(
        # Zweite Legende (von p_b) ausblenden — sie hat showlegend = TRUE doppelt
        showlegend = TRUE,
        legend = list(
          orientation = "h",          # horizontal
          x           = 0,            # linksbündig
          y           = 1.08,         # über dem Plot
          xanchor     = "left",
          yanchor     = "bottom",
          tracegroupgap = 0,
          itemsizing  = "constant",
          font        = list(size = 11)
        )
      ) %>%
      # Doppelte Legendeneinträge entfernen: Traces von p_b auf showlegend=FALSE setzen
      style(showlegend = FALSE, traces = 3:4)   # p_b hat Trace-Index 3 und 4
  })
  
  
  output$plot_growth_fractions <- renderPlotly({
    
    # ── Panel (b): logistic merchantability q_p(t) ─────────────────────────────
    ages_c <- seq(1, 20, by = 0.1)
    df_c <- tibble(
      age                    = ages_c,
      `Stem (d > 15 cm)`    = f_stem_merch(ages_c),
      `Branch (d > 7 cm)`   = f_branch_merch(ages_c)
    ) |>
      pivot_longer(-age, names_to = "compartment", values_to = "q") |>
      mutate(compartment = factor(compartment,
                                  levels = c("Stem (d > 15 cm)", "Branch (d > 7 cm)")))
    
    t50_vals <- tibble(
      compartment = factor(c("Stem (d > 15 cm)", "Branch (d > 7 cm)"),
                           levels = c("Stem (d > 15 cm)", "Branch (d > 7 cm)")),
      t50 = c(8.19, 7.55)
    )
    
    merch_cols <- c("Stem (d > 15 cm)"  = "#1b9e77",
                    "Branch (d > 7 cm)" = "#d95f02")
    
    p_c <- ggplot(df_c, aes(x = age, y = q, colour = compartment)) +
      geom_vline(data = t50_vals,
                 aes(xintercept = t50, colour = compartment),
                 linetype = "dashed", linewidth = 0.45, show.legend = FALSE) +
      geom_hline(yintercept = 0.5,
                 linetype = "dotted", colour = "grey60", linewidth = 0.4) +
      geom_line(linewidth = 0.8) +
      annotate("text", x = 8.19 + 0.3, y = 0.08,
               label = "t[50,s] = 8.2 years",
               size = 2.4, hjust = 0, colour = "#1b9e77") +
      annotate("text", x = 7.55 - 0.3, y = 0.20,
               label = "t[50,6] = 7.6 years",
               size = 2.4, hjust = 1, colour = "#d95f02") +
      scale_colour_manual(name = "Compartment", values = merch_cols) +
      scale_x_continuous(breaks = seq(0, 20, 2)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         limits = c(0, max(df_c$q)*1.1)) +
      labs(x = "Stand age (yr)",
           y = "q(p,t)  [merchantable fraction]") +
      theme_minimal(base_size = 8) +
      theme(legend.position  = c(0.15, 0.75),
            legend.key.size  = unit(0.4, "cm"),
            legend.text      = element_text(size = 7),
            legend.title     = element_text(size = 7, face = "bold"),
            legend.background = element_rect(fill = "white", colour = NA),
            panel.grid.minor = element_blank())
    
    # ── Assemble with patchwork ─────────────────────────────────────────────────
    ggplotly(p_c)
  })
  
  # --------------------------------------------------------------------------
  # 7) MAP OUTPUT — Network & Map Tab
  # --------------------------------------------------------------------------
  
  output$map_network <- renderLeaflet({
    
    req(rv$afs_workspace, rv$sites_sf, rv$sites, rv$storages, rv$consumers)
    
    # ── Farb-Konstanten ──────────────────────────────────────────────────────
    COL_HUB <- "#c05000"
    COL_P1  <- "#6a0dad"
    COL_P2  <- "#08519c"
    COL_P3  <- "#a50026"
    
    # ── Bounding Box ─────────────────────────────────────────────────────────
    bbox_wgs84 <- st_as_sfc(
      st_bbox(c(xmin = 10.6, xmax = 13.2, ymin = 50.9, ymax = 52.8),
              crs = st_crs(4326))
    )

    # ── AFS-Sites: Cluster-Zuweisung + Opportunitätskosten ───────────────────
    cluster_assig <- rv$afs_workspace$site_cluster_assig %>%
      select(site_id, hac_cluster) %>%
      mutate(hac_cluster = as.factor(hac_cluster))
    
    sites_sf_clust <- rv$sites_sf %>%
      mutate(site_id = seq_len(nrow(rv$sites_sf))) %>%
      left_join(cluster_assig, by = "site_id") %>%
      st_transform(4326) %>%
      st_intersection(bbox_wgs84)
    
    # Opportunitätskosten: aus milp_instance (falls vorhanden), sonst rv$sites
    # NEU — robust gegen fehlende C_opp:
    sites_for_opp <- if (!is.null(rv$milp_instance) &&
                         "C_opp" %in% names(rv$milp_instance$sites)) {
      # Nach Optimierung: aus milp_instance
      rv$milp_instance$sites %>%
        select(site_id, C_opp) %>%
        rename(hac_cluster = site_id) %>%
        mutate(hac_cluster = as.factor(hac_cluster))
      
    } else if ("C_opp" %in% names(rv$sites)) {
      # rv$sites hat zufällig schon C_opp
      rv$sites %>%
        select(site_id, C_opp) %>%
        rename(hac_cluster = site_id) %>%
        mutate(hac_cluster = as.factor(hac_cluster))
      
    } else {
      # Vor Optimierung: C_opp mit input$opp_mean synthetisch erzeugen
      set.seed(123578)
      c.opp.vec <- rnorm(nrow(rv$sites), 1, 0.3)
      rv$sites %>%
        mutate(C_opp = input$opp_mean * c.opp.vec) %>%
        select(site_id, C_opp) %>%
        rename(hac_cluster = site_id) %>%
        mutate(hac_cluster = as.factor(hac_cluster))
    }
    
    sites_sf_clust <- sites_sf_clust %>%
      left_join(sites_for_opp, by = "hac_cluster")
    
    # ── Aktive Sites markieren (nur wenn Lösung vorhanden) ───────────────────
    if (!is.null(rv$ext)) {
      active_ids <- unique(rv$ext$Xij$site_id)
      sites_sf_clust <- sites_sf_clust %>%
        mutate(active = hac_cluster %in% as.factor(active_ids))
    } else {
      sites_sf_clust <- sites_sf_clust %>%
        mutate(active = FALSE)
    }
    
    # Farbpalette Opportunitätskosten
    pal_opp <- colorNumeric(
      palette  = c("#1a9641", "#ffffbf", "#d7191c"),
      domain   = sites_sf_clust$C_opp,
      na.color = "transparent"
    )
    
    # Inaktive Sites nach Optimierung in Grau dimmen
    sites_sf_clust <- sites_sf_clust %>%
      mutate(
        fill_color   = if_else(active | !any(active), pal_opp(C_opp), "#bdbdbd"),
        fill_opacity = if_else(active | !any(active), 0.85, 0.35)
      )
    
    # ── Hubs / Storages ───────────────────────────────────────────────────────
    stor_sf <- st_as_sf(
      rv$storages %>% mutate(ptsize = 8),
      coords = c("lng", "lat"), crs = 4326
    ) %>% st_intersection(bbox_wgs84)
    
    # ── Consumers ─────────────────────────────────────────────────────────────
    consumers_plt <- rv$consumers %>%
      mutate(
        total_demand = demand_P1 + demand_P2 + demand_P3,
        kategorie = case_when(
          demand_P1 >= demand_P2 & demand_P1 >= demand_P3 & demand_P1 > 0 ~ "Chemical / Pulp (P1)",
          demand_P2 >= demand_P1 & demand_P2 >= demand_P3 & demand_P2 > 0 ~ "Pulp / Paper (P2)",
          demand_P3 > 0 ~ "Energy / Biogas (P3)",
          TRUE          ~ "Other"
        ),
        ptsize = pmin(14, pmax(5, 5 + log1p(total_demand) / 1.8))
      )
    
    cons_sf <- st_as_sf(
      consumers_plt, coords = c("lng", "lat"), crs = 4326
    ) %>% st_intersection(bbox_wgs84)
    
    pal_cons <- colorFactor(
      palette = c(
        "Chemical / Pulp (P1)" = COL_P1,
        "Pulp / Paper (P2)"    = COL_P2,
        "Energy / Biogas (P3)" = COL_P3,
        "Other"                = "grey60"
      ),
      domain = cons_sf$kategorie
    )
    
    # ── Leaflet aufbauen ──────────────────────────────────────────────────────
    leaflet(options = leafletOptions(zoomControl = TRUE), width = "100%") %>%
      
      #addProviderTiles(providers$CartoDB.Positron, group = "Basemap") %>%
      
      addProviderTiles(providers$Esri.WorldGrayCanvas) %>% 
      
      addPolygons(
        data         = sites_sf_clust,
        stroke       = TRUE, color = "grey30", weight = 0.5,
        fillColor    = ~fill_color,
        fillOpacity  = ~fill_opacity,
        popup        = ~paste0(
          "<strong>Cluster-ID:</strong> ", hac_cluster,
          "<br><strong>Opp. Kosten:</strong> ", round(C_opp, 1), " €/ha",
          "<br><strong>Status:</strong> ",
          if_else(active, "&#10003; Aktiv", "&#8212; Inaktiv")
        ),
        group        = "AFS Sites"
      ) %>%
      
      addCircleMarkers(
        data        = stor_sf,
        radius      = ~ptsize,
        color       = COL_HUB, stroke = TRUE, weight = 2,
        fillColor   = COL_HUB, fillOpacity = 0.95,
        popup       = ~"<strong>Hub / Storage</strong>",
        group       = "Hubs"
      ) %>%
      
      addCircleMarkers(
        data        = cons_sf,
        radius      = ~ptsize,
        color       = "white", weight = 1,
        fillColor   = ~pal_cons(kategorie), fillOpacity = 0.9,
        popup       = ~paste0(
          "<strong>Consumer-Typ:</strong> ", kategorie,
          "<br><strong>Gesamtnachfrage:</strong> ",
          round(total_demand, 1), " kt"
        ),
        group       = "Consumers"
      ) %>%
      
      addLegend(
        position = "bottomright",
        pal      = pal_opp,
        values   = sites_sf_clust$C_opp,
        title    = "Opp. Kosten (€/ha)",
        opacity  = 0.8
      ) %>%
      
      addLegend(
        position = "topright",
        pal      = pal_cons,
        values   = cons_sf$kategorie,
        title    = "Consumer-Typ",
        opacity  = 0.95
      ) %>%
      
      addLayersControl(
        overlayGroups = c("AFS Sites", "Hubs", "Consumers"),
        options       = layersControlOptions(collapsed = FALSE)
      ) %>%
      
      fitBounds(lng1 = 10.6, lat1 = 50.9, lng2 = 13.2, lat2 = 52.8)
  })
  
  # In server.R
  output$plot_biomass_time <- renderPlotly({
    
    req(rv$ext)
    #browser()
    prod_labels <- c("1" = "Stem (P1)", "2" = "Branches (P2)", "3" = "Residues (P3)")
    prod_colors <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3")
    
    # Deliveries aggregated by period × del_product
    delivery_ts <- rv$ext$Xjk %>%
      group_by(period, product = factor(del_product)) %>%
      summarise(volume = sum(value), .groups = "drop")
    
    # Raw biomass shipped (Xij) aggregated by period — shown as dashed reference
    shipped_ts <- rv$ext$Xij %>%
      group_by(period) %>%
      summarise(shipped = sum(value), .groups = "drop")
    
    
    # ── Total demand per period ─────────────────────────────────────────────────
    demand_max <- rv$milp_instance$demand %>%
      group_by(product = factor(product, levels = c("3", "2", "1"))) %>%
      summarise(max_demand = sum(D_max)/rv$milp_instance$n_periods, .groups = "drop") %>% 
      mutate(max_demand = cumsum(max_demand))
    
    
    p <- ggplot(delivery_ts, aes(x = period, y = volume, fill = product)) +
      geom_area(position = "stack", alpha = 0.80, colour = "white", linewidth = 0.2) +
      # Raw shipped biomass reference (dashed grey)
      geom_line(data = shipped_ts,
                aes(x = period, y = shipped),
                inherit.aes = FALSE,
                colour = "#636363", linetype = "dashed", linewidth = 0.7) +
      # Max consumer demand per product (dotted, colour-matched to product)
      geom_hline(data = demand_max,
                 aes(yintercept = max_demand, colour = product),
                 linetype = "dashed", linewidth = 1.0, inherit.aes = FALSE) +
      scale_fill_manual(
        name   = "Product",
        values = prod_colors,
        labels = prod_labels
      ) +
      scale_colour_manual(
        name   = "Max. demand",
        values = prod_colors,
        labels = prod_labels
      ) +
      scale_x_continuous(name = "Planning period (year)",
                         breaks = seq(1, rv$milp_instance$n_periods, by = 2)) +
      scale_y_continuous(name = "Biomass quantity (t fresh biomass)",
                         labels = scales::comma) +
      annotate("text", x = max(shipped_ts$period) - 1,
               y = max(shipped_ts$shipped) * 0.95,
               label = "Raw shipped\n(site-hub)", size = 2.5,
               colour = "#636363", hjust = 1) +
      theme_minimal(base_size = 9) +
      theme(legend.position = "right",
            panel.grid.minor = element_blank())
    
    # ── Konvertierung zu Plotly ─────────────────────────────────────────────────
    ggplotly(p, tooltip = c("x", "y", "fill", "colour")) %>%
      layout(
        legend = list(orientation = "h", x = 0, y = -0.2),
        hovermode = "x unified"
      )
  })
  
  
  # --------------------------------------------------------------------------
  # SANKEY: Biomassflüsse Sites → Hubs → Consumers nach Produktart
  # In server.R, ergänzen nach dem plot_biomass_time-Block
  # --------------------------------------------------------------------------
  
  output$plot_sankey <- renderPlotly({
    req(rv$ext, rv$milp_instance, rv$storages, rv$consumers)
    
    
    
    prod_labels <- c("1" = "Stem (P1)", "2" = "Branches (P2)", "3" = "Residues (P3)")
    prod_colors <- c("1" = "#1b9e77",   "2" = "#d95f02",        "3" = "#7570b3")
    
    # ── 1) Labels für Knoten vorbereiten ─────────────────────────────────────
    # Sites (Cluster-IDs aus milp_instance)
    site_ids  <- sort(unique(rv$ext$Xij$site_id))
    site_lbls <- paste0("Site ", site_ids)
    
    # Produkt-Knoten (P1/P2/P3 als Zwischenstufe zwischen Site und Hub)
    prod_ids  <- c("1", "2", "3")
    prod_lbls <- prod_labels[prod_ids]
    
    # Hubs
    hub_ids   <- sort(unique(rv$storages$storage_id))
    hub_lbls  <- paste0("Hub ", hub_ids)
    
    # Consumers
    cons_ids  <- sort(unique(rv$consumers$consumer_id))
    cons_lbls <- ifelse(
      !is.na(rv$consumers$name[match(cons_ids, rv$consumers$consumer_id)]),
      rv$consumers$name[match(cons_ids, rv$consumers$consumer_id)],
      paste0("Consumer ", cons_ids)
    )
    
    # ── Knoten-Index-Tabelle (0-basiert für Plotly) ───────────────────────────
    nodes <- data.frame(
      id    = c(paste0("site_",  site_ids),
                paste0("prod_",  prod_ids),
                paste0("hub_",   hub_ids),
                paste0("cons_",  cons_ids)),
      label = c(site_lbls, prod_lbls, hub_lbls, cons_lbls),
      stringsAsFactors = FALSE
    )
    nodes$idx <- seq_len(nrow(nodes)) - 1L   # 0-basiert
    
    node_idx <- function(key) nodes$idx[match(key, nodes$id)]
    
    #browser()
    
    # ── 2) Kanten: Sites → Produkt-Knoten (aggregiert über alle Perioden) ──
    xij_agg <- rv$ext$Xij %>%
      rename(product = p) %>% 
      mutate(product = as.character(product)) %>%
      group_by(site_id, product) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      filter(value > 0) %>%
      mutate(
        source = node_idx(paste0("site_", site_id)),
        target = node_idx(paste0("prod_", product)),
        color  = prod_colors[product]
      )
    
    # ── 3) Kanten: Produkt-Knoten → Hubs ─────────────────────────────────────
    xij_hub_agg <- rv$ext$Xij %>%
      rename(product = p) %>% 
      mutate(product = as.character(product)) %>%
      group_by(product, hub_id) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      filter(value > 0) %>%
      mutate(
        source = node_idx(paste0("prod_", product)),
        target = node_idx(paste0("hub_",  hub_id)),
        color  = prod_colors[product]
      )
    
    # ── 4) Kanten: Hubs → Consumers (nach del_product) ───────────────────────
    xjk_agg <- rv$ext$Xjk %>%
      mutate(del_product = as.character(del_product)) %>%
      group_by(hub_id, consumer_id, del_product) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      filter(value > 0) %>%
      mutate(
        source = node_idx(paste0("hub_",  hub_id)),
        target = node_idx(paste0("cons_", consumer_id)),
        color  = prod_colors[del_product]
      )
    
    # ── 5) Alle Kanten zusammenführen ─────────────────────────────────────────
    all_links <- bind_rows(
      xij_agg     %>% select(source, target, value, color),
      xij_hub_agg %>% select(source, target, value, color),
      xjk_agg     %>% select(source, target, value, color)
    ) %>%
      mutate(
        value     = round(value / 1000, 2),   # t → kt
        link_color = paste0(
          gsub("rgb", "rgba",
               col2rgb(color) %>%
                 apply(2, function(x) paste0("rgba(", paste(x, collapse=","), ",0.45)"))
          )
        )
      )
    
    # Hilfsfunktion: hex → rgba-String (Plotly-kompatibel)
    hex2rgba <- function(hex, alpha = 0.45) {
      rgb_vals <- col2rgb(hex)
      paste0("rgba(", rgb_vals[1], ",", rgb_vals[2], ",", rgb_vals[3], ",", alpha, ")")
    }
    
    all_links <- all_links %>%
      mutate(link_color = sapply(color, hex2rgba, alpha = 0.45))
    
    # ── 6) Knotenfarben ───────────────────────────────────────────────────────
    node_color <- case_when(
      startsWith(nodes$id, "site_") ~ "rgba(100,160,100,0.85)",
      nodes$id == "prod_1"          ~ "rgba(27,158,119,0.9)",
      nodes$id == "prod_2"          ~ "rgba(217,95,2,0.9)",
      nodes$id == "prod_3"          ~ "rgba(117,112,179,0.9)",
      startsWith(nodes$id, "hub_")  ~ "rgba(192,80,0,0.85)",
      startsWith(nodes$id, "cons_") ~ "rgba(50,80,160,0.85)",
      TRUE ~ "rgba(150,150,150,0.7)"
    )
    
    # ── 7) Plotly Sankey ──────────────────────────────────────────────────────
    plot_ly(
      type = "sankey",
      orientation = "h",
      arrangement = "snap",
      
      node = list(
        label = nodes$label,
        color = node_color,
        pad   = 15,
        thickness = 20,
        line  = list(color = "white", width = 0.5)
      ),
      
      link = list(
        source = all_links$source,
        target = all_links$target,
        value  = all_links$value,
        color  = all_links$link_color,
        label  = paste0(
          round(all_links$value, 1), " kt",
          " [", nodes$label[all_links$source + 1L],
          " → ", nodes$label[all_links$target + 1L], "]"
        )
      )
    ) %>%
      layout(
        title = list(
          text = "Biomass Flows: Sites → Products → Hubs → Consumers (kt total)",
          font = list(size = 13)
        ),
        font = list(size = 10),
        margin = list(l = 5, r = 5, t = 40, b = 5)
      )
  
    })
  
}
