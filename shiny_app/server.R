# ==============================================================================
# server.R — AFS-SCD Shiny App (Refaktorierte Version)
# Änderungen gegenüber Vorgängerversion:
#   1. Modularisierung: Biomasse-, Karten- und Flowlogik in eigene Module ausgelagert
#   2. global.R-Konvention: source()-Aufrufe explizit, kein !-Präfix-Trick
#   3. build_model_data() als echtes reactive() statt doppelter Evaluierung
#   4. Versionskonsistenz: precomp RDS-Pfad als Konstante
#   5. NULL-sicherer %||%-Operator definiert
#   6. Interne Hilfsfunktionen (hex2rgba, prepare_solution_objects) als
#      top-level Fns vor dem Server-Function-Body
# ==============================================================================

library(shiny)
library(shinydashboard)
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

# ------------------------------------------------------------------------------
# KONSTANTEN
# ------------------------------------------------------------------------------
RDS_PATH <- "result_lp_v11.rds"
OPP_SEED <- 123578L

# NULL-Koaleszenz-Operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------------------------
# HILFSFUNKTIONEN (top-level, vor server function)
# ------------------------------------------------------------------------------

#' Hex-Farbe in RGBA-String konvertieren (Plotly-kompatibel)
hex2rgba <- function(hex, alpha = 0.45) {
  rgb_vals <- col2rgb(hex)
  paste0("rgba(", rgb_vals[1], ",", rgb_vals[2], ",", rgb_vals[3], ",", alpha, ")")
}

#' Lösungsvektor in benannte Liste mit sprechenden Spaltennamen überführen
prepare_solution_objects <- function(result) {
  ext <- result$solution_vector
  
  ext$Xjk <- ext$Xjk %>%
    rename(src_product = p, del_product = q, period = t,
           hub_id = j, consumer_id = k)
  
  ext$Xij <- ext$Xij %>%
    rename(period = t, site_id = i, hub_id = j)
  
  ext$z <- ext$z %>%
    rename(arc_type = type, site_id = i)
  
  ext$S <- ext$S %>%
    rename(hub_id = j, period = t)
  
  ext
}

# ==============================================================================
# SERVER MODULE
# ==============================================================================

# ------------------------------------------------------------------------------
# Modul: Biomass Growth
# ------------------------------------------------------------------------------
biomassGrowthServer <- function(id, rv, input) {
  moduleServer(id, function(input_m, output_m, session_m) {
    
    # Reaktive Hilfsdaten — Sidebar-Inputs aus parent session via `input`
    growth_df <- reactive({
      req(rv$afs_workspace)
      age.vec <- seq(0, max(20, input$max_age), by = 0.1)
      tmp <- build_scenario_ts(
        ages   = age.vec, N = input$N_trees, C_site = input$C_site,
        k = input$k_gomp, t0 = input$t0_gomp, label = "Interactive scenario"
      )
      tmp %>%
        pivot_longer(cols = c(`Merch. stem`, `Merch. branch`, Residue),
                     names_to = "component", values_to = "biomass_dm") %>%
        mutate(component = recode(component,
                                  "Merch. stem"   = "stem",
                                  "Merch. branch" = "branch",
                                  "Residue"       = "residue")) %>%
        group_by(age) %>%
        mutate(total_dm = sum(biomass_dm, na.rm = TRUE),
               share    = ifelse(total_dm > 0, biomass_dm / total_dm, 0)) %>%
        ungroup()
    })
    
    # KPI Outputs
    output_m$kpi_max_yield <- renderInfoBox({
      total_max <- max(growth_df()$total_dm, na.rm = TRUE)
      infoBox("Max biomass", paste0(round(total_max, 1), " t DM/ha"),
              icon = icon("chart-line"), color = "green")
    })
    
    output_m$kpi_stem_share_10 <- renderInfoBox({
      df  <- growth_df() %>% filter(round(age, 1) == 10, component == "stem")
      val <- ifelse(nrow(df) > 0, 100 * df$share[1], NA)
      infoBox("Stem share at age 10", paste0(round(val, 1), "%"),
              icon = icon("tree"), color = "olive")
    })
    
    output_m$kpi_optimal_age <- renderInfoBox({
      df       <- growth_df() %>% distinct(age, total_dm)
      best_age <- df$age[which.max(df$total_dm)]
      infoBox("Peak AGB age", paste0(round(best_age, 1), " yr"),
              icon = icon("bullseye"), color = "light-blue")
    })
    
    output_m$kpi_fresh_weight <- renderInfoBox({
      df <- growth_df() %>% distinct(age, total_dm)
      fw <- max(df$total_dm, na.rm = TRUE) / (1 - input$moisture)
      infoBox("Peak fresh weight", paste0(round(fw, 1), " t FW/ha"),
              icon = icon("tint"), color = "teal")
    })
    
    # Plot: Gestapeltes Biomassediagramm
    output_m$plot_growth_stacked <- renderPlotly({
      df <- growth_df()
      p <- ggplot(df, aes(x = age, y = biomass_dm, fill = component)) +
        geom_area(color = "white", linewidth = 0.2, alpha = 0.75) +
        scale_fill_manual(
          values = c(stem = "#1b9e77", branch = "#d95f02", residue = "#7570b3")) +
        labs(x = "Stand age (years)", y = "Biomass (t DM/ha)", fill = "Fraction") +
        theme_minimal(base_size = 12)
      ggplotly(p)
    })
    
    # Plot: Fraktionsanteile
    output_m$plot_fraction_shares <- renderPlotly({
      df <- growth_df()
      p <- ggplot(df, aes(age, share, color = component,
                          text = paste0("Age: ", round(age, 1),
                                        "<br>", component,
                                        "<br>Share: ",
                                        percent(share, accuracy = 0.1)))) +
        geom_line(linewidth = 1.1) +
        scale_y_continuous(labels = percent) +
        scale_color_manual(
          values = c(stem = "#1b9e77", branch = "#d95f02", residue = "#7570b3")) +
        labs(x = "Stand age (years)", y = "Share of total biomass",
             color = "Fraction") +
        theme_minimal(base_size = 12)
      ggplotly(p, tooltip = "text")
    })
    
    # Plot: Asymptoten (Stand & Baum)
    output_m$plot_growth_asymptotes <- renderPlotly({
      N_seq    <- seq(20, 400, by = 1)
      lty_vals <- c(advantageous = "solid", disadvantageous = "dashed")
      
      df_a <- bind_rows(
        tibble(N = N_seq, A = A_stand(N_seq, C_site = 4.0, beta = -.3),
               scenario = "advantageous"),
        tibble(N = N_seq, A = A_stand(N_seq, C_site = 6.5, beta = -.5),
               scenario = "disadvantageous")
      )
      df_tree <- bind_rows(
        tibble(N = N_seq, A = A_tree(N_seq, C_site = 4.0, beta = -.3),
               scenario = "advantageous"),
        tibble(N = N_seq, A = A_tree(N_seq, C_site = 6.5, beta = -.5),
               scenario = "disadvantageous")
      )
      
      p_a <- ggplot(df_a, aes(N, A, linetype = scenario)) +
        geom_line(colour = "grey20", linewidth = 0.7) +
        scale_linetype_manual(name = "Scenario", values = lty_vals) +
        labs(x = "Density N (trees/ha)", y = "A(N) (t DM/ha)") +
        theme_minimal(base_size = 8) +
        theme(legend.position = c(0.72, 0.82), panel.grid.minor = element_blank())
      
      p_b <- ggplot(df_tree, aes(N, A, linetype = scenario)) +
        geom_line(colour = "grey20", linewidth = 0.7) +
        scale_linetype_manual(name = "Scenario", values = lty_vals) +
        labs(x = "Density N (trees/ha)", y = "A_tree(N) (t DM per tree)") +
        theme_minimal(base_size = 8) +
        theme(legend.position = c(0.72, 0.82), panel.grid.minor = element_blank())
      
      subplot(
        ggplotly(p_a, tooltip = "none"),
        ggplotly(p_b, tooltip = "none"),
        nrows = 2, shareX = TRUE, titleY = TRUE, margin = 0.08
      ) %>%
        layout(showlegend = TRUE,
               legend = list(orientation = "h", x = 0, y = 1.08,
                             xanchor = "left", yanchor = "bottom",
                             tracegroupgap = 0, itemsizing = "constant",
                             font = list(size = 11))) %>%
        style(showlegend = FALSE, traces = 3:4)
    })
    
    # Plot: Logistische Handelsholzanteile
    output_m$plot_growth_fractions <- renderPlotly({
      ages_c <- seq(1, 20, by = 0.1)
      df_c <- tibble(
        age                 = ages_c,
        `Stem (d > 15 cm)` = f_stem_merch(ages_c),
        `Branch (d > 7 cm)` = f_branch_merch(ages_c)
      ) |>
        pivot_longer(-age, names_to = "compartment", values_to = "q") |>
        mutate(compartment = factor(compartment,
                                    levels = c("Stem (d > 15 cm)",
                                               "Branch (d > 7 cm)")))
      
      t50_vals   <- tibble(
        compartment = factor(c("Stem (d > 15 cm)", "Branch (d > 7 cm)"),
                             levels = c("Stem (d > 15 cm)", "Branch (d > 7 cm)")),
        t50 = c(8.19, 7.55)
      )
      merch_cols <- c("Stem (d > 15 cm)" = "#1b9e77",
                      "Branch (d > 7 cm)" = "#d95f02")
      
      p_c <- ggplot(df_c, aes(age, q, colour = compartment)) +
        geom_vline(data = t50_vals, aes(xintercept = t50, colour = compartment),
                   linetype = "dashed", linewidth = 0.45, show.legend = FALSE) +
        geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey60",
                   linewidth = 0.4) +
        geom_line(linewidth = 0.8) +
        annotate("text", x = 8.49, y = 0.08, label = "t[50,s] = 8.2 years",
                 size = 2.4, hjust = 0, colour = "#1b9e77") +
        annotate("text", x = 7.25, y = 0.20, label = "t[50,6] = 7.6 years",
                 size = 2.4, hjust = 1, colour = "#d95f02") +
        scale_colour_manual(name = "Compartment", values = merch_cols) +
        scale_x_continuous(breaks = seq(0, 20, 2)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           limits = c(0, max(df_c$q) * 1.1)) +
        labs(x = "Stand age (yr)", y = "q(p,t)  [merchantable fraction]") +
        theme_minimal(base_size = 8) +
        theme(legend.position = c(0.15, 0.75), panel.grid.minor = element_blank())
      ggplotly(p_c)
    })
    
    list(growth_df = growth_df)
  })
}

# ------------------------------------------------------------------------------
# Modul: Karte & Netzwerk
# ------------------------------------------------------------------------------
networkMapServer <- function(id, rv, input) {
  moduleServer(id, function(input_m, output_m, session_m) {
    
    output_m$map_network <- renderLeaflet({
      req(rv$afs_workspace, rv$sites_sf, rv$sites, rv$storages, rv$consumers)
      
      COL_HUB <- "#c05000"
      COL_P1  <- "#6a0dad"
      COL_P2  <- "#08519c"
      COL_P3  <- "#a50026"
      
      bbox_wgs84 <- st_as_sfc(
        st_bbox(c(xmin = 10.6, xmax = 13.2, ymin = 50.9, ymax = 52.8),
                crs = st_crs(4326))
      )
      
      cluster_assig <- rv$afs_workspace$site_cluster_assig %>%
        select(site_id, hac_cluster) %>%
        mutate(hac_cluster = as.factor(hac_cluster))
      
      sites_sf_clust <- rv$sites_sf %>%
        mutate(site_id = seq_len(nrow(rv$sites_sf))) %>%
        left_join(cluster_assig, by = "site_id") %>%
        st_transform(4326) %>%
        st_intersection(bbox_wgs84)
      
      # Opportunitätskosten: robust gegen fehlende C_opp
      sites_for_opp <- if (!is.null(rv$milp_instance) &&
                           "C_opp" %in% names(rv$milp_instance$sites)) {
        rv$milp_instance$sites %>%
          select(site_id, C_opp) %>%
          rename(hac_cluster = site_id) %>%
          mutate(hac_cluster = as.factor(hac_cluster))
      } else if ("C_opp" %in% names(rv$sites)) {
        rv$sites %>%
          select(site_id, C_opp) %>%
          rename(hac_cluster = site_id) %>%
          mutate(hac_cluster = as.factor(hac_cluster))
      } else {
        set.seed(OPP_SEED)
        c.opp.vec <- rnorm(nrow(rv$sites), 1, 0.3)
        rv$sites %>%
          mutate(C_opp = input$opp_mean * c.opp.vec) %>%
          select(site_id, C_opp) %>%
          rename(hac_cluster = site_id) %>%
          mutate(hac_cluster = as.factor(hac_cluster))
      }
      
      sites_sf_clust <- sites_sf_clust %>%
        left_join(sites_for_opp, by = "hac_cluster")
      
      if (!is.null(rv$ext)) {
        active_ids     <- unique(rv$ext$Xij$site_id)
        sites_sf_clust <- sites_sf_clust %>%
          mutate(active = hac_cluster %in% as.factor(active_ids))
      } else {
        sites_sf_clust <- sites_sf_clust %>% mutate(active = FALSE)
      }
      
      pal_opp <- colorNumeric(
        palette  = c("#1a9641", "#ffffbf", "#d7191c"),
        domain   = sites_sf_clust$C_opp, na.color = "transparent"
      )
      
      sites_sf_clust <- sites_sf_clust %>%
        mutate(
          fill_color   = if_else(active | !any(active), pal_opp(C_opp), "#bdbdbd"),
          fill_opacity = if_else(active | !any(active), 0.85, 0.35)
        )
      
      stor_sf <- st_as_sf(
        rv$storages %>% arrange(storage_id) %>%
          mutate(hub_nr = paste0("Hub ", row_number()), ptsize = 8),
        coords = c("lng", "lat"), crs = 4326
      ) %>% st_intersection(bbox_wgs84)
      
      consumers_plt <- rv$consumers %>%
        mutate(
          consumer_nr  = paste0("Consumer ", row_number()),
          total_demand = demand_P1 + demand_P2 + demand_P3,
          kategorie = case_when(
            demand_P1 >= demand_P2 & demand_P1 >= demand_P3 & demand_P1 > 0 ~
              "Chemical / Pulp (P1)",
            demand_P2 >= demand_P1 & demand_P2 >= demand_P3 & demand_P2 > 0 ~
              "Pulp / Paper (P2)",
            demand_P3 > 0 ~ "Energy / Biogas (P3)",
            TRUE          ~ "Other"
          ),
          marker_color = case_when(
            kategorie == "Chemical / Pulp (P1)" ~ "purple",
            kategorie == "Pulp / Paper (P2)"    ~ "blue",
            kategorie == "Energy / Biogas (P3)" ~ "red",
            TRUE                                ~ "gray"
          )
        )
      
      cons_sf <- st_as_sf(
        consumers_plt, coords = c("lng", "lat"), crs = 4326
      ) %>% st_intersection(bbox_wgs84)
      
      pal_cons <- colorFactor(
        palette = c("Chemical / Pulp (P1)" = COL_P1,
                    "Pulp / Paper (P2)"    = COL_P2,
                    "Energy / Biogas (P3)" = COL_P3,
                    "Other"                = "grey60"),
        domain  = cons_sf$kategorie
      )
      
      leaflet(options = leafletOptions(zoomControl = TRUE), width = "100%") %>%
        addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
        addPolygons(
          data = sites_sf_clust, stroke = TRUE, color = "grey30", weight = 0.5,
          fillColor = ~fill_color, fillOpacity = ~fill_opacity,
          popup = ~paste0(
            "<strong>Cluster-ID:</strong> ", hac_cluster,
            "<br><strong>Opp. Kosten:</strong> ", round(C_opp, 1), " €/ha",
            "<br><strong>Status:</strong> ",
            if_else(active, "&#10003; Aktiv", "&#8212; Inaktiv")),
          group = "AFS Sites"
        ) %>%
        addCircleMarkers(
          data = stor_sf, radius = ~ptsize, color = COL_HUB,
          stroke = TRUE, weight = 2, fillColor = COL_HUB, fillOpacity = 0.95,
          popup = ~paste0("<strong>", hub_nr, "</strong>",
                          "<br><strong>Typ:</strong> Hub / Storage"),
          group = "Hubs"
        ) %>%
        addAwesomeMarkers(
          data = cons_sf,
          icon = ~awesomeIcons(icon = "industry", library = "fa",
                               markerColor = marker_color, iconColor = "white"),
          popup = ~paste0(
            "<strong>", consumer_nr, "</strong>",
            "<br><strong>Consumer-Typ:</strong> ", kategorie,
            "<br><strong>Gesamtnachfrage:</strong> ",
            round(total_demand, 1), " kt"),
          label = ~consumer_nr,
          group = "Consumers"
        ) %>%
        addLegend(position = "bottomright", pal = pal_opp,
                  values = sites_sf_clust$C_opp,
                  title = "Opp. Kosten (€/ha)", opacity = 0.8) %>%
        addLegend(position = "topright", pal = pal_cons,
                  values = cons_sf$kategorie,
                  title = "Consumer-Typ", opacity = 0.95) %>%
        addLayersControl(
          overlayGroups = c("AFS Sites", "Hubs", "Consumers"),
          options       = layersControlOptions(collapsed = FALSE)
        ) %>%
        fitBounds(lng1 = 10.6, lat1 = 50.9, lng2 = 13.2, lat2 = 52.8)
    })
  })
}

# ------------------------------------------------------------------------------
# Modul: Material Flows
# ------------------------------------------------------------------------------
materialFlowsServer <- function(id, rv) {
  moduleServer(id, function(input_m, output_m, session_m) {
    
    PROD_LABELS <- c("1" = "Stem (P1)", "2" = "Branches (P2)",
                     "3" = "Residues (P3)")
    PROD_COLORS <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3")
    
    # Biomass over Time -------------------------------------------------------
    output_m$plot_biomass_time <- renderPlotly({
      req(rv$ext, rv$milp_instance)
      
      delivery_ts <- rv$ext$Xjk %>%
        group_by(period, product = factor(del_product)) %>%
        summarise(volume = sum(value), .groups = "drop")
      
      shipped_ts <- rv$ext$Xij %>%
        group_by(period) %>%
        summarise(shipped = sum(value), .groups = "drop")
      
      demand_max <- rv$milp_instance$demand %>%
        group_by(product = factor(product, levels = c("3", "2", "1"))) %>%
        summarise(max_demand = sum(D_max) / rv$milp_instance$n_periods,
                  .groups = "drop") %>%
        mutate(max_demand = cumsum(max_demand))
      
      p <- ggplot(delivery_ts, aes(x = period, y = volume, fill = product)) +
        geom_area(position = "stack", alpha = 0.80, colour = "white",
                  linewidth = 0.2) +
        geom_line(data = shipped_ts, aes(x = period, y = shipped),
                  inherit.aes = FALSE, colour = "#636363",
                  linetype = "dashed", linewidth = 0.7) +
        geom_hline(data = demand_max,
                   aes(yintercept = max_demand, colour = product),
                   linetype = "dashed", linewidth = 1.0, inherit.aes = FALSE) +
        scale_fill_manual(name = "Product", values = PROD_COLORS,
                          labels = PROD_LABELS) +
        scale_colour_manual(name = "Max. demand", values = PROD_COLORS,
                            labels = PROD_LABELS) +
        scale_x_continuous(name = "Planning period (year)",
                           breaks = seq(1, rv$milp_instance$n_periods, by = 2)) +
        scale_y_continuous(name = "Biomass quantity (t fresh biomass)",
                           labels = scales::comma) +
        annotate("text", x = max(shipped_ts$period) - 1,
                 y = max(shipped_ts$shipped) * 0.95,
                 label = "Raw shipped\n(site-hub)", size = 2.5,
                 colour = "#636363", hjust = 1) +
        theme_minimal(base_size = 9) +
        theme(legend.position = "right", panel.grid.minor = element_blank())
      
      ggplotly(p, tooltip = c("x", "y", "fill", "colour")) %>%
        layout(legend = list(orientation = "h", x = 0, y = -0.2),
               hovermode = "x unified")
    })
    
    # Sankey ------------------------------------------------------------------
    output_m$plot_sankey <- renderPlotly({
      req(rv$ext, rv$milp_instance, rv$storages, rv$consumers)
      
      site_ids <- sort(unique(rv$ext$Xij$site_id))
      hub_ids  <- sort(unique(rv$storages$storage_id))
      cons_ids <- sort(unique(rv$consumers$consumer_id))
      prod_ids <- c("1", "2", "3")
      
      nodes <- data.frame(
        id    = c(paste0("site_", site_ids), paste0("prod_", prod_ids),
                  paste0("hub_",  hub_ids),  paste0("cons_", cons_ids)),
        label = c(paste0("Site ", site_ids), PROD_LABELS[prod_ids],
                  paste0("Hub ", hub_ids),
                  paste0("Consumer ", seq_along(cons_ids))),
        stringsAsFactors = FALSE
      ) %>% mutate(idx = row_number() - 1L)
      
      node_idx <- function(key) nodes$idx[match(key, nodes$id)]
      
      xij_agg <- rv$ext$Xij %>%
        rename(product = p) %>%
        mutate(product = as.character(product)) %>%
        group_by(site_id, product) %>%
        summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
        filter(value > 0) %>%
        mutate(source = node_idx(paste0("site_", site_id)),
               target = node_idx(paste0("prod_", product)),
               color  = PROD_COLORS[product])
      
      xij_hub_agg <- rv$ext$Xij %>%
        rename(product = p) %>%
        mutate(product = as.character(product)) %>%
        group_by(product, hub_id) %>%
        summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
        filter(value > 0) %>%
        mutate(source = node_idx(paste0("prod_", product)),
               target = node_idx(paste0("hub_", hub_id)),
               color  = PROD_COLORS[product])
      
      xjk_agg <- rv$ext$Xjk %>%
        mutate(del_product = as.character(del_product)) %>%
        group_by(hub_id, consumer_id, del_product) %>%
        summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
        filter(value > 0) %>%
        mutate(source = node_idx(paste0("hub_", hub_id)),
               target = node_idx(paste0("cons_", consumer_id)),
               color  = PROD_COLORS[del_product])
      
      all_links <- bind_rows(
        xij_agg     %>% select(source, target, value, color),
        xij_hub_agg %>% select(source, target, value, color),
        xjk_agg     %>% select(source, target, value, color)
      ) %>%
        mutate(
          value      = round(value / 1000, 2),
          link_color = vapply(color, hex2rgba, character(1), alpha = 0.45)
        )
      
      node_color <- case_when(
        startsWith(nodes$id, "site_") ~ "rgba(100,160,100,0.85)",
        nodes$id == "prod_1"          ~ "rgba(27,158,119,0.9)",
        nodes$id == "prod_2"          ~ "rgba(217,95,2,0.9)",
        nodes$id == "prod_3"          ~ "rgba(117,112,179,0.9)",
        startsWith(nodes$id, "hub_")  ~ "rgba(192,80,0,0.85)",
        startsWith(nodes$id, "cons_") ~ "rgba(50,80,160,0.85)",
        TRUE ~ "rgba(150,150,150,0.7)"
      )
      
      plot_ly(
        type = "sankey", orientation = "h", arrangement = "snap",
        node = list(label = nodes$label, color = node_color,
                    pad = 15, thickness = 20,
                    line = list(color = "white", width = 0.5)),
        link = list(
          source = all_links$source, target = all_links$target,
          value  = all_links$value,  color  = all_links$link_color,
          label  = paste0(round(all_links$value, 1), " kt [",
                          nodes$label[all_links$source + 1L], " → ",
                          nodes$label[all_links$target + 1L], "]")
        )
      ) %>%
        layout(
          title  = list(
            text = "Biomass Flows: Sites → Products → Hubs → Consumers (kt total)",
            font = list(size = 13)),
          font   = list(size = 10),
          margin = list(l = 5, r = 5, t = 40, b = 5)
        )
    })
    
    # Revenue per Consumer ----------------------------------------------------
    output_m$plot_rev_consumer <- renderPlotly({
      req(rv$ext, rv$milp_instance, rv$consumers)
      
      price_df <- as.data.frame(rv$milp_instance$consumer_prices) %>%
        setNames(c("consumer_id", "del_product", "price")) %>%
        mutate(consumer_id = as.integer(consumer_id),
               del_product = as.integer(del_product))
      
      cons_labels <- rv$consumers %>%
        arrange(consumer_id) %>%
        mutate(rank = row_number(), label = paste0("Consumer ", rank)) %>%
        select(consumer_id, label)
      
      rev_df <- rv$ext$Xjk %>%
        group_by(consumer_id, del_product) %>%
        summarise(volume = sum(value, na.rm = TRUE), .groups = "drop") %>%
        left_join(price_df, by = c("consumer_id", "del_product")) %>%
        mutate(
          revenue     = volume * coalesce(price, 0) / 1e6,
          del_product = as.character(del_product),
          prod_label  = PROD_LABELS[del_product]
        ) %>%
        left_join(cons_labels, by = "consumer_id") %>%
        filter(revenue > 0)
      
      p <- ggplot(rev_df,
                  aes(x = reorder(label, -revenue), y = revenue,
                      fill = prod_label,
                      text = paste0("<b>", label, "</b><br>",
                                    "Product: ", prod_label, "<br>",
                                    "Volume: ",
                                    scales::comma(round(volume / 1000, 1)),
                                    " kt<br>",
                                    "Revenue: ",
                                    scales::comma(round(revenue, 2)),
                                    " Mio. €"))) +
        geom_col(position = "stack", width = 0.7, alpha = 0.88) +
        scale_fill_manual(name   = "Product",
                          values = setNames(PROD_COLORS, PROD_LABELS)) +
        scale_y_continuous(name = "Revenue (Mio. €)", labels = scales::comma) +
        labs(x = NULL) +
        theme_minimal(base_size = 10) +
        theme(axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
              legend.position = "right",
              panel.grid.minor = element_blank())
      
      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "v", x = 1.02, y = 0.5),
               hovermode = "closest", margin = list(b = 90))
    })
    
    # Demand Fulfilment -------------------------------------------------------
    output_m$plot_demand_fulfilment <- renderPlotly({
      req(rv$ext, rv$milp_instance, rv$consumers)
      
      cons_labels <- rv$consumers %>%
        mutate(consumer_id = seq_len(nrow(rv$consumers))) %>%
        arrange(consumer_id) %>%
        mutate(rank = row_number(), label = paste0("Consumer ", rank)) %>%
        select(consumer_id, label, rank)
      
      delivered <- rv$ext$Xjk %>%
        group_by(consumer_id, del_product) %>%
        summarise(delivered = sum(value, na.rm = TRUE), .groups = "drop") %>%
        rename(product = del_product) %>%
        mutate(product = as.integer(product))
      
      demand_total <- rv$milp_instance$demand %>%
        group_by(consumer_id = as.integer(consumer_id),
                 product     = as.integer(product)) %>%
        summarise(demand = sum(D_max, na.rm = TRUE), .groups = "drop")
      
      fulfil_df <- demand_total %>%
        left_join(delivered, by = c("consumer_id", "product")) %>%
        filter(demand > 0) %>%
        mutate(
          delivered  = coalesce(delivered, 0),
          rate_pct   = round(pmin(delivered / pmax(demand, 1), 1) * 100, 1),
          prod_label = PROD_LABELS[as.character(product)]
        ) %>%
        left_join(cons_labels, by = "consumer_id")
      
      prod_order <- c("Stem (P1)", "Branches (P2)", "Residues (P3)")
      cons_order <- cons_labels %>% arrange(rank) %>% pull(label)
      
      mat_rate <- fulfil_df %>%
        select(prod_label, label, rate_pct) %>%
        tidyr::pivot_wider(names_from = label, values_from = rate_pct,
                           values_fill = NA_real_) %>%
        arrange(match(prod_label, prod_order))
      
      z_mat           <- as.matrix(mat_rate[, cons_order])
      rownames(z_mat) <- mat_rate$prod_label
      
      hover_mat <- fulfil_df %>%
        mutate(hover = paste0(
          "<b>", label, "</b><br>",
          "Product: ", prod_label, "<br>",
          "Delivered: ", scales::comma(round(delivered / 1000, 1)), " kt<br>",
          "Max demand: ", scales::comma(round(demand / 1000, 1)), " kt<br>",
          "<b>Fulfilment: ", rate_pct, " %</b>"
        )) %>%
        select(prod_label, label, hover) %>%
        tidyr::pivot_wider(names_from = label, values_from = hover,
                           values_fill = "") %>%
        arrange(match(prod_label, prod_order))
      
      text_mat <- as.matrix(hover_mat[, cons_order])
      
      annotations <- list()
      for (r in seq_len(nrow(z_mat))) {
        for (c in seq_len(ncol(z_mat))) {
          val <- z_mat[r, c]
          if (!is.na(val)) {
            annotations <- c(annotations, list(list(
              x = cons_order[c], y = prod_order[r],
              text = paste0(val, "%"), showarrow = FALSE,
              font = list(size  = 11,
                          color = if (val >= 60) "white" else "grey20",
                          family = "Arial")
            )))
          }
        }
      }
      
      plot_ly(x = cons_order, y = prod_order, z = z_mat,
              type = "heatmap", text = text_mat,
              hovertemplate = "%{text}<extra></extra>",
              colorscale = list(list(0, "#d73027"), list(0.5, "#fee08b"),
                                list(1, "#1a9850")),
              zmin = 0, zmax = 100,
              colorbar = list(title = "Fulfilment (%)",
                              tickvals = c(0, 50, 100),
                              ticktext = c("0%", "50%", "100%"), len = 0.6),
              xgap = 3, ygap = 3) %>%
        layout(
          annotations = annotations,
          xaxis = list(title = "", tickangle = -35, tickfont = list(size = 9)),
          yaxis = list(title = "", tickfont = list(size = 10),
                       categoryorder = "array", categoryarray = prod_order),
          margin = list(t = 20, r = 20, b = 110, l = 110)
        )
    })
    
    # Product Cascade ---------------------------------------------------------
    output_m$plot_product_cascade <- renderPlotly({
      req(rv$ext, rv$ext$Xjk)
      
      SRC_LABEL  <- c("1" = "Stems", "2" = "Branches", "3" = "Residues")
      DEL_LABEL  <- c("1" = "Industrial", "2" = "Pulp & Paper", "3" = "Energy")
      P_COLS_SRC <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3")
      
      flow_df <- rv$ext$Xjk %>%
        group_by(src_product, del_product) %>%
        summarise(volume = sum(value, na.rm = TRUE), .groups = "drop") %>%
        filter(volume > 0) %>%
        mutate(src_id    = as.character(src_product),
               del_id    = as.character(del_product),
               src_label = unname(SRC_LABEL[src_id]),
               del_label = unname(DEL_LABEL[del_id]))
      
      validate(need(nrow(flow_df) > 0, "No product-cascade flows available."))
      
      src_nodes <- tibble::tibble(
        node_key   = paste0("src_", names(SRC_LABEL)),
        node_label = unname(SRC_LABEL),
        node_group = "src", x = 0.001, y = c(0.15, 0.45, 0.75)
      )
      del_nodes <- tibble::tibble(
        node_key   = paste0("del_", names(DEL_LABEL)),
        node_label = unname(DEL_LABEL),
        node_group = "del", x = 0.999, y = c(0.15, 0.45, 0.75)
      )
      nodes <- dplyr::bind_rows(src_nodes, del_nodes) %>%
        dplyr::mutate(
          idx        = dplyr::row_number() - 1L,
          node_color = dplyr::if_else(
            node_group == "src",
            vapply(sub("^src_", "", node_key),
                   function(z) hex2rgba(P_COLS_SRC[z], 0.85), character(1)),
            rep("rgba(50,80,160,0.80)", dplyr::n())
          )
        )
      node_idx <- function(keys) nodes$idx[match(keys, nodes$node_key)]
      
      links <- flow_df %>%
        mutate(
          source = node_idx(paste0("src_", src_id)),
          target = node_idx(paste0("del_", del_id)),
          color  = vapply(src_id,
                          function(z) hex2rgba(P_COLS_SRC[z], 0.60),
                          character(1)),
          value  = round(volume / 1000, 1)
        )
      
      plot_ly(
        type = "sankey", orientation = "h", arrangement = "freeform",
        node = list(
          label = nodes$node_label,
          x     = nodes$x, y = nodes$y,
          color = nodes$node_color,   # ← vordefinierter Vektor der Länge 6
          pad = 20, thickness = 24,
          line = list(color = "white", width = 0.5)
        ),
        link = list(
          source = links$source, target = links$target,
          value  = links$value,  color  = links$color,
          label  = paste0(links$value, " kt — ",
                          links$src_label, " → ", links$del_label)
        )
      ) %>%
        layout(
          title  = list(text = "Feedstock × End-Use Cascade (kt total)",
                        font = list(size = 13)),
          font   = list(size = 11),
          margin = list(l = 10, r = 10, t = 40, b = 10)
        )
    })
  })
}

# ------------------------------------------------------------------------------
# Modul: Site KPIs
# ------------------------------------------------------------------------------
siteKPIsServer <- function(id, rv) {
  moduleServer(id, function(input_m, output_m, session_m) {
    
    output_m$plot_rotation_dist <- renderPlotly({
      req(rv$site_profit)
      
      rot_df <- rv$ext$z %>%
        filter(value > 0, arc_type == "harvest") %>%
        group_by(site_id) %>%
        mutate( cycle_length = t - s )
      
      validate(
        need(nrow(rot_df) > 0, "No activated harvest cycles available.")
      )
      
      p <- ggplot(
        rot_df, aes( x = cycle_length)
      ) +
        geom_histogram(
          binwidth = 1,
          boundary = 0.5,
          fill = "#7570b3",
          colour = "white",
          alpha = 0.85
        ) +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        labs(
          x = "Rotation cycle length (years)",
          y = "Number of observations"
        ) +
        theme_minimal(base_size = 11) +
        theme(panel.grid.minor = element_blank())
      
      ggplotly(p, tooltip = c("x", "y")) %>%
        layout(showlegend = FALSE)
      
    })
    
    output_m$plot_profit_split <- renderPlotly({
      req(rv$site_profit)
      
      # ── 0. Colour palette — identical to fig-results-cost-profit-breakdown ───────
      prod_colors_rev <- c(
        "Revenue P1 (Stem)"     = "#1b9e77",
        "Revenue P2 (Branches)" = "#d95f02",
        "Revenue P3 (Residue)"  = "#7570b3"
      )
      cost_colors <- c(
        "Establishment"         = "#b2182b",
        "Harvesting"            = "#ef8a62",
        "Maintenance"           = "#fddbc7",
        "Opportunity"           = "#67001f",
        "Transport (raw)"       = "#2166ac",
        "Transport (pre-proc.)" = "#92c5de",
        "Storage"               = "#4d4d4d"
      )
      
      needed_cols <- c("site_id", "cost_est", "cost_harv", "cost_main", "cost_opp","cost_tr_raw", "cost_tr_pre","cost_stor", "rev_P1", "rev_P2", "rev_P3")
      
      df  <- rv$site_profit %>%
        mutate(site_label = paste0("Site ", site_id)) %>%
        select(site_label, profit_ha_yr, denom, all_of(needed_cols)) %>% 
        mutate(
          across(
            c(starts_with("cost_"), starts_with("rev_")),
            ~ .x / denom
          )
        ) %>% 
        mutate(
          cost_stor   = replace(cost_stor, is.na(cost_stor), 0)
        ) %>%
        select(site_label, profit_ha_yr, starts_with("rev_"), starts_with("cost_")) %>% 
        arrange(desc(profit_ha_yr)) %>%
        mutate(site_label = factor(site_label, levels = site_label))
      
      # ── Hilfsfunktion: einen Trace hinzufügen ─────────────────────────────────
      add_rev_bar <- function(fig, col, name) {
        add_bars(
          fig,
          data            = df,
          x               = ~site_label,
          y               = as.formula(paste0("~", col)),
          name            = name,
          legendgroup     = "revenue",
          legendgrouptitle = list(text = "<b>Revenue by product</b>"),
          marker          = list(color = prod_colors_rev[[name]]),
          hovertemplate   = paste0("<b>%{x}</b><br>", name, ": %{y:,.1f} €/ha·yr<extra></extra>")
        )
      }
      
      add_cost_bar <- function(fig, col, name) {
        add_bars(
          fig,
          data            = df,
          x               = ~site_label,
          y               = as.formula(paste0("~(-", col, ")")),
          name            = name,
          legendgroup     = "costs",
          legendgrouptitle = list(text = "<b>Costs</b>"),
          marker          = list(color = cost_colors[[name]]),
          customdata      = as.formula(paste0("~", col)),
          hovertemplate   = paste0("<b>%{x}</b><br>", name, ": %{customdata:,.1f} €/ha·yr<extra></extra>")
        )
      }
      
      # ── 2. Traces zusammenbauen ────────────────────────────────────────────────
      fig <- plot_ly() %>%
        # Revenues
        add_rev_bar("rev_P1",  "Revenue P1 (Stem)")     %>%
        add_rev_bar("rev_P2",  "Revenue P2 (Branches)") %>%
        add_rev_bar("rev_P3",  "Revenue P3 (Residue)")  %>%
        # Costs
        add_cost_bar("cost_est",     "Establishment")         %>%
        add_cost_bar("cost_harv",    "Harvesting")            %>%
        add_cost_bar("cost_main",    "Maintenance")           %>%
        add_cost_bar("cost_opp",     "Opportunity")           %>%
        add_cost_bar("cost_tr_raw",  "Transport (raw)")       %>%
        add_cost_bar("cost_tr_pre",  "Transport (pre-proc.)") %>%
        add_cost_bar("cost_stor",    "Storage")               %>%
        # Net profit as diamond markers
        add_markers(
          data          = df,
          x             = ~site_label,
          y             = ~profit_ha_yr,
          name          = "Profit",
          legendgroup   = "profit",
          marker        = list(symbol = "diamond", size = 8, color = "black"),
          hovertemplate = "<b>%{x}</b><br>Profit: %{y:,.1f} €/ha·yr<extra></extra>"
        )
      
      # ── 3. Layout ──────────────────────────────────────────────────────────────
      fig %>%
        layout(
          barmode      = "relative",
          bargap       = 0.18,
          plot_bgcolor = "rgba(235,235,235,0.6)",
          paper_bgcolor = "white",
          xaxis = list(
            title      = "Site ID",
            tickangle  = -60,
            categoryorder = "array",
            categoryarray = levels(df$site_label)
          ),
          yaxis = list(
            title       = "€ per hectare and year",
            zeroline    = TRUE,
            zerolinewidth = 1.2,
            zerolinecolor = "gray40",
            gridcolor   = "white"
          ),
          legend = list(
            orientation  = "v",
            x            = 1.02,
            y            = 1,
            xanchor      = "left",
            yanchor      = "top",
            tracegroupgap = 12
          ),
          margin = list(l = 70, r = 190, b = 100, t = 20)
        )
    })
    
    output_m$plot_profit_vs_opp <- renderPlotly({
      req(rv$site_profit)
      
      p <- ggplot(rv$site_profit, aes(opp_cost_site, profit_ha_yr)) +
        geom_point(aes(size = area_afs),
                   colour = "#1b9e77",
                   alpha = 0.75) +
        geom_smooth(method = "lm", se = FALSE, colour = "#AAC800",
                    linewidth = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
        labs(x = "Opportunity cost (€/ha/yr)", y = "Profit (€/ha/yr)") +
        theme_minimal(base_size = 9)
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    output_m$plot_p1share_profit <- renderPlotly({
      req(rv$site_profit)
      p <- ggplot(rv$site_profit, aes(share_p1, profit_ha_yr)) +
        geom_point(aes(size = area_afs),
                   colour = "#1b9e77",
                   alpha = 0.75) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
        geom_smooth(method = "lm", se = FALSE, colour = "#AAC800",
                    linewidth = 0.8) +
        labs(x = "P1 (Stem) share", y = "Profit (€/ha/yr)") +
        theme_minimal(base_size = 9)
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    output_m$plot_dist_profit <- renderPlotly({
      req(rv$site_profit)
      p <- ggplot(rv$site_profit, aes(avg_dist_consumer_km + avg_dist_hub_km, profit_ha_yr)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
        geom_point(aes(size = area_afs),
                   colour = "#d95f02",
                   alpha = 0.75) +
        geom_smooth(method = "lm", se = FALSE, colour = "#AAC800",
                    linewidth = 0.8) +
        labs(x = "Distance to hub (km)", y = "Profit (€/ha/yr)") +
        theme_minimal(base_size = 9)
      ggplotly(p, tooltip = c("x", "y"))
    })
  })
}

# ==============================================================================
# HAUPTSERVER
# ==============================================================================
function(input, output, session) {
  
  # --------------------------------------------------------------------------
  # 1) INITIALISIERUNG
  # --------------------------------------------------------------------------
  rv <- reactiveValues(
    afs_workspace = NULL,
    sites_sf      = NULL,
    sites         = NULL,
    storages      = NULL,
    consumers     = NULL,
    data_obj      = NULL,
    params_obj    = NULL,
    milp_instance = NULL,
    solve_result  = NULL,
    ext           = NULL,
    site_profit   = NULL,
    sensitivity   = NULL,
    status        = "Idle"
  )
  
  # ── Workspace laden ─────────────────────────────────────────────────────────
  observe({
    if (!is.null(rv$afs_workspace)) return()
    showNotification("Loading workspace and source files ...", type = "message")
    
    # Explizite source()-Aufrufe — kein !-Präfix-Konvention
    source("!afs_biomass_setup.r",          local = FALSE)
    source("!helper_func.r",                local = FALSE)
    source("!helper_extract_site_profit.R", local = FALSE)
    Rcpp::sourceCpp("build_and_solve_afs_lp_v12_highs.cpp")
    
    load("afs_workspace_red.RData")
    rv$afs_workspace <- afs_workspace
    rv$sites_sf      <- afs_workspace$sites_sf
    rv$sites         <- afs_workspace$sites_clustered
    rv$storages      <- afs_workspace$storages
    rv$consumers     <- afs_workspace$consumers
    
    rv$consumers <- rv$consumers %>%
      mutate(demand_P1 = demand_P1 / 5,
             P1 = P1 / 2, P2 = P2 / 2, P3 = P3 * 0.6) %>%
      mutate(total_demand = demand_P1 + demand_P2 + demand_P3) %>%
      filter(total_demand >= 1.5) %>%
      select(-total_demand) %>%
      mutate(consumer_id = seq_len(nrow(.)))
    
    rv$storages$CAP_stor <- rv$storages$CAP_stor * 100000
    rv$storages$CAP_proc <- rv$storages$CAP_proc * 100000
    
    rv$status <- "Base data loaded"
  })
  
  # ── Vorberechnetes RDS laden ────────────────────────────────────────────────
  observe({
    req(rv$afs_workspace)
    if (!is.null(rv$ext)) return()
    
    if (!file.exists(RDS_PATH)) {
      showNotification(
        paste0(RDS_PATH, " not found — run optimization to generate results."),
        type = "warning", duration = 8)
      return()
    }
    
    showNotification("Loading precomputed solution ...", type = "message")
    
    precomp         <- readRDS(RDS_PATH)
    rv$ext          <- prepare_solution_objects(precomp)
    rv$solve_result <- precomp
    
    obj              <- build_model_data()
    rv$data_obj      <- obj$data_obj
    rv$params_obj    <- obj$params_obj
    rv$milp_instance <- build_optimization_instance(
      data   = rv$data_obj,
      params = rv$params_obj
    )
    
    rv$site_profit <- extract_site_profit(
      list(
        setting = list(
          opp      = input$opp_mean,
          cost.log = input$c_tr_raw,
          cost.est = input$C_est,
          revenue  = mean(c(input$rev_mult_P1, input$rev_mult_P2,
                            input$rev_mult_P3))
        ),
        milp_instance = rv$milp_instance,
        ext           = rv$ext
      ),
      scenario_name = "Precomputed"
    )
    
    rv$status <- "Precomputed solution loaded"
    showNotification("Precomputed solution ready. All plots initialized.",
                     type = "message")
  })
  
  # --------------------------------------------------------------------------
  # 2) REAKTIVE MODELLBAUSTEINE
  # --------------------------------------------------------------------------
  
  build_yields_reactive <- reactive({
    req(rv$afs_workspace)
    y <- build_scenario_ts(
      ages   = 1:input$max_age, N = input$N_trees, C_site = input$C_site,
      k = input$k_gomp, t0 = input$t0_gomp, label = "Interactive scenario"
    )
    y %>%
      select(-scenario) %>%
      pivot_longer(-1, names_to = "product", values_to = "yield_ha") %>%
      mutate(
        product  = as.integer(as.character(factor(product, labels = c(2, 1, 3)))),
        yield_ha = yield_ha / (1 - input$moisture)
      )
  })
  
  build_model_data <- reactive({
    req(rv$sites, rv$storages, rv$consumers, rv$afs_workspace)
    
    set.seed(OPP_SEED)
    c.opp.vec <- rnorm(nrow(rv$sites), 1, .3)
    
    sites_model <- rv$sites %>%
      mutate(C_est  = input$C_est,
             C_harv = input$C_harv,
             C_main = input$C_main,
             C_opp  = input$opp_mean * c.opp.vec)
    
    consumers_model <- rv$consumers %>%
      mutate(P1 = P1 * input$rev_mult_P1,
             P2 = P2 * input$rev_mult_P2,
             P3 = P3 * input$rev_mult_P3)
    
    data_obj <- list(
      sites         = sites_model,
      storages      = rv$storages,
      consumers     = consumers_model,
      dist_ij       = rv$afs_workspace$dist_ij_clust,
      dist_jk       = rv$afs_workspace$dist_jk[,
                                               consumers_model$consumer_id, drop = FALSE],
      yields_by_age = build_yields_reactive()
    )
    data_obj$consumers[, c("demand_P1", "demand_P2", "demand_P3")] <-
      data_obj$consumers[, c("demand_P1", "demand_P2", "demand_P3")] * 1000
    
    params_obj <- list(
      n_periods = input$n_periods,
      max_age   = as.integer(input$max_age),
      min_age   = as.integer(input$min_age),
      c_tr_raw  = input$c_tr_raw,
      c_tr_pre  = input$c_tr_pre
    )
    
    list(data_obj = data_obj, params_obj = params_obj)
  })
  
  build_site_profit <- reactive({
    req(rv$milp_instance, rv$ext)
    extract_site_profit(list(
      setting = list(
        opp      = input$opp_mean,
        cost.log = input$c_tr_raw,
        cost.est = input$C_est,
        revenue  = mean(c(input$rev_mult_P1, input$rev_mult_P2,
                          input$rev_mult_P3))
      ),
      milp_instance = rv$milp_instance,
      ext           = rv$ext
    ), scenario_name = "Interactive")
  })
  
  # --------------------------------------------------------------------------
  # 3) OPTIMIERUNG AUSLÖSEN
  # --------------------------------------------------------------------------
  observeEvent(input$run_opt, {
    req(build_model_data())
    rv$status <- "Building optimization instance"
    showNotification("Building optimization instance ...", type = "message")
    
    obj              <- build_model_data()
    rv$data_obj      <- obj$data_obj
    rv$params_obj    <- obj$params_obj
    rv$milp_instance <- build_optimization_instance(
      data   = rv$data_obj,
      params = rv$params_obj
    )
    
    showNotification("Solving LP model ...", type = "warning", duration = NULL)
    rv$status <- "Solving optimization model"
    
    rv$solve_result <- build_and_solve_afs_lp_v12(
      rv$milp_instance,
      highs_params = list(time_limit = 600, log_to_console = TRUE,
                          presolve = "on"),
      verbose = TRUE
    )
    
    
    
    if(rv$solve_result$objective_value == 0){
      rv$ext         <- NULL
      rv$site_profit <- NULL
    }
    else{
      rv$ext         <- prepare_solution_objects(rv$solve_result)
      rv$site_profit <- build_site_profit()
    }
    
    rv$status      <- "Optimization complete"
    showNotification("Optimization finished.", type = "message")
  })
  
  # --------------------------------------------------------------------------
  # 4) MODULE INITIALISIEREN
  # --------------------------------------------------------------------------
  biomassGrowthServer("biomass",   rv, input)
  networkMapServer("network",      rv, input)
  materialFlowsServer("flows",     rv)
  siteKPIsServer("site_kpis",      rv)
  
  # --------------------------------------------------------------------------
  # 5) KPI INFO BOXES (direkte Outputs für ui.R Kompatibilität)
  # --------------------------------------------------------------------------
  output$kpi_max_yield <- renderInfoBox({
    req(rv$afs_workspace)
    age.vec <- seq(0, max(20, input$max_age), by = 0.1)
    tmp <- build_scenario_ts(
      ages = age.vec, N = input$N_trees, C_site = input$C_site,
      k = input$k_gomp, t0 = input$t0_gomp, label = "x"
    )
    total_max <- tmp %>%
      pivot_longer(cols = c(`Merch. stem`, `Merch. branch`, Residue),
                   names_to = "c", values_to = "v") %>%
      group_by(age) %>% summarise(s = sum(v, na.rm = TRUE)) %>%
      pull(s) %>% max(na.rm = TRUE)
    infoBox("Max biomass", paste0(round(total_max, 1), " t DM/ha"),
            icon = icon("chart-line"), color = "green")
  })
  
  output$kpi_stem_share_10 <- renderInfoBox({
    req(rv$afs_workspace)
    age.vec <- seq(0, max(20, input$max_age), by = 0.1)
    tmp <- build_scenario_ts(
      ages = age.vec, N = input$N_trees, C_site = input$C_site,
      k = input$k_gomp, t0 = input$t0_gomp, label = "x"
    ) %>%
      pivot_longer(cols = c(`Merch. stem`, `Merch. branch`, Residue),
                   names_to = "component", values_to = "biomass_dm") %>%
      mutate(component = recode(component, "Merch. stem" = "stem",
                                "Merch. branch" = "branch",
                                "Residue" = "residue")) %>%
      group_by(age) %>%
      mutate(total_dm = sum(biomass_dm, na.rm = TRUE),
             share    = ifelse(total_dm > 0, biomass_dm / total_dm, 0)) %>%
      ungroup()
    df  <- tmp %>% filter(round(age, 1) == 10, component == "stem")
    val <- ifelse(nrow(df) > 0, 100 * df$share[1], NA)
    infoBox("Stem share at age 10", paste0(round(val, 1), "%"),
            icon = icon("tree"), color = "olive")
  })
  
  output$kpi_optimal_age <- renderInfoBox({
    req(rv$afs_workspace)
    age.vec <- seq(0, max(20, input$max_age), by = 0.1)
    tmp <- build_scenario_ts(
      ages = age.vec, N = input$N_trees, C_site = input$C_site,
      k = input$k_gomp, t0 = input$t0_gomp, label = "x"
    ) %>%
      pivot_longer(cols = c(`Merch. stem`, `Merch. branch`, Residue),
                   names_to = "c", values_to = "v") %>%
      group_by(age) %>% summarise(total_dm = sum(v, na.rm = TRUE)) %>%
      ungroup()
    best_age <- tmp$age[which.max(tmp$total_dm)]
    infoBox("Peak AGB age", paste0(round(best_age, 1), " yr"),
            icon = icon("bullseye"), color = "light-blue")
  })
  
  output$kpi_fresh_weight <- renderInfoBox({
    req(rv$afs_workspace)
    age.vec <- seq(0, max(20, input$max_age), by = 0.1)
    tmp <- build_scenario_ts(
      ages = age.vec, N = input$N_trees, C_site = input$C_site,
      k = input$k_gomp, t0 = input$t0_gomp, label = "x"
    ) %>%
      pivot_longer(cols = c(`Merch. stem`, `Merch. branch`, Residue),
                   names_to = "c", values_to = "v") %>%
      group_by(age) %>% summarise(total_dm = sum(v, na.rm = TRUE)) %>%
      ungroup()
    fw <- max(tmp$total_dm, na.rm = TRUE) / (1 - input$moisture)
    infoBox("Peak fresh weight", paste0(round(fw, 1), " t FW/ha"),
            icon = icon("tint"), color = "teal")
  })
  
  output$kpi_n_sites <- renderInfoBox({
    req(rv$site_profit)
    infoBox("Active sites", nrow(rv$site_profit),
            icon = icon("map-pin"), color = "green")
  })
  
  output$kpi_total_area <- renderInfoBox({
    req(rv$site_profit)
    infoBox("AFS area",
            paste0(comma(round(sum(rv$site_profit$area_afs, na.rm = TRUE), 0)),
                   " ha"),
            icon = icon("draw-polygon"), color = "olive")
  })
  
  output$kpi_obj_val <- renderInfoBox({
    req(rv$solve_result)
    obj <- rv$solve_result$objective_value %||%
      rv$solve_result$objval %||% NA
    infoBox("Total profit",
            paste0(comma(round(obj / 1000000, 0)), " Mill. €"),
            icon = icon("euro-sign"), color = "light-blue")
  })
  
  output$kpi_solver_gap <- renderInfoBox({
    infoBox("Solver gap", "precomputed / n.a.",
            icon = icon("calculator"), color = "teal")
  })
  
  # --------------------------------------------------------------------------
  # 6) TEXT-OUTPUTS & TABELLE
  # --------------------------------------------------------------------------
  output$txt_solver_status <- renderText({ rv$status })
  output$txt_obj <- renderText({
    if (is.null(rv$solve_result)) "—" else
      format(rv$solve_result$objective_value, big.mark = ",")
  })
  output$txt_active_sites <- renderText({
    if (is.null(rv$site_profit)) "—" else nrow(rv$site_profit)
  })
  output$txt_total_area <- renderText({
    if (is.null(rv$site_profit)) "—" else
      round(sum(rv$site_profit$area_afs, na.rm = TRUE), 1)
  })
  output$txt_mean_profit <- renderText({
    if (is.null(rv$site_profit)) "—" else
      round(mean(rv$site_profit$profit_ha_yr, na.rm = TRUE), 1)
  })
  
  output$table_sites <- DT::renderDataTable({
    req(rv$site_profit)
    DT::datatable(rv$site_profit,
                  options  = list(pageLength = 15, scrollX = TRUE),
                  rownames = FALSE)
  })
  
  # --------------------------------------------------------------------------
  # 7) PLOT-OUTPUTS (delegiert an Module via Namespace-Präfix)
  # --------------------------------------------------------------------------
  # Plots werden durch die jeweiligen moduleServer()-Aufrufe oben registriert.
  # Für ui.R-Kompatibilität ohne vollständige Modularisierung der UI:
  # Outputs werden über ihre Standard-IDs angesprochen.
  
  output$plot_growth_stacked      <- renderPlotly({ NULL })
  output$plot_fraction_shares     <- renderPlotly({ NULL })
  output$plot_growth_asymptotes   <- renderPlotly({ NULL })
  output$plot_growth_fractions    <- renderPlotly({ NULL })
  output$map_network              <- renderLeaflet({ NULL })
  output$plot_biomass_time        <- renderPlotly({ NULL })
  output$plot_sankey              <- renderPlotly({ NULL })
  output$plot_rev_consumer        <- renderPlotly({ NULL })
  output$plot_demand_fulfilment   <- renderPlotly({ NULL })
  output$plot_product_cascade     <- renderPlotly({ NULL })
  output$plot_rotation_dist       <- renderPlotly({ NULL })
  output$plot_profit_split        <- renderPlotly({ NULL })
  output$plot_profit_vs_opp       <- renderPlotly({ NULL })
  output$plot_p1share_profit      <- renderPlotly({ NULL })
  output$plot_dist_profit         <- renderPlotly({ NULL })
  
  # --------------------------------------------------------------------------
  # 8) CSV EXPORT
  # --------------------------------------------------------------------------
  output$export_csv <- downloadHandler(
    filename = function() paste0("afs_results_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(rv$site_profit)
      write.csv(rv$site_profit, file, row.names = FALSE)
    }
  )
}