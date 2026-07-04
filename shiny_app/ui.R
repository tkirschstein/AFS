
library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)
library(leaflet)
library(DT)

dashboardPage(
  skin = "green",
  
  # --------------------------------------------------------------------------
  # HEADER
  # --------------------------------------------------------------------------
  dashboardHeader(
    title = "AFS-SCD: Case Study Explorer",
    titleWidth = 300
  ),
  
  # --------------------------------------------------------------------------
  # SIDEBAR — Eingabeparameter des Optimierungsmodells
  # --------------------------------------------------------------------------
  dashboardSidebar(
    width = 300,
    useShinyjs(),
    
    # -- Wachstumsparameter --------------------------------------------------
    h4("Biomass Growth", style = "color:#a8d5a2; margin:12px 12px 4px;"),
    sliderInput("N_trees",   "Tree density (trees/ha)",
                min = 100, max = 500, value = 250, step = 25),
    sliderInput("C_site",    "Site quality C_site",
                min = 3.0, max = 9.0, value = 6.5, step = 0.5),
    sliderInput("k_gomp",    "Gompertz k (growth rate)",
                min = 0.10, max = 0.35, value = 0.194, step = 0.005),
    sliderInput("t0_gomp",   "Gompertz t0 (inflection yr)",
                min = 6, max = 14, value = 9.7, step = 0.1),
    sliderInput("moisture",  "Moisture content at felling",
                min = 0.30, max = 0.65, value = 0.50, step = 0.05),
    
    hr(),
    
    # -- Planungshorizont & Rotation -----------------------------------------
    h4("Planning & Rotation", style = "color:#a8d5a2; margin:12px 12px 4px;"),
    sliderInput("n_periods", "Planning horizon (years)",
                min = 20, max = 60, value = 40, step = 5),
    sliderInput("min_age",   "Min rotation age (years)",
                min = 1, max = 7, value = 3, step = 1),
    sliderInput("max_age",   "Max rotation age (years)",
                min = 8, max = 25, value = 20, step = 1),
    
    hr(),
    
    # -- Standortkosten -------------------------------------------------------
    h4("Site Economics", style = "color:#a8d5a2; margin:12px 12px 4px;"),
    numericInput("C_est",    "Establishment cost (€/ha)",
                 value = 2500, min = 500, max = 8000, step = 100),
    numericInput("C_harv",   "Harvest cost (€/ha)",
                 value = 150, min = 20, max = 600, step = 10),
    numericInput("C_main",   "Maintenance cost (€/ha/yr)",
                 value = 10, min = 0, max = 200, step = 5),
    numericInput("opp_mean", "Mean opportunity cost (€/ha/yr)",
                 value = 500, min = 0, max = 1200, step = 25),
    
    hr(),
    
    # -- Logistikkosten -------------------------------------------------------
    h4("Logistics", style = "color:#a8d5a2; margin:12px 12px 4px;"),
    numericInput("c_tr_raw", "Transport cost raw (€/t·km)",
                 value = 0.08, min = 0.01, max = 0.50, step = 0.01),
    numericInput("c_tr_pre", "Transport cost chips (€/t·km)",
                 value = 0.06, min = 0.01, max = 0.50, step = 0.01),
    
    hr(),
    
    # -- Sensitivitätsszenarien -----------------------------------------------
    h4("Revenue Multipliers", style = "color:#a8d5a2; margin:12px 12px 4px;"),
    sliderInput("rev_mult_P1",  "Revenue multiplier P1",
                min = 0, max = 10.0, value = 1.0, step = 0.25),
    sliderInput("rev_mult_P2",  "Revenue multiplier P2",
                min = 0, max = 10.0, value = 1.0, step = 0.25),
    sliderInput("rev_mult_P3",  "Revenue multiplier P3",
                min = 0, max = 10.0, value = 1.0, step = 0.25),
    
    hr(),
    
    # -- Steuerung ------------------------------------------------------------
    actionButton("run_opt", "Run Optimization",
                 icon = icon("play"),
                 style = "background:#2c7a2c; color:white; width:90%; margin:6px 5%;"),
    br(),
    downloadButton("export_csv", "Export Results",
                   style = "width:90%; margin:6px 5%;"),
    br(),
    tags$small(
      style = "color:#aaa; margin:8px 12px; display:block;",
      "Results load from precomputed RDS when available."
    )
  ),
  
  # --------------------------------------------------------------------------
  # BODY
  # --------------------------------------------------------------------------
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .box-title { font-weight: bold; }
        .info-box-number { font-size: 1.4em; }
        .nav-tabs-custom .nav-tabs li.active a { border-top-color: #2c7a2c; }
        .shiny-notification { font-size: 0.9em; }
      "))
    ),
    
    tabsetPanel(
      id = "main_tabs",
      
      # ======================================================================
      # TAB 1: Wachstumsmodell
      # ======================================================================
      tabPanel(
        title = tagList(icon("seedling"), "Biomass Growth"),
        value = "tab_growth",
        
        fluidRow(
          # KPI-Boxen
          infoBoxOutput("kpi_max_yield", width = 3),
          infoBoxOutput("kpi_stem_share_10", width = 3),
          infoBoxOutput("kpi_optimal_age", width = 3),
          infoBoxOutput("kpi_fresh_weight", width = 3)
        ),
        
        fluidRow(
          box(
            title = "Stand-level Biomass Components",
            status = "success", solidHeader = TRUE, width = 12,
            helpText("Stacked biomass fractions: stem (P1), branch (P2), residue (P3) over stand age."),
            plotlyOutput("plot_growth_stacked", height = 380)
          )
        ),
        
        fluidRow(
          box(
            title = "Tree and stand-level asymptotes",
            status = "warning", solidHeader = TRUE, width = 6,
            collapsible = TRUE, collapsed = TRUE,
            plotlyOutput("plot_growth_asymptotes", height = 300)
          ),
          box(
            title = "Biomass fraction growth",
            status = "primary", solidHeader = TRUE, width = 6,
            collapsible = TRUE, collapsed = TRUE,
            plotlyOutput("plot_growth_fractions", height = 300)
          )
        )
      ),
      
      # ======================================================================
      # TAB 2: Netzwerk & Karte  (Fig. 4 / Fig. 5 equivalent)
      # ======================================================================
      tabPanel(
        title = tagList(icon("map-marked-alt"), "Network & Map"),
        value = "tab_network",
        
        fluidRow(
          infoBoxOutput("kpi_n_sites",   width = 3),
          infoBoxOutput("kpi_total_area", width = 3),
          infoBoxOutput("kpi_obj_val",    width = 3),
          infoBoxOutput("kpi_solver_gap", width = 3)
        ),
        
        fluidRow(
          box(
            title = "Supply Chain Network — Active Sites & Flows",
            status = "primary", solidHeader = TRUE, width = 8,
            helpText("Green markers = active AFS sites; orange = hubs/storages;
                     purple = industrial consumers. Line width ~ flow volume."),
            leafletOutput("map_network", height = 560)
          ),
          box(
            title = "Result Summary",
            status = "success", solidHeader = TRUE, width = 4,
            h5("Solver Status"), textOutput("txt_solver_status"), br(),
            h5("Objective Value (€)"), textOutput("txt_obj"), br(),
            h5("Active Sites"), textOutput("txt_active_sites"), br(),
            h5("Total AFS Area (ha)"), textOutput("txt_total_area"), br(),
            h5("Mean Profit (€/ha/yr)"), textOutput("txt_mean_profit"), br(),
            hr(),
            h5("Legend"),
            tags$ul(
              tags$li(HTML("<span style='color:darkgreen;'>&#11044;</span> Active AFS site")),
              tags$li(HTML("<span style='color:orange;'>&#9650;</span> Hub / Storage")),
              tags$li(HTML("<span style='color:purple;'>&#9632;</span> Industrial consumer")),
              tags$li(HTML("<span style='color:steelblue;'>&#9135;</span> Raw biomass flow")),
              tags$li(HTML("<span style='color:darkblue;'>&#9135;</span> Processed flow"))
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Site Details Table",
            status = "primary", solidHeader = TRUE, width = 12,
            collapsible = TRUE, collapsed = TRUE,
            DT::dataTableOutput("table_sites")
          )
        )
      ),
      
      # ======================================================================
      # TAB 3: Biomassströme & Ergebnisse (Fig. 6 equivalent)
      # ======================================================================
      tabPanel(
        title = tagList(icon("chart-area"), "Material Flows"),
        value = "tab_flows",
        
        fluidRow(
          box(
            title = "Harvested Biomass by Product over Time",
            status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_biomass_time", height = 320)
          ),
          box(
            title = "Alluvial: Biomass Cascade (Site → Hub → Consumer)",
            status = "info", solidHeader = TRUE, width = 6,
            helpText("Sankey diagram of total flows aggregated over planning horizon."),
            plotlyOutput("plot_sankey", height = 320)
          )
        ),
        
        fluidRow(
          box(
            title = "Revenue by Product & Consumer",
            status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_rev_consumer", height = 320)
          ),
          box(
            title = "Storage Levels over Time",
            status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_storage_time", height = 320)
          )
        ),
        
        fluidRow(
          box(
            title = "Demand Fulfilment Rate by Consumer & Product",
            status = "primary", solidHeader = TRUE, width = 12,
            plotlyOutput("plot_demand_fulfilment", height = 300)
          )
        )
      ),
      
      # ======================================================================
      # TAB 4: Site-level KPIs  (Fig. 7 equivalent)
      # ======================================================================
      tabPanel(
        title = tagList(icon("map-pin"), "Site KPIs"),
        value = "tab_site_kpis",
        
        fluidRow(
          box(
            title = "Profit per ha·yr by Opportunity Cost",
            status = "primary", solidHeader = TRUE, width = 6,
            helpText("Each dot = one AFS site. X-axis = site-specific opportunity cost (€/ha/yr).
                     Dashed line = break-even."),
            plotlyOutput("plot_profit_vs_opp", height = 380)
          ),
          box(
            title = "Site Activation vs. Opportunity Cost",
            status = "info", solidHeader = TRUE, width = 6,
            helpText("Share of activated sites as function of opportunity cost threshold."),
            plotlyOutput("plot_activation_rate", height = 380)
          )
        ),
        
        fluidRow(
          box(
            title = "Rotation Length Distribution",
            status = "warning", solidHeader = TRUE, width = 4,
            plotlyOutput("plot_rotation_dist", height = 300)
          ),
          box(
            title = "P1 (Stem) Share vs. Profit",
            status = "success", solidHeader = TRUE, width = 4,
            helpText("Sites with higher stem share earn more per ha·yr."),
            plotlyOutput("plot_p1share_profit", height = 300)
          ),
          box(
            title = "Distance to Hub vs. Profit",
            status = "primary", solidHeader = TRUE, width = 4,
            plotlyOutput("plot_dist_profit", height = 300)
          )
        )
      )
    ) # end tabsetPanel
  ) # end dashboardBody
) # end dashboardPage

