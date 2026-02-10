# ============================================================================
# Agroforestry Supply Chain Optimization Shiny App - Version 10.0
# ============================================================================

library(shiny)
library(leaflet)
library(sf)
library(slam)
library(dplyr)
library(ompr)
library(ompr.roi)
library(ROI)
library(ROI.plugin.glpk)
library(ROI.plugin.gurobi)
library(ROI.plugin.highs)
library(DT)
library(Matrix)
library(shinydashboard)
library(shinyjs)
library(tidyr)
library(purrr)
library(tibble)
library(rhandsontable)
library(ggplot2)
library(plotly)
library(data.table)

# Load helper functions and optimization model
source("!helper_func.r")
source("build_agroforestry_lp_v2.R")
source("build_agroforestry_lp_v10.r")
source("!helper_instance_builder_v8a.R")

# Load spatial data
load("!base_data.rdata")

# drop small feldblocks
data_fb_filtered <- data_fb_filtered %>% filter(area >= 10)
fb_zentroide_wgs84 <- fb_zentroide_wgs84 %>% filter(area >= 10)

# ============================================================================
# CONSTANTS & CONFIGURATION
# ============================================================================

FLOW_VOLUME_THRESHOLD <- 0.01
HARVEST_THRESHOLD     <- 0.01
DEFAULT_MAP_CENTER    <- list(lng = 12.0, lat = 51.8)
DEFAULT_MAP_ZOOM      <- 8
DEFAULT_PRODUCT_NAMES <- c("Chemical", "Pulp", "Energy")
DEFAULT_PRODUCT_COLORS <- c("darkgreen", "darkblue", "darkorange")

# ============================================================================
# ALLOMETRIC YIELD FUNCTIONS FOR 3 PRODUCTS
# ============================================================================

allometric_product1 <- function(age, a_max = 30, k = 0.25) {
  a_max * (1 - exp(-k * age))
}

allometric_product2 <- function(age, a_max = 25, k = 0.4, inflection = 5) {
  a_max / (1 + exp(-k * (age - inflection)))
}

allometric_product3 <- function(age, a_max = 28, b = 3, p = 1.5) {
  a_max * (age / (age + b))^p
}

# ============================================================================
# BUILD YIELDS TABLE - NEW FORMAT: rows=products, columns=ages
# ============================================================================

build_yields_by_age_v9 <- function(max_age,
                                   a_max_p1 = 30, k_p1 = 0.25,
                                   a_max_p2 = 25, k_p2 = 0.4, inflection_p2 = 5,
                                   a_max_p3 = 28, b_p3 = 3, p_p3 = 1.5) {
  ages <- 1:max_age
  
  yields_p1 <- sapply(ages, function(a) allometric_product1(a, a_max_p1, k_p1))
  yields_p2 <- sapply(ages, function(a) allometric_product2(a, a_max_p2, k_p2, inflection_p2))
  yields_p3 <- sapply(ages, function(a) allometric_product3(a, a_max_p3, b_p3, p_p3))
  
  # Create matrix: rows = products, columns = ages
  yields_matrix <- rbind(yields_p1, yields_p2, yields_p3) 
  
  # Convert to data frame with product names as row labels
  df <- as.data.frame(yields_matrix)
  colnames(df) <- paste0("Age_", ages)
  df <- cbind(Product = c("Chemical", "Pulp", "Energy"), df)
  
  return(df)
}

# ============================================================================
# CONVERT WIDE FORMAT TO LONG FORMAT FOR OPTIMIZATION
# ============================================================================

yields_wide_to_long <- function(yields_wide) {
  # yields_wide: Product | Age_1 | Age_2 | ... | Age_max_age
  # Returns: product | age | yield_ha
  
  product_names <- yields_wide$Product
  n_products <- length(product_names)
  
  # Extract numeric columns (Age_1, Age_2, ...)
  age_cols <- colnames(yields_wide)[grepl("^Age_", colnames(yields_wide))]
  ages <- as.numeric(gsub("Age_", "", age_cols))
  
  # Reshape to long format
  yields_long <- data.frame(
    age = rep(ages, each = n_products),
    product = rep(1:n_products, by = length(ages)),
    yield_ha = as.numeric(unlist(yields_wide[, age_cols]))
  )
  
  return(yields_long)
}

# ============================================================================
# HELPER: Aggregate solution for time-series analysis

# ============================================================================

build_time_series_outputs <- function(opt_result, instance) {
  Y   <- opt_result$solution$Y
  S   <- opt_result$solution$S
  Xjk <- opt_result$solution$Xjk
  Xij <- opt_result$solution$Xij
  Dkpt <- opt_result$solution$Dkpt
  
  # Harvest per period & product in kt
  biomass_prod <- Y %>%
    group_by(t, p) %>%
    summarise(volume = sum(value)/1000, .groups = "drop") %>%
    rename(period = t, product = p)
  
  # Storage level per period & product
  stock_levels <- S %>%
    group_by(t, p) %>%
    summarise(stock = sum(value)/1000, .groups = "drop") %>%
    rename(period = t, product = p)
  
  # Demand satisfied per consumer, period, product
  demand_sat <- Xjk %>%
    group_by(k, t, pp) %>%
    summarise(volume = sum(value), .groups = "drop") %>%
    left_join(instance$demand, by = c("k" = "consumer_id", "t" = "period","pp" = "product")) %>%
    rename(period = t, product = pp, consumer_id = k) %>% 
    mutate(share = ifelse(D_max > 0, volume / D_max, 0)) %>% 
    left_join(instance$consumer_prices, by = c("consumer_id", "product") ) %>% #prices
    mutate(revenue = volume * price) %>% 
    # scale D_max and volume to kt
    mutate(D_max = D_max / 1000,
           value = volume / 1000)
  
  
  # wide table biomass produced periods x products
  biomass_wide <- biomass_prod %>%
    pivot_wider(names_from = product, values_from = volume, 
                names_prefix = "Product_",
                values_fill = 0) %>% 
    select(-period) %>% 
    t()
  
  # demand volume satisfied wide
  demand_sat_wide <- demand_sat %>%
    group_by(period, product) %>%
    summarise(total_volume = sum(volume), .groups = "drop") %>%
    pivot_wider(names_from = product, values_from = total_volume,
                names_prefix = "Product_",
                values_fill = 0) %>%
    select(-period) %>%
    t()
  
  # revenues wide
  revenues_wide <- demand_sat %>%
    group_by(period, product) %>%
    summarise(total_revenue = sum(revenue), .groups = "drop") %>%
    pivot_wider(names_from = product, values_from = total_revenue,
                names_prefix = "Product_",
                values_fill = 0) %>%
    select(-period) %>%
    t()
  
  
  list(
    biomass_prod = biomass_prod,
    stock_levels = stock_levels,
    demand_sat   = demand_sat,
    biomass_wide = biomass_wide,
    demand_sat_wide = demand_sat_wide,
    revenues_wide = revenues_wide
  )
}

# ============================================================================
# DATA PREPARATION FUNCTION
# ============================================================================

prepare_data <- function(cutoff_area) {
  sites_sf      <- data_fb_filtered
  sites_coords  <- fb_zentroide_wgs84
  
  if ("HBN" %in% colnames(sites_sf)) {
    if ("AL" %in% unique(sites_sf$HBN)) {
      sites_sf <- sites_sf %>% filter(HBN == "AL", area >= cutoff_area)
      sites_coords <- sites_coords %>%
        filter(HBN == "AL", area >= cutoff_area) %>%
        sf::st_coordinates()
    }
  }
  
  sites <- data.frame(
    site_id   = 1:nrow(sites_coords),
    name      = paste0("Site ", 1:nrow(sites_coords)),
    lat       = sites_coords[, "Y"],
    lng       = sites_coords[, "X"],
    area_ha   = sites_sf$area,
    C_est     = 5000,
    C_harv    = 50,
    stringsAsFactors = FALSE
  )
  
  list(
    sites_sf = sites_sf,
    sites    = sites
  )
}

# ============================================================================
# HELPER: Safe ggplot to plotly conversion
# ============================================================================

safe_ggplotly <- function(gg_plot, tooltip_str = "auto") {
  tryCatch(
    {
      if (tooltip_str == "auto") {
        ggplotly(gg_plot, tooltip = "text")
      } else {
        ggplotly(gg_plot, tooltip = tooltip_str)
      }
    },
    error = function(e) {
      message("ggplotly conversion failed: ", conditionMessage(e))
      return(NULL)
    }
  )
}

# ============================================================================
# SHINY UI
# ============================================================================

ui <- dashboardPage(
  dashboardHeader(title = "Agroforestry Supply Chain Optimization v9.0"),
  dashboardSidebar(
    width = 320,
    h4("Growth Parameters", style = "color: #2c5f2d;"),
    sliderInput("n_periods",
                "Planning Horizon (Periods)",
                min = 3, max = 50, value = 25, step = 1
    ),
    sliderInput("min_age",
                "Min Plantation Age (Years)",
                min = 1, max = 8, value = 3, step = 1
    ),
    sliderInput("max_age",
                "Max Plantation Age (Years)",
                min = 3, max = 15, value = 8, step = 1
    ),
    hr(),
    h4("Economic Parameters", style = "color: #2c5f2d;"),
    numericInput("c_est_area",
                 "AFS Establishment cost (€/ha)",
                 value = 100, min = 00, max = 500000, step = 10
    ),
    numericInput("c_opp_area",
                 "AFS opportunity cost (€/ha/year)",
                 value = 10, min = 0, max = 500000, step = 1
    ),
    numericInput("c_harv_area",
                 "AFS harvesting cost (€/ha)",
                 value = 10, min = 0, max = 500000, step = 1
    ),
    numericInput("c_tr_raw",
                 "Transport Cost Raw (€/km/t)",
                 value = 0.05, min = 0.01, max = 5, step = 0.01
    ),
    numericInput("c_tr_pre",
                 "Transport Cost Processed (€/km/t)",
                 value = 0.1, min = 0.01, max = 5, step = 0.01
    ),
    hr(),
    h4("Optimzation Parameters", style = "color: #2c5f2d;"),
    numericInput("time_limit",
                 "Time Limit (seconds)",
                 value = 300, min = 30, max = 3600, step = 30
    ),
    numericInput("mip_gap",
                 "MIP Gap (%)",
                 value = 5, min = 0.1, max = 20, step = 0.1
    ),
    hr(),
    h4("Visualization", style = "color: #2c5f2d;"),
    checkboxInput("show_ac_sites", "Show agricultural sites", value = TRUE),
    sliderInput("show_ac_area",
                "Minimal site size (ha)",
                min = 1, max = 1000, value = 300, step = 1
    ),
    checkboxInput("show_flows", "Show Supply Chain Flows", value = TRUE),
    checkboxInput("show_centroids", "Show Site Centroids", value = TRUE),
    checkboxInput("show_storages", "Show Pretreatment Facilities", value = TRUE),
    checkboxInput("show_consumers", "Show Consumer Facilities", value = TRUE),
    hr(),
    h4("Actions", style = "color: #2c5f2d;"),
    actionButton("update_data",
                 "Update Data",
                 class = "btn-warning",
                 width = "100%",
                 style = "background: #ff9800; color: white; margin-bottom: 10px;"
    ),
    br(),
    actionButton("run_optimization",
                 "Run Optimization",
                 class = "btn-primary",
                 width = "100%",
                 style = "background: #2c5f2d; color: white; margin-bottom: 10px;"
    ),
    br(),
    actionButton("reset_view",
                 "Reset View",
                 width = "100%",
                 style = "margin-bottom: 10px;"
    ),
    br(),
    downloadButton("export_results",
                   "Export Results",
                   width = "100%"
    )
  ),
  dashboardBody(
    useShinyjs(),
    tabBox(
      width = 12,
      
      # --- TAB 1: Network & Summary (WITH DATA INPUTS) ---
      tabPanel(
        title = "Network & Summary",
        fluidRow(
          column(
            width = 9,
            box(
              width = 12,
              title = "Supply Chain Network Visualization",
              status = "primary",
              solidHeader = TRUE,
              leafletOutput("optimization_map", height = 600)
            )
          ),
          column(
            width = 3,
            box(
              width = 12,
              title = "Optimization Results",
              status = "success",
              solidHeader = TRUE,
              h5("Status"),
              textOutput("status_text"),
              br(),
              h5("Objective Value"),
              textOutput("objective_value"),
              br(),
              h5("Total Flow Volume"),
              textOutput("total_volume"),
              br(),
              h5("Active Sites"),
              textOutput("active_sites"),
              br(),
              h5("Solver Status"),
              textOutput("solver_status")
            ),
            box(
              width = 12,
              title = "Legend",
              p("🟢 Production Sites"),
              p("🟠 Pretreatment Facilities"),
              p("🟣 Consumer Facilities"),
              p("↦ Supply Chain Flows")
            )
          )
        ),
        
        # --- FACILITIES DATA (MOVED FROM TAB 2) ---
        fluidRow(
          column(
            width = 6,
            box(
              width = 12,
              title = "Pretreatment Facilities (Sawmills)",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              collapsed = TRUE,
              helpText("Details of pretreatment facilities. Processing and storage capacities in kt of biomass per year. Storage cost are given in € per ton and year.
                       You can add/remove rows and edit coordinates. Click 'Update Data' to recalculate distances and update map."),
              rHandsontableOutput("pre_fac")
            )
          ),
          column(
            width = 6,
            box(
              width = 12,
              title = "Consumer Facilities (Biochem Sites)",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              collapsed = TRUE,
              helpText("Add details of consumers. Total demand in kt of biomass per year  and product. Prices for each product with P1 = Chemical (€/t), P2 = Pulp (€/t), P3 = Energy (€/t). You can add/remove rows and edit coordinates."),
              rHandsontableOutput("cons_fac")
            )
          )
        ),
        
        # --- GROWTH & YIELD DATA (MOVED FROM TAB 3) ---
        fluidRow(
          box(
            width = 12,
            title = "Allometric Yield Functions (3 Products)",
            status = "warning",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            h5("Product 1 (Chemical): y(a) = a_max · (1 - exp(-k·a))"),
            h5("Product 2 (Pulp): y(a) = a_max / (1 + exp(-k·(a-inflection)))"),
            h5("Product 3 (Energy): y(a) = a_max · (a / (a+b))^p")
          )
        ),
        fluidRow(
          box(
            width = 4,
            title = "Product 1 Parameters",
            status = "warning",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            numericInput("a_max_p1", "a_max (t/ha)", value = 30, min = 1, max = 100),
            numericInput("k_p1", "k (growth rate)", value = 0.25, min = 0.01, max = 2, step = 0.01)
          ),
          box(
            width = 4,
            title = "Product 2 Parameters",
            status = "warning",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            numericInput("a_max_p2", "a_max (t/ha)", value = 25, min = 1, max = 100),
            numericInput("k_p2", "k (growth rate)", value = 0.4, min = 0.01, max = 2, step = 0.01),
            numericInput("inflection_p2", "Inflection (years)", value = 5, min = 1, max = 10)
          ),
          box(
            width = 4,
            title = "Product 3 Parameters",
            status = "warning",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            numericInput("a_max_p3", "a_max (t/ha)", value = 28, min = 1, max = 100),
            numericInput("b_p3", "b (shape)", value = 3, min = 0.1, max = 10, step = 0.1),
            numericInput("p_p3", "p (exponent)", value = 1.5, min = 0.5, max = 3, step = 0.1)
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Yield Curves for All 3 Products",
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            plotlyOutput("plot_yield_curves_all", height = 400)
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Underlying Yield Table (η_{p,a}) - Edit to Customize",
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            helpText("NEW FORMAT: Rows = Products, Columns = Ages. Edit yield values directly. Click 'Update Data' to apply. Table updates dynamically with max_age parameter."),
            rHandsontableOutput("table_yields_by_age")
          )
        )
      ),
      
      # --- TAB 2: Result Analysis (WITH DETAILED RESULTS TABLE) ---
      tabPanel(
        title = "Result Analysis",
        fluidRow(
          box(
            width = 6,
            title = "Produced Biomass over Time",
            status = "primary",
            solidHeader = TRUE,
            plotlyOutput("plot_biomass_prod", height = 320)
          ),
          box(
            width = 6,
            title = "Storage Stock Levels over Time",
            status = "primary",
            solidHeader = TRUE,
            plotlyOutput("plot_stock_levels", height = 320)
          )
        ),
        fluidRow(
          box(
            width = 6,
            title = "Demand Satisfied over Time",
            status = "success",
            solidHeader = TRUE,
            plotlyOutput("plot_demand_sat", height = 320)
          ),
          
          # place for an plot
        ),
        fluidRow(
          box(
            width = 6,
            title = "Revenues by Product and Site",
            status = "primary",
            solidHeader = TRUE,
            plotlyOutput("plot_revenues_sites", height = 320)
          ),
          box(
            width = 6,
            title = "Revenues by Product and Period",
            status = "success",
            solidHeader = TRUE,
            plotlyOutput("plot_revenues", height = 320)
          )
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "Harvested Biomass by Production Site over Time",
            status = "primary",
            solidHeader = TRUE,
            helpText("Aggregated harvested biomass per production site and period (all products combined)."),
            plotlyOutput("plot_harvest_by_site", height = 380)
          )
        ),
        
        # --- DETAILED RESULTS TABLE (MOVED FROM TAB 1) ---
        fluidRow(
          column(
            width = 12,
            box(
              width = 12,
              title = "Detailed Results",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              collapsed = FALSE,
              DT::dataTableOutput("flows_table")
            )
          )
        )
      )
    )
  )
)

# ============================================================================
# SHINY SERVER
# ============================================================================

server <- function(input, output, session) {
  
  # Initialize data lazily
  base_data <- reactive({
    showNotification("Loading base data...", type = "message", duration = 2)
    res <- prepare_data(cutoff_area = input$show_ac_area)
    res
  })
  
  # Reactive values for user-editable data
  user_data <- reactiveValues(
    storages = NULL,
    consumers = NULL,
    yields_wide = NULL,
    dist_ij = NULL,
    dist_jk = NULL,
    data_updated = FALSE
  )
  
  # Initialize facilities data from spatial data
  observeEvent(base_data(), {
    if (is.null(user_data$storages)) {
      sawmills_coords <- sf::st_coordinates(sawmills_sf)
      user_data$storages <- data.frame(
        storage_id = 1:nrow(sawmills_coords),
        name = sawmills_sf$name,
        lat = sawmills_coords[, "Y"],
        lng = sawmills_coords[, "X"],
        CAP_stor = 500,
        CAP_proc = c(700,500,250,40,150,250,90),
        c_stor = 1,
        stringsAsFactors = FALSE
      )
    }
    
    if (is.null(user_data$consumers)) {
      biochem_coords <- sf::st_coordinates(biochem_sites_sf)
      user_data$consumers <- data.frame(
        consumer_id = 1:nrow(biochem_coords),
        name = biochem_sites_sf$name,
        lat = biochem_coords[, "Y"],
        lng = biochem_coords[, "X"],
        demand_P1 = c(150,0  ,100,0  ,0),
        demand_P2 = c(0  ,150,70 ,0  ,0),
        demand_P3 = c(0  ,50 ,0  ,120,40),
        P1 = c(150,0,150,0,0),
        P2 = c(0,120,70,0,0),
        P3 = c(0,60,0,50,40),
        stringsAsFactors = FALSE
      )
    }
  }, once = TRUE)
  
  # Yields table: parameter-derived version (wide format) - DYNAMIC WITH max_age
  yields_wide_from_params <- reactive({
    build_yields_by_age_v9(
      max_age = input$max_age,
      a_max_p1 = input$a_max_p1, k_p1 = input$k_p1,
      a_max_p2 = input$a_max_p2, k_p2 = input$k_p2, inflection_p2 = input$inflection_p2,
      a_max_p3 = input$a_max_p3, b_p3 = input$b_p3, p_p3 = input$p_p3
    )
  })
  
  # Update yields_wide when max_age or parameters change
  observeEvent(list(input$max_age, input$a_max_p1, input$k_p1, 
                    input$a_max_p2, input$k_p2, input$inflection_p2,
                    input$a_max_p3, input$b_p3, input$p_p3), {
                      user_data$yields_wide <- yields_wide_from_params()
                    }, ignoreNULL = TRUE)
  
  # Get current yields (either from params or user edits)
  yields_wide_current <- reactive({
    if (!is.null(user_data$yields_wide)) {
      user_data$yields_wide
    } else {
      yields_wide_from_params()
    }
  })
  
  # ============================================================================
  # rHandsontable outputs
  # ============================================================================
  
  output$pre_fac <- renderRHandsontable({
    req(user_data$storages)
    
    
    rhandsontable(user_data$storages, rowHeaders = NULL, stretchH = "all") %>%
      hot_col("storage_id", readOnly = TRUE) %>%
      hot_col("name", readOnly = FALSE) %>%
      hot_col("lat", readOnly = FALSE) %>%
      hot_col("lng", readOnly = FALSE) %>%
      hot_col("CAP_stor", readOnly = FALSE) %>%
      hot_col("CAP_proc", readOnly = FALSE) %>%
      hot_col("c_stor", readOnly = FALSE)
  })
  
  output$cons_fac <- renderRHandsontable({
    req(user_data$consumers)
    rhandsontable(user_data$consumers, rowHeaders = NULL, stretchH = "all") %>%
      hot_col("consumer_id", readOnly = TRUE) %>%
      hot_col("name", readOnly = FALSE) %>%
      hot_col("lat", readOnly = FALSE) %>%
      hot_col("lng", readOnly = FALSE) %>%
      hot_col("demand_P1", readOnly = FALSE) %>%
      hot_col("demand_P2", readOnly = FALSE) %>%
      hot_col("demand_P3", readOnly = FALSE) %>%
      hot_col("P1", readOnly = FALSE) %>%
      hot_col("P2", readOnly = FALSE) %>%
      hot_col("P3", readOnly = FALSE)
  })
  
  output$table_yields_by_age <- renderRHandsontable({
    req(yields_wide_current())
    df <- yields_wide_current()
    
    
    rhandsontable(df, rowHeaders = NULL, stretchH = "all") %>%
      hot_col("Product", readOnly = TRUE)
  })
  
  # ============================================================================
  # UPDATE DATA BUTTON - Recalculate distances and update map
  # ============================================================================
  
  observeEvent(input$update_data, {
    showNotification("Updating data...", type = "message", duration = 2)
    
    tryCatch({
      # Read current table states
      if (!is.null(input$pre_fac)) {
        updated_storages <- hot_to_r(input$pre_fac)
        if (!is.null(updated_storages) && nrow(updated_storages) > 0) {
          # Reassign storage_id sequentially
          updated_storages$storage_id <- 1:nrow(updated_storages)
          user_data$storages <- updated_storages
        }
      }
      
      if (!is.null(input$cons_fac)) {
        updated_consumers <- hot_to_r(input$cons_fac)
        if (!is.null(updated_consumers) && nrow(updated_consumers) > 0) {
          # Reassign consumer_id sequentially
          updated_consumers$consumer_id <- 1:nrow(updated_consumers)
          user_data$consumers <- updated_consumers
        }
      }
      
      if (!is.null(input$table_yields_by_age)) {
        user_data$yields_wide <- hot_to_r(input$table_yields_by_age)
      }
      
      
      # Recalculate distance matrices
      current_sites <- base_data()$sites
      
      cat("\n=== Recalculating distance matrices ===\n")
      cat("Sites:", nrow(current_sites), "\n")
      cat("Storages:", nrow(user_data$storages), "\n")
      cat("Consumers:", nrow(user_data$consumers), "\n")
      
      user_data$dist_ij <- calculate_distance_matrix_osm(
        starts = current_sites,
        destinations = user_data$storages
      )$distance_matrix_km
      
      user_data$dist_jk <- calculate_distance_matrix_osm(
        starts = user_data$storages,
        destinations = user_data$consumers
      )$distance_matrix_km
      
      user_data$data_updated <- TRUE
      
      showNotification("Data updated successfully! Map will refresh.", type = "message", duration = 3)
      
    }, error = function(e) {
      showNotification(paste("Error updating data:", conditionMessage(e)), type = "error", duration = 5)
      cat("Update error:", conditionMessage(e), "\n")
    })
  })
  
  # Results storage
  results <- reactiveValues(
    optimization = NULL,
    instance = NULL,
    solution = NULL,
    status = "Ready",
    objective = NA,
    total_volume = NA,
    active_sites = 0,
    solver_status = "—",
    timeseries = NULL
  )
  
  # Plot yield curves (plotly)
  output$plot_yield_curves_all <- renderPlotly({
    yields_wide <- yields_wide_current()
    yields_long <- yields_wide_to_long(yields_wide)
    
    p <- ggplot(yields_long, aes(x = age, y = yield_ha,
                                 color = factor(product),
                                 linetype = factor(product))) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2) +
      labs(
        x = "Age (years)",
        y = "Yield (t/ha)",
        title = "Yield Curves for All 3 Products"
      ) +
      scale_color_manual(
        values = DEFAULT_PRODUCT_COLORS,
        labels = DEFAULT_PRODUCT_NAMES
      ) +
      scale_linetype_manual(
        values = c("1" = "solid", "2" = "dashed", "3" = "dotted"),
        labels = DEFAULT_PRODUCT_NAMES
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 11),
        legend.position = "bottom"
      )
    safe_ggplotly(p)
  })
  
  # ============================================================================
  # RUN OPTIMIZATION
  # ============================================================================
  
  observeEvent(input$run_optimization, {
    showNotification("Starting optimization...", type = "message", duration = 2)
    
    #tryCatch({
    # Get current data
    current_sites <- isolate(base_data()$sites)
    current_storages <- isolate(user_data$storages)
    current_consumers <- isolate(user_data$consumers)
    
    # update harvesting and establishment cost
    current_sites$C_est <- input$c_est_area * current_sites$area_ha
    current_sites$C_harv <- input$c_harv_area * current_sites$area_ha
    
    # scale capacity and demand from kt to tons
    current_storages$CAP_stor <- current_storages$CAP_stor * 1000
    current_storages$CAP_proc <- current_storages$CAP_proc * 1000
    current_consumers$demand_P1 <- current_consumers$demand_P1 * 1000
    current_consumers$demand_P2 <- current_consumers$demand_P2 * 1000
    current_consumers$demand_P3 <- current_consumers$demand_P3 * 1000
    
    # Use updated distance matrices if available, otherwise calculate
    if (!is.null(user_data$dist_ij) && !is.null(user_data$dist_jk)) {
      dist_ij <- user_data$dist_ij
      dist_jk <- user_data$dist_jk
    } else {
      cat("\n=== Calculating distance matrices ===\n")
      dist_ij <- calculate_distance_matrix_osm(
        starts = current_sites,
        destinations = current_storages
      )$distance_matrix_km
      
      dist_jk <- calculate_distance_matrix_osm(
        starts = current_storages,
        destinations = current_consumers
      )$distance_matrix_km
      
      user_data$dist_ij <- dist_ij
      user_data$dist_jk <- dist_jk
    }
    
    # Prepare data object for optimization
    opt_data <- list(
      sites = current_sites,
      storages = current_storages,
      consumers = current_consumers,
      dist_ij = dist_ij,
      dist_jk = dist_jk
    )
    
    # Build parameters
    params <- list(
      n_periods = input$n_periods,
      max_age = input$max_age,
      min_age = input$min_age,
      c_tr_raw = input$c_tr_raw,
      c_tr_pre = input$c_tr_pre,
      c_opp = input$c_opp_area,
      n_products = 3
    )
    
    # Convert yields from wide to long format
    yields_wide <- isolate(yields_wide_current())
    
    # set yields for ages between 1 and minimal age to 0
    yields_wide[, 2:input$min_age] <- 0
    
    custom_yields <- yields_wide_to_long(yields_wide)
    
    # Build optimization instance
    instance <- generate_instance_v8_final(
      data = opt_data,
      params = params,
      yields_by_age = custom_yields
    )
    
    cat("\n=== Solving optimization problem ===\n")
    
    # Solve
    opt_result <- solve_agroforestry_sparse(
      instance,
      time_limit = input$time_limit,
      mip_gap = input$mip_gap / 100,
      verbose = TRUE,
      solver = "gurobi"
    )
    
    # Extract and process results
    result_ext <- extract_result(opt_result)
    ts_out <- build_time_series_outputs(opt_result, instance)
    
    # Update results
    results$optimization <- opt_result
    results$instance <- instance
    results$solution <- result_ext
    results$timeseries <- ts_out
    results$status <- "Complete"
    results$objective <- opt_result$objective
    results$active_sites <- length(result_ext$sites_est)
    results$total_volume <- sum(result_ext$Y$value)
    results$solver_status <- as.character(opt_result$status$code)
    
    #cat("\n=== Optimization complete ===\n")
    #cat("Objective:", results$objective, "\n")
    #cat("Active sites:", results$active_sites, "\n")
    
    showNotification("Optimization complete!", type = "message", duration = 3)
    
    # }, error = function(e) {
    #  results$status <- paste("Error:", conditionMessage(e))
    #  results$solver_status <- "Failed"
    #  showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = 5)
    #    cat("\n=== Optimization ERROR ===\n")
    #    cat("Error message:", conditionMessage(e), "\n")
    #  })
  })
  
  # ============================================================================
  # RENDER MAP - Updates when data changes
  # ============================================================================
  
  output$optimization_map <- renderLeaflet({
    # Trigger re-render when data is updated
    user_data$data_updated
    
    m <- leaflet() %>%
      addTiles() %>%
      setView(lng = DEFAULT_MAP_CENTER$lng, lat = DEFAULT_MAP_CENTER$lat, zoom = DEFAULT_MAP_ZOOM) %>%
      addScaleBar(position = "bottomleft")
    
    current_data <- base_data()
    
    if (input$show_ac_sites && !is.null(current_data$sites_sf)) {
      m <- m %>%
        addPolygons(
          data = current_data$sites_sf,
          color = "darkgreen",
          weight = 1,
          fillOpacity = 0.5,
          popup = paste("Area (ha):", round(current_data$sites_sf$area, 0))
        )
    }
    
    if (input$show_centroids && !is.null(current_data$sites)) {
      m <- m %>% addCircleMarkers(
        data = current_data$sites,
        lng = ~lng, lat = ~lat,
        radius = 2,
        color = "darkgreen",
        fillColor = "lightgreen",
        weight = 2,
        opacity = 0.8,
        fillOpacity = 0.6,
        popup = ~paste0(name, " | Area: ", area_ha, " ha"),
        group = "Production Sites"
      )
    }
    
    if (input$show_storages && !is.null(user_data$storages)) {
      m <- m %>% addCircleMarkers(
        data = user_data$storages,
        lng = ~lng, lat = ~lat,
        radius = 4,
        color = "darkorange",
        fillColor = "orange",
        weight = 2,
        opacity = 0.9,
        fillOpacity = 0.8,
        popup = ~paste0(name, " | Storage: ", CAP_stor, " t"),
        group = "Pretreatment Facilities"
      )
    }
    
    if (input$show_consumers && !is.null(user_data$consumers)) {
      m <- m %>% addCircleMarkers(
        data = user_data$consumers,
        lng = ~lng, lat = ~lat,
        radius = 4,
        color = "darkviolet",
        fillColor = "plum",
        weight = 2,
        opacity = 0.9,
        fillOpacity = 0.8,
        popup = ~paste0(name, " | Total demand: ", demand_P1+demand_P2+demand_P3, " kt"),
        group = "Consumer Facilities"
      )
    }
    
    # Add optimization flows if available
    if (input$show_flows && !is.null(results$solution)) {
      
      flows_df_ij <- results$solution$Xij
      flows_df_jk <- results$solution$Xjk
      
      # aggregate flows over time
      flows_df_ij <- flows_df_ij %>%
        group_by(i, j, p) %>%
        summarise(volume = sum(value), .groups = "drop")
      flows_df_jk <- flows_df_jk %>%
        group_by(j, k, p) %>%
        summarise(volume = sum(value), .groups = "drop")
      
      # Get flow coordinates raw biomass
      flows_geom_ij <- flows_df_ij %>%
        left_join(
          results$instance$sites %>%
            select(site_id, lat, lng) %>%
            rename(from_lat = lat, from_lng = lng),
          by = c("i" = "site_id")
        ) %>%
        left_join(
          results$instance$storages %>%
            select(storage_id, lat, lng) %>%
            rename(to_lat = lat, to_lng = lng),
          by = c("j" = "storage_id")
        )
      
      # flows from pretreatment
      flows_geom_jk <- flows_df_jk %>%
        left_join(
          results$instance$storages %>%
            select(storage_id, lat, lng) %>%
            rename(from_lat = lat, from_lng = lng),
          by = c("j" = "storage_id")
        ) %>%
        left_join(
          results$instance$consumers %>%
            select(consumer_id, lat, lng) %>%
            rename(to_lat = lat, to_lng = lng),
          by = c("k" = "consumer_id")
        )
      
      # raw biomass flows
      for (idx in seq_len(nrow(flows_geom_ij))) {
        flow <- flows_geom_ij[idx, ]
        if (!is.na(flow$from_lat) && !is.na(flow$to_lat)) {
          color  <- c("red", "blue", "green")[flow$p]
          weight <- 1
          route  <- calculate_route_info(
            src_lat = flow$from_lat,
            src_lon = flow$from_lng,
            dst_lat = flow$to_lat,
            dst_lon = flow$to_lng
          )
          
          m <- m %>% addPolylines(
            data   = route$route,
            color  = color,
            weight = weight,
            opacity = 0.8,
            popup  = paste0(
              "Raw Biomass Flow\n",
              "Volume: ", round(flow$volume, 0), " t\n",
              "Product: ", flow$p
            ),
            group = "Raw Biomass Flows"
          )
        }
      }
      
      # pre-treated biomass flows
      for (idx in seq_len(nrow(flows_geom_jk))) {
        flow <- flows_geom_jk[idx, ]
        if (!is.na(flow$from_lat) && !is.na(flow$to_lat)) {
          color  <- c("darkred", "darkblue", "darkgreen")[flow$p]
          weight <- 1
          route  <- calculate_route_info(
            src_lat = flow$from_lat,
            src_lon = flow$from_lng,
            dst_lat = flow$to_lat,
            dst_lon = flow$to_lng
          )
          
          m <- m %>% addPolylines(
            data   = route$route,
            color  = color,
            weight = weight,
            opacity = 0.8,
            popup  = paste0(
              "Processed Biomass Flow\n",
              "Volume: ", round(flow$volume, 0), " t\n",
              "Product: ", flow$p
            ),
            group = "Processed Biomass Flows"
          )
        } 
      }
      m <- m %>% addLayersControl(
        overlayGroups = c(
          "Production Sites", "Pretreatment Facilities", "Consumer Facilities",
          "Raw Biomass Flows", "Processed Biomass Flows"
        ),
        options = layersControlOptions(collapsed = FALSE)
      )  
      
    }
    m
  })
  
  # Output text displays
  output$status_text <- renderText({ results$status })
  
  output$objective_value <- renderText({
    if (is.na(results$objective)) "—"
    else paste0("€", format(round(results$objective, 0), big.mark = ","))
  })
  
  output$total_volume <- renderText({
    if (is.na(results$total_volume)) "—"
    else paste0(round(results$total_volume/1000, 0), " kt")
  })
  
  output$active_sites <- renderText({ results$active_sites })
  output$solver_status <- renderText({ results$solver_status })
  
  # Flows table
  output$flows_table <- DT::renderDataTable({
    if (is.null(results$solution) || nrow(results$solution$Xij) == 0) {
      data.frame(Message = "Run optimization to see detailed results")
    } else {
      results$solution$Xij %>%
        filter(value > FLOW_VOLUME_THRESHOLD) %>%
        select(i, j, p, t, value) %>%
        rename(
          "From Site" = i,
          "To Facility" = j,
          "Product" = p,
          "Period" = t,
          "Volume (t)" = value
        ) %>%
        DT::datatable(options = list(pageLength = 10, scrollX = TRUE))
    }
  })
  
  # --- Plotly outputs ---
  
  output$plot_biomass_prod <- renderPlotly({
    req(results$timeseries)
    df <- results$timeseries$biomass_prod
    p <- ggplot(df, aes(x = period, y = volume, color = factor(product))) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2) +
      labs(x = "Period", y = "Biomass (kt)", title = "Produced Biomass over Time", color = "Product") +
      scale_color_manual(values = DEFAULT_PRODUCT_COLORS, labels = DEFAULT_PRODUCT_NAMES) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
    safe_ggplotly(p)
  })
  
  output$plot_stock_levels <- renderPlotly({
    req(results$timeseries)
    df <- results$timeseries$stock_levels
    p <- ggplot(df, aes(x = period, y = stock, color = factor(product))) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2) +
      labs(x = "Period", y = "Stock (kt)", title = "Storage Stock Levels over Time", color = "Product") +
      scale_color_manual(values = DEFAULT_PRODUCT_COLORS, labels = DEFAULT_PRODUCT_NAMES) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
    safe_ggplotly(p)
  })
  
  output$plot_demand_sat <- renderPlotly({
    req(results$timeseries)
    df <- results$timeseries$demand_sat %>%
      group_by(period, product) %>%
      summarise(volume = sum(volume), .groups = "drop")
    p <- ggplot(df, aes(x = period, y = volume, color = factor(product))) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2) +
      labs(x = "Period", y = "Volume (kt)", color = "Product") +
      scale_color_manual(values = DEFAULT_PRODUCT_COLORS, labels = DEFAULT_PRODUCT_NAMES) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
    safe_ggplotly(p)
  })
  
  output$plot_revenues <- renderPlotly({
    req(results$timeseries)
    df <- results$timeseries$demand_sat
    p <- ggplot(df, aes(x = period, y = revenue, fill = factor(product))) +
      geom_col(position = "stack") +
      labs(x = "Period", y = "Revenue (€)", title = "Revenues by Product and Period", fill = "Product") +
      scale_fill_manual(values = DEFAULT_PRODUCT_COLORS, labels = DEFAULT_PRODUCT_NAMES) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
    safe_ggplotly(p)
  })
  
  
  output$plot_revenues_sites <- renderPlotly({
    req(results$timeseries)
    
    df <- results$timeseries$demand_sat
    # replace consumer_id by site name
    df <- df %>%
      left_join(
        results$instance$consumers %>% select(consumer_id, name),
        by = "consumer_id"
      ) %>%
      rename(consumer_name = name)
    
    p <- ggplot(df, aes(x = consumer_name, y = revenue, fill = factor(product))) +
      geom_col(position = "stack") +
      labs(x = "Period", y = "Revenue (€)", fill = "Product") +
      scale_fill_manual(values = DEFAULT_PRODUCT_COLORS, labels = DEFAULT_PRODUCT_NAMES) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    
    safe_ggplotly(p)
  })
  
  
  output$plot_harvest_by_site <- renderPlotly({
    req(results$solution)
    Y <- results$solution$Y
    df <- Y %>%
      filter(value > HARVEST_THRESHOLD) %>%
      group_by(i, t) %>%
      summarise(volume = sum(value), .groups = "drop") %>%
      rename(site_id = i, period = t) %>%
      left_join(
        results$instance$sites %>% select(site_id, name),
        by = "site_id"
      )
    
    if (nrow(df) == 0) {
      return(NULL)
    }
    
    p <- ggplot(df, aes(x = period, y = volume, color = factor(site_id))) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 1.8) +
      labs(
        x = "Period",
        y = "Biomass Harvested (t)",
        title = "Harvested Biomass by Production Site over Time",
        color = "Site"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
      )
    safe_ggplotly(p)
  })
  
  # Reset view
  observeEvent(input$reset_view, {
    leafletProxy("optimization_map") %>%
      setView(lng = DEFAULT_MAP_CENTER$lng, lat = DEFAULT_MAP_CENTER$lat, zoom = DEFAULT_MAP_ZOOM)
    results$status <- "Ready"
    results$objective <- NA
    results$total_volume <- NA
    results$active_sites <- 0
    results$solver_status <- "—"
    results$solution <- NULL
    results$timeseries <- NULL
  })
  
  # Export results
  output$export_results <- downloadHandler(
    filename = function() {
      paste0("agroforestry_results_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      if (is.null(results$solution) || nrow(results$solution$Xij) == 0) {
        write.csv(data.frame(Message = "No optimization results available"), file, row.names = FALSE)
      } else {
        write.csv(
          results$solution$Xij %>% filter(value > FLOW_VOLUME_THRESHOLD),
          file,
          row.names = FALSE
        )
      }
    }
  )
}

# ============================================================================
# Run the application
# ============================================================================

shinyApp(ui = ui, server = server)
