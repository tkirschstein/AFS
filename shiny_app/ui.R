
library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)
library(leaflet)
library(DT)

# ==============================================================================
# DiP Sachsen-Anhalt Farbschema
#   Primärgrün   : #00B400  (DiP Hauptgrün - Navigationshintergrund)
#   Dunkelgrün   : #005A46  (DiP Akzentgrün - Logo-Hintergrund, Hover)
#   Hellgrün     : #AAC800  (DiP Gelbgrün   - Highlights, Labels)
#   Sidebar-BG   : #1c2b1c  (Dunkelgrüner Sidebar-Hintergrund)
#   Sekundärblau : #0095C3  (Info-Boxen, Karte)
# ==============================================================================

dip_css <- "
/* ── Header ─────────────────────────────────────────────────────────────── */
.skin-green .main-header .logo,
.skin-green .main-header .logo:hover {
  background-color: #005A46 !important;
  color: #ffffff !important;
  font-weight: 700;
}
.skin-green .main-header .navbar {
  background-color: #00B400 !important;
}
.skin-green .main-header .navbar .sidebar-toggle,
.skin-green .main-header .navbar .sidebar-toggle:hover {
  background-color: rgba(0,0,0,0.12) !important;
}
/* ── Sidebar ─────────────────────────────────────────────────────────────── */
.skin-green .main-sidebar,
.skin-green .left-side {
  background-color: #1c2b1c !important;
}
.skin-green .sidebar-menu > li > a {
  color: #c8dfc8 !important;
}
.skin-green .sidebar-menu > li:hover > a,
.skin-green .sidebar-menu > li.active > a {
  background-color: #005A46 !important;
  color: #ffffff !important;
  border-left-color: #AAC800 !important;
}
/* ── Slider ──────────────────────────────────────────────────────────────── */
.irs-bar, .irs-bar-edge {
  background: #00B400 !important;
  border-color: #005A46 !important;
}
.irs-slider {
  background: #00B400 !important;
  border: 2px solid #005A46 !important;
}
.irs-from, .irs-to, .irs-single {
  background: #005A46 !important;
  color: #fff !important;
}
/* ── Buttons ─────────────────────────────────────────────────────────────── */
.btn-success {
  background-color: #00B400 !important;
  border-color: #005A46 !important;
  color: #fff !important;
}
.btn-success:hover, .btn-success:focus {
  background-color: #009200 !important;
  border-color: #005A46 !important;
}
.btn-primary {
  background-color: #0095C3 !important;
  border-color: #006da3 !important;
}
/* ── Box-Header (status-Farben) ───────────────────────────────────────────── */
.box.box-success > .box-header { background: #005A46 !important; color: #fff !important; }
.box.box-primary > .box-header { background: #00873e !important; color: #fff !important; }
.box.box-info    > .box-header { background: #0095C3 !important; color: #fff !important; }
.box.box-warning > .box-header { background: #7a9a00 !important; color: #fff !important; }
/* ── Info Boxes ──────────────────────────────────────────────────────────── */
.info-box-number { font-size: 1.4em; }
.bg-green, .bg-green > .info-box-icon  { background-color: #00B400 !important; }
.bg-aqua,  .bg-aqua  > .info-box-icon  { background-color: #0095C3 !important; }
.bg-yellow,.bg-yellow > .info-box-icon { background-color: #7a9a00 !important; }
.bg-red,   .bg-red   > .info-box-icon  { background-color: #8B0000 !important; }
/* ── Aktiver Tab ─────────────────────────────────────────────────────────── */
.nav-tabs-custom > .nav-tabs > li.active { border-top-color: #00B400 !important; }
.nav-tabs > li.active > a,
.nav-tabs > li.active > a:hover,
.nav-tabs > li.active > a:focus { color: #005A46 !important; font-weight: 600; }
/* ── Content-Bereich ─────────────────────────────────────────────────────── */
.content-wrapper, .right-side { background-color: #f4f6f2 !important; }
/* ── Notifications ───────────────────────────────────────────────────────── */
.shiny-notification { font-size: 0.9em; border-left: 4px solid #00B400; }
/* ── Logo-Container Sidebar (unten) ──────────────────────────────────────── */
#sidebar-logos {
  padding: 10px 8px 8px 8px;
  display: flex;
  align-items: center;
  justify-content: space-around;
  gap: 8px;
  background: rgba(255,255,255,0.05);
  border-top: 1px solid rgba(255,255,255,0.12);
  margin-top: 6px;
}
#sidebar-logos img {
  max-height: 46px;
  max-width: 118px;
  width: auto;
  object-fit: contain;
  filter: brightness(1.15) saturate(0.9);
}
/* ── Box-Titel ───────────────────────────────────────────────────────────── */
.box-title { font-weight: 700; }
/* ── HR in Sidebar ───────────────────────────────────────────────────────── */
.main-sidebar hr { border-color: rgba(255,255,255,0.12) !important; margin: 6px 12px; }
"

dashboardPage(
  skin = "green",

  # --------------------------------------------------------------------------
  # HEADER — Logos per tags$li in die Navbar eingebettet
  # --------------------------------------------------------------------------
  dashboardHeader(
    title = tagList(
      tags$span("AFS-SCD", style = "font-weight:700; color:#AAC800;"),
      tags$span(": Case Study Explorer", style = "font-size:13px; color:#d0ead0;")
    ),
    titleWidth = 300,

    # Logos rechts in der Headerleiste
    tags$li(
      class = "dropdown",
      style = "padding: 8px 12px; display:flex; align-items:center; gap:10px;",
      tags$a(
        href = "https://www.dip-sachsen-anhalt.de", target = "_blank",
        tags$img(
          src    = "DiP_Foerderlogos_SA-300x236.png",
          height = "34px",
          alt    = "DiP Sachsen-Anhalt Foerderlogo",
          title  = "DiP Sachsen-Anhalt – Modellregion der Biooekonomie",
          style  = "max-width:130px; object-fit:contain; filter:brightness(10);"
        )
      ),
      tags$a(
        href = "https://www.bmbf.de", target = "_blank",
        tags$img(
          src    = "Logo_BFTR-1.svg",
          height = "34px",
          alt    = "BFTR Foerderlogo",
          title  = "Bundesministerium fuer Forschung, Technologie und Raumfahrt",
          style  = "max-width:130px; object-fit:contain; filter:brightness(10);"
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # SIDEBAR
  # --------------------------------------------------------------------------
  dashboardSidebar(
    width = 300,
    useShinyjs(),

    # -- Wachstumsparameter --------------------------------------------------
    h4("Biomass Growth",
       style = "color:#AAC800; margin:12px 12px 4px; font-size:12px; text-transform:uppercase; letter-spacing:0.06em;"),
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
    h4("Planning & Rotation",
       style = "color:#AAC800; margin:12px 12px 4px; font-size:12px; text-transform:uppercase; letter-spacing:0.06em;"),
    sliderInput("n_periods", "Planning horizon (years)",
                min = 20, max = 60, value = 40, step = 5),
    sliderInput("min_age",   "Min rotation age (years)",
                min = 1, max = 7, value = 3, step = 1),
    sliderInput("max_age",   "Max rotation age (years)",
                min = 8, max = 25, value = 20, step = 1),

    hr(),

    # -- Standortkosten -------------------------------------------------------
    h4("Site Economics",
       style = "color:#AAC800; margin:12px 12px 4px; font-size:12px; text-transform:uppercase; letter-spacing:0.06em;"),
    numericInput("C_est",    "Establishment cost (\u20ac/ha)",
                 value = 2500, min = 500, max = 8000, step = 100),
    numericInput("C_harv",   "Harvest cost (\u20ac/ha)",
                 value = 150, min = 20, max = 600, step = 10),
    numericInput("C_main",   "Maintenance cost (\u20ac/ha/yr)",
                 value = 10, min = 0, max = 200, step = 5),
    numericInput("opp_mean", "Mean opportunity cost (\u20ac/ha/yr)",
                 value = 500, min = 0, max = 1200, step = 25),

    hr(),

    # -- Logistikkosten -------------------------------------------------------
    h4("Logistics",
       style = "color:#AAC800; margin:12px 12px 4px; font-size:12px; text-transform:uppercase; letter-spacing:0.06em;"),
    numericInput("c_tr_raw", "Transport cost raw (\u20ac/t\u00b7km)",
                 value = 0.08, min = 0.01, max = 0.50, step = 0.01),
    numericInput("c_tr_pre", "Transport cost chips (\u20ac/t\u00b7km)",
                 value = 0.06, min = 0.01, max = 0.50, step = 0.01),

    hr(),

    # -- Revenue Multipliers -------------------------------------------------
    h4("Revenue Multipliers",
       style = "color:#AAC800; margin:12px 12px 4px; font-size:12px; text-transform:uppercase; letter-spacing:0.06em;"),
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
                 style = paste0("background:#00B400; color:white; width:90%; margin:6px 5%;",
                                "font-weight:600; border:none; border-radius:4px;")),
    br(),
    downloadButton("export_csv", "Export Results",
                   style = paste0("width:90%; margin:6px 5%; background:#005A46;",
                                  "color:white; border:none; border-radius:4px;")),
    br(),
    tags$small(
      style = "color:#88aa88; margin:8px 12px; display:block;",
      "Results load from precomputed RDS when available."
    ),

    # -- Logos am Ende der Sidebar --------------------------------------------
    tags$div(
      id = "sidebar-logos",
      tags$a(
        href = "https://www.dip-sachsen-anhalt.de", target = "_blank",
        tags$img(
          src   = "DiP_Foerderlogos_SA-300x236.png",
          alt   = "DiP Sachsen-Anhalt",
          title = "DiP Sachsen-Anhalt – Modellregion der Biooekonomie"
        )
      ),
      tags$a(
        href = "https://www.bmbf.de", target = "_blank",
        tags$img(
          src   = "Logo_BFTR-1.svg",
          alt   = "BFTR",
          title = "Bundesministerium fuer Forschung, Technologie und Raumfahrt"
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # BODY
  # --------------------------------------------------------------------------
  dashboardBody(
    tags$head(
      tags$style(HTML(dip_css))
    ),

    tabsetPanel(
      id = "main_tabs",

      # ====================================================================
      # TAB 1: Wachstumsmodell
      # ====================================================================
      tabPanel(
        title = tagList(icon("seedling"), "Biomass Growth"),
        value = "tab_growth",

        fluidRow(
          infoBoxOutput("kpi_max_yield",     width = 3),
          infoBoxOutput("kpi_stem_share_10", width = 3),
          infoBoxOutput("kpi_optimal_age",   width = 3),
          infoBoxOutput("kpi_fresh_weight",  width = 3)
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
            collapsible = TRUE, collapsed = FALSE,
            plotlyOutput("plot_growth_asymptotes", height = 300)
          ),
          box(
            title = "Biomass fraction growth",
            status = "primary", solidHeader = TRUE, width = 6,
            collapsible = TRUE, collapsed = FALSE,
            plotlyOutput("plot_growth_fractions", height = 300)
          )
        )
      ),

      # ====================================================================
      # TAB 2: Netzwerk & Karte
      # ====================================================================
      tabPanel(
        title = tagList(icon("map-marked-alt"), "Network & Map"),
        value = "tab_network",

        fluidRow(
          infoBoxOutput("kpi_n_sites",    width = 3),
          infoBoxOutput("kpi_total_area", width = 3),
          infoBoxOutput("kpi_obj_val",    width = 3),
          infoBoxOutput("kpi_solver_gap", width = 3)
        ),

        fluidRow(
          box(
            title = "Supply Chain Network \u2014 Active Sites & Flows",
            status = "primary", solidHeader = TRUE, width = 8,
            helpText("Green markers = active AFS sites; orange = hubs/storages;
                     purple = industrial consumers. Line width ~ flow volume."),
            leafletOutput("map_network", height = 560)
          ),
          box(
            title = "Result Summary",
            status = "success", solidHeader = TRUE, width = 4,
            h5("Solver Status"),             textOutput("txt_solver_status"), br(),
            h5("Objective Value (\u20ac)"),  textOutput("txt_obj"),           br(),
            h5("Active Sites"),              textOutput("txt_active_sites"),   br(),
            h5("Total AFS Area (ha)"),       textOutput("txt_total_area"),     br(),
            h5("Mean Profit (\u20ac/ha/yr)"),textOutput("txt_mean_profit"),    br(),
            hr(),
            h5("Legend"),
            tags$ul(
              tags$li(HTML("<span style='color:#00B400;'>&#11044;</span> Active AFS site")),
              tags$li(HTML("<span style='color:orange;'>&#9650;</span> Hub / Storage")),
              tags$li(HTML("<span style='color:purple;'>&#9632;</span> Industrial consumer")),
              tags$li(HTML("<span style='color:#0095C3;'>&#9135;</span> Raw biomass flow")),
              tags$li(HTML("<span style='color:#005A46;'>&#9135;</span> Processed flow"))
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

      # ====================================================================
      # TAB 3: Material Flows
      # ====================================================================
      tabPanel(
        title = tagList(icon("chart-area"), "Material Flows"),
        value = "tab_flows",

        fluidRow(
          box(
            title = "Alluvial: Biomass Cascade (Site \u2192 Hub \u2192 Consumer)",
            status = "info", solidHeader = TRUE, width = 12,
            helpText("Sankey diagram of total flows aggregated over planning horizon."),
            plotlyOutput("plot_sankey", height = 320)
          )
        ),

        fluidRow(
          box(
            title = "Harvested Biomass by Product over Time",
            status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_biomass_time", height = 320)
          ),
          box(
            title = "Demand Fulfilment Rate by Consumer & Product",
            status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_demand_fulfilment", height = 300)
          )
        ),

        fluidRow(
          box(
            title = "Revenue by Product & Consumer",
            status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_rev_consumer", height = 320)
          ),
          box(
            title = "Aggregated product cascade",
            status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_product_cascade", height = 320)
          )
        )
      ),

      # ====================================================================
      # TAB 4: Site KPIs
      # ====================================================================
      tabPanel(
        title = tagList(icon("map-pin"), "Site KPIs"),
        value = "tab_site_kpis",

        fluidRow(
          box(
            title = "Profit per ha\u00b7yr by Opportunity Cost",
            status = "info", solidHeader = TRUE, width = 6,
            helpText("Each dot = one AFS site. X-axis = site-specific opportunity cost (\u20ac/ha/yr).
                     Dashed line = break-even."),
            plotlyOutput("plot_rotation_dist", height = 380)
          ),
          box(
            title = "Profit split per site",
            status = "success", solidHeader = TRUE, width = 6,
            helpText("Comparison of revenues and cost split by sub-types."),
            plotlyOutput("plot_profit_split", height = 380)
          )
        ),

        fluidRow(
          box(
            title = "Opportunity cost vs. Profit",
            status = "success", solidHeader = TRUE, width = 4,
            plotlyOutput("plot_profit_vs_opp", height = 300)
          ),
          box(
            title = "P1 (Stem) Share vs. Profit",
            status = "success", solidHeader = TRUE, width = 4,
            plotlyOutput("plot_p1share_profit", height = 300)
          ),
          box(
            title = "Distance vs. Profit",
            status = "success", solidHeader = TRUE, width = 4,
            plotlyOutput("plot_dist_profit", height = 300)
          )
        )
      )

    ) # end tabsetPanel
  )   # end dashboardBody
)     # end dashboardPage
