# =============================================================================
# AFS Poplar Biomass Model — Setup chunk for paper.Rmd
# Sources: Jha 2018, Civitarese 2019, Niemczyk 2021, Huber 2017, Graves 2007
# =============================================================================

# ── 1. Gompertz growth model ─────────────────────────────────────────────────
#' Total dry biomass per ha of AFS system (t DM/ha)
#' @param t  Stand age (years)
#' @param A  Asymptote (t DM/ha); site- and density-specific
#' @param k  Growth rate coefficient (yr^-1); default 0.194 (central EU)
#' @param t0 Inflection age (yr); default 9.7
gomp <- function(t, A, k = 0.194, t0 = 9.7) {
  A * exp(-exp(-k * (t - t0)))
}

# ── 2. Density–asymptote relationship (power law) ────────────────────────────
#' Per-tree Gompertz asymptote as function of planting density
#' Fitted exponent beta = 0.380 from Niemczyk 2021 & calibration data
#' @param N      Planting density (trees/ha)
#' @param C_site Site constant: 4475 (conservative), 6471 (good site)
A_tree <- function(N, C_site = 4475) C_site * N^(-0.380)

#' Stand-level asymptote (t DM/ha) = A_tree * N / 1000
A_stand <- function(N, C_site = 4475) A_tree(N, C_site) * N / 1000

# Scenario presets
SCENARIOS <- list(
  conservative = list(C_site = 4475, k = 0.194, t0 = 9.7,
                      label = "Conservative (central EU, ~113 trees/ha)"),
  good_site    = list(C_site = 6471, k = 0.173, t0 = 10.1,
                      label = "Good site (~139 trees/ha, fertile soils)")
)

# ── 3. Aboveground / belowground split ───────────────────────────────────────
AGB_FRAC <- 0.89   # aboveground fraction of total DM (Jha 2018)

# ── 4. Biomass compartment fractions (age-dependent, linear) ─────────────────
#' Returns named list of AGB fractions (stem, merch_branch, residue)
#' Anchored to Civitarese 2019 (t=3,6,9) and Jha 2018 (t=13)
#' Valid for t in [3, 20] years; winter harvest assumed (leaves shed)
fractions_agb <- function(t) {
  s <- -0.00524 * t + 0.8396   # stem fraction (raw)
  m <-  0.00643 * t + 0.0562   # merchantable branch fraction (raw)
  r <- -0.00119 * t + 0.1042   # residue fraction (raw)
  tot <- s + m + r
  list(stem         = s / tot,
       merch_branch = m / tot,
       residue      = r / tot)
}

# ── 5. Diameter-class merchantable fractions (logistic) ──────────────────────
#' Fraction of stem biomass with d > 15 cm (merchantable timber)
#' Fitted to Civitarese 2019 / Niemczyk 2021; r=0.467, t50=8.19 yr
f_stem_merch <- function(t) 1 / (1 + exp(-0.467 * (t - 8.19)))

#' Fraction of branch biomass with insertion d > 7 cm (merchantable)
#' Fitted to Civitarese 2019 / Jha 2018; r=0.698, t50=7.55 yr
f_branch_merch <- function(t) 1 / (1 + exp(-0.698 * (t - 7.55)))

# ── 6. Fresh weight conversion ───────────────────────────────────────────────
MC_FELLING <- 0.55   # moisture content at felling (Civitarese 2019)
#' Convert dry matter to fresh weight
#' @param dm_t  Dry matter (t)
dm_to_fresh <- function(dm_t, mc = MC_FELLING) dm_t / (1 - mc)

# ── 7. Master harvest function ───────────────────────────────────────────────
#' Compute full harvest breakdown for an AFS stand
#'
#' @param t         Harvest age (years)
#' @param N         Planting density (trees/ha)
#' @param area_ha   Total stand area (ha)
#' @param C_site    Site constant (4475 = conservative, 6471 = good site)
#' @param k,t0      Gompertz shape parameters
#' @return data.frame with DM and fresh weight for each compartment (total stand)
harvest_biomass <- function(t, N, area_ha = 1,
                            C_site = 4475, k = 0.194, t0 = 9.7) {
  A        <- A_stand(N, C_site)
  total_dm <- gomp(t, A, k, t0) * area_ha       # t DM total (incl. roots)
  agb_dm   <- total_dm * AGB_FRAC               # t DM aboveground
  
  fr       <- fractions_agb(t)
  stem_dm   <- agb_dm * fr$stem
  branch_dm <- agb_dm * fr$merch_branch
  resid_dm  <- agb_dm * fr$residue
  
  # Apply diameter-class merchantability logistic
  stem_merch_dm   <- stem_dm   * f_stem_merch(t)
  stem_resid_dm   <- stem_dm   * (1 - f_stem_merch(t))
  branch_merch_dm <- branch_dm * f_branch_merch(t)
  branch_resid_dm <- branch_dm * (1 - f_branch_merch(t))
  total_resid_dm  <- resid_dm + stem_resid_dm + branch_resid_dm
  
  data.frame(
    compartment = c("Stem (d > 15 cm)",
                    "Merch. branches (d > 7 cm)",
                    "Residue / fine material"),
    dm_t    = round(c(stem_merch_dm, branch_merch_dm, total_resid_dm), 1),
    fresh_t = round(dm_to_fresh(
      c(stem_merch_dm, branch_merch_dm, total_resid_dm)), 1)
  )
}

# ── 8. Quick-reference table (Table 1 in paper) ──────────────────────────────
#' Reproduce standard Table 1: 100 ha, t=10 yr, both scenarios
table_harvest_100ha <- function(t = 10, area_ha = 100, N = 113) {
  cons <- harvest_biomass(t, N, area_ha,
                          C_site = 4475, k = 0.194, t0 = 9.7)
  good <- harvest_biomass(t, N, area_ha,
                          C_site = 6471, k = 0.173, t0 = 10.1)
  cons$scenario <- "Conservative"
  good$scenario <- "Good site"
  rbind(cons, good)[, c("scenario","compartment","dm_t","fresh_t")]
}

# ── 9. Parameter overview (for inline reporting in text) ─────────────────────
MODEL_PARAMS <- list(
  A_cons      = A_stand(113, 4475),  # t DM/ha asymptote, conservative
  A_good      = A_stand(139, 6471),  # t DM/ha asymptote, good site
  beta_density = 0.380,              # competition exponent
  k_cons      = 0.194,               # Gompertz k, conservative
  k_good      = 0.173,               # Gompertz k, good site
  t0_cons     = 9.7,                 # inflection yr, conservative
  t0_good     = 10.1,                # inflection yr, good site
  r_stem      = 0.467,               # logistic r, stem d > 15 cm
  t50_stem    = 8.19,                # logistic t50, stem d > 15 cm
  r_branch    = 0.698,               # logistic r, branch d > 7 cm
  t50_branch  = 7.55,                # logistic t50, branch d > 7 cm
  agb_frac    = AGB_FRAC,
  mc_felling  = MC_FELLING
)