# ============================================================================
# TEST FILE 3: COMPLEX OBSERVATION MODELS COMPARISON
# ============================================================================
# Purpose: Test complex and more realistic monitoring scenarios with irregular patterns
# Focus Areas:
# 1. Area-specific observation patterns (some indices only in Area 1 vs Area 2)
# 2. Irregular temporal coverage (gaps, sporadic sampling)
# 3. Fleet-specific catch reporting
# 4. Mixed observation errors
# 5. Realistic survey timing variations

 rm(list=ls())
 devtools::load_all()
 library(ggplot2)
 library(dplyr)
 library(tidyr)

# ============================================================================
# SHARED SETUP FOR ALL TESTS
# ============================================================================

cat("=====================================\n")
cat("FISHSIMGTG MULTIFLEET VALIDATION - FILE 3\n")
cat("Complex Spatial-Temporal Observation Models\n")
cat("=====================================\n")

# Life History - consistent across all tests
lh_obj <- new("LifeHistory")
lh_obj@title <- "Kole"
lh_obj@speciesName <- "Ctenochaetus strigosus"
lh_obj@Linf <- 17.7
lh_obj@K <- 0.423
lh_obj@t0 <- -0.51
lh_obj@L50 <- 8.4
lh_obj@L95delta <- 1.2
lh_obj@M <- 0.08
lh_obj@L_type<-"FL"
lh_obj@L_units<-"cm"
lh_obj@LW_A <- 0.046
lh_obj@LW_B <- 2.85
lh_obj@Steep <- 0.54
lh_obj@recSD <- 0  # Deterministic for comparison
lh_obj@recRho <- 0
lh_obj@isHermaph<-FALSE
lh_obj@R0<-10000

# Time Area
ta <- new("TimeArea")
ta@title <- "Validation Test"
ta@gtg <- 13
ta@areas <- 2
ta@recArea <- c(0.99, 0.01)
ta@iterations <- 3
ta@historicalYears <- 12
ta@historicalBio <- 0.5
ta@historicalBioType <- "relB"
ta@move <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = FALSE)

# Historical effort
ta@historicalEffort <- matrix(c(
  # Area 1: declining then stabilizing
  1.8, 1.6, 1.4, 1.2, 1.0, 0.9, 0.8, 0.8, 0.7, 0.7, 0.6, 0.6,
  # Area 2: variable pattern
  1.2, 1.5, 1.1, 1.3, 0.9, 1.1, 0.8, 1.0, 0.7, 0.9, 0.6, 0.8
), nrow = 12, ncol = 2, byrow = FALSE)

# Stochastic
stochastic_obj <- new("Stochastic")
stochastic_obj@historicalBio <- c(0.45, 0.65)
stochastic_obj@Steep <- c(0.60, 0.80)


# ============================================================================
# DEFINE FISHERY SELECTIVITIES
# ============================================================================

# Fleet 1
fleet1_sel_hist <- new("Fishery")
fleet1_sel_hist@title <- "Fleet 1 - Large Fish"
fleet1_sel_hist@vulType <- "logistic"
fleet1_sel_hist@vulParams <- c(12.0, 0.08)
fleet1_sel_hist@retType <- "full"
fleet1_sel_hist@retMax <- 1
fleet1_sel_hist@Dmort <- 0

fleet1_sel_proj <- new("Fishery")
fleet1_sel_proj@title <- "Fleet 1 - Large Fish"
fleet1_sel_proj@vulType <- "logistic"
fleet1_sel_proj@vulParams <- c(12.5, 0.08)
fleet1_sel_proj@retType <- "full"
fleet1_sel_proj@retMax <- 1
fleet1_sel_proj@Dmort <- 0

# Fleet 2
fleet2_sel_hist <- new("Fishery")
fleet2_sel_hist@title <- "Fleet 2 - Medium Fish"
fleet2_sel_hist@vulType <- "logistic"
fleet2_sel_hist@vulParams <- c(9.5, 0.12)
fleet2_sel_hist@retType <- "full"
fleet2_sel_hist@retMax <- 1
fleet2_sel_hist@Dmort <- 0

fleet2_sel_proj <- new("Fishery")
fleet2_sel_proj@title <- "Fleet 2 - Medium Fish"
fleet2_sel_proj@vulType <- "logistic"
fleet2_sel_proj@vulParams <- c(9.8, 0.12)
fleet2_sel_proj@retType <- "full"
fleet2_sel_proj@retMax <- 1
fleet2_sel_proj@Dmort <- 0

# Fleet 3
fleet3_sel_hist <- new("Fishery")
fleet3_sel_hist@title <- "Fleet 3 - Small Fish"
fleet3_sel_hist@vulType <- "logistic"
fleet3_sel_hist@vulParams <- c(7.5, 0.15)
fleet3_sel_hist@retType <- "full"
fleet3_sel_hist@retMax <- 1
fleet3_sel_hist@Dmort <- 0

fleet3_sel_proj <- new("Fishery")
fleet3_sel_proj@title <- "Fleet 3 - Small Fish"
fleet3_sel_proj@vulType <- "logistic"
fleet3_sel_proj@vulParams <- c(7.8, 0.15)
fleet3_sel_proj@retType <- "full"
fleet3_sel_proj@retMax <- 1
fleet3_sel_proj@Dmort <- 0


# ============================================================================
# DEFINE SURVEY SELECTIVITY - 4 INDEPENDENT SURVEYS
# ============================================================================

# Survey A
surveyA_sel_hist <- new("Fishery")
surveyA_sel_hist@title <- "Survey A - Juveniles"
surveyA_sel_hist@vulType <- "logistic"
surveyA_sel_hist@vulParams <- c(6.0, 0.20)
surveyA_sel_hist@retType <- "full"
surveyA_sel_hist@retMax <- 1
surveyA_sel_hist@Dmort <- 0

surveyA_sel_proj <- new("Fishery")
surveyA_sel_proj@title <- "Survey A - Juveniles"
surveyA_sel_proj@vulType <- "logistic"
surveyA_sel_proj@vulParams <- c(6.2, 0.20)
surveyA_sel_proj@retType <- "full"
surveyA_sel_proj@retMax <- 1
surveyA_sel_proj@Dmort <- 0

# Survey B
surveyB_sel_hist <- new("Fishery")
surveyB_sel_hist@title <- "Survey B - Adults"
surveyB_sel_hist@vulType <- "logistic"
surveyB_sel_hist@vulParams <- c(11.0, 0.10)
surveyB_sel_hist@retType <- "full"
surveyB_sel_hist@retMax <- 1
surveyB_sel_hist@Dmort <- 0

surveyB_sel_proj <- new("Fishery")
surveyB_sel_proj@title <- "Survey B - Adults"
surveyB_sel_proj@vulType <- "logistic"
surveyB_sel_proj@vulParams <- c(11.3, 0.10)
surveyB_sel_proj@retType <- "full"
surveyB_sel_proj@retMax <- 1
surveyB_sel_proj@Dmort <- 0

# Survey C
surveyC_sel_hist <- new("Fishery")
surveyC_sel_hist@title <- "Survey C - Broad"
surveyC_sel_hist@vulType <- "logistic"
surveyC_sel_hist@vulParams <- c(9.0, 0.25)
surveyC_sel_hist@retType <- "full"
surveyC_sel_hist@retMax <- 1
surveyC_sel_hist@Dmort <- 0

surveyC_sel_proj <- new("Fishery")
surveyC_sel_proj@title <- "Survey C - Broad"
surveyC_sel_proj@vulType <- "logistic"
surveyC_sel_proj@vulParams <- c(9.2, 0.25)
surveyC_sel_proj@retType <- "full"
surveyC_sel_proj@retMax <- 1
surveyC_sel_proj@Dmort <- 0

# Survey D
surveyD_sel_hist <- new("Fishery")
surveyD_sel_hist@title <- "Survey D - Research"
surveyD_sel_hist@vulType <- "logistic"
surveyD_sel_hist@vulParams <- c(8.5, 0.18)
surveyD_sel_hist@retType <- "full"
surveyD_sel_hist@retMax <- 1
surveyD_sel_hist@Dmort <- 0

surveyD_sel_proj <- new("Fishery")
surveyD_sel_proj@title <- "Survey D - Research"
surveyD_sel_proj@vulType <- "logistic"
surveyD_sel_proj@vulParams <- c(8.7, 0.18)
surveyD_sel_proj@retType <- "full"
surveyD_sel_proj@retMax <- 1
surveyD_sel_proj@Dmort <- 0

# Shared random seed
test_seed <- 12345


# ============================================================================
# CREATE MULTIFLEET OBJECT - 3 FLEETS WITH DIFFERENT SPATIAL PATTERNS
# ============================================================================

multifleet_complex <- new("Multifleet")
multifleet_complex@nfleets <- 3
multifleet_complex@fleet_proportions <- c(0.5, 0.3, 0.2)  # Fleet 1 is the dominant
multifleet_complex@allocation_type <- "catch"
multifleet_complex@fleet_selectivity_hist_list <- list(fleet1_sel_hist, fleet2_sel_hist, fleet3_sel_hist)
#multifleet_complex@fleet_selectivity_proj_list <- list(fleet1_sel_proj, fleet2_sel_proj, fleet3_sel_proj)

# New structure: For 2 areas, 2 fleets (ask bill)
 multifleet_complex@fleet_selectivity_proj_list <- list(
  # Area 1
  list(fleet1_sel_proj, fleet2_sel_proj, fleet3_sel_proj),
  # Area 2
  list(fleet1_sel_proj, fleet2_sel_proj, fleet3_sel_proj)
 )

#multifleet_complex@fleet_selectivity_proj_list <- list(fleet1_sel_proj, fleet2_sel_proj, fleet3_sel_proj)




#adding the array of fleet historical eefort
multifleet_complex@fleet_historicalEffort <- array(dim = c(ta@historicalYears, ta@areas, 3))
multifleet_complex@fleet_historicalEffort[,,1] <- ta@historicalEffort
multifleet_complex@fleet_historicalEffort[,,2] <- ta@historicalEffort
multifleet_complex@fleet_historicalEffort[,,3] <- ta@historicalEffort

# ============================================================================
# MANAGEMENT STRATEGIES - MULTI FLEET WITH OBS MODELS
# ============================================================================

complexMP <- function(phase, dataObject) {
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase == 1) {
    combined_data <- list()

    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject)
      for(col_name in names(index_result)) {
        combined_data[[col_name]] <- index_result[[col_name]]
      }
    }

    if(!is.null(CatchObsObj)) {
      catch_result <- calculate_single_CatchObs(dataObject)
      for(col_name in names(catch_result)) {
        combined_data[[col_name]] <- catch_result[[col_name]]
      }
    }

    if(!is.null(LengthCompObj)) {
      lc_result <- calculate_single_LengthComp(dataObject)
      for(col_name in names(lc_result)) {
        combined_data[[col_name]] <- lc_result[[col_name]]
      }
    }

    return(combined_data)
  }

  if(phase == 2) return(list())
  if(phase == 3) {

    Flocal <- data.frame()
    #creating one row per area-fleet combination
    for(m in 1:areas) {
      for(f in 1:nfleets) {
        #row contain: [year, iteration, area, fleet, F_value]
        Flocal <- rbind(Flocal, c(j, k, m, f, 0.05))  # Conservative F = 0.05
      }
    }



    #results in 4 rows: (1,1), (1,2), (2,1), (2,2) for 2 areas × 2 fleets
    #                   (A1,F1) (A1,F2)

    return(list(
      year = Flocal[,1],           # [j, j, j, j]
      iteration = Flocal[,2],      # [k, k, k, k]
      area = Flocal[,3],           # [1, 1, 2, 2]
      fleet = Flocal[,4],          # [1, 2, 1, 2]
      Flocal = Flocal[,5]          # [0.05, 0.05, 0.05, 0.05]
    ))
  }
}


strategy_complex <- new("Strategy")
strategy_complex@title <- "Complex Spatial-Temporal Strategy"
strategy_complex@projectionYears <- 7
strategy_complex@projectionName <- "complexMP"
strategy_complex@projectionParams <- list()

# ============================================================================
# COMPREHENSIVE COMPLEX INDEX OBJECT
# ============================================================================

complex_indices <- new("Index")
complex_indices@indexID <- "Complex_Spatial_Temporal"
complex_indices@title <- "Complex Spatial-Temporal Indices"
complex_indices@useWeight <- TRUE

#All survey selectivities
complex_indices@selectivity_hist_list <- list(surveyA_sel_hist, surveyB_sel_hist,
                                              surveyC_sel_hist, surveyD_sel_hist)
complex_indices@selectivity_proj_list <- list(surveyA_sel_proj, surveyB_sel_proj,
                                              surveyC_sel_proj, surveyD_sel_proj)

total_years <- ta@historicalYears + strategy_complex@projectionYears  # 12 + 7 = 19 years

complex_indices@survey_design <- list(

  # 1. Survey A (Juvenile): FI, Area 1 only, irregular sampling (years)
  list(
    indextype = "FI",
    areas = c(1),
    indexYears = c(1, 3, 5, 8, 11, 14, 17),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.2,  # occur early in year
    q_hist_bounds = c(0.0003, 0.0008),
    q_proj_bounds = c(0.0004, 0.0010),
    hyperstability_hist_bounds = c(1.0, 1.0),
    hyperstability_proj_bounds = c(1.0, 1.0),
    obsError_CV_hist_bounds = c(0.15, 0.25),
    obsError_CV_proj_bounds = c(0.10, 0.20)
  ),

  # 2. Survey B (Adult): FI, Area 2 only, every 3 years
  list(
    indextype = "FI",
    areas = c(2),
    indexYears = seq(2, total_years, 3),  # years 2, 5, 8, 11, 14, 17
    selectivity_hist_idx = 2,
    selectivity_proj_idx = 2,
    survey_timing = 0.5,  # occur mid-year
    q_hist_bounds = c(0.0001, 0.0004),
    q_proj_bounds = c(0.0002, 0.0005),
    hyperstability_hist_bounds = c(0.9, 1.1),
    hyperstability_proj_bounds = c(0.9, 1.1),
    obsError_CV_hist_bounds = c(0.08, 0.15),
    obsError_CV_proj_bounds = c(0.08, 0.15)
  ),

  # 3. Survey C (Broad): FI, Both areas, but with gaps in years
  list(
    indextype = "FI",
    areas = c(1, 2),
    indexYears = c(1, 2, 3, 6, 7, 8, 11, 15, 16, 19),
    selectivity_hist_idx = 3,
    selectivity_proj_idx = 3,
    survey_timing = 0.8,  # occur late in year
    q_hist_bounds = c(0.0002, 0.0006),
    q_proj_bounds = c(0.0003, 0.0008),
    hyperstability_hist_bounds = c(0.85, 1.15),
    hyperstability_proj_bounds = c(0.85, 1.15),
    obsError_CV_hist_bounds = c(0.12, 0.22),
    obsError_CV_proj_bounds = c(0.12, 0.22)
  ),

  # 4. Survey D (Research survey): FI, Area 1 and Area 2
  list(
    indextype = "FI",
    areas = c(1, 2),  # occur in both areas
    indexYears = c(4, 6, 9, 12, 15, 18),  # research years
    selectivity_hist_idx = 4,
    selectivity_proj_idx = 4,
    survey_timing = 0.6,  # after mid year
    q_hist_bounds = c(0.00015, 0.0005),
    q_proj_bounds = c(0.0002, 0.0007),
    hyperstability_hist_bounds = c(0.95, 1.05),
    hyperstability_proj_bounds = c(0.95, 1.05),
    obsError_CV_hist_bounds = c(0.10, 0.18),
    obsError_CV_proj_bounds = c(0.10, 0.18)
  ),

  # 5. Fleet 1 CPUE: FD, Area 1 annual
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    indexYears = 1:total_years,  # complete annual coverage
    q_hist_bounds = c(0.0001, 0.0003),
    q_proj_bounds = c(0.00015, 0.0004),
    hyperstability_hist_bounds = c(0.8, 1.2),
    hyperstability_proj_bounds = c(0.8, 1.2),
    obsError_CV_hist_bounds = c(0.20, 0.35),
    obsError_CV_proj_bounds = c(0.15, 0.30)
  ),

  # 6. Fleet 1 CPUE: FD, Area 2, sporadic cpue data
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(2),
    indexYears = c(1, 4, 7, 9, 13, 16, 19),
    q_hist_bounds = c(0.00008, 0.0002),
    q_proj_bounds = c(0.0001, 0.0003),
    hyperstability_hist_bounds = c(0.7, 1.3),
    hyperstability_proj_bounds = c(0.7, 1.3),
    obsError_CV_hist_bounds = c(0.25, 0.40),
    obsError_CV_proj_bounds = c(0.20, 0.35)
  ),

  # 7. Fleet 2 CPUE: FD, Both areas, every other year cpue data
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(1, 2),
    indexYears = seq(1, total_years, 2),
    q_hist_bounds = c(0.00012, 0.0004),
    q_proj_bounds = c(0.0002, 0.0005),
    hyperstability_hist_bounds = c(0.75, 1.25),
    hyperstability_proj_bounds = c(0.75, 1.25),
    obsError_CV_hist_bounds = c(0.18, 0.32),
    obsError_CV_proj_bounds = c(0.18, 0.32)
  ),

  # 8. Fleet 3 CPUE: FD, Area 2, poor projection sampling
  list(
    indextype = "FD",
    fleet_id = 3,
    areas = c(2),
    indexYears = c(1:12, 15, 18),
    q_hist_bounds = c(0.0002, 0.0007),
    q_proj_bounds = c(0.0003, 0.0009),
    hyperstability_hist_bounds = c(0.6, 1.4),
    hyperstability_proj_bounds = c(0.6, 1.4),
    obsError_CV_hist_bounds = c(0.22, 0.38),
    obsError_CV_proj_bounds = c(0.30, 0.50)
  ),

  # 9. Fleet 3 CPUE: FD, Area 1 (expansion fishery testing)
  list(
    indextype = "FD",
    fleet_id = 3,
    areas = c(1),
    indexYears = c(6, 10, 14, 17), # only few years
    q_hist_bounds = c(0.00015, 0.0005),
    q_proj_bounds = c(0.0002, 0.0006),
    hyperstability_hist_bounds = c(0.5, 1.5),
    hyperstability_proj_bounds = c(0.5, 1.5),
    obsError_CV_hist_bounds = c(0.30, 0.50),
    obsError_CV_proj_bounds = c(0.35, 0.55)
  )
)

cat("\nCreated complex index object with 9 indices:\n")
cat("  - 4 FI Surveys with different spatial-temporal patterns\n")
cat("  - 5 Fleet-specific CPUE indices with typical constraints\n")
cat("  - Irregular temporal coverage simulating real-world scenarios\n")

# ============================================================================
# COMPLEX FLEET-SPECIFIC CATCH OBSERVATIONS
# ============================================================================

complex_catch <- new("CatchObs")
complex_catch@catchID <- "Complex_Fleet_Catch"
complex_catch@title <- "Complex Fleet-Specific Catch"

complex_catch@fleet_configs <- list(
  # Fleet 1: catch data in from area 1
  list(
    fleet_id = 1,
    areas = c(1),
    catchYears = 1:total_years,

    reporting_rates = c(
      rep(0.95, ta@historicalYears),  # Historical
      rep(0.98, strategy_complex@projectionYears)  # Projection improvement
    ),
    obs_CVs = matrix(
      c(c(rep(0.08, ta@historicalYears), rep(0.06, strategy_complex@projectionYears)),  # Low bounds
        c(rep(0.15, ta@historicalYears), rep(0.12, strategy_complex@projectionYears))), # Upper bounds
      ncol = 2
    )
  ),

  # Fleet 2: catch data in from area 1 and 2
  list(
    fleet_id = 2,
    areas = c(1, 2),
    catchYears = 1:total_years,
    reporting_rates = c(
      rep(0.80, ta@historicalYears),
      rep(0.85, strategy_complex@projectionYears)
    ),
    obs_CVs = matrix(
      c(c(rep(0.15, ta@historicalYears), rep(0.12, strategy_complex@projectionYears)),
        c(rep(0.28, ta@historicalYears), rep(0.25, strategy_complex@projectionYears))),
      ncol = 2
    )
  ),

  # Fleet 3: catch data in from area 1 and 2
  list(
    fleet_id = 3,
    areas = c(1, 2),
    catchYears = c(1:total_years),
    reporting_rates = c(
      seq(0.50, 0.70, length.out = ta@historicalYears),  #gradual historical improvement
      seq(0.75, 0.90, length.out = strategy_complex@projectionYears)  #gradual projection improvement
    ),
    obs_CVs = matrix(
      c(c(seq(0.35, 0.25, length.out = ta@historicalYears),
          seq(0.22, 0.15, length.out = strategy_complex@projectionYears)),
        c(seq(0.50, 0.40, length.out = ta@historicalYears),
          seq(0.35, 0.25, length.out = strategy_complex@projectionYears))),
      ncol = 2
    )
  )
)

cat("\nCreated complex catch observations:\n")

# ============================================================================
# COMPLEX LENGTH COMPOSITION OBSERVATIONS
# ============================================================================

complex_lcomp <- new("LCompObs")
complex_lcomp@indexID <- "Complex_LComp"
complex_lcomp@title <- "Complex Length Composition"
complex_lcomp@length_bin_width <- 1
complex_lcomp@selectivity_hist_list <- list(surveyA_sel_hist, surveyB_sel_hist,
                                            surveyC_sel_hist, surveyD_sel_hist)

complex_lcomp@selectivity_proj_list <- list(surveyA_sel_proj, surveyB_sel_proj,
                                            surveyC_sel_proj, surveyD_sel_proj)

complex_lcomp@survey_design <- list(
  # 1. Fleet 1 LC: Area 1 annual
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    years = 1:total_years,
    sample_sizes = c(rep(150, ta@historicalYears), rep(200, strategy_complex@projectionYears))
  ),

  # 2. Fleet 1 LC: Area 2 sporadic
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(2),
    years = c(2, 5, 8, 11, 14, 17),
    sample_sizes = c(80, 90, 100, 110, 120, 130)
  ),

  # 3. Fleet 2 LC: Both areas, biennial
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(1, 2),
    years = seq(1, total_years, 2),
    sample_sizes = seq(100, 160, length.out = 10)
  ),

  # 4. Fleet 3 LC: Area 2, irregular sampling
  list(
    indextype = "FD",
    fleet_id = 3,
    areas = c(2),
    years = c(1, 3, 6, 8, 10, 13, 15, 18),
    sample_sizes = c(200, 180, 220, 190, 210, 230, 250, 270)
  ),

  # 5. Survey A LC: Area 1, irregular sampling
  list(
    indextype = "FI",
    areas = c(1),
    years = c(2, 5, 9, 13, 16),
    sample_sizes = c(120, 130, 140, 150, 160),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.2
  ),

  # 6. Survey B LC: Area 2, irregular sampling
  list(
    indextype = "FI",
    areas = c(2),
    years = c(3, 7, 12, 17),
    sample_sizes = c(180, 200, 220, 240),
    selectivity_hist_idx = 2,
    selectivity_proj_idx = 2,
    survey_timing = 0.5
  ),

  # 7. Survey C LC: Both areas, irregular sampling
  list(
    indextype = "FI",
    areas = c(1, 2),
    years = c(4, 8, 15, 19),
    sample_sizes = c(300, 350, 400, 450),
    selectivity_hist_idx = 3,
    selectivity_proj_idx = 3,
    survey_timing = 0.8
  ),

  # 8. Survey D LC: Research survey, specific years, large sample size
  list(
    indextype = "FI",
    areas = c(1, 2),
    years = c(6, 12, 18),
    sample_sizes = c(500, 550, 600),
    selectivity_hist_idx = 4,
    selectivity_proj_idx = 4,
    survey_timing = 0.6
  )
)

cat("\nCreated complex length composition with 8 programs:\n")
cat("  - 4 FD with different spatial patterns\n")
cat("  - 4 FI with different temporal coverage\n")
cat("  - Sample sizes ranging from 80 to 600\n")

# ============================================================================
# RUN COMPLEX SPATIAL-TEMPORAL SIMULATION
# ============================================================================
result_complex <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StrategyObj = strategy_complex,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_complex,
  IndexObj = complex_indices,
  CatchObsObj = complex_catch,
  LengthCompObj = complex_lcomp,
  wd = getwd(),
  fileName = "test3_complex_spatial_temporal",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "complexMP"
)

cat("Complex simulation completed successfully\n")


# ============================================================================
# VALIDATION
# ============================================================================


result_complex <- readProjection(getwd(), "test3_complex_spatial_temporal")

#Population dynamics
cat("population pynamics:\n")
cat("  sim years:", dim(result_complex$dynamics$SB)[1], "\n")
cat("  areas:", dim(result_complex$dynamics$SB)[3], "\n")
cat("  iter:", dim(result_complex$dynamics$SB)[2], "\n")

if(!is.null(result_complex$dynamics$multifleet)) {
  mf <- result_complex$dynamics$multifleet
  cat("  Multifleet: TRUE with", mf$nfleets, "fleets\n")
  cat("  fleet proportions achieved:", paste(round(mf$actual_catch_proportions, 3), collapse = ", "), "\n")
}


#outputs
result_complex$dynamics$SB
# result_complex$dynamics$VB
# result_complex$dynamics$Ftotal
result_complex$dynamics$SPR
result_complex$dynamics$recN
result_complex$dynamics$multifleet$Ftotal_by_fleet


#Plot simulationr esults:
#Life history
lhOut<-LHwrapper(lh_obj, ta, doPlot = TRUE)

#Fishery sel
selOut_F1_hist<-selWrapper(lh = lhOut, ta, FisheryObj = fleet1_sel_hist, doPlot = TRUE)
selOut_F1_proj<-selWrapper(lh = lhOut, ta, FisheryObj = fleet1_sel_proj, doPlot = TRUE)

selOut_F2_hist<-selWrapper(lh = lhOut, ta, FisheryObj = fleet2_sel_hist, doPlot = TRUE)
selOut_F2_proj<-selWrapper(lh = lhOut, ta, FisheryObj = fleet2_sel_hist, doPlot = TRUE)

selOut_F2_hist<-selWrapper(lh = lhOut, ta, FisheryObj = fleet3_sel_hist, doPlot = TRUE)
selOut_F2_proj<-selWrapper(lh = lhOut, ta, FisheryObj = fleet3_sel_proj, doPlot = TRUE)

#Survey sel
selOut_A_hist<-selWrapper(lh = lhOut, ta, FisheryObj = surveyA_sel_hist, doPlot = TRUE)
selOut_A_proj<-selWrapper(lh = lhOut, ta, FisheryObj = surveyA_sel_proj, doPlot = TRUE)

selOut_B_hist<-selWrapper(lh = lhOut, ta, FisheryObj = surveyB_sel_hist, doPlot = TRUE)
selOut_B_proj<-selWrapper(lh = lhOut, ta, FisheryObj = surveyB_sel_proj, doPlot = TRUE)

selOut_C_hist<-selWrapper(lh = lhOut, ta, FisheryObj = surveyC_sel_hist, doPlot = TRUE)
selOut_C_proj<-selWrapper(lh = lhOut, ta, FisheryObj = surveyC_sel_proj, doPlot = TRUE)

selOut_D_hist<-selWrapper(lh = lhOut, ta, FisheryObj = surveyD_sel_hist, doPlot = TRUE)
selOut_D_proj<-selWrapper(lh = lhOut, ta, FisheryObj = surveyD_sel_proj, doPlot = TRUE)

# loading plot fucntion for toher plots

# source("fishSimGTG-plot-functions.R")

#obs: need to add units
plot_SB(result_complex)
plot_SB(result_complex, areas=1)
plot_SB(result_complex, areas=2)
plot_SB(result_complex, areas=c(1,2))


plot_catchB(result_complex)
plot_catchB(result_complex,areas=1)
plot_catchB(result_complex,areas=2)

plot_catchN(result_complex)
plot_catchN(result_complex, areas=1)
plot_catchN(result_complex, areas=2)

plot_discB(result_complex)
plot_discB(result_complex,areas=1)
plot_discB(result_complex,areas=2)

plot_discN(result_complex)
plot_discN(result_complex,areas=1)
plot_discN(result_complex,areas=2)

plot_catchB_multi(result_complex)
plot_catchB_multi(result_complex, areas=1)
plot_catchB_multi(result_complex, areas=2)

plot_catchN_multi(result_complex,show_individual = TRUE)
plot_catchN_multi(result_complex, areas=1)
plot_catchN_multi(result_complex, areas=2)

plot_SPR(result_complex)
plot_recN(result_complex)

result_complex$dynamics$multifleet$Ftotal_by_fleet

plot_Ftotal_multi(result_complex,areas=c(1,2))
plot_Ftotal_multi(result_complex,areas=1)
plot_Ftotal_multi(result_complex,areas=2)


#plot obs models (indices)
plot_survey_indices(result_complex)
plot_cpue_indices(result_complex)

plot_all_indices(result_complex)

#plot individual indices
plot_indices(result_complex,
             index_pattern = "IDX_CPUE.*Fleet_1",
             show_individual = TRUE,
             title = "Fleet 1 CPUE Only")

plot_indices(result_complex,
             index_pattern = "IDX_CPUE.*Fleet_2",
             show_individual = TRUE,
             title = "Fleet 2 CPUE Only")

plot_indices(result_complex,
             index_pattern = "IDX_CPUE.*Fleet_3",
             show_individual = TRUE,
             title = "Fleet 3 CPUE Only")



plot_catch_observations_both(result_complex,show_individual = TRUE)
plot_catch_observations_multifleet(result_complex,show_individual = TRUE)


# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_complex,show_individual = TRUE)
plot_fishery_length_comp(result_complex, areas=1,show_individual = TRUE)
plot_fishery_length_comp(result_complex, areas=2,show_individual = TRUE)

plot_survey_length_comp(result_complex,show_individual = TRUE)
plot_survey_length_comp(result_complex, areas=1,show_individual = TRUE)
plot_survey_length_comp(result_complex, areas=2,show_individual = TRUE)



# NEW: Custom filtering for fleets and areas
plot_length_composition_by_area(result_complex,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets



plot_length_composition_by_area(result_complex,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets


plot_length_composition_by_area(result_complex,
                                program_pattern = "LC_Survey",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets






