# ============================================================================
# TEST FILE 2: BASIC OBSERVATION MODELS COMPARISON
# ============================================================================
# Purpose: Test observation models in single fleet vs multifleet contexts
# Tests:
# 1. Single fleet with all observation types
# 2. Multifleet (3 fleets) with comprehensive observation models
# 3. Mixed temporal coverage (some indices annual, others periodic)
# 4. Different spatial coverage (some indices single area, others multi-area)
# 5. Length composition with multiple sampling programs
# 6. Validation of selectivity independence between FI and FD

rm(list=ls())
devtools::load_all()
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# SHARED SETUP FOR ALL TESTS
# ============================================================================

cat("=====================================\n")
cat("FISHSIMGTG MULTIFLEET VALIDATION - FILE 2\n")
cat("Comprehensive Observation Models Testing\n")
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
ta@iterations <- 20
ta@historicalYears <- 10
ta@historicalBio <- 0.5
ta@historicalBioType <- "relB"
ta@move <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = FALSE)

# Historical effort - declining trend
ta@historicalEffort <- matrix(c(1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6,
                                1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6),
                              nrow = 10, ncol = 2, byrow = FALSE)

# Stochastic
stochastic_obj <- new("Stochastic")
stochastic_obj@historicalBio <- c(0.55, 0.65)
stochastic_obj@Steep <- c(0.65, 0.75)


# ============================================================================
# DEFINE FISHERY SELECTIVITIES
# ============================================================================

# Historical fishery (for the single fleet appraoch)
hist_fishery <- new("Fishery")
hist_fishery@title<-"Historical Fishery"
hist_fishery@vulType<-"logistic"
hist_fishery@vulParams<-c(10.2, 0.1)
hist_fishery@retType<-"full"
hist_fishery@retMax <- 1
hist_fishery@Dmort <- 0

# Projection fisheries (for the single fleet approach)
proj_fishery_area1 <- new("Fishery")
proj_fishery_area1@title<-"Proj Area 1"
proj_fishery_area1@vulType<-"logistic"
proj_fishery_area1@vulParams<-c(10.2, 0.1)
proj_fishery_area1@retType<-"full"
proj_fishery_area1@retMax <- 1
proj_fishery_area1@Dmort <- 0

proj_fishery_area2 <- new("Fishery")
proj_fishery_area2@title<-"Proj Area 2"
proj_fishery_area2@vulType<-"logistic"
proj_fishery_area2@vulParams<-c(10.2, 0.1)
proj_fishery_area2@retType<-"full"
proj_fishery_area2@retMax <- 1
proj_fishery_area2@Dmort <- 0

proj_fishery_list <- list(proj_fishery_area1, proj_fishery_area2)


# Survey selectivities (multiple independent surveys)
survey1_sel_hist <- new("Fishery")
survey1_sel_hist@title <- "Research Survey 1 Historical"
survey1_sel_hist@vulType <- "logistic"
survey1_sel_hist@vulParams <- c(7.5, 0.15)   # Small fish survey
survey1_sel_hist@retType <- "full"
survey1_sel_hist@retMax <- 1
survey1_sel_hist@Dmort <- 0

survey2_sel_hist <- new("Fishery")
survey2_sel_hist@title <- "Research Survey 2 Historical"
survey2_sel_hist@vulType <- "logistic"
survey2_sel_hist@vulParams <- c(9.5, 0.18)
survey2_sel_hist@retType <- "full"
survey2_sel_hist@retMax <- 1
survey2_sel_hist@Dmort <- 0

survey1_sel_proj <- new("Fishery")
survey1_sel_proj@title <- "Research Survey 1 Projection"
survey1_sel_proj@vulType <- "logistic"
survey1_sel_proj@vulParams <- c(7.8, 0.15)
survey1_sel_proj@retType <- "full"
survey1_sel_proj@retMax <- 1
survey1_sel_proj@Dmort <- 0

survey2_sel_proj <- new("Fishery")
survey2_sel_proj@title <- "Research Survey 2 Projection"
survey2_sel_proj@vulType <- "logistic"
survey2_sel_proj@vulParams <- c(9.8, 0.18)
survey2_sel_proj@retType <- "full"
survey2_sel_proj@retMax <- 1
survey2_sel_proj@Dmort <- 0

# Shared random seed
test_seed <- 12345


# ============================================================================
# TEST 1: SINGLE FLEET WITH COMPREHENSIVE OBSERVATION MODELS
# ============================================================================

# ============================================================================
# MANAGEMENT STRATEGIES - SINGLE FLEET
# ============================================================================

# Single fleet
singleCompMP <- function(phase, dataObject) {
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
    year <- rep(j, areas)
    iteration <- rep(k, areas)
    area <- 1:areas
    Flocal <- rep(0.10, areas)  # Conservative F

    return(list(year=year, iteration=iteration, area=area, Flocal=Flocal))
  }
}

strategy_single_comp <- new("Strategy")
strategy_single_comp@title <- "Single Fleet Comprehensive"
strategy_single_comp@projectionYears <- 5
strategy_single_comp@projectionName <- "singleCompMP"
strategy_single_comp@projectionParams <- list()

# ============================================================================
# COMPREHENSIVE OBSERVATION MODELS FOR SINGLE FLEET
# ============================================================================
# Comprehensive index object for single fleet
single_comprehensive_index <- new("Index")
single_comprehensive_index@indexID <- "SingleFleet_Comprehensive"
single_comprehensive_index@title <- "Single Fleet Comprehensive Indices"
single_comprehensive_index@useWeight <- TRUE

# Survey selectivities
single_comprehensive_index@selectivity_hist_list <- list(survey1_sel_hist, survey2_sel_hist)
single_comprehensive_index@selectivity_proj_list <- list(survey1_sel_proj, survey2_sel_proj)


single_comprehensive_index@survey_design <- list(
  # 1. FD CPUE - Data covers both areas, annual
  list(
    indextype = "FD",
    areas = c(1, 2),
    indexYears = 1:15,  # All years
    q_hist_bounds = c(0.0001, 0.0002),
    q_proj_bounds = c(0.00015, 0.00025),
    hyperstability_hist_bounds = c(0.9, 1.1),
    hyperstability_proj_bounds = c(0.9, 1.1),
    obsError_CV_hist_bounds = c(0.15, 0.25),
    obsError_CV_proj_bounds = c(0.15, 0.25)
  ),

  # 2. FD CPUE - Area 1 only, biennial
  list(
    indextype = "FD",
    areas = c(1),
    indexYears = seq(1, 15, 2),  # Every other year
    q_hist_bounds = c(0.00015, 0.0003),
    q_proj_bounds = c(0.0002, 0.0004),
    hyperstability_hist_bounds = c(0.8, 1.2),
    hyperstability_proj_bounds = c(0.8, 1.2),
    obsError_CV_hist_bounds = c(0.20, 0.35),
    obsError_CV_proj_bounds = c(0.20, 0.35)
  ),

  # 3. FI Survey 1 - Both areas, every 3 years
  list(
    indextype = "FI",
    areas = c(1, 2),
    indexYears = seq(1, 15, 3),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.3,  # early in year
    q_hist_bounds = c(0.0002, 0.0005),
    q_proj_bounds = c(0.0003, 0.0006),
    hyperstability_hist_bounds = c(1.0, 1.0),
    hyperstability_proj_bounds = c(1.0, 1.0),
    obsError_CV_hist_bounds = c(0.10, 0.20),
    obsError_CV_proj_bounds = c(0.10, 0.20)
  ),

  # 4. FI Survey 2 - Area 2 only
  list(
    indextype = "FI",
    areas = c(2),
    indexYears = c(4, 8, 12, 15),  # Specific years
    selectivity_hist_idx = 2,
    selectivity_proj_idx = 2,
    survey_timing = 0.7,  # Late in year
    q_hist_bounds = c(0.00025, 0.0008),
    q_proj_bounds = c(0.0004, 0.001),
    hyperstability_hist_bounds = c(0.95, 1.05),
    hyperstability_proj_bounds = c(0.95, 1.05),
    obsError_CV_hist_bounds = c(0.15, 0.30),
    obsError_CV_proj_bounds = c(0.15, 0.30)
  )
)


# Add this verification code right after defining single_comprehensive_index@survey_design
cat("=== SURVEY DESIGN VERIFICATION ===\n")
for(i in 1:length(single_comprehensive_index@survey_design)) {
  design <- single_comprehensive_index@survey_design[[i]]
  cat("Survey", i, ":\n")
  cat("  indextype:", design$indextype, "\n")
  cat("  areas:", paste(design$areas, collapse = "_"), "\n")
  if(design$indextype == "FI") {
    cat("  selectivity_hist_idx:", design$selectivity_hist_idx, "\n")
    cat("  survey_timing:", design$survey_timing, "\n")
  }
  cat("\n")
}
cat("=====================================\n")

# Comprehensive catch observations for single fleet
single_comprehensive_catch <- new("CatchObs")
single_comprehensive_catch@catchID <- "SingleFleet_Comprehensive_Catch"
single_comprehensive_catch@title <- "Single Fleet Comprehensive Catch"
single_comprehensive_catch@areas <- c(1, 2)
single_comprehensive_catch@catchYears <- 1:15
# Variable reporting rates - improving over time
single_comprehensive_catch@reporting_rates <- c(
  seq(0.8, 1.0, length.out = 10),  # Historical improvement
  rep(1.0, 5)  # Perfect in projection
)
# Variable observation error - decreasing over time
single_comprehensive_catch@obs_CVs <- matrix(
  c(c(seq(0.3, 0.2, length.out = 10), rep(0.15, 5)),   # Lower bounds
    c(seq(0.45, 0.3, length.out = 10), rep(0.25, 5))), # Upper bounds
  ncol = 2
)

# Comprehensive length composition for single fleet
single_comprehensive_lcomp <- new("LCompObs")
single_comprehensive_lcomp@indexID <- "SingleFleet_Comprehensive_LComp"
single_comprehensive_lcomp@title <- "Single Fleet Comprehensive Length Comp"
single_comprehensive_lcomp@length_bin_width <- 1
single_comprehensive_lcomp@selectivity_hist_list <- list(survey1_sel_hist, survey2_sel_hist)
single_comprehensive_lcomp@selectivity_proj_list <- list(survey1_sel_proj, survey2_sel_proj)

# Calculate years for length composition sampling
fishery_lc_years <- seq(2, 15, 2)  # every other year fishery sampling
survey1_lc_years <- seq(1, 15, 3)  # every 3 years survey 1
survey2_lc_years <- c(4, 8, 12, 15) # Irregular survey 2

single_comprehensive_lcomp@survey_design <- list(
  # 1. Fishery length composition
  list(
    indextype = "FD",
    areas = c(1, 2),
    years = fishery_lc_years,
    sample_sizes = seq(80, 200, length.out = length(fishery_lc_years))
  ),

  # 2. Survey 1 length composition
  list(
    indextype = "FI",
    areas = c(1, 2),
    years = survey1_lc_years,
    sample_sizes = c(150, 180,100, 195, 200),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.3
  ),

  # 3. Survey 2 length composition - Area 1 only
  list(
    indextype = "FI",
    areas = c(2),
    years = survey2_lc_years,
    sample_sizes = c(120, 140, 160, 180),
    selectivity_hist_idx = 2,
    selectivity_proj_idx = 2,
    survey_timing = 0.7
  )
)


cat("\n=== TEST 1: SINGLE FLEET COMPREHENSIVE OBSERVATION MODELS ===\n")

result_single_comp <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  ProFisheryObj_list = proj_fishery_list,
  StrategyObj = strategy_single_comp,
  StochasticObj = stochastic_obj,
  MultifleetObj = NULL,
  IndexObj = single_comprehensive_index,
  CatchObsObj = single_comprehensive_catch,
  LengthCompObj = single_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test1_single_with_obs",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "singleCompMP"
)

cat("Single fleet comprehensive simulation with obs models completed\n")

result_single_comp <- readProjection(getwd(), "test1_single_with_obs")

# Population outputs

result_single_comp$dynamics$SB

result_single_comp$dynamics$VB

result_single_comp$dynamics$Ftotal

result_single_comp$dynamics$recN

result_single_comp$dynamics$SPR


# obs model outputs
result_single_comp$HCR$decisionData$IDX_CPUE_1
result_single_comp$HCR$decisionData$IDX_CPUE_1_areas
result_single_comp$HCR$decisionData$IDX_CPUE_2
result_single_comp$HCR$decisionData$IDX_CPUE_2_areas
result_single_comp$HCR$decisionData$IDX_Survey_3
result_single_comp$HCR$decisionData$IDX_Survey_3_areas
result_single_comp$HCR$decisionData$IDX_Survey_4
result_single_comp$HCR$decisionData$IDX_Survey_4_areas

#Exploring length obs
result_single_comp$HCR$decisionData$LC_Fishery_1_indextype
result_single_comp$HCR$decisionData$LC_Fishery_1_areas

result_single_comp$HCR$decisionData$LC_Survey_2_indextype
result_single_comp$HCR$decisionData$LC_Survey_2_areas

result_single_comp$HCR$decisionData$LC_Survey_3_indextype
result_single_comp$HCR$decisionData$LC_Survey_3_areas


#exploring cacth obs
result_single_comp$HCR$decisionData$true_catch
result_single_comp$HCR$decisionData$observed_catch

result_single_comp$HCR$decisionData$observed_catch_area_1
result_single_comp$HCR$decisionData$observed_catch_area_2


# ============================================================================
# TEST 2: MULTIFLEET 2 FLEETS WITH COMPREHENSIVE OBSERVATION MODELS
# ============================================================================

#multifleet selectivity definitions
# Fleet selectivities for multifleet tests (historical)
fleet1_sel_hist <- new("Fishery")
fleet1_sel_hist@title <- "Fleet 1"
fleet1_sel_hist@vulType <- "logistic"
fleet1_sel_hist@vulParams <- c(10.2, 0.1)  # Same as base single fleet for 1-fleet test
fleet1_sel_hist@retType <- "full"
fleet1_sel_hist@retMax <- 1
fleet1_sel_hist@Dmort <- 0

fleet2_sel_hist <- new("Fishery")
fleet2_sel_hist@title <- "Fleet 2"
fleet2_sel_hist@vulType <- "logistic"
fleet2_sel_hist@vulParams <- c(9, 0.1)  # Different selectivity
fleet2_sel_hist@retType <- "full"
fleet2_sel_hist@retMax <- 1
fleet2_sel_hist@Dmort <- 0

fleet3_sel_hist <- new("Fishery")
fleet3_sel_hist@title <- "Fleet 3"
fleet3_sel_hist@vulType <- "logistic"
fleet3_sel_hist@vulParams <- c(11.2, 0.1)  # Different selectivity
fleet3_sel_hist@retType <- "full"
fleet3_sel_hist@retMax <- 1
fleet3_sel_hist@Dmort <- 0

# Fleet selectivities for multifleet tests (projection)
fleet1_sel_proj <- new("Fishery")
fleet1_sel_proj@title <- "Fleet 1"
fleet1_sel_proj@vulType <- "logistic"
fleet1_sel_proj@vulParams <- c(10.2, 0.1)  # Same as base single fleet for 1-fleet test
fleet1_sel_proj@retType <- "full"
fleet1_sel_proj@retMax <- 1
fleet1_sel_proj@Dmort <- 0

fleet2_sel_proj <- new("Fishery")
fleet2_sel_proj@title <- "Fleet 2"
fleet2_sel_proj@vulType <- "logistic"
fleet2_sel_proj@vulParams <- c(10, 0.1)  # Different selectivity
fleet2_sel_proj@retType <- "full"
fleet2_sel_proj@retMax <- 1
fleet2_sel_proj@Dmort <- 0

fleet3_sel_proj <- new("Fishery")
fleet3_sel_proj@title <- "Fleet 3"
fleet3_sel_proj@vulType <- "logistic"
fleet3_sel_proj@vulParams <- c(11.2, 0.1)  # Different selectivity
fleet3_sel_proj@retType <- "full"
fleet3_sel_proj@retMax <- 1
fleet3_sel_proj@Dmort <- 0

# ============================================================================
# MANAGEMENT STRATEGIES - MULTIFLEET FLEET
# ============================================================================

# Multifleet strategy
multiCompMP <- function(phase, dataObject) {
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

  #need to modify phase 3
  if(phase == 3) {

    #return fleet-specific F values for each area-fleet combination
    result_data <- data.frame()

    for(area in 1:areas) {
      for(fleet in 1:nfleets) {
        result_data <- rbind(result_data, data.frame(
          year = j,
          iteration = k,
          area = area,
          fleet = fleet,
          Flocal = 0.05  # conservative F for testing
        ))
      }
    }



    return(list(year=result_data$year, iteration=result_data$iteration, area=result_data$area,
                fleet=result_data$fleet, Flocal=result_data$Flocal))
  }
}

# Strategy objects
strategy_multi_comp <- new("Strategy")
strategy_multi_comp@title <- "Multifleet Comprehensive"
strategy_multi_comp@projectionYears <- 5
strategy_multi_comp@projectionName <- "multiCompMP"
strategy_multi_comp@projectionParams <- list()


# ============================================================================
# MULTIFLEET (2 FLEETS) OBJECT
# ============================================================================

multifleet_2fleet <- new("Multifleet")
multifleet_2fleet@nfleets <- 2
multifleet_2fleet@fleet_proportions <- c(0.6, 0.4)
multifleet_2fleet@allocation_type <- "catch"
multifleet_2fleet@fleet_selectivity_hist_list <- list(fleet1_sel_hist, fleet2_sel_hist)
#multifleet_2fleet@fleet_selectivity_proj_list <- list(fleet1_sel_proj, fleet2_sel_proj)

#changing prj sel with the nested structure
multifleet_2fleet@fleet_selectivity_proj_list <- list(
  # Area 1 projections - both fleets
  list(fleet1_sel_proj, fleet2_sel_proj),
  # Area 2 projections - both fleets
  list(fleet1_sel_proj, fleet2_sel_proj)
)


#adding the array of fleet historical eefort
multifleet_2fleet@fleet_historicalEffort <- array(dim = c(ta@historicalYears, ta@areas, 2))
multifleet_2fleet@fleet_historicalEffort[,,1] <- ta@historicalEffort
multifleet_2fleet@fleet_historicalEffort[,,2] <- ta@historicalEffort


# ============================================================================
# COMPREHENSIVE OBSERVATION MODELS FOR MULTIFLEET (2 FLEETS)
# ============================================================================

# Comprehensive index object for multifleet
multi_comprehensive_index <- new("Index")
multi_comprehensive_index@indexID <- "Multifleet_Comprehensive"
multi_comprehensive_index@title <- "Multifleet Comprehensive Indices"
multi_comprehensive_index@useWeight <- TRUE

# Survey selectivities (same as single fleet)
multi_comprehensive_index@selectivity_hist_list <- list(survey1_sel_hist, survey2_sel_hist)
multi_comprehensive_index@selectivity_proj_list <- list(survey1_sel_proj, survey2_sel_proj)

multi_comprehensive_index@survey_design <- list(
  # 1. FI Survey 1 - Both areas
  list(
    indextype = "FI",
    areas = c(1, 2),
    indexYears = seq(1, 15, 3),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.3,
    q_hist_bounds = c(0.0002, 0.0005),
    q_proj_bounds = c(0.0003, 0.0006),
    hyperstability_hist_bounds = c(1.0, 1.0),
    hyperstability_proj_bounds = c(1.0, 1.0),
    obsError_CV_hist_bounds = c(0.10, 0.20),
    obsError_CV_proj_bounds = c(0.10, 0.20)
  ),

  # 2. FI Survey 2 - Area 2 only
  list(
    indextype = "FI",
    areas = c(2),
    indexYears = c(4, 8, 12, 15),
    selectivity_hist_idx = 2,
    selectivity_proj_idx = 2,
    survey_timing = 0.7,
    q_hist_bounds = c(0.00025, 0.0008),
    q_proj_bounds = c(0.0004, 0.001),
    hyperstability_hist_bounds = c(0.95, 1.05),
    hyperstability_proj_bounds = c(0.95, 1.05),
    obsError_CV_hist_bounds = c(0.15, 0.30),
    obsError_CV_proj_bounds = c(0.15, 0.30)
  ),

  # Fishery dependent indices (FD) with fleet ID
  # 3. Fleet 1 CPUE - Area 1 only, annual
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    indexYears = 1:15,
    q_hist_bounds = c(0.0001, 0.0003),
    q_proj_bounds = c(0.00015, 0.0004),
    hyperstability_hist_bounds = c(0.9, 1.1),
    hyperstability_proj_bounds = c(0.9, 1.1),
    obsError_CV_hist_bounds = c(0.15, 0.25),
    obsError_CV_proj_bounds = c(0.15, 0.25)
  ),

  # 4. Fleet 1 CPUE - Area 2 only
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(2),
    indexYears = seq(2, 15, 2),
    q_hist_bounds = c(0.00015, 0.0004),
    q_proj_bounds = c(0.0002, 0.0005),
    hyperstability_hist_bounds = c(0.8, 1.2),
    hyperstability_proj_bounds = c(0.8, 1.2),
    obsError_CV_hist_bounds = c(0.20, 0.35),
    obsError_CV_proj_bounds = c(0.20, 0.35)
  ),

  # 5. Fleet 2 CPUE - Both areas, annual
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(1),
    indexYears = 1:15,
    q_hist_bounds = c(0.0002, 0.0006),
    q_proj_bounds = c(0.0003, 0.0008),
    hyperstability_hist_bounds = c(0.7, 1.3),
    hyperstability_proj_bounds = c(0.7, 1.3),
    obsError_CV_hist_bounds = c(0.25, 0.40),
    obsError_CV_proj_bounds = c(0.25, 0.40)
  ),

  # 6. Fleet 2 CPUE - Area 2 only, triennial
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(2),
    indexYears = seq(1, 15, 2),
    q_hist_bounds = c(0.00025, 0.0007),
    q_proj_bounds = c(0.0004, 0.0009),
    hyperstability_hist_bounds = c(0.6, 1.4),
    hyperstability_proj_bounds = c(0.6, 1.4),
    obsError_CV_hist_bounds = c(0.30, 0.45),
    obsError_CV_proj_bounds = c(0.30, 0.45)
  )
)

# Comprehensive catch observations for multifleet (2 fleets)
multi_comprehensive_catch <- new("CatchObs")
multi_comprehensive_catch@catchID <- "Multifleet_Comprehensive_Catch"
multi_comprehensive_catch@title <- "Fleet-Specific Comprehensive Catch"

multi_comprehensive_catch@fleet_configs <- list(
  # Fleet 1 configuration
  list(
    fleet_id = 1,
    areas = c(1, 2),
    catchYears = 1:15,
    reporting_rates = c(seq(0.9, 1.0, length.out = 10), rep(1.0, 5)),
    obs_CVs = matrix(
      c(c(seq(0.15, 0.10, length.out = 10), rep(0.08, 5)),
        c(seq(0.25, 0.20, length.out = 10), rep(0.15, 5))),
      ncol = 2
    )
  ),

  # Fleet 2 configuration
  list(
    fleet_id = 2,
    areas = c(1, 2),
    catchYears = seq(1, 15, 1),
    reporting_rates = c(seq(0.7, 0.9, length.out = 10), rep(0.95, 5)),
    obs_CVs = matrix(
      c(c(seq(0.25, 0.20, length.out = 10), rep(0.18, 5)),
        c(seq(0.40, 0.35, length.out = 10), rep(0.30, 5))),
      ncol = 2
    )
  )
)

# Comprehensive length composition for multifleet (2 fleets)
multi_comprehensive_lcomp <- new("LCompObs")
multi_comprehensive_lcomp@indexID <- "Multifleet_Comprehensive_LComp"
multi_comprehensive_lcomp@title <- "Fleet-Specific Comprehensive Length Comp"
multi_comprehensive_lcomp@length_bin_width <- 1
multi_comprehensive_lcomp@selectivity_hist_list <- list(survey1_sel_hist, survey2_sel_hist)
multi_comprehensive_lcomp@selectivity_proj_list <- list(survey1_sel_proj, survey2_sel_proj)

multi_comprehensive_lcomp@survey_design <- list(
  # 1. Fleet 1 fishery length composition
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    years = seq(1, 15, 2),
    sample_sizes = seq(100, 250, length.out = 8)
  ),

  # 2. Fleet 2 fishery length composition
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(2),
    years = seq(1, 15, 1),
    sample_sizes = seq(80, 180, length.out = 15)
  ),

  # 3. Survey 1 length composition - Area 1 only
  list(
    indextype = "FI",
    areas = c(1,2),
    years = seq(1,15,3),
    sample_sizes = c(300,250,300,300,300),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.3
  ),

  # 4. Survey 2 length composition - Area 2 only
  list(
    indextype = "FI",
    areas = c(2),
    years = c(4, 8, 12, 15),
    sample_sizes = c(150, 180, 200, 220),
    selectivity_hist_idx = 2,
    selectivity_proj_idx = 2,
    survey_timing = 0.7
  )

)

cat("\n=== TEST 2: MULTIFLEET COMPREHENSIVE ===\n")


result_multi_comp <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StrategyObj = strategy_multi_comp,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_2fleet,
  IndexObj = multi_comprehensive_index,
  CatchObsObj = multi_comprehensive_catch,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test2_multi_comprehensive",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "multiCompMP"
)
cat("Multifleet (2 fleeets) comprehensive simulation completed\n")

result_multi_comp  <- readProjection(getwd(), "test2_multi_comprehensive")

# Population outputs

result_multi_comp$dynamics$SB
result_multi_comp$dynamics$VB
result_multi_comp$dynamics$Ftotal
result_multi_comp$dynamics$recN
result_multi_comp$dynamics$SPR

# obs model outputs
result_multi_comp$HCR$decisionData$IDX_Survey_1
result_multi_comp$HCR$decisionData$IDX_Survey_1_areas

result_multi_comp$HCR$decisionData$IDX_Survey_2
result_multi_comp$HCR$decisionData$IDX_Survey_2_areas


result_multi_comp$HCR$decisionData$IDX_CPUE_3_Fleet_1
result_multi_comp$HCR$decisionData$IDX_CPUE_3_Fleet_1_areas

result_multi_comp$HCR$decisionData$IDX_CPUE_4_Fleet_1
result_multi_comp$HCR$decisionData$IDX_CPUE_4_Fleet_1_areas

result_multi_comp$HCR$decisionData$IDX_CPUE_5_Fleet_2
result_multi_comp$HCR$decisionData$IDX_CPUE_5_Fleet_2_areas

result_multi_comp$HCR$decisionData$IDX_CPUE_6_Fleet_2
result_multi_comp$HCR$decisionData$IDX_CPUE_6_Fleet_2_areas


#LC obs
result_multi_comp$HCR$decisionData$LC_Fishery_1_Fleet_1_indextype
result_multi_comp$HCR$decisionData$LC_Fishery_1_Fleet_1_areas

result_multi_comp$HCR$decisionData$LC_Fishery_2_Fleet_2_indextype
result_multi_comp$HCR$decisionData$LC_Fishery_2_Fleet_2_areas

result_multi_comp$HCR$decisionData$LC_Survey_3_indextype
result_multi_comp$HCR$decisionData$LC_Survey_3_areas

result_multi_comp$HCR$decisionData$LC_Survey_4_indextype
result_multi_comp$HCR$decisionData$LC_Survey_4_areas


#catch obs
result_multi_comp$HCR$decisionData$total_true_catch
result_multi_comp$HCR$decisionData$total_observed_catch

result_multi_comp$HCR$decisionData$fleet_1_true_catch
result_multi_comp$HCR$decisionData$fleet_1_observed_catch

result_multi_comp$HCR$decisionData$fleet_2_true_catch
result_multi_comp$HCR$decisionData$fleet_2_observed_catch




# # Clean up intermediate files
# file.remove("test1_single_fleet.rds")
# file.remove("test2_multifleet_1_fleetE.rds")
# file.remove("test2_multifleet_1_fleetC.rds")
# file.remove("test3_multifleet_2_fleets.rds")
# file.remove("test4_multifleet_2_fleets.rds")
# file.remove("test5_multifleet_3_fleets.rds")

