# ============================================================================
# TEST: constant high eff and nothing else changing
# ============================================================================
#1. high constant effort and nothing else changing
#2. Test deterministic scenarios with comprehensive observation models


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
ta@iterations <- 3
ta@historicalYears <- 15
ta@historicalBio <- 0.5
ta@historicalBioType <- "relB"
ta@move <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = FALSE)

# Historical effort - declining trend
# ta@historicalEffort <- matrix(
#   c(seq(2.0, 1.0, length.out = 15),  #area 1: linear decline
#     seq(1.8, 0.8, length.out = 15)), #area 2: linear decline
#   nrow = 15, ncol = 2, byrow = FALSE)

ta@historicalEffort <- matrix(6.5, nrow = 15, ncol = 2)

# Stochastic
stochastic_obj <- NULL

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
survey1_sel_hist@vulParams <- c(10.2, 0.1)
survey1_sel_hist@retType <- "full"
survey1_sel_hist@retMax <- 1
survey1_sel_hist@Dmort <- 0

survey2_sel_hist <- new("Fishery")
survey2_sel_hist@title <- "Research Survey 2 Historical"
survey2_sel_hist@vulType <- "logistic"
survey2_sel_hist@vulParams <- c(10.2, 0.1)
survey2_sel_hist@retType <- "full"
survey2_sel_hist@retMax <- 1
survey2_sel_hist@Dmort <- 0

survey1_sel_proj <- new("Fishery")
survey1_sel_proj@title <- "Research Survey 1 Projection"
survey1_sel_proj@vulType <- "logistic"
survey1_sel_proj@vulParams <- c(10.2, 0.1)
survey1_sel_proj@retType <- "full"
survey1_sel_proj@retMax <- 1
survey1_sel_proj@Dmort <- 0

survey2_sel_proj <- new("Fishery")
survey2_sel_proj@title <- "Research Survey 2 Projection"
survey2_sel_proj@vulType <- "logistic"
survey2_sel_proj@vulParams <- c(10.2, 0.1)
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

  #phase 1 called once per year/iter combination
  if(phase == 1) {
    combined_data <- list()

    #each observation model returns ONE ROW for this j,k combination
    #the model calls phase 1 repeatedly and builds the complete dataset
    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject) #returns 1-row tibble
      #add each column from the tibble to combined_data
      for(col_name in names(index_result)) {
        combined_data[[col_name]] <- index_result[[col_name]]
      }
    }

    if(!is.null(CatchObsObj)) {
      catch_result <- calculate_single_CatchObs(dataObject) #returns 1-row tibble
      # add catch observation columns
      for(col_name in names(catch_result)) {
        combined_data[[col_name]] <- catch_result[[col_name]]
      }
    }
    #same
    if(!is.null(LengthCompObj)) {
      lc_result <- calculate_single_LengthComp(dataObject)
      for(col_name in names(lc_result)) {
        combined_data[[col_name]] <- lc_result[[col_name]]
      }
    }

    return(combined_data) #returns single row of obsrvtation
  }

  if(phase == 2) return(list())
  #phase 3 called once per year/iter
  if(phase == 3) {
    #create vectors for each area
    year <- rep(j, areas)
    iteration <- rep(k, areas)
    area <- 1:areas
    Flocal <- rep(0.3, areas)  # Conservative F

    return(list(year=year, iteration=iteration, area=area, Flocal=Flocal))
  }
}

#Single fleet Phase 3: Returns 2 F values (one per area)

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
    q_hist_bounds = c(0.0002, 0.0002),
    q_proj_bounds = c(0.0002, 0.0002),
    hyperstability_hist_bounds = c(1.0, 1.0),
    hyperstability_proj_bounds = c(1.0, 1.0),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
  ),

  # 2. FD CPUE - Area 1 only, biennial
  list(
    indextype = "FD",
    areas = c(1),
    indexYears = seq(1, 15, 2),  # Every other year
    q_hist_bounds = c(0.0003, 0.0003),
    q_proj_bounds = c(0.0003, 0.0003),
    hyperstability_hist_bounds = c(1.0, 1.0),
    hyperstability_proj_bounds = c(1.0, 1.0),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
  ),

  # 3. FI Survey 1 - Both areas, every 3 years
  list(
    indextype = "FI",
    areas = c(1, 2),
    indexYears = seq(1, 15, 3),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.3,  # early in year
    q_hist_bounds = c(0.0004, 0.0004),
    q_proj_bounds = c(0.0004, 0.0004),
    hyperstability_hist_bounds = c(1.0, 1.0),
    hyperstability_proj_bounds = c(1.0, 1.0),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
  ),

  # 4. FI Survey 2 - Area 2 only
  list(
    indextype = "FI",
    areas = c(2),
    indexYears = c(4, 8, 12, 15),  # Specific years
    selectivity_hist_idx = 2,
    selectivity_proj_idx = 2,
    survey_timing = 0.7,  # Late in year
    q_hist_bounds = c(0.0004, 0.0004),
    q_proj_bounds = c(0.0004, 0.0004),
    hyperstability_hist_bounds = c(1, 1),
    hyperstability_proj_bounds = c(1, 1),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
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
single_comprehensive_catch@reporting_rates <- rep(1.0, 20)
# Variable observation error - decreasing over time
single_comprehensive_catch@obs_CVs <- matrix(rep(c(0, 0), each = 20), ncol = 2)

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
    sample_sizes = rep(500, length(fishery_lc_years))
  ),

  # 2. Survey 1 length composition
  list(
    indextype = "FI",
    areas = c(1, 2),
    years = survey1_lc_years,
    sample_sizes = c(500, 500,500, 500, 500),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.3
  ),

  # 3. Survey 2 length composition - Area 1 only
  list(
    indextype = "FI",
    areas = c(2),
    years = survey2_lc_years,
    sample_sizes = c(500, 500, 500, 500),
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
  StochasticObj = NULL,
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


# obs model outputs (data frame now)
str(result_single_comp$HCR$decisionData)
result_single_comp$HCR$decisionData$IDX_CPUE_1


#plots
# loading plot fucntion for toher plots

# #obs: need to add units
plot_SB(result_single_comp)
plot_SB(result_single_comp, areas=1)
plot_SB(result_single_comp, areas=2)
plot_SB(result_single_comp, areas=c(1,2))

plot_VB(result_single_comp)
plot_VB(result_single_comp, areas=1)
plot_VB(result_single_comp, areas=2)

plot_Ftotal(result_single_comp)
plot_Ftotal(result_single_comp, areas=1)
plot_Ftotal(result_single_comp, areas=2)


plot_catchB(result_single_comp)
plot_catchB(result_single_comp,areas=1)
plot_catchB(result_single_comp,areas=2)

plot_catchN(result_single_comp)
plot_catchN(result_single_comp, areas=1)
plot_catchN(result_single_comp, areas=2)

plot_discB(result_single_comp)
plot_discB(result_single_comp,areas=1)
plot_discB(result_single_comp,areas=2)

plot_discN(result_single_comp)
plot_discN(result_single_comp,areas=1)
plot_discN(result_single_comp,areas=2)

plot_catchB_multi(result_single_comp)
plot_catchB_multi(result_single_comp, areas=1)
plot_catchB_multi(result_single_comp, areas=2)

plot_catchN_multi(result_single_comp,show_individual = TRUE)
plot_catchN_multi(result_single_comp, areas=1)
plot_catchN_multi(result_single_comp, areas=2)

plot_SPR(result_single_comp)
plot_recN(result_single_comp)


# #plot obs models (indices)
plot_survey_indices(result_single_comp)
plot_cpue_indices(result_single_comp)

plot_all_indices(result_single_comp)

plot_catch_observations_both(result_single_comp,show_individual = FALSE)
plot_catch_observations_both(result_single_comp,show_individual = TRUE)


# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_single_comp,show_individual = TRUE)

plot_survey_length_comp(result_single_comp,show_individual = TRUE)
plot_survey_length_comp(result_single_comp, areas=1,show_individual = TRUE)
plot_survey_length_comp(result_single_comp, areas=2,show_individual = TRUE)



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
fleet2_sel_hist@vulParams <- c(10.2, 0.1)
fleet2_sel_hist@retType <- "full"
fleet2_sel_hist@retMax <- 1
fleet2_sel_hist@Dmort <- 0

fleet3_sel_hist <- new("Fishery")
fleet3_sel_hist@title <- "Fleet 3"
fleet3_sel_hist@vulType <- "logistic"
fleet3_sel_hist@vulParams <- c(10.2, 0.1)
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
fleet2_sel_proj@vulParams <- c(10.2, 0.1)
fleet2_sel_proj@retType <- "full"
fleet2_sel_proj@retMax <- 1
fleet2_sel_proj@Dmort <- 0

fleet3_sel_proj <- new("Fishery")
fleet3_sel_proj@title <- "Fleet 3"
fleet3_sel_proj@vulType <- "logistic"
fleet3_sel_proj@vulParams <- c(10.2, 0.1)
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

    #debugging
    if(phase == 1 && j == 2 && k == 1) {
      cat("=== ARRAY AVAILABILITY CHECK ===\n")
      cat("N exists:", exists("N"), "\n")
      cat("Z exists:", exists("Z"), "\n")
      cat("catchNage exists:", exists("catchNage"), "\n")
      cat("VB exists:", exists("VB"), "\n")
      cat("RB exists:", exists("RB"), "\n")
      cat("catchB exists:", exists("catchB"), "\n")
      cat("is_multifleet:", is_multifleet, "\n")
      if(is_multifleet) {
        cat("RB_by_fleet exists:", exists("RB_by_fleet"), "\n")
        cat("catchB_by_fleet exists:", exists("catchB_by_fleet"), "\n")
      }
      cat("================================\n")
    }
    #if indexObj exist
    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject) #returns 1-row tibble
      #add each column from the tibble to combined_data
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

  #Vania edit's to match Bill's edits
  if(phase == 3) {

    Flocal <- data.frame()

    #creating one row per area-fleet combination
    for(m in 1:areas) {
      for(f in 1:nfleets) {
        #row contain: [year, iteration, area, fleet, F_value]
        Flocal <- rbind(Flocal, c(j, k, m, f, 0.12))  # Conservative F = 0.05
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

#Multi Fleet Phase 3: Returns 4 F values (one per area-fleet combination)

# For a simulation with 15 years × 6 iterations = 90
# Phase 1 gets called 90 times
# Each call returns one row of observation data
# the model rbinds these together into decisionData

#Now, each obs column becomes a vector of length 90
#For example:
#IDX_CPUE_1
#Year 1, Iter 1: value_1
#Year 1, Iter 2: value_2

#This creates: [value_1, value_2, value_3, value_4, ...] - a long vector of length 90




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

# New structure: For 2 areas, 2 fleets
multifleet_2fleet@fleet_selectivity_proj_list <- list(
  # Area 1
  list(fleet1_sel_proj, fleet2_sel_proj),
  # Area 2
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
    q_hist_bounds = c(0.0002, 0.0002),
    q_proj_bounds = c(0.0002, 0.0002),
    hyperstability_hist_bounds = c(1.0, 1.0),
    hyperstability_proj_bounds = c(1.0, 1.0),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
  ),

  # 2. FI Survey 2 - Area 2 only
  list(
    indextype = "FI",
    areas = c(2),
    indexYears = c(4, 8, 12, 15),
    selectivity_hist_idx = 2,
    selectivity_proj_idx = 2,
    survey_timing = 0.7,
    q_hist_bounds = c(0.0004, 0.0004),
    q_proj_bounds = c(0.0004, 0.0004),
    hyperstability_hist_bounds = c(1, 1),
    hyperstability_proj_bounds = c(1, 1),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
  ),

  # Fishery dependent indices (FD) with fleet ID
  # 3. Fleet 1 CPUE - Area 1 only, annual
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    indexYears = 1:15,
    q_hist_bounds = c(0.0003, 0.0003),
    q_proj_bounds = c(0.0003, 0.0003),
    hyperstability_hist_bounds = c(1, 1),
    hyperstability_proj_bounds = c(1, 1),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
  ),

  # 4. Fleet 1 CPUE - Area 2 only
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(2),
    indexYears = seq(2, 15, 2),
    q_hist_bounds = c(0.0004, 0.0004),
    q_proj_bounds = c(0.0004, 0.0004),
    hyperstability_hist_bounds = c(1, 1),
    hyperstability_proj_bounds = c(1, 1),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
  ),

  # 5. Fleet 2 CPUE - Both areas, annual
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(1),
    indexYears = 1:15,
    q_hist_bounds = c(0.0006, 0.0006),
    q_proj_bounds = c(0.0006, 0.0006),
    hyperstability_hist_bounds = c(1, 1),
    hyperstability_proj_bounds = c(1, 1),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
  ),

  # 6. Fleet 2 CPUE - Area 2 only, triennial
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(2),
    indexYears = seq(1, 15, 2),
    q_hist_bounds = c(0.0009, 0.0009),
    q_proj_bounds = c(0.0009, 0.0009),
    hyperstability_hist_bounds = c(1, 1),
    hyperstability_proj_bounds = c(1, 1),
    obsError_CV_hist_bounds = c(0, 0),
    obsError_CV_proj_bounds = c(0, 0)
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
    catchYears = 1:20,
    reporting_rates = rep(1.0, 20),

    obs_CVs = matrix(rep(c(0, 0), each = 20), ncol = 2)
  ),

  # Fleet 2 configuration
  list(
    fleet_id = 2,
    areas = c(1, 2),
    catchYears = seq(1, 20, 1),
    reporting_rates = rep(1.0, 20),
    obs_CVs = matrix(rep(c(0, 0), each = 20), ncol = 2)
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
    sample_sizes = rep(500, length = 8)
  ),

  # 2. Fleet 2 fishery length composition
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(2),
    years = seq(1, 15, 1),
    sample_sizes = rep(500, length = 15)
  ),

  # 3. Survey 1 length composition - Area 1 only
  list(
    indextype = "FI",
    areas = c(1,2),
    years = seq(1,15,3),
    sample_sizes = c(500,500,500,500,500),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.3
  ),

  # 4. Survey 2 length composition - Area 2 only
  list(
    indextype = "FI",
    areas = c(2),
    years = c(4, 8, 12, 15),
    sample_sizes = c(500, 500, 500, 500),
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
  StochasticObj = NULL,
  MultifleetObj = multifleet_2fleet,
  IndexObj = multi_comprehensive_index,
  CatchObsObj = multi_comprehensive_catch,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test2_multi_comprehensive",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
  #customToCluster = "multiCompMP"
)
cat("Multifleet (2 fleeets) comprehensive simulation completed\n")

result_multi_comp  <- readProjection(getwd(), "test2_multi_comprehensive")

# Population outputs

result_multi_comp$dynamics$SB
#result_multi_comp$dynamics$VB
#result_multi_comp$dynamics$Ftotal
result_multi_comp$dynamics$recN
result_multi_comp$dynamics$SPR
result_multi_comp$dynamics$multifleet$Ftotal_by_fleet

# obs model outputs (new structure data frame)
result_multi_comp$HCR$decisionData
result_multi_comp$HCR$decisionData$IDX_Survey_1


# # Clean up intermediate files
# file.remove("test1_single_fleet.rds")
# file.remove("test2_multifleet_1_fleetE.rds")
# file.remove("test2_multifleet_1_fleetC.rds")
# file.remove("test3_multifleet_2_fleets.rds")
# file.remove("test4_multifleet_2_fleets.rds")
# file.remove("test5_multifleet_3_fleets.rds")


#obs: need to add units
plot_SB(result_multi_comp)
plot_SB(result_multi_comp, areas=1)
plot_SB(result_multi_comp, areas=2)
plot_SB(result_multi_comp, areas=c(1,2))


plot_catchB(result_multi_comp)
plot_catchB(result_multi_comp,areas=1)
plot_catchB(result_multi_comp,areas=2)

plot_catchN(result_multi_comp)
plot_catchN(result_multi_comp, areas=1)
plot_catchN(result_multi_comp, areas=2)

# plot_discB(result_multi_comp)
# plot_discB(result_multi_comp,areas=1)
# plot_discB(result_multi_comp,areas=2)

plot_discN(result_multi_comp)
plot_discN(result_multi_comp,areas=1)
plot_discN(result_multi_comp,areas=2)

plot_catchB_multi(result_multi_comp)
plot_catchB_multi(result_multi_comp, areas=1)
plot_catchB_multi(result_multi_comp, areas=2)

plot_catchN_multi(result_multi_comp,show_individual = TRUE)
plot_catchN_multi(result_multi_comp, areas=1)
plot_catchN_multi(result_multi_comp, areas=2)

plot_SPR(result_multi_comp)
plot_recN(result_multi_comp)


#plot obs models (indices)
plot_survey_indices(result_multi_comp)
plot_cpue_indices(result_multi_comp)

plot_all_indices(result_multi_comp)

#plot individual indices
plot_indices(result_multi_comp,
             index_pattern = "IDX_CPUE.*Fleet_1",
             show_individual = TRUE,
             title = "Fleet 1 CPUE Only")

plot_indices(result_multi_comp,
             index_pattern = "IDX_CPUE.*Fleet_2",
             show_individual = TRUE,
             title = "Fleet 2 CPUE Only")





plot_catch_observations_both(result_multi_comp,show_individual = TRUE)
plot_catch_observations_multifleet(result_multi_comp,show_individual = TRUE)


# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_multi_comp,show_individual = TRUE)
plot_fishery_length_comp(result_multi_comp, areas=1,show_individual = TRUE)
plot_fishery_length_comp(result_multi_comp, areas=2,show_individual = TRUE)

plot_survey_length_comp(result_multi_comp,show_individual = TRUE)
plot_survey_length_comp(result_multi_comp, areas=1,show_individual = TRUE)
plot_survey_length_comp(result_multi_comp, areas=2,show_individual = TRUE)



# NEW: Custom filtering for fleets and areas
plot_length_composition_by_area(result_multi_comp,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets



plot_length_composition_by_area(result_multi_comp,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets


plot_length_composition_by_area(result_multi_comp,
                                program_pattern = "LC_Survey",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets






