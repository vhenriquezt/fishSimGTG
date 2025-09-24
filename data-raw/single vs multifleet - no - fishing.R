# ============================================================================
# TEST FILE 5: Shutting off fishing
# ============================================================================
# Purpose: Test observation models in single fleet vs multifleet contexts
# Tests:
# 1. Removing fishing


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
ta@historicalYears <- 10
ta@historicalBio <- 0.95
ta@historicalBioType <- "relB"
ta@move <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = FALSE)

# Historical effort - declining trend
# ta@historicalEffort <- matrix(c(1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6,
#                                 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6),
#                               nrow = 10, ncol = 2, byrow = FALSE)


ta@historicalEffort <- matrix(0.001, nrow = 10, ncol = 2, byrow = FALSE)

# Stochastic (not used)
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
    Flocal <- rep(0.0, areas)  # Conservative F

    return(list(year=year, iteration=iteration, area=area, Flocal=Flocal))
  }
}

#Single fleet Phase 3: Returns 2 F values (one per area)

strategy_single_comp <- new("Strategy")
strategy_single_comp@title <- "Single Fleet Comprehensive"
strategy_single_comp@projectionYears <- 40
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

# No CPUE because there is no fishing
single_comprehensive_index@survey_design <- list(

  # 1. FI Survey 1 - Both areas, every 3 years
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

  # 1. FI Survey 2 - Area 2 only
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


cat("\n=== TEST NO FISHING: SINGLE FLEET COMPREHENSIVE OBSERVATION MODELS ===\n")

result_single_comp <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  ProFisheryObj_list = proj_fishery_list,
  StrategyObj = strategy_single_comp,
  StochasticObj = NULL,
  MultifleetObj = NULL,
  IndexObj = single_comprehensive_index,
  CatchObsObj = NULL,
  LengthCompObj = single_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test_single_no_fishing",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "singleCompMP"
)

cat("Single fleet comprehensive simulation with obs models completed\n")

test_single_no_fishing <- readProjection(getwd(), "test_single_no_fishing")

# Population outputs

test_single_no_fishing$dynamics$SB

test_single_no_fishing$dynamics$VB

test_single_no_fishing$dynamics$Ftotal

test_single_no_fishing$dynamics$recN

test_single_no_fishing$dynamics$SPR


# obs model outputs (data frame now)
test_single_no_fishing$HCR$decisionData



#plots
# loading plot fucntion for toher plots

# #obs: need to add units
plot_SB(test_single_no_fishing)
plot_SB(test_single_no_fishing, areas=1)
plot_SB(test_single_no_fishing, areas=2)
plot_SB(test_single_no_fishing, areas=c(1,2))

plot_VB(test_single_no_fishing)
plot_VB(test_single_no_fishing, areas=1)
plot_VB(test_single_no_fishing, areas=2)

plot_Ftotal(test_single_no_fishing)
plot_Ftotal(test_single_no_fishing, areas=1)
plot_Ftotal(test_single_no_fishing, areas=2)


plot_catchB(test_single_no_fishing)
plot_catchB(test_single_no_fishing,areas=1)
plot_catchB(test_single_no_fishing,areas=2)

plot_catchN(test_single_no_fishing)
plot_catchN(test_single_no_fishing, areas=1)
plot_catchN(test_single_no_fishing, areas=2)

plot_discB(test_single_no_fishing)
plot_discB(test_single_no_fishing,areas=1)
plot_discB(test_single_no_fishing,areas=2)

plot_discN(test_single_no_fishing)
plot_discN(test_single_no_fishing,areas=1)
plot_discN(test_single_no_fishing,areas=2)

plot_catchB_multi(test_single_no_fishing)
plot_catchB_multi(test_single_no_fishing, areas=1)
plot_catchB_multi(test_single_no_fishing, areas=2)

plot_catchN_multi(test_single_no_fishing,show_individual = TRUE)
plot_catchN_multi(test_single_no_fishing, areas=1)
plot_catchN_multi(test_single_no_fishing, areas=2)

plot_SPR(test_single_no_fishing)
plot_recN(test_single_no_fishing)


# #plot obs models (indices)
plot_survey_indices(test_single_no_fishing)

plot_all_indices(test_single_no_fishing)


# plot LC obs models
plot_survey_length_comp(test_single_no_fishing,show_individual = TRUE)
plot_survey_length_comp(test_single_no_fishing, areas=1,show_individual = TRUE)
plot_survey_length_comp(test_single_no_fishing, areas=2,show_individual = TRUE)



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


    #if indexObj exist
    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject) #returns 1-row tibble
      #add each column from the tibble to combined_data
      for(col_name in names(index_result)) {
        combined_data[[col_name]] <- index_result[[col_name]]
      }
    }

    # if(!is.null(CatchObsObj)) {
    #   catch_result <- calculate_single_CatchObs(dataObject)
    #   for(col_name in names(catch_result)) {
    #     combined_data[[col_name]] <- catch_result[[col_name]]
    #   }
    # }

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
        Flocal <- rbind(Flocal, c(j, k, m, f, 0.01))  # Conservative F = 0.05
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
strategy_multi_comp@projectionYears <- 40
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

  # 1. Survey 1 length composition - Area 1 only
  list(
    indextype = "FI",
    areas = c(1,2),
    years = seq(1,15,3),
    sample_sizes = c(300,250,300,300,300),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.3
  ),

  # 2. Survey 2 length composition - Area 2 only
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
  #StochasticObj = NULL,
  MultifleetObj = multifleet_2fleet,
  IndexObj = multi_comprehensive_index,
  #CatchObsObj = NULL,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test_multifleet_no_fishing",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
  #customToCluster = "multiCompMP"
)
cat("Multifleet (2 fleeets) comprehensive simulation completed\n")


#there is a problem here:
#No-fishing example fails because F = 0 creates a validation issue in the multifleet selectivity structure


test_multifleet_no_fishing  <- readProjection(getwd(), "test_multifleet_no_fishing")

# Population outputs

test_multifleet_no_fishing$dynamics$SB
test_multifleet_no_fishing$dynamics$recN
test_multifleet_no_fishing$dynamics$SPR

# obs model outputs (new structure data frame)
test_multifleet_no_fishing$HCR$decisionData
test_multifleet_no_fishing$HCR$decisionData$IDX_Survey_1


# # Clean up intermediate files
# file.remove("test1_single_fleet.rds")
# file.remove("test2_multifleet_1_fleetE.rds")
# file.remove("test2_multifleet_1_fleetC.rds")
# file.remove("test3_multifleet_2_fleets.rds")
# file.remove("test4_multifleet_2_fleets.rds")
# file.remove("test5_multifleet_3_fleets.rds")


#obs: need to add units
plot_SB(test_multifleet_no_fishing)
plot_SB(test_multifleet_no_fishing, areas=1)
plot_SB(test_multifleet_no_fishing, areas=2)
plot_SB(test_multifleet_no_fishing, areas=c(1,2))


plot_catchB(test_multifleet_no_fishing)
plot_catchB(test_multifleet_no_fishing,areas=1)
plot_catchB(test_multifleet_no_fishing,areas=2)

plot_catchN(test_multifleet_no_fishing)
plot_catchN(test_multifleet_no_fishing, areas=1)
plot_catchN(test_multifleet_no_fishing, areas=2)

# plot_discB(test_multifleet_no_fishing)
# plot_discB(test_multifleet_no_fishing,areas=1)
# plot_discB(test_multifleet_no_fishing,areas=2)

plot_discN(test_multifleet_no_fishing)
plot_discN(test_multifleet_no_fishing,areas=1)
plot_discN(test_multifleet_no_fishing,areas=2)

plot_catchB_multi(test_multifleet_no_fishing)
plot_catchB_multi(test_multifleet_no_fishing, areas=1)
plot_catchB_multi(test_multifleet_no_fishing, areas=2)

plot_catchN_multi(test_multifleet_no_fishing,show_individual = TRUE)
plot_catchN_multi(test_multifleet_no_fishing, areas=1)
plot_catchN_multi(test_multifleet_no_fishing, areas=2)

plot_SPR(test_multifleet_no_fishing)
plot_recN(test_multifleet_no_fishing)


#plot obs models (indices)
plot_survey_indices(test_multifleet_no_fishing)
plot_cpue_indices(test_multifleet_no_fishing)

plot_all_indices(test_multifleet_no_fishing)

#plot individual indices
plot_indices(test_multifleet_no_fishing,
             index_pattern = "IDX_CPUE.*Fleet_1",
             show_individual = TRUE,
             title = "Fleet 1 CPUE Only")

plot_indices(test_multifleet_no_fishing,
             index_pattern = "IDX_CPUE.*Fleet_2",
             show_individual = TRUE,
             title = "Fleet 2 CPUE Only")





plot_catch_observations_both(test_multifleet_no_fishing,show_individual = TRUE)
plot_catch_observations_multifleet(test_multifleet_no_fishing,show_individual = TRUE)


# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(test_multifleet_no_fishing,show_individual = TRUE)
plot_fishery_length_comp(test_multifleet_no_fishing, areas=1,show_individual = TRUE)
plot_fishery_length_comp(test_multifleet_no_fishing, areas=2,show_individual = TRUE)

plot_survey_length_comp(test_multifleet_no_fishing,show_individual = TRUE)
plot_survey_length_comp(test_multifleet_no_fishing, areas=1,show_individual = TRUE)
plot_survey_length_comp(test_multifleet_no_fishing, areas=2,show_individual = TRUE)



# NEW: Custom filtering for fleets and areas
plot_length_composition_by_area(test_multifleet_no_fishing,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets



plot_length_composition_by_area(test_multifleet_no_fishing,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets


plot_length_composition_by_area(test_multifleet_no_fishing,
                                program_pattern = "LC_Survey",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets






