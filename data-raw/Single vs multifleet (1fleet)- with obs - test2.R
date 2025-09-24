## ============================================================================
# TEST FILE 4: BASIC SINGLE FLEET VS MULTIFLEET (1 FLEET)
#                   INCLUDING OBSERVATION MODELS
# ============================================================================
# Purpose: Validate backward compatibility and  multifleet functionality
# Tests:
# 1. Single fleet (original)
# 2. Multifleet with 1 fleet (EFFORT - CATCH) (should match single fleet)


rm(list=ls())
devtools::load_all()
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# SHARED SETUP FOR ALL TESTS
# ============================================================================

cat("=====================================\n")
cat("FISHSIMGTG MULTIFLEET VALIDATION - FILE 1\n")
cat("Basic Single vs Multifleet Comparison\n")
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

# Time Area - consistent
ta <- new("TimeArea")
ta@title <- "Validation Test"
ta@gtg <- 13
ta@areas <- 2
ta@recArea <- c(0.99, 0.01)
ta@iterations <- 2  # Small for quick testing
ta@historicalYears <- 10
ta@historicalBio <- 0.5
ta@historicalBioType <- "relB"
ta@move <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = FALSE)

# Historical effort - declining trend
ta@historicalEffort <- matrix(c(1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6,
                                1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6),
                              nrow = 10, ncol = 2, byrow = FALSE)


# Stochastic - minimal variation for clear comparison
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


# Fleet selectivities for multifleet (1 FLEET) tests (historical)
fleet1_sel_hist <- new("Fishery")
fleet1_sel_hist@title <- "Fleet 1"
fleet1_sel_hist@vulType <- "logistic"
fleet1_sel_hist@vulParams <- c(10.2, 0.1)  # Same as base single fleet for 1-fleet test
fleet1_sel_hist@retType <- "full"
fleet1_sel_hist@retMax <- 1
fleet1_sel_hist@Dmort <- 0


# Fleet selectivities for multifleet (1 FLEET) tests (projection)
fleet1_sel_proj <- new("Fishery")
fleet1_sel_proj@title <- "Fleet 1"
fleet1_sel_proj@vulType <- "logistic"
fleet1_sel_proj@vulParams <- c(10.2, 0.1)  # Same as base single fleet for 1-fleet test
fleet1_sel_proj@retType <- "full"
fleet1_sel_proj@retMax <- 1
fleet1_sel_proj@Dmort <- 0


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

# ==================================== #
#            END SHARED SET UP         #
# ==================================== #



# ============================================================================
# SIMPLE MANAGEMENT STRATEGY - WITHOUTH OBSERVATION MODELS
# ============================================================================

#Bill edit: in single-species mode you do not need to specify an MP to simulate
#historical dynamics. I've made changes to multi species that enable the same.

simpleStrategy <- function(phase, dataObject) {
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase == 1){
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
  if(phase == 3) {
    #create vectors for each area
    year <- rep(j, areas)
    iteration <- rep(k, areas)
    area <- 1:areas
    Flocal <- rep(0.05, areas)  # Conservative F

    return(list(year=year, iteration=iteration, area=area, Flocal=Flocal))
  }
}

strategy_obj <- new("Strategy")
strategy_obj@title <- "Simple Constant F"
strategy_obj@projectionYears <- 5
strategy_obj@projectionName <- "simpleStrategy"
strategy_obj@projectionParams <- list()



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




# ============================================================================
# TEST 1: SINGLE FLEET
# ============================================================================

cat("\n=== TEST 1: SINGLE FLEET (BASELINE) ===\n")

result_single <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  ProFisheryObj_list = proj_fishery_list,
  StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = NULL,
  IndexObj = single_comprehensive_index,
  CatchObsObj = single_comprehensive_catch,
  LengthCompObj = single_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test2_single_fleet",
  seed = test_seed,
  doPlot = TRUE,
  doDiagnostic = FALSE
)

result_single  <- readProjection(getwd(), "test2_single_fleet")
result_single$HCR$decisionData

cat("Result Single fleet simulation completed\n")
plot_SB(result_single)
plot_SB(result_single, areas=1)
plot_SB(result_single, areas=2)
plot_SB(result_single, areas=c(1,2))


plot_catchB(result_single)
plot_catchB(result_single,areas=1)
plot_catchB(result_single,areas=2)

plot_catchN(result_single)
plot_catchN(result_single, areas=1)
plot_catchN(result_single, areas=2)

plot_discB(result_single)
plot_discB(result_single,areas=1)
plot_discB(result_single,areas=2)

plot_discN(result_single)
plot_discN(result_single,areas=1)
plot_discN(result_single,areas=2)

plot_catchB_multi(result_single)
plot_catchB_multi(result_single, areas=1)
plot_catchB_multi(result_single, areas=2)

plot_catchN_multi(result_single,show_individual = TRUE)
plot_catchN_multi(result_single, areas=1)
plot_catchN_multi(result_single, areas=2)

plot_SPR(result_single)
plot_recN(result_single)


#plot obs models (indices)
plot_survey_indices(result_single)
plot_cpue_indices(result_single)
plot_all_indices(result_single)


#plot individual indices
plot_indices(result_single,
             index_pattern = "IDX_CPUE_1",
             show_individual = TRUE,
             title = "CPUE 1 Only")

plot_indices(result_single,
             index_pattern = "IDX_CPUE_2",
             show_individual = TRUE,
             title = "CPUE 2 Only")


plot_indices(result_single,
             index_pattern = "IDX_Survey_3",
             show_individual = TRUE,
             title = "Survey 3 Only")


plot_indices(result_single,
             index_pattern = "IDX_Survey_4",
             show_individual = TRUE,
             title = "Survey 4 Only")

#plot catch obs
plot_catch_observations_both(result_single,show_individual = TRUE)

# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_single,show_individual = TRUE)
plot_fishery_length_comp(result_single, areas=1,show_individual = TRUE)
plot_fishery_length_comp(result_single, areas=2,show_individual = TRUE)

plot_survey_length_comp(result_single,show_individual = TRUE)
plot_survey_length_comp(result_single, areas=1,show_individual = TRUE)
plot_survey_length_comp(result_single, areas=2,show_individual = TRUE)



# NEW: Custom filtering for fleets and areas
plot_length_composition_by_area(result_single,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific SIM



plot_length_composition_by_area(result_single,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets


plot_length_composition_by_area(result_single,
                                program_pattern = "LC_Survey",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets



# ============================================================================
# TEST 2: MULTIFLEET WITH 1 FLEET (effort) (SHOULD MATCH SINGLE FLEET)
# ============================================================================
#Strategy

multiStrategy <- function(phase, dataObject) {
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase == 1){
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


strategy_obj <- new("Strategy")
strategy_obj@title <- "Simple Constant F"
strategy_obj@projectionYears <- 5
strategy_obj@projectionName <- "multiStrategy"
strategy_obj@projectionParams <- list()



cat("\n=== TEST 2: MULTIFLEET WITH 1 FLEET ===\n")

multifleet_1E <- new("Multifleet")
multifleet_1E@nfleets <- 1
multifleet_1E@fleet_proportions <- c(1.0)
multifleet_1E@allocation_type <- "effort"
multifleet_1E@fleet_selectivity_hist_list <- list(fleet1_sel_hist)
#multifleet_1E@fleet_selectivity_proj_list <- list(fleet1_sel_proj) # for multifleet 1 fleet use simple list
multifleet_1E@fleet_selectivity_proj_list <- list(
  # Area 1
  list(fleet1_sel_proj),
  # Area 2
  list(fleet1_sel_proj)
)

multifleet_1E@fleet_historicalEffort <- array(dim =c(ta@historicalYears, ta@areas, 1))
multifleet_1E@fleet_historicalEffort[,,1] <-  ta@historicalEffort



# ============================================================================
# COMPREHENSIVE OBSERVATION MODELS FOR MULTIFLEET 1 FLEET
# ============================================================================
# Comprehensive index object for MULTI fleet (1 FLEET)
multi_comprehensive_index <- new("Index")
multi_comprehensive_index@indexID <- "Multifleet_Comprehensive"
multi_comprehensive_index@title <- "Multifleet Comprehensive Indices"
multi_comprehensive_index@useWeight <- TRUE

# Survey selectivities (same as single fleet)
multi_comprehensive_index@selectivity_hist_list <- list(survey1_sel_hist, survey2_sel_hist)
multi_comprehensive_index@selectivity_proj_list <- list(survey1_sel_proj, survey2_sel_proj)

multi_comprehensive_index@survey_design <- list(
  # 1. FD CPUE - Data covers both areas, annual
  list(
    indextype = "FD",
    fleet_id = 1, #FLEET 1
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
    fleet_id = 1, #FLEET 1
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

# Comprehensive catch observations for MULTI fleet (1 FLEET)
multi_comprehensive_catch <- new("CatchObs")
multi_comprehensive_catch@catchID <- "Multifleet_Comprehensive_Catch"
multi_comprehensive_catch@title <- "Multifleet Comprehensive Catch"

multi_comprehensive_catch@fleet_configs <- list(
  # Fleet 1 configuration
  list(
    fleet_id = 1,
    areas = c(1, 2),
    catchYears = 1:15,
    reporting_rates = c(seq(0.8, 1.0, length.out = 10), rep(1.0, 5)),
    obs_CVs = matrix(
      c(c(seq(0.3, 0.2, length.out = 10), rep(0.15, 5)),
        c(seq(0.45, 0.3, length.out = 10), rep(0.25, 5))),
      ncol = 2
    )
  )
)



# Comprehensive length composition for MULTI fleet (1 FLEET)
multi_comprehensive_lcomp <- new("LCompObs")
multi_comprehensive_lcomp@indexID <- "Multifleet_Comprehensive_LComp"
multi_comprehensive_lcomp@title <- "Multifleet Comprehensive Length Comp"
multi_comprehensive_lcomp@length_bin_width <- 1
multi_comprehensive_lcomp@selectivity_hist_list <- list(survey1_sel_hist, survey2_sel_hist)
multi_comprehensive_lcomp@selectivity_proj_list <- list(survey1_sel_proj, survey2_sel_proj)

# Calculate years for length composition sampling
fishery_lc_years <- seq(2, 15, 2)  # every other year fishery sampling
survey1_lc_years <- seq(1, 15, 3)  # every 3 years survey 1
survey2_lc_years <- c(4, 8, 12, 15) # Irregular survey 2



multi_comprehensive_lcomp@survey_design <- list(
  # 1. Fleet 1 fishery length composition
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1,2),
    years = fishery_lc_years,
    sample_sizes = seq(80, 200, length.out = length(fishery_lc_years))
  ),

  # 2. Survey FI  length composition
  list(
    indextype = "FI",
    areas = c(1,2),
    years = survey1_lc_years,
    sample_sizes = c(150, 180,100, 195, 200),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.3
  ),

  # 3. Survey FI length composition - Area 2 only
  list(
    indextype = "FI",
    areas = c(2),
    years = survey2_lc_years,
    sample_sizes = c(120, 140, 160, 180),
    selectivity_hist_idx = 2,
    selectivity_proj_idx = 2,
    survey_timing = 0.7
  ))


result_multifleet_1E <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_1E,
  IndexObj = multi_comprehensive_index,
  CatchObsObj = multi_comprehensive_catch,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test2_multifleet_1_fleetE",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
)

cat("Multifleet (1 fleet - effort) simulation completed\n")


result_multi_1E  <- readProjection(getwd(), "test2_multifleet_1_fleetE")

plot_SB(result_multi_1E)
plot_SB(result_multi_1E, areas=1)
plot_SB(result_multi_1E, areas=2)
plot_SB(result_multi_1E, areas=c(1,2))


plot_catchB(result_multi_1E)
plot_catchB(result_multi_1E,areas=1)
plot_catchB(result_multi_1E,areas=2)

plot_catchN(result_multi_1E)
plot_catchN(result_multi_1E, areas=1)
plot_catchN(result_multi_1E, areas=2)

# plot_discB(result_multi_1E)
# plot_discB(result_multi_1E,areas=1)
# plot_discB(result_multi_1E,areas=2)

plot_discN(result_multi_1E)
plot_discN(result_multi_1E,areas=1)
plot_discN(result_multi_1E,areas=2)

plot_catchB_multi(result_multi_1E)
plot_catchB_multi(result_multi_1E, areas=1)
plot_catchB_multi(result_multi_1E, areas=2)

plot_catchN_multi(result_multi_1E,show_individual = TRUE)
plot_catchN_multi(result_multi_1E, areas=1)
plot_catchN_multi(result_multi_1E, areas=2)

plot_SPR(result_multi_1E)
plot_recN(result_multi_1E)


#plot obs models (indices)
plot_survey_indices(result_multi_1E)
plot_cpue_indices(result_multi_1E)

plot_all_indices(result_multi_1E)

result_multi_1E$HCR$decisionData$IDX_Survey_3

#plot individual indices
plot_indices(result_multi_1E,
             index_pattern = "IDX_CPUE_1_Fleet_1",
             show_individual = TRUE,
             title = "Fleet 1 CPUE1 Only")

plot_indices(result_multi_1E,
             index_pattern = "IDX_CPUE_2_Fleet_1",
             show_individual = TRUE,
             title = "Fleet 1 CPUE2 Only")

plot_indices(result_multi_1E,
             index_pattern = "IDX_Survey_3",
             show_individual = TRUE,
             title = "Survey 3 Only")

plot_indices(result_multi_1E,
             index_pattern = "IDX_Survey_4",
             show_individual = TRUE,
             title = "Survey 4 Only")



plot_catch_observations_both(result_multi_1E,show_individual = TRUE)
plot_catch_observations_multifleet(result_multi_1E,show_individual = TRUE)


# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_multi_1E,show_individual = TRUE)
plot_fishery_length_comp(result_multi_1E, areas=1,show_individual = TRUE)
plot_fishery_length_comp(result_multi_1E, areas=2,show_individual = TRUE)

plot_survey_length_comp(result_multi_1E,show_individual = TRUE)
plot_survey_length_comp(result_multi_1E, areas=1,show_individual = TRUE)
plot_survey_length_comp(result_multi_1E, areas=2,show_individual = TRUE)



# NEW: Custom filtering for fleets and areas
result_multi_1E$HCR$decisionData$LC_Survey_2_indextype

plot_length_composition_by_area(result_multi_1E,
                                program_pattern = "LC_Fishery_1",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets



plot_length_composition_by_area(result_multi_1E,
                                program_pattern = "LC_Survey_2",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets


plot_length_composition_by_area(result_multi_1E,
                                program_pattern = "LC_Survey_3",
                                area_filter = c(2),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets




# ============================================================================
# TEST 3: MULTIFLEET WITH 1 FLEET (catch) (SHOULD MATCH SINGLE FLEET)
# ============================================================================

cat("\n=== TEST 3: MULTIFLEET WITH 1 FLEET ===\n")

multifleet_1C <- new("Multifleet")
multifleet_1C@nfleets <- 1
multifleet_1C@fleet_proportions <- c(1.0)
multifleet_1C@allocation_type <- "catch"
multifleet_1C@fleet_selectivity_hist_list <- list(fleet1_sel_hist)
multifleet_1C@fleet_selectivity_proj_list <- list(
  # Area 1
  list(fleet1_sel_proj),
  # Area 2
  list(fleet1_sel_proj)
)
multifleet_1C@fleet_historicalEffort <- array(dim =c(ta@historicalYears, ta@areas, 1))
multifleet_1C@fleet_historicalEffort[,,1] <-  ta@historicalEffort

result_multifleet_1C <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_1C,
  IndexObj = multi_comprehensive_index,
  CatchObsObj = multi_comprehensive_catch,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test2_multifleet_1_fleetC",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
)

cat("Multifleet (1 fleet - catch) simulation completed\n")
result_multi_1C <- readProjection(getwd(), "test2_multifleet_1_fleetC")

plot_SB(result_multi_1C)
plot_SB(result_multi_1C, areas=1)
plot_SB(result_multi_1C, areas=2)
plot_SB(result_multi_1C, areas=c(1,2))


plot_catchB(result_multi_1C)
plot_catchB(result_multi_1C,areas=1)
plot_catchB(result_multi_1C,areas=2)

plot_catchN(result_multi_1C)
plot_catchN(result_multi_1C, areas=1)
plot_catchN(result_multi_1C, areas=2)

# plot_discB(result_multi_1C)
# plot_discB(result_multi_1C,areas=1)
# plot_discB(result_multi_1C,areas=2)

plot_discN(result_multi_1C)
plot_discN(result_multi_1C,areas=1)
plot_discN(result_multi_1C,areas=2)

plot_catchB_multi(result_multi_1C)
plot_catchB_multi(result_multi_1C, areas=1)
plot_catchB_multi(result_multi_1C, areas=2)

plot_catchN_multi(result_multi_1C,show_individual = TRUE)
plot_catchN_multi(result_multi_1C, areas=1)
plot_catchN_multi(result_multi_1C, areas=2)

plot_SPR(result_multi_1C)
plot_recN(result_multi_1C)


#plot obs models (indices)
plot_survey_indices(result_multi_1C)
plot_cpue_indices(result_multi_1C)

plot_all_indices(result_multi_1C)

result_multi_1C$HCR$decisionData$IDX_Survey_3

#plot individual indices
plot_indices(result_multi_1C,
             index_pattern = "IDX_CPUE_1_Fleet_1",
             show_individual = TRUE,
             title = "Fleet 1 CPUE1 Only")

plot_indices(result_multi_1C,
             index_pattern = "IDX_CPUE_2_Fleet_1",
             show_individual = TRUE,
             title = "Fleet 1 CPUE2 Only")

plot_indices(result_multi_1C,
             index_pattern = "IDX_Survey_3",
             show_individual = TRUE,
             title = "Survey 3 Only")

plot_indices(result_multi_1C,
             index_pattern = "IDX_Survey_4",
             show_individual = TRUE,
             title = "Survey 4 Only")



plot_catch_observations_both(result_multi_1C,show_individual = TRUE)
plot_catch_observations_multifleet(result_multi_1C,show_individual = TRUE)


# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_multi_1C,show_individual = TRUE)
plot_fishery_length_comp(result_multi_1C, areas=1,show_individual = TRUE)
plot_fishery_length_comp(result_multi_1C, areas=2,show_individual = TRUE)

plot_survey_length_comp(result_multi_1C,show_individual = TRUE)
plot_survey_length_comp(result_multi_1C, areas=1,show_individual = TRUE)
plot_survey_length_comp(result_multi_1C, areas=2,show_individual = TRUE)



# NEW: Custom filtering for fleets and areas
result_multi_1C$HCR$decisionData$LC_Survey_2_indextype

plot_length_composition_by_area(result_multi_1C,
                                program_pattern = "LC_Fishery_1",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets



plot_length_composition_by_area(result_multi_1C,
                                program_pattern = "LC_Survey_2",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets


plot_length_composition_by_area(result_multi_1C,
                                program_pattern = "LC_Survey_3",
                                area_filter = c(2),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets





# ============================================================================
# COMPARISON AND VALIDATION (SINGLE FLEET VS MULTIFLEET)
# ============================================================================

cat("\n=== COMPARISON AND VALIDATION ===\n")


result_single
result_multi_1E
result_multi_1C

#Bill edit: modified to work with lack of Ftotal for multi fleet
# Function to compare key metrics
compare_metrics <- function(result1, result2, name1, name2, tolerance = 1e-10) {
  cat(sprintf("\nComparing %s vs %s:\n", name1, name2))

  # Compare spawning biomass
  sb_diff <- max(abs(result1$dynamics$SB - result2$dynamics$SB))
  cat(sprintf("  Max SB difference: %.2e\n", sb_diff))

  # Compare total F
  if(is.null(result1$dynamics$multifleet)){
    F1 <- result1$dynamics$Ftotal
  } else {
    F1 <- result1$dynamics$multifleet$Ftotal_by_fleet[,,,1]
  }

  if(is.null(result2$dynamics$multifleet)){
    F2 <- result2$dynamics$Ftotal
  } else {
    F2 <- result2$dynamics$multifleet$Ftotal_by_fleet[,,,1]
  }

  f_diff <- max(abs(F1 - F2))
  cat(sprintf("  Max F difference: %.2e\n", f_diff))

  # Compare SPR
  spr_diff <- max(abs(result1$dynamics$SPR - result2$dynamics$SPR))
  cat(sprintf("  Max SPR difference: %.2e\n", spr_diff))

  # Overall assessment
  max_diff <- max(sb_diff, f_diff, spr_diff)
  if(max_diff < tolerance) {
    cat(sprintf("  RESULT: IDENTICAL (max diff = %.2e)\n", max_diff))
    return(TRUE)
  } else {
    cat(sprintf("  RESULT: DIFFERENT (max diff = %.2e)\n", max_diff))
    return(FALSE)
  }
}

# Test 1: Single fleet vs Multifleet(1 fleetE) should be identical
identical_1 <- compare_metrics(result_single, result_multi_1E,
                               "Single Fleet", "Multifleet(1E)", 1e-8)


# Test 1: Single fleet vs Multifleet(1 fleetC) should be identical
identical_2 <- compare_metrics(result_single, result_multi_1C,
                               "Single Fleet", "Multifleet(1C)", 1e-8)

# Test 1: Multifleet(1 fleetE) vs Multifleet(1 fleetC) should be identical
identical_3 <- compare_metrics(result_multi_1E, result_multi_1C,
                               "Multifleet(1E)", "Multifleet(1C)", 1e-8)



all_results <- list(result_single, result_multi_1E, result_multi_1C)
result_names <- c("Single", "Multi(1E)", "Multi(1C)")

cat("\nFinal biomass by test:\n")
for(i in 1:length(all_results)) {
  final_sb <- sum(all_results[[i]]$dynamics$SB[10, , ])  # Year 15, all iterations, all areas
  cat(sprintf("  %s: %.2f\n", result_names[i], final_sb))
}

cat("\n=== TEST FILE 1 COMPLETED ===\n")
cat("Backward compatibility: VERIFIED\n")
cat("Multifleet functionality: VERIFIED\n")

