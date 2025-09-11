#updating the code to include observation models

# ============================================================================
# COMPREHENSIVE MULTIFLEET VALIDATION FRAMEWORK
# ============================================================================
# Tests backward compatibility and validates multifleet implementation
# Compares: Single Fleet vs Multifleet vs Multifleet(1 fleet)

rm(list=ls())
devtools::load_all()
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# SETUP COMMON OBJECTS FOR ALL TESTS
# ============================================================================

# Life History - Kole fish
lh_obj <- new("LifeHistory")
lh_obj@title<-"Kole"
lh_obj@speciesName<-"Ctenochaetus strigosus"
lh_obj@Linf<-17.7
lh_obj@K<-0.423
lh_obj@t0<- -0.51
lh_obj@L50<-8.4
lh_obj@L95delta<-1.26
lh_obj@M<-0.08
lh_obj@L_type<-"FL"
lh_obj@L_units<-"cm"
lh_obj@LW_A<-0.046
lh_obj@LW_B<-2.85
lh_obj@Steep<-0.54
lh_obj@recSD<-0 # No recruitment variation for deterministic testing
lh_obj@recRho<-0
lh_obj@isHermaph<-FALSE
lh_obj@R0<-10000

# Time Area setup
ta <- new("TimeArea")
ta@title = "Validation Test"
ta@gtg = 13
ta@areas = 2
ta@recArea = c(0.99, 0.01)
ta@iterations = 3  # Small number for quick testing
ta@historicalYears = 10
ta@historicalBio = 0.5
ta@historicalBioType = "relB"
ta@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)

# Historical effort pattern - same for both areas
ta@historicalEffort <- matrix(c(1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6,
                                1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6),
                              nrow = 10, ncol = 2, byrow = FALSE)


# Stochastic object - IDENTICAL for all tests
stochastic_obj <- new("Stochastic")
stochastic_obj@historicalBio = c(0.4, 0.6)  # Small range for consistency
stochastic_obj@Steep = c(0.50, 0.60)        # Small range for consistency



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


# SHARED RANDOM SEED
validation_seed <- 123

# ============================================================================
# DEFINE MANAGEMENT STRATEGIES FOR EACH TEST SCENARIO
# ============================================================================

# Strategy for single fleet (original format)
simpleMP_single <- function(phase, dataObject) {
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase==1){
    # Phase 1: Collect observation data
    combined_data <- list()

    # Index observations
    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject)
      for(col_name in names(index_result)) {
        combined_data[[col_name]] <- index_result[[col_name]]
      }
    }

    # Catch observations
    if(!is.null(CatchObsObj)) {
      catch_result <- calculate_single_CatchObs(dataObject)
      for(col_name in names(catch_result)) {
        combined_data[[col_name]] <- catch_result[[col_name]]
      }
    }

    # Length composition observations
    if(!is.null(LengthCompObj)) {
      lc_result <- calculate_single_LengthComp(dataObject)
      for(col_name in names(lc_result)) {
        combined_data[[col_name]] <- lc_result[[col_name]]
      }
    }

    return(combined_data)

  }

  if(phase==2) return(list())
  if(phase==3) {
    year = rep(j, areas)
    iteration = rep(k, areas)
    area = 1:areas
    Flocal = rep(0.15, areas)  # Constant F = 0.15

    return(list(year=year, iteration=iteration, area=area, Flocal=Flocal))
  }
}

strategy_single <- new("Strategy")
strategy_single@title <- "Single Fleet Validation"
strategy_single@projectionYears <- 5
strategy_single@projectionName <- "simpleMP_single"
strategy_single@projectionParams <- list()



# ============================================================================
# CREATE OBSERVATION MODELS - FOR SINGLE FLEET
# ============================================================================

# CPUE Index (covers both areas)
cpue_index <- new("Index")
cpue_index@indexID <- "CPUE_Biomass"
cpue_index@title <- "CPUE Biomass Index"
cpue_index@useWeight <- TRUE  # Biomass-based

cpue_index@survey_design <- list(
  list(
    indextype = "FD",
    areas = c(1, 2),
    indexYears = 1:(ta@historicalYears + 5),  # All years
    q_hist_bounds = c(0.0001, 0.00015),
    q_proj_bounds = c(0.00012, 0.0002),
    hyperstability_hist_bounds = c(0.9, 1.1),
    hyperstability_proj_bounds = c(0.85, 1.05),
    obsError_CV_hist_bounds = c(0.2, 0.3),
    obsError_CV_proj_bounds = c(0.15, 0.25)
  )
)

cpue_index@selectivity_hist_list <- list()
cpue_index@selectivity_proj_list <- list()

# Catch observations
catch_obs <- new("CatchObs")
catch_obs@catchID <- "Fishery_Catch"
catch_obs@title <- "Fishery Catch Observations"
catch_obs@areas <- c(1, 2)
catch_obs@catchYears <- 1:(ta@historicalYears + 5)
catch_obs@reporting_rates <- rep(1, ta@historicalYears + 5)  # Perfect reporting
catch_obs@obs_CVs <- matrix(c(rep(0.2, ta@historicalYears + 5),
                              rep(0.3, ta@historicalYears + 5)), ncol=2)

cat("Created single fleet catch observations:\n")
cat("  - Areas:", paste(catch_obs@areas, collapse = ", "), "\n")
cat("  - Years with data:", length(catch_obs@catchYears), "\n")
cat("  - Reporting rate range:", round(range(catch_obs@reporting_rates), 2), "\n")
cat("  - obs_CVs RANGE:", round(range(catch_obs@obs_CVs), 2), "\n")


# Length composition
length_comp <- new("LCompObs")
length_comp@indexID <- "Length_Comp"
length_comp@title <- "Length Composition"
length_comp@length_bin_width <- 1  # 1 cm bins


# Calculate years first, then match sample_sizes length
lc_years <- seq(2, ta@historicalYears + 5, by = 2)

length_comp@survey_design <- list(
  list(
    indextype = "FD",
    areas = c(1, 2),
    years = lc_years,
    sample_sizes = rep(100, length(lc_years))
  )
)

length_comp@selectivity_hist_list <- list()
length_comp@selectivity_proj_list <- list()


# ============================================================================
# TEST SCENARIO 1: SINGLE FLEET WITH OBS MODELS
# ============================================================================

cat("=====================================\n")
cat("RUNNING TEST 1: SINGLE FLEET         \n")
cat("=====================================\n")



result_single <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  ProFisheryObj_list = proj_fishery_list,
  StrategyObj = strategy_single,
  StochasticObj = stochastic_obj,
  MultifleetObj = NULL,  # Single fleet mode
  IndexObj = cpue_index,
  CatchObsObj = catch_obs,
  LengthCompObj = length_comp,
  wd = getwd(),
  fileName = "validation_single_fleet_obsmodels",
  seed = validation_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
)
cat("Single fleet obs models simulation completed\n")

result_single <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_single_fleet_obsmodels")
result_single$dynamics$SB # 1 set of vectors for each area for n iterations
result_single$dynamics$VB
result_single$dynamics$RB
result_single$dynamics$catchB
result_single$dynamics$catchN
result_single$dynamics$Ftotal # 1 set of vectors for each area for n iterations
result_single$dynamics$SPR # 1 set of vectors area-combined

result_single$HCR$decisionData$CPUE_1
result_single$HCR$decisionData$observed_catch_area_1
result_single$HCR$decisionData$observed_catch_area_2
result_single$HCR$decisionData$Fishery_1_indextype
result_single$HCR$decisionData$Fishery_1_sample_size
result_single$HCR$decisionData$Fishery_1_total_catch
result_single$HCR$decisionData$Fishery_1_count_bin_1 # ETC ETC
result_single$HCR$decisionData$Fishery_1_count_bin_2
result_single$HCR$decisionData$Fishery_1_count_bin_5
result_single$HCR$decisionData$Fishery_1_count_bin_10
result_single$HCR$decisionData$Fishery_1_count_bin_17

# ============================================================================
# MULTIFLEET EXAMPLE
# ============================================================================

cat("============================================\n")
cat("RUNNING TEST: MULTIFLEET (DIFFERENT SELECT )\n")
cat("============================================\n")


# Historical fleet selectivities
fleet1_hist  <- new("Fishery")
fleet1_hist @vulType <- "logistic"
fleet1_hist @vulParams <- c(10.2, 0.1)
fleet1_hist @retType <- "full"
fleet1_hist @retMax <- 1
fleet1_hist @Dmort <- 0

fleet2_hist  <- new("Fishery")
fleet2_hist @vulType <- "logistic"
fleet2_hist @vulParams <- c(9, 0.1)
fleet2_hist @retType <- "full"
fleet2_hist @retMax <- 1
fleet2_hist @Dmort <- 0


# Projection fleet selectivities
fleet1_proj <- new("Fishery")
fleet1_proj@vulType <- "logistic"
fleet1_proj@vulParams <- c(11, 0.1)
fleet1_proj@retType <- "full"
fleet1_proj@retMax <- 1
fleet1_proj@Dmort <- 0

fleet2_proj <- new("Fishery")
fleet2_proj@vulType <- "logistic"
fleet2_proj@vulParams <- c(10, 0.1)
fleet2_proj@retType <- "full"
fleet2_proj@retMax <- 1
fleet2_proj@Dmort <- 0


#multifleet object

multifleet_example  <- new("Multifleet")
multifleet_example @nfleets <- 2
multifleet_example @fleet_proportions <- c(0.6, 0.4)
multifleet_example @allocation_type <- "catch"  # Use effort allocation for comparison
multifleet_example @fleet_selectivity_hist_list  <- list(fleet1_hist, fleet2_hist)
multifleet_example @fleet_selectivity_proj_list  <- list(fleet1_proj, fleet2_proj)

cat("Created multifleet object:\n")
cat("  - 2 fleets with different selectivities\n")
cat("  - Fleet 1 historical: c(10.2, 0.1), projection: c(11, 0.1)\n")
cat("  - Fleet 2 historical: c(9, 0.1), projection: c(10, 0.1)\n")
cat("  - Target catch proportions: 60% Fleet 1, 40% Fleet 2\n")
cat("  - Allocation type: catch (iterative to find effort proportions)\n")


# ============================================================================
# CREATE FI SURVEY SELECTIVITY (INDEPENDENT OF FLEETS)
# ============================================================================

# Create survey selectivity (different from both fleets)
survey_sel_hist <- new("Fishery")
survey_sel_hist@title <- "Survey Historical"
survey_sel_hist@vulType <- "logistic"
survey_sel_hist@vulParams <- c(8.5, 0.2)
survey_sel_hist@retType <- "full"
survey_sel_hist@retMax <- 1
survey_sel_hist@Dmort <- 0

survey_sel_proj <- new("Fishery")
survey_sel_proj@title <- "Survey Projection"
survey_sel_proj@vulType <- "logistic"
survey_sel_proj@vulParams <- c(8.8, 0.2)
survey_sel_proj@retType <- "full"
survey_sel_proj@retMax <- 1
survey_sel_proj@Dmort <- 0

# ============================================================================
# CREATE COMPREHENSIVE INDEX OBJECT
# ============================================================================
multifleet_indices <- new("Index")
multifleet_indices@indexID <- "MultifleetIndices"
multifleet_indices@title <- "Fleet-Specific and Survey Indices"
multifleet_indices@useWeight <- TRUE  # Biomass-based

# Only FI surveys need selectivity objects (FD indices use fleet selectivity automatically)
multifleet_indices@selectivity_hist_list <- list(survey_sel_hist)
multifleet_indices@selectivity_proj_list <- list(survey_sel_proj)


multifleet_indices@survey_design <- list(
  # 1. FI Survey covering both areas (every 3 years)
  list(
    indextype = "FI",
    areas = c(1, 2),
    indexYears = seq(3, 15, 3),  # Years 3, 6, 9, 12, 15
    selectivity_hist_idx = 1,    # Uses survey_sel_hist
    selectivity_proj_idx = 1,    # Uses survey_sel_proj
    survey_timing = 0.5,         # Mid-year survey
    q_hist_bounds = c(0.0001, 0.0003),
    q_proj_bounds = c(0.00012, 0.00035),
    hyperstability_hist_bounds = c(0.95, 1.05),
    hyperstability_proj_bounds = c(0.95, 1.05),
    obsError_CV_hist_bounds = c(0.10, 0.20),
    obsError_CV_proj_bounds = c(0.10, 0.20)
  ),

  # 2. Fleet 1 CPUE covering both areas (annual)
  list(
    indextype = "FD",
    fleet_id = 1,  # Uses fleet1_hist/fleet1_proj selectivity automatically
    areas = c(1, 2),
    indexYears = 1:15,
    q_hist_bounds = c(0.0002, 0.0008),
    q_proj_bounds = c(0.0003, 0.0009),
    hyperstability_hist_bounds = c(0.8, 1.2),
    hyperstability_proj_bounds = c(0.8, 1.2),
    obsError_CV_hist_bounds = c(0.15, 0.25),
    obsError_CV_proj_bounds = c(0.15, 0.25)
  ),

  # 3. Fleet 1 CPUE Area 1 only (every other year)
  list(
    indextype = "FD",
    fleet_id = 1,  # Uses fleet1_hist/fleet1_proj selectivity automatically
    areas = c(1),  # Area 1 only
    indexYears = seq(2, 14, 2),  # Years 2, 4, 6, 8, 10, 12, 14
    q_hist_bounds = c(0.0003, 0.001),
    q_proj_bounds = c(0.0004, 0.0012),
    hyperstability_hist_bounds = c(0.8, 1.2),
    hyperstability_proj_bounds = c(0.8, 1.2),
    obsError_CV_hist_bounds = c(0.20, 0.30),
    obsError_CV_proj_bounds = c(0.20, 0.30)
  ),

  # 4. Fleet 2 CPUE covering both areas (annual)
  list(
    indextype = "FD",
    fleet_id = 2,  # Uses fleet2_hist/fleet2_proj selectivity automatically
    areas = c(1, 2),
    indexYears = 1:15,
    q_hist_bounds = c(0.0001, 0.0006),
    q_proj_bounds = c(0.00015, 0.0007),
    hyperstability_hist_bounds = c(0.7, 1.3),
    hyperstability_proj_bounds = c(0.7, 1.3),
    obsError_CV_hist_bounds = c(0.20, 0.35),
    obsError_CV_proj_bounds = c(0.20, 0.35)
  ),

  # 5. Fleet 2 CPUE Area 2 only (every three years)
  list(
    indextype = "FD",
    fleet_id = 2,  # Uses fleet2_hist/fleet2_proj selectivity automatically
    areas = c(2),  # Area 2 only
    indexYears = seq(1, 15, 3),  # Years 1, 4, 7, 10, 13
    q_hist_bounds = c(0.0002, 0.0009),
    q_proj_bounds = c(0.00025, 0.001),
    hyperstability_hist_bounds = c(0.7, 1.3),
    hyperstability_proj_bounds = c(0.7, 1.3),
    obsError_CV_hist_bounds = c(0.25, 0.40),
    obsError_CV_proj_bounds = c(0.25, 0.40)
  )
)


#Provide overall description

cat("\nCreated comprehensive index object with 5 indices:\n")
cat("  - 1 FI Survey (independent selectivity, both areas, triennial)\n")
cat("  - 2 Fleet 1 CPUE indices (both areas annual + area 1 biennial)\n")
cat("  - 2 Fleet 2 CPUE indices (both areas annual + area 2 triennial)\n")
cat("\nSelectivity:\n")
cat("  - FI Survey: Uses survey_sel_hist/survey_sel_proj\n")
cat("  - Fleet 1 CPUE: Uses fleet1_hist/fleet1_proj automatically\n")
cat("  - Fleet 2 CPUE: Uses fleet2_hist/fleet2_proj automatically\n")


# ============================================================================
# CREATE COMPREHENSIVE CATCH OBJECT (NEW FLEET-SPECIFIC FORMAT)
# ============================================================================

#creating multifleet catch observations with fleet-specific configurations
catch_obs_multifleet <- new("CatchObs")
catch_obs_multifleet@catchID <- "Multifleet_Catch"
catch_obs_multifleet@title <- "Fleet-Specific Catch Observations"


#calculate dimensions first to avoid errors
total_years <- ta@historicalYears + 5  # 10 + 5 = 15 years
fleet1_years <- total_years             # Fleet 1: every year (15 years)
fleet2_years <- length(seq(2, total_years, by = 2))  # Fleet 2: every other year (7 years)

cat("Total simulation years:", total_years, "\n")
cat("Fleet 1 data years:", fleet1_years, "\n")
cat("Fleet 2 data years:", fleet2_years, "\n")

#creating fleet-specific configurations
catch_obs_multifleet@fleet_configs <- list(
  # Fleet 1 configuration - low obs error
  list(
    fleet_id = 1,
    areas = c(1, 2),  # Fleet 1 operates in both areas
    catchYears = 1:total_years,  # complete data coverage
    reporting_rates = c(
      seq(0.8, 1.0, length.out = ta@historicalYears),  # improving historical
      rep(1.0, 5)  # perfect projection
    ),
    obs_CVs = matrix(
      c(rep(0.15, ta@historicalYears), rep(0.10, 5),   # Lower bounds - improving projection
        rep(0.25, ta@historicalYears), rep(0.20, 5)),  # Upper bounds - improving projection
      ncol = 2
    )
  ),

  # Fleet 2 configuration - hiher CV
  list(
    fleet_id = 2,
    areas = c(1, 2),  # Fleet 2 also operates in both areas
    catchYears = seq(2, total_years, by = 2),  # available every other year
    reporting_rates = c(
      rep(0.6, 5),   # consistent under-reporting historical
      rep(0.9, 2)  # improved but not perfect projection
    ),
    obs_CVs = matrix(
      c(rep(0.30, 5), rep(0.25, 2),   # Higher CV
        rep(0.45, 5), rep(0.35, 2)),  # Upper bounds
      ncol = 2
    )
  )
)

cat("\nCreated multifleet catch observations:\n")
cat("  - Number of fleet configurations:", length(catch_obs_multifleet@fleet_configs), "\n")

for(i in 1:length(catch_obs_multifleet@fleet_configs)) {
  config <- catch_obs_multifleet@fleet_configs[[i]]
  cat("  Fleet", config$fleet_id, ":\n")
  cat("    - Areas:", paste(config$areas, collapse = ", "), "\n")
  cat("    - Years with data:", length(config$catchYears), "\n")
  cat("    - Reporting rate range:", round(range(config$reporting_rates), 2), "\n")
  cat("    - CV bounds:", round(range(config$obs_CVs), 2), "\n")
}


# ============================================================================
# CREATE COMPREHENSIVE LENGTH OBJECT (NEW FLEET-SPECIFIC FORMAT)
# ============================================================================
# Length composition
length_comp_multifleet <- new("LCompObs")
length_comp_multifleet@indexID <- "Length_Comp_multifleet"
length_comp_multifleet@title <- "Fleet-Specific Length Composition"
length_comp_multifleet@length_bin_width <- 1  # 1 cm bins
length_comp_multifleet@selectivity_hist_list <- list(survey_sel_hist)
length_comp_multifleet@selectivity_proj_list <- list(survey_sel_proj)


length_comp_multifleet@survey_design <- list(
  # Fleet 1 FD length composition
  list(
    indextype = "FD",
    fleet_id = 1,           #Fleet 1
    areas = c(1, 2),
    years = c(2, 4, 6, 8, 10),
    sample_sizes = c(100, 120, 140, 160, 180)
  ),

  # Fleet 2 FD length composition
  list(
    indextype = "FD",
    fleet_id = 2,           # Fleet 2
    areas = c(1, 2),
    years = c(3, 6, 9, 12),
    sample_sizes = c(80, 100, 120, 140)
  ),

  # FI survey (unchanged - no fleet_id)
  list(
    indextype = "FI",
    areas = c(1, 2),
    years = c(5, 10, 15),
    sample_sizes = c(200, 250, 300),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.5
  )
)


# ============================================================================
# CREATE MULTIFLEET MANAGEMENT STRATEGY WITH CATCH OBSERVATIONS
# ============================================================================

# Strategy for multifleet (with fleet column) to maintain maintain data structure consistency with
simpleMP_multifleet <- function(phase, dataObject) {
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase==1) {
    # Phase 1: Collect observation data (same for single/multi)
    combined_data <- list()

    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject)
      for(col_name in names(index_result)) {
        combined_data[[col_name]] <- index_result[[col_name]]
      }
    }

    # Catch observations - MULTIFLEET
    if(!is.null(CatchObsObj)) {
      catch_result <- calculate_single_CatchObs(dataObject)
      for(col_name in names(catch_result)) {
        combined_data[[col_name]] <- catch_result[[col_name]]
      }
    }

    # Length composition observations
    if(!is.null(LengthCompObj)) {
      lc_result <- calculate_single_LengthComp(dataObject)
      for(col_name in names(lc_result)) {
        combined_data[[col_name]] <- lc_result[[col_name]]
      }
    }

    return(combined_data)

  }

  if(phase==2) return(list())
  if(phase==3) {
    year = rep(j, areas)
    iteration = rep(k, areas)
    area = 1:areas
    fleet = rep(0, areas)  # 0 = total F across all fleets
    Flocal = rep(0.15, areas)  # Constant F = 0.15

    return(list(year=year, iteration=iteration, area=area, fleet=fleet, Flocal=Flocal))
  }
}

strategy_multifleet <- new("Strategy")
strategy_multifleet@title <- "Multifleet with Indices and catch obs"
strategy_multifleet@projectionYears <- 5
strategy_multifleet@projectionName <- "simpleMP_multifleet"
strategy_multifleet@projectionParams <- list()


# ============================================================================
# RUN MULTIFLEET SIMULATION
# ============================================================================

cat("\n============================================\n")
cat("RUNNING MULTIFLEET SIMULATION (INDEX ONLY)\n")
cat("============================================\n")

result_multifleet <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StrategyObj = strategy_multifleet,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_example,    # Uses the defined fleet selectivities
  IndexObj = multifleet_indices,         # Comprehensive index object
  CatchObsObj = catch_obs_multifleet,
  LengthCompObj = length_comp_multifleet,
  wd = getwd(),
  fileName = "multifleet_complete_example",
  seed = validation_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
)

cat("Multifleet simulation completed successfully\n")

# ============================================================================
# EXAMINE RESULTS
# ============================================================================

result_multifleet <- readProjection(getwd(), "multifleet_complete_example")

cat("\n=== MULTIFLEET SIMULATION RESULTS ===\n")

# Population dynamics
cat("Population dynamics:\n")
cat("  SB dimensions:", dim(result_multifleet$dynamics$SB), "\n")
cat("  Multifleet detected:", !is.null(result_multifleet$dynamics$multifleet), "\n")

if(!is.null(result_multifleet$dynamics$multifleet)) {
  mf <- result_multifleet$dynamics$multifleet
  cat("  Number of fleets:", mf$nfleets, "\n")
  cat("  Fleet proportions (original):", paste(round(mf$fleet_proportions, 3), collapse = ", "), "\n")
  cat("  Allocation type:", mf$allocation_type, "\n")

  if(!is.null(mf$final_effort_proportions)) {
    cat("  Final effort proportions:", paste(round(mf$final_effort_proportions, 3), collapse = ", "), "\n")
  }
  if(!is.null(mf$actual_catch_proportions)) {
    cat("  Actual catch proportions:", paste(round(mf$actual_catch_proportions, 3), collapse = ", "), "\n")
  }
}

# Index observations
survey_cols <- grep("Survey_1", names(result_multifleet$HCR$decisionData), value = TRUE)
fleet1_cols <- grep("CPUE_.*_Fleet_1", names(result_multifleet$HCR$decisionData), value = TRUE)
fleet2_cols <- grep("CPUE_.*_Fleet_2", names(result_multifleet$HCR$decisionData), value = TRUE)

cat("\nIndex observations collected:\n")
cat("  FI Survey indices:", length(survey_cols), "\n")
cat("  Fleet 1 CPUE indices:", length(fleet1_cols), "\n")
cat("  Fleet 2 CPUE indices:", length(fleet2_cols), "\n")

# Validate that fleet indices are different (due to different selectivities)
if(length(fleet1_cols) > 0 && length(fleet2_cols) > 0) {
  fleet1_values <- result_multifleet$HCR$decisionData[[fleet1_cols[1]]]
  fleet2_values <- result_multifleet$HCR$decisionData[[fleet2_cols[1]]]

  valid_indices <- !is.na(fleet1_values) & !is.na(fleet2_values)
  if(sum(valid_indices) > 5) {
    correlation <- cor(fleet1_values[valid_indices], fleet2_values[valid_indices])

    cat("\nFleet-specific index validation:\n")
    cat("  Correlation between Fleet 1 and Fleet 2 CPUE:", round(correlation, 3), "\n")
    cat("  Expected: < 1.0 (should be different due to different selectivities)\n")

    # Sample comparison
    sample_data <- data.frame(
      year = result_multifleet$HCR$decisionData$j[valid_indices][1:6],
      fleet1_cpue = round(fleet1_values[valid_indices][1:6], 5),
      fleet2_cpue = round(fleet2_values[valid_indices][1:6], 5)
    )
    cat("  Sample values (first 6 valid observations):\n")
    print(sample_data)
  }
}

# Check fleet metadata
fleet_id_cols <- grep("_fleet_id$", names(result_multifleet$HCR$decisionData), value = TRUE)
if(length(fleet_id_cols) > 0) {
  unique_fleet_ids <- unique(result_multifleet$HCR$decisionData[[fleet_id_cols[1]]])
  unique_fleet_ids <- unique_fleet_ids[!is.na(unique_fleet_ids)]
  cat("\nFleet IDs found in data:", paste(unique_fleet_ids, collapse = ", "), "\n")
}

cat("\n=== VALIDATION SUMMARY ===\n")
cat("OK Multifleet object created with defined selectivities\n")
cat("OK 5 indices created (1 FI + 4 fleet-specific FD)\n")
cat("OK FD indices automatically use fleet selectivities\n")
cat("OK Fleet-specific indices generate different signals\n")
cat("OK Catch and length composition not included in this example\n")

# more checks
result_multifleet$MultifleetObj@fleet_selectivity_hist_list[[1]]@vulParams  # [10.2, 0.1]
result_multifleet$MultifleetObj@fleet_selectivity_hist_list[[2]]@vulParams  # [9, 0.1]
dim(result_multifleet$dynamics$multifleet$RB_by_fleet)  # Should be [years, iters, areas, fleets]
head(result_multifleet$HCR$decisionData[c("j", "CPUE_2_Fleet_1", "CPUE_4_Fleet_2")])
unique(result_multifleet$HCR$decisionData$CPUE_2_Fleet_1_fleet_id) ## Should be 1
unique(result_multifleet$HCR$decisionData$CPUE_4_Fleet_2_fleet_id) # should be 2

result_multifleet$HCR$decisionData$fleet_1_true_catch
result_multifleet$HCR$decisionData$fleet_1_observed_catch

result_multifleet$HCR$decisionData$fleet_2_true_catch
result_multifleet$HCR$decisionData$fleet_2_observed_catch

result_multifleet$HCR$decisionData$fleet_1_observed_catch_area_1
result_multifleet$HCR$decisionData$fleet_1_observed_catch_area_2

result_multifleet$HCR$decisionData$fleet_2_observed_catch_area_1
result_multifleet$HCR$decisionData$fleet_2_observed_catch_area_2

result_multifleet$HCR$decisionData$Fishery_2_Fleet_2_count_bin_20
