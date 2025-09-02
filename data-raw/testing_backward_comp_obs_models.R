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

# Stochastic object - IDENTICAL for all tests
stochastic_obj <- new("Stochastic")
stochastic_obj@historicalBio = c(0.4, 0.6)  # Small range for consistency
stochastic_obj@Steep = c(0.50, 0.60)        # Small range for consistency

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

# Strategy for multifleet (with fleet column) to maintain maintain data structure consistency with
simpleMP_multi <- function(phase, dataObject) {
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

strategy_single <- new("Strategy")
strategy_single@title <- "Single Fleet Validation"
strategy_single@projectionYears <- 5
strategy_single@projectionName <- "simpleMP_single"
strategy_single@projectionParams <- list()

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
result_single$HCR$decisionData$Fishery_1_total_catch


# ============================================================================
# CREATE OBSERVATION MODELS - FOR MULTIFLEET (2 FLEETS)
# ============================================================================

# -----------------------------
# FLEET-SPECIFIC CPUE (FD)
# -----------------------------
# CPUE Fleet 1 (covers both areas)
cpue_fleet1 <- new("Index")
cpue_fleet1@indexID <- "CPUE_Fleet1"
cpue_fleet1@title <- "Fleet 1 CPUE Biomass"
cpue_fleet1@useWeight <- TRUE  # Biomass-based

cpue_fleet1@survey_design <- list(
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

cpue_fleet1@selectivity_hist_list <- list()
cpue_fleet1@selectivity_proj_list <- list()

# CPUE Fleet 2 (covers both areas)
cpue_fleet2 <- new("Index")
cpue_fleet2@indexID <- "CPUE_Fleet2"
cpue_fleet2@title <- "Fleet 2 CPUE Biomass"
cpue_fleet2@useWeight <- TRUE  # Biomass-based

cpue_fleet2@survey_design <- list(
  list(
    indextype = "FD",
    areas = c(1, 2),
    indexYears = 1:(ta@historicalYears + 5),  # All years
    q_hist_bounds = c(0.00008, 0.00012),  # Slightly different catchability
    q_proj_bounds = c(0.0001, 0.00018),
    hyperstability_hist_bounds = c(0.95, 1.15),
    hyperstability_proj_bounds = c(0.9, 1.1),
    obsError_CV_hist_bounds = c(0.25, 0.35),  # Higher observation error
    obsError_CV_proj_bounds = c(0.2, 0.3)
  )
)

cpue_fleet2@selectivity_hist_list <- list()
cpue_fleet2@selectivity_proj_list <- list()

# -----------------------------
# FISHERY-INDEPENDENT SURVEY (AREA 1 ONLY)
# -----------------------------

# Create survey selectivity (different from fishing)
survey_selectivity <- new("Fishery")
survey_selectivity@vulType <- "logistic"
survey_selectivity@vulParams <- c(7.0, 2.0)  # Smaller sizes than fishery
survey_selectivity@retType <- "full"
survey_selectivity@retMax <- 1
survey_selectivity@Dmort <- 0

fi_survey <- new("Index")
fi_survey@indexID <- "FI_Survey"
fi_survey@title <- "Fishery-Independent Survey Area 1"
fi_survey@useWeight <- TRUE  # Biomass-based

fi_survey@survey_design <- list(
  list(
    indextype = "FI",
    areas = c(1),  # Area 1 only
    indexYears = seq(3, ta@historicalYears + 5, by = 2),  # Every other year starting year 3
    q_hist_bounds = c(0.0002, 0.0003),
    q_proj_bounds = c(0.00025, 0.00035),
    hyperstability_hist_bounds = c(1.0, 1.0),  # No hyperstability for survey
    hyperstability_proj_bounds = c(1.0, 1.0),
    obsError_CV_hist_bounds = c(0.15, 0.2),
    obsError_CV_proj_bounds = c(0.1, 0.15),
    selectivity_hist_idx = 1,  # Index for selectivity object
    selectivity_proj_idx = 1,
    survey_timing = 0.5  # Mid-year survey
  )
)

fi_survey@selectivity_hist_list <- list(survey_selectivity)
fi_survey@selectivity_proj_list <- list(survey_selectivity)

# -----------------------------
# FLEET-SPECIFIC CATCH OBSERVATIONS
# -----------------------------

# Catch observations Fleet 1
catch_obs_fleet1 <- new("CatchObs")
catch_obs_fleet1@catchID <- "Catch_Fleet1"
catch_obs_fleet1@title <- "Fleet 1 Catch Observations"
catch_obs_fleet1@areas <- c(1, 2)
catch_obs_fleet1@catchYears <- 1:(ta@historicalYears + 5)
catch_obs_fleet1@reporting_rates <- rep(0.95, ta@historicalYears + 5)  # 95% reporting
catch_obs_fleet1@obs_CVs <- matrix(c(rep(0.15, ta@historicalYears + 5),
                                     rep(0.25, ta@historicalYears + 5)), ncol=2)

# Catch observations Fleet 2
catch_obs_fleet2 <- new("CatchObs")
catch_obs_fleet2@catchID <- "Catch_Fleet2"
catch_obs_fleet2@title <- "Fleet 2 Catch Observations"
catch_obs_fleet2@areas <- c(1, 2)
catch_obs_fleet2@catchYears <- 1:(ta@historicalYears + 5)
catch_obs_fleet2@reporting_rates <- rep(1.05, ta@historicalYears + 5)  # 105% reporting (over-reporting)
catch_obs_fleet2@obs_CVs <- matrix(c(rep(0.2, ta@historicalYears + 5),
                                     rep(0.3, ta@historicalYears + 5)), ncol=2)

# -----------------------------
# FLEET-SPECIFIC LENGTH COMPOSITION (FD)
# -----------------------------

# Length composition Fleet 1
lc_years_fleet1 <- seq(2, ta@historicalYears + 5, by = 2)  # Every other year

length_comp_fleet1 <- new("LCompObs")
length_comp_fleet1@indexID <- "LC_Fleet1"
length_comp_fleet1@title <- "Fleet 1 Length Composition"
length_comp_fleet1@length_bin_width <- 1  # 1 cm bins

length_comp_fleet1@survey_design <- list(
  list(
    indextype = "FD",
    areas = c(1, 2),
    years = lc_years_fleet1,
    sample_sizes = rep(150, length(lc_years_fleet1))  # 150 samples per event
  )
)

length_comp_fleet1@selectivity_hist_list <- list()
length_comp_fleet1@selectivity_proj_list <- list()

# Length composition Fleet 2
lc_years_fleet2 <- seq(1, ta@historicalYears + 5, by = 3)  # Every third year

length_comp_fleet2 <- new("LCompObs")
length_comp_fleet2@indexID <- "LC_Fleet2"
length_comp_fleet2@title <- "Fleet 2 Length Composition"
length_comp_fleet2@length_bin_width <- 1  # 1 cm bins

length_comp_fleet2@survey_design <- list(
  list(
    indextype = "FD",
    areas = c(1, 2),
    years = lc_years_fleet2,
    sample_sizes = rep(100, length(lc_years_fleet2))  # 100 samples per event
  )
)

length_comp_fleet2@selectivity_hist_list <- list()
length_comp_fleet2@selectivity_proj_list <- list()


# ============================================================================
# TEST SCENARIO 2: MULTIFLEET WITH 2 FLEETS (IDENTICAL SELECTIVITY)
# ============================================================================

cat("============================================\n")
cat("RUNNING TEST 2: MULTIFLEET (2 IDENTICAL SEL)\n")
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
fleet2_hist @vulParams <- c(10.2, 0.1)
fleet2_hist @retType <- "full"
fleet2_hist @retMax <- 1
fleet2_hist @Dmort <- 0


# Projection fleet selectivities
fleet1_proj <- new("Fishery")
fleet1_proj@vulType <- "logistic"
fleet1_proj@vulParams <- c(10.2, 0.1)
fleet1_proj@retType <- "full"
fleet1_proj@retMax <- 1
fleet1_proj@Dmort <- 0

fleet2_proj <- new("Fishery")
fleet2_proj@vulType <- "logistic"
fleet2_proj@vulParams <- c(10.2, 0.1)
fleet2_proj@retType <- "full"
fleet2_proj@retMax <- 1
fleet2_proj@Dmort <- 0


#multifleet object

multifleet_identical <- new("Multifleet")
multifleet_identical@nfleets <- 2
multifleet_identical@fleet_proportions <- c(0.6, 0.4)
multifleet_identical@allocation_type <- "catch"  # Use effort allocation for comparison
multifleet_identical@fleet_selectivity_hist_list  <- list(fleet1_hist, fleet2_hist)
multifleet_identical@fleet_selectivity_proj_list  <- list(fleet1_proj, fleet2_proj)


strategy_multi2 <- new("Strategy")
strategy_multi2@title <- "Multifleet 2 Identical Validation"
strategy_multi2@projectionYears <- 5
strategy_multi2@projectionName <- "simpleMP_multi"
strategy_multi2@projectionParams <- list()

result_multi2 <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  ProFisheryObj_list = proj_fishery_list,
  StrategyObj = strategy_multi2,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_identical,
  wd = getwd(),
  fileName = "validation_multifleet_2identical",
  seed = validation_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
)

cat("Multifleet (2 identical) simulation completed\n")

result_multi2 <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_multifleet_2identical")
result_multi2$dynamics$SB
result_multi2$dynamics$VB

