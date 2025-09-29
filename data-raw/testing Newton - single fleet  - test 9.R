#testing Newton - single fleet (WIP)


rm(list=ls())
#devtools::document()  #to generate NAMESPACE and .Rd files
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
ta@iterations <- 6
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

  #For a simple testing I am setting a TAC based on average recent catch
  #This runs once per year/iteration combination, after all Phase 1 data collection is complete
  if(phase == 2) {

  #calculate the last year of the historical period
  yrHist <- TimeAreaObj@historicalYears + 1  # End of historical period
  n_years_avg <- 3
  start_year <- max(2, j - n_years_avg) #ensures we never go before year 2 (prevent to going year 1 or earlier)

  #calculate TAC for each area
  TAC_by_area <- numeric(areas) #empty vector to store TAC by area

  for(m in 1:areas) {
    #extract recent observed catches for this area
    #create a col named "observed_catch_area_1" or "observed_catch_area_2"
    #decisionData: data frame built by Phase 1, with one row per year/iteration/area
    #for example:
    #If j=12, k=1, m=1, start_year=9:
    #extracts observed catches from Area 1, iteration 1, years 9, 10, 11
    #these are the 3 most recent years of catch data available


    recent_catches <- decisionData[[paste0("observed_catch_area_", m)]][
      which(decisionData$k == k &
              decisionData$j >= start_year &
              decisionData$j < j) #excludes current year (cannot use data because data have not been collected yet!)
    ]

    #if catch data is available use it, otherwise fix a value
    if(length(recent_catches) > 0 && !all(is.na(recent_catches))) {
      TAC_by_area[m] <- mean(recent_catches, na.rm = TRUE) #average, ignoring any NAs
    } else {
      #use a small fraction of current retained biomass at currecnt  year j, iteration k, area m
      #first projection year if catch observations were not collected in years 9-11
      TAC_by_area[m] <- 0.10 * RB[j, k, m]
    }

    #apply a 10% buffer  (precautionary appraoch) because catch obs have errors
    TAC_by_area[m] <- TAC_by_area[m] * 0.9
  }

   #return the decision for this year/iteration/area combination
  return(list(
    year = rep(j, areas),    #current year for both areas
    iteration = rep(k, areas), #current iteration for both areas
    area = 1:areas,            #area identifiers
    TAC = TAC_by_area          # TAC for Area 1 and Area 2
  ))
  }

  #All of this will get stored in decisionAnnual data frame and will be available in Phase 3

  #phase 3 called once per year/iter
  if(phase == 3) {
    #initialize F vector for each area
    F_by_area <- numeric(areas)

    #extract TAC decisions for this year/iteration made in phase 2
    #decisionAnnual is a data frame that accumulates all Phase 2 decisions

    #TAC_decisions dataframe:
    #year iteration area  TAC
    #12         1    1 42.0
    #12         1    2 38.5



    TAC_decisions <- decisionAnnual[decisionAnnual$year == j &
                                      decisionAnnual$iteration == k, ]

    #process each area independently
    for(m in 1:areas) {
      #get TAC for this specific area
      TAC_target <- TAC_decisions$TAC[TAC_decisions$area == m]

      #prepare inputs for solveTAC_to_F_fishSimGTG
      #for single fleet: TAC_targets is a single value for the solver
      TAC_targets <- c(TAC_target)

      #no effort-managed fleets in this case
      effort_F_by_fleet <- numeric(0)

      #call Newton-Raphson solver
      F_result <- solveTAC_to_F_fishSimGTG(
        j = j,
        k = k,
        TAC_targets = TAC_targets,
        N = N,
        lh = lh,
        selGroup = selGroup[[m]],  #area-specific selectivity
        M_rate = lh$LifeHistory@M,
        effort_F_by_fleet = effort_F_by_fleet,
        is_multifleet = FALSE,
        areas = 1,  #process one area at a time
        nfleets = 1,
        TAC_type = "keep",  #use retained catch
        control = list(maxiterF = 300, tolF = 1e-4)
      )

      #extract F value (single fleet returns scalar)
      F_by_area[m] <- F_result

      #check convergence
      if(!is.null(attr(F_result, "converged")) && !attr(F_result, "converged")) {
        warning(paste("Newton-Raphson did not converge for area", m,
                      "year", j, "iteration", k))
      }
    }

    #return F values for all areas
    return(list(
      year = rep(j, areas),
      iteration = rep(k, areas),
      area = 1:areas,
      Flocal = F_by_area
    ))
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

