#testing Newton - single fleet
#validation


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
ta@historicalEffort <- matrix(c(1.5, 3.4, 3.3, 3.8, 3.1, 3.0, 2.9, 1.8, 1.7, 1.6,
                                1.2, 2.4, 2.3, 2.9, 2.1, 2.0, 1.9, 0.8, 0.7, 0.6),
                              nrow = 10, ncol = 2, byrow = FALSE)

plot(ta@historicalEffort[,1])
plot(ta@historicalEffort[,2])

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
proj_fishery_area2@vulParams<-c(8.5, 0.1)
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

#adding a global variable to explore the outputs of Newton-Rapson iterations
nr_diagnostics_storage <<- data.frame()

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
  yrHist <- TimeAreaObj@historicalYears + 1  # the end of historical period, e.g., 11
  n_years_avg <- 3
  # e.g., if j=12, start_year=9
  start_year <- max(2, j - n_years_avg) #ensures we never go before year 2 (prevent to going year 1 or earlier)

  #calculate TAC for each area
  TAC_by_area <- numeric(areas) #empty vector to store TAC by area

  #define area-specific multipliers
  area_multipliers <- c(0.8, 1.6)  # Area 1: reduce 20%, Area 2: increase 20%

  #diagonstic
  cat(sprintf("\n=== PHASE 2 DEBUG: Year %d, Iteration %d ===\n", j-1, k))


  for(m in 1:areas) {
    #extract recent observed catches for this area
    recent_catches <- decisionData[[paste0("observed_catch_area_", m)]][
      which(decisionData$k == k &
              decisionData$j >= start_year &        #Years 9, 10, ...
              decisionData$j < j) #excludes current year (cannot use data because data have not been collected yet!)
    ]

    #example for Area 1, iteration 1, year 12:
    #recent_catches = [catch_year9_area1, catch_year10_area1, catch_year11_area1]

    # Why j < j (not j <= j)?
    # year j=12 is currently being calculated
    # we do not have catch observations for year 12 yet (they happen in Phase 1 of next year)
    #So we can only use years 9, 10, 11



    #show what data we are susing
    cat(sprintf("Area %d recent catches: %s\n",
                m, paste(round(recent_catches, 2), collapse = ", ")))


    #Calculate base TAC
    if(length(recent_catches) > 0 && !all(is.na(recent_catches))) {
      base_TAC <- mean(recent_catches, na.rm = TRUE)
      cat(sprintf("Area %d: Using observed catches, base_TAC = %.2f\n", m, base_TAC))
    } else {
      base_TAC <- 0.10 * RB[j, k, m]
      cat(sprintf("Area %d: No observed catches, using fallback = %.2f\n", m, base_TAC))
    }

    # Apply area-specific multiplier
    TAC_by_area[m] <- base_TAC * area_multipliers[m]
    cat(sprintf("Area %d: Final TAC (after multiplier) = %.2f\n\n", m, TAC_by_area[m]))
  }

  return(list(
    year = rep(j, areas),                # e.g., [12, 12]
    iteration = rep(k, areas),           # e.g., [1, 1]
    area = 1:areas,                      # [1, 2]
    TAC = TAC_by_area                    # e.g., [80, 120]
  ))
  }#closing phase 2


  #All of this will get stored in decisionAnnual data frame and will be available in Phase 3

  #phase 3 called once per year/iter (runs once per year per iteration during the projection period, after Phase 2)
  if(phase == 3) {
    # Initialize F vector for each area
    F_by_area <- numeric(areas)           # will store [F_area1, F_area2]


    #extract TAC decisions for this year/iteration (from Phase 2)
    TAC_decisions <- decisionAnnual[decisionAnnual$year == j &
                                      decisionAnnual$iteration == k, ]

    # Process each area independently (process area 1, then area 2)
    for(m in 1:areas) {
      # Get TAC for this area
      TAC_target <- TAC_decisions$TAC[TAC_decisions$area == m]

      # For m=1: TAC_target = 80
      # For m=2: TAC_target = 120




      # Prepare area-specific data
      # (create temporary N arrays with ONLY this area data)
      N_temp <- lapply(1:lh$gtg, function(gtg_idx) {
        temp_array <- array(0, dim = c(dim(N[[gtg_idx]])[1],
                                       dim(N[[gtg_idx]])[2],
                                       1))                   #note: only 1 area!
        #ocpy data from area m into position 1
        temp_array[, , 1] <- N[[gtg_idx]][, , m]
        # If m=1: copy area 1 data
        # If m=2: copy area 2 data
        temp_array
      })

      #create temporary selectivity with ONLY this area selectivity
      selGroup_temp <- list(selGroup[[m]])
      # If m=1: selGroup_temp[[1]] = selectivity for area 1
      # If m=2: selGroup_temp[[1]] = selectivity for area 2


      #why this coping?
      # - the solver expects data structure: N[[gtg]][age, year, area]
      # - but... we want to solve for ONE area at a time
      # - so I create a "fake" 1-area structure containing just the current area's data
      # - the solver sees areas = 1 and processes area = 1 (which is actually area m)



      # #debug
      # cat(sprintf("\n=== DEBUG Area %d ===\n", m))
      # cat("N_temp abundance check (GTG 1, Age 10, Year j, Area 1):", N_temp[[1]][10, j, 1], "\n")
      # cat("Original N (GTG 1, Age 10, Year j, Area", m, "):", N[[1]][10, j, m], "\n")
      # cat("Match?", N_temp[[1]][10, j, 1] == N[[1]][10, j, m], "\n")
      # cat("Selectivity check (keep, GTG 1, Age 10):", selGroup_temp[[1]]$keep[[1]][10], "\n")
      # cat("Original selectivity (keep, GTG 1, Age 10, Area", m, "):", selGroup[[m]]$keep[[1]][10], "\n")
      # cat("==================\n")



      # Call Newton-Raphson solver with area-specific data
      F_result <- solveTAC_to_F_fishSimGTG(
        j = j,                                # Year 12
        k = k,                                # Iteration 1
        TAC_targets = c(TAC_target),          # e.g., for m=1: [80], for m=2: [120]
        N = N_temp,                           # abundance for THIS area only
        lh = lh,                              # LH
        selGroup = selGroup_temp,             # selectivity for THIS area only
        M_rate = lh$LifeHistory@M,
        effort_F_by_fleet = numeric(0),
        is_multifleet = FALSE,
        areas = 1,                            # solver think there is  only 1 area (the current area m)
        nfleets = 1,                          # working with single fleet mode
        TAC_type = "keep",
        control = list(maxiterF = 300, tolF = 1e-4)
      )

      cat(sprintf("PHASE3 DEBUG: Area %d received F = %.6f from solver\n", m, F_result))

      #extract F value (F_result is a scalar: e.g., 0.15 for area 1, 0.22 for area 2)
      F_by_area[m] <- F_result



      # What is happening inside the newton-rapson solver (for Area 1, TAC=80):
      # initial guess: F ≈ 80 / vulnerable_biomass_area1 ≈ 0.12
      # iIteration 1: Predict catch with F=0.12 > get 75 > too low
      # derivative says: increase F by 0.03
      # iteration 2: Predict catch with F=0.15 > get 80.1 → close enough to converge
      # the fucntion return F=0.15



      #store Newton-Raphson diagnostics (all from solver attributes)
      nr_diagnostics_storage <<- rbind(nr_diagnostics_storage, data.frame(
        year = j,
        iteration = k,
        area = m,
        TAC = TAC_target,
        initial_F_guess = ifelse(is.null(attr(F_result, "initial_guess")),
                                 NA_real_, as.numeric(attr(F_result, "initial_guess"))), # e.g., 0.12
        final_F = as.numeric(F_result),                                                  # e.g., 0.15
        converged = ifelse(is.null(attr(F_result, "converged")),
                           FALSE, attr(F_result, "converged")),
        nr_iterations = ifelse(is.null(attr(F_result, "iterations")),
                               NA_integer_, as.integer(attr(F_result, "iterations"))),
        nr_error = ifelse(is.null(attr(F_result, "final_error")),
                          NA_real_, as.numeric(attr(F_result, "final_error")[1])),
        predicted_catch = ifelse(is.null(attr(F_result, "predicted_catch")),
                                 NA_real_, as.numeric(attr(F_result, "predicted_catch")[1])),
        target_catch = TAC_target,
        stringsAsFactors = FALSE
      ))
    }

    #return F values for all areas (standard structure that match fixedStrategy)
    #result stored in "decisionLocal"
    return(list(
      year = rep(j, areas),
      iteration = rep(k, areas),
      area = 1:areas,
      Flocal = F_by_area

    ))
  } #closing phase 3
} #closing startegy

#Single fleet Phase 3: Returns 2 F values (one per area)



# Each year:
# Phase 1: Observe catches from previous year
# Phase 2: Calculate TAC based on rolling 3-year average
# Phase 3: Convert TAC to F using Newton-Raphson
# Apply F: Population changes, produces new catches

# Key structure to keep in mind:
# TAC is forward-looking: TAC for year j is set based on catches from years j-3, j-2, j-1
# Solver works area-by-area: Each area gets its own independent F calculation
# Newton-Raphson iterates: Starts with a guess, refines until predicted catch matches target
# Baranov equation: Accounts for both fishing and natural mortality simultaneously
# Data structures matter in NR: so we need to copy area-specific data into a "1-area structure" for the solver
# Diagnostics track convergence: Each solver run records initial guess, final F, iterations, error

#======================The End==================================#


# Create new TAC-based strategy
strategy_single_tac <- new("Strategy")
strategy_single_tac@title <- "Single Fleet TAC with Newton-Raphson"
strategy_single_tac@projectionYears <- 5
strategy_single_tac@projectionName <- "singleCompMP"
strategy_single_tac@projectionParams <- list()  # No additional params needed

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
  StrategyObj = strategy_single_tac, #updated
  StochasticObj = stochastic_obj,
  MultifleetObj = NULL,
  IndexObj = single_comprehensive_index,
  CatchObsObj = single_comprehensive_catch,
  LengthCompObj = single_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test1_single_with_tac",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "singleCompMP"
)

cat("Single fleet comprehensive simulation with obs models completed\n")

result_single_comp <- readProjection(getwd(), "test1_single_with_tac")

# Population outputs
result_single_comp$dynamics$SB
result_single_comp$dynamics$VB
result_single_comp$dynamics$Ftotal
result_single_comp$dynamics$catchB
result_single_comp$dynamics$recN
result_single_comp$dynamics$SPR

result_single_comp$HCR$decisionData$observed_catch_area_1
result_single_comp$HCR$decisionData$observed_catch_area_2

result_single_comp$HCR$decisionAnnual$TAC
result_single_comp$HCR$decisionAnnual

result_single_comp$dynamics$Ftotal
result_single_comp$HCR$decisionLocal


nr_diag <- nr_diagnostics_storage
#write.csv(nr_diag, "nr_diag.csv")

#check for failures
failures <- nr_diag[!nr_diag$converged, ]
if(nrow(failures) > 0) {
  cat("WARNING: Found", nrow(failures), "convergence failures:\n")
  print(failures)
} else {
  cat("All solver calls converged successfully!\n")
}
View(nr_diag)


# obs model outputs (data frame now)
str(result_single_comp$HCR$decisionData)
result_single_comp$HCR$decisionData$IDX_CPUE_1


#plots

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


#===========================================================================================#
#========================== Exploring NR outputs and performance ===========================#
#===========================================================================================#

#check SB
#Area 1 (TAC reduced by 20%): SB should increase
plot_SB(result_single_comp,areas=1)
#Area 2 (TAC increased by 60%): SB should decline
plot_SB(result_single_comp,areas=2)

#check F
plot_Ftotal(result_single_comp)

#cacth
plot_catchB_multi(result_single_comp,areas=1)
plot_catchB_multi(result_single_comp,areas=2)

#TAC
plot_TAC_by_area(result_single_comp, areas = 1)
plot_TAC_by_area(result_single_comp, areas = 2)
plot_TAC(result_single_comp, areas = "all")  # All areas in one plot

# Total catch across all areas (original behavior)
plot_catch_observations_both(result_single_comp)

# Area 1 only
plot_catch_observations_both(result_single_comp, areas = 1)
# Area 2 only
plot_catch_observations_both(result_single_comp, areas = 2)

# Both areas with faceting
plot_catch_observations_both(result_single_comp, areas = c(1, 2))

#combined catch obs:
plot_catch_observations_both(result_single_comp,show_individual = TRUE)


#calculate relative error
nr_diag$relative_error <- abs(nr_diag$predicted_catch - nr_diag$target_catch) / nr_diag$target_catch
ggplot(nr_diag, aes(x = target_catch, y = predicted_catch, color = factor(area))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Target vs Predicted Catch - should be on 1:1 line")


#check for extreme values
extreme_F <- nr_diag[nr_diag$final_F > 2 | nr_diag$final_F < 0.001, ]
if(nrow(extreme_F) > 0) {
  cat("WARNING: Extreme F values detected\n")
  print(extreme_F)
}

#realized catch vs target TAC
catchB <- result_single_comp$dynamics$catchB
TAC_decisions <- result_single_comp$HCR$decisionAnnual


#compare for each year/iteration/area
proj_start <- 12
for(yr in proj_start:(proj_start+2)) {
  for(iter in 1:3) {
    for(area in 1:2) {
      TAC <- TAC_decisions$TAC[TAC_decisions$year == yr &
                                 TAC_decisions$iteration == iter &
                                 TAC_decisions$area == area]
      realized <- catchB[yr, iter, area]

      error <- abs(realized - TAC) / TAC * 100

      cat(sprintf("Yr %d, Iter %d, Area %d: TAC=%.2f, Realized=%.2f, Error=%.1f%%\n",
                  yr-1, iter, area, TAC, realized, error))

      if(error > 5) {
        cat("  WARNING: Error >5%\n")
      }
    }
  }
}



#=============================================================================#
#===================      EFFORT BASED STRATEGy ==============================#
#=============================================================================#
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
ta@historicalEffort <- matrix(c(1.5, 3.4, 3.3, 3.8, 3.1, 3.0, 2.9, 1.8, 1.7, 1.6,
                                1.2, 2.4, 2.3, 2.9, 2.1, 2.0, 1.9, 0.8, 0.7, 0.6),
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
proj_fishery_area2@vulParams<-c(8.5, 0.1)
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
# MANAGEMENT STRATEGIES - SINGLE FLEET (SIMPLE EFFORT BASED)
# ============================================================================

# Single fleet
singleCompMP_Effort <- function(phase, dataObject) {
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

  # Simple effort decision rule

   if(phase == 2) {
     effort_multipliers <- c(0.8, 1.1) # define fixed effort multipliers for each area c(reduce, increase)

     # Return effort decisions
     return(list(
       year = rep(j, areas),
       iteration = rep(k, areas),
       area = 1:areas,
       effort_multiplier = effort_multipliers
     ))
   }


  if(phase == 3) {

    yrHist <- TimeAreaObj@historicalYears + 1

    #get effort decisions from Phase 2
    effort_decisions <- decisionAnnual[decisionAnnual$year == j &
                                         decisionAnnual$iteration == k, ]

    F_by_area <- numeric(areas) #that would be effort by area

    for(m in 1:areas) {
      effort_mult <- effort_decisions$effort_multiplier[effort_decisions$area == m]
      F_by_area[m] <- Ftotal[yrHist, k, m] * effort_mult
    }


    return(list(
      year = rep(j, areas),
      iteration = rep(k, areas),
      area = 1:areas,
      Flocal = F_by_area
    ))
  }
}

# very basic strategy effort example:
# - Phase 2 just returns fixed multipliers: [0.8, 1.1]
# - no extra calculations
# - each year, every iteration: same multipliers (for now)
# - Area 1 always gets 80% of historical effort 9decrease)
# - Area 2 always gets 110% of historical effort (increase)



# Create new TAC-based strategy
strategy_single_eff <- new("Strategy")
strategy_single_eff@title <- "Single Fleet Effort Strategy"
strategy_single_eff@projectionYears <- 5
strategy_single_eff@projectionName <- "singleCompMP_Effort"
strategy_single_eff@projectionParams <- list()  # No additional params needed

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

result_single_comp_effort <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  ProFisheryObj_list = proj_fishery_list,
  StrategyObj = strategy_single_eff, #updated
  StochasticObj = stochastic_obj,
  MultifleetObj = NULL,
  IndexObj = single_comprehensive_index,
  CatchObsObj = single_comprehensive_catch,
  LengthCompObj = single_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test1_single_with_effort",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "singleCompMP_Effort"
)

cat("Single fleet comprehensive simulation with obs models completed\n")

result_single_comp_effort <- readProjection(getwd(), "test1_single_with_effort")

# Population outputs

result_single_comp_effort$dynamics$SB

result_single_comp_effort$dynamics$VB

result_single_comp_effort$dynamics$Ftotal
result_single_comp_effort$dynamics$catchB


result_single_comp_effort$dynamics$recN

result_single_comp_effort$dynamics$SPR

result_single_comp_effort$HCR$decisionData$observed_catch_area_1
result_single_comp_effort$HCR$decisionData$observed_catch_area_2

result_single_comp_effort$HCR$decisionAnnual$TAC
result_single_comp_effort$HCR$decisionAnnual

result_single_comp_effort$dynamics$Ftotal

result_single_comp_effort$HCR$decisionLocal


# obs model outputs (data frame now)
str(result_single_comp_effort$HCR$decisionData)
result_single_comp_effort$HCR$decisionData$IDX_CPUE_1


#plots
# loading plot fucntion for toher plots

# #obs: need to add units
plot_SB(result_single_comp_effort)
plot_SB(result_single_comp_effort, areas=1)
plot_SB(result_single_comp_effort, areas=2)
plot_SB(result_single_comp_effort, areas=c(1,2))

plot_VB(result_single_comp_effort)
plot_VB(result_single_comp_effort, areas=1)
plot_VB(result_single_comp_effort, areas=2)

plot_Ftotal(result_single_comp_effort)
plot_Ftotal(result_single_comp_effort, areas=1)
plot_Ftotal(result_single_comp_effort, areas=2)


plot_catchB(result_single_comp_effort)
plot_catchB(result_single_comp_effort,areas=1)
plot_catchB(result_single_comp_effort,areas=2)

plot_catchN(result_single_comp_effort)
plot_catchN(result_single_comp_effort, areas=1)
plot_catchN(result_single_comp_effort, areas=2)

plot_discB(result_single_comp_effort)
plot_discB(result_single_comp_effort,areas=1)
plot_discB(result_single_comp_effort,areas=2)

plot_discN(result_single_comp_effort)
plot_discN(result_single_comp_effort,areas=1)
plot_discN(result_single_comp_effort,areas=2)

plot_catchB_multi(result_single_comp_effort)
plot_catchB_multi(result_single_comp_effort, areas=1)
plot_catchB_multi(result_single_comp_effort, areas=2)

plot_catchN_multi(result_single_comp_effort,show_individual = TRUE)
plot_catchN_multi(result_single_comp_effort, areas=1)
plot_catchN_multi(result_single_comp_effort, areas=2)

plot_SPR(result_single_comp_effort)
plot_recN(result_single_comp_effort)


# #plot obs models (indices)
plot_survey_indices(result_single_comp_effort)
plot_cpue_indices(result_single_comp_effort)

plot_all_indices(result_single_comp_effort)

plot_catch_observations_both(result_single_comp_effort,show_individual = FALSE)
plot_catch_observations_both(result_single_comp_effort,show_individual = TRUE)


# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_single_comp_effort,show_individual = TRUE)

plot_survey_length_comp(result_single_comp_effort,show_individual = TRUE)
plot_survey_length_comp(result_single_comp_effort, areas=1,show_individual = TRUE)
plot_survey_length_comp(result_single_comp_effort, areas=2,show_individual = TRUE)

