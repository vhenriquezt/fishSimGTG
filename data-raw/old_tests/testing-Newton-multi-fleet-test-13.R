#Testing Newton-Rhapson with 3 fleets

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
  # Area 1: increase then managed reduction
  c(1.0, 1.2, 1.4, 1.6, 1.8, 1.7, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0),
  # Area 2: slow development
  c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.2, 1.1, 1.0, 0.9, 0.9, 0.8)
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
fleet1_sel_hist@vulParams <- c(11.0, 0.08)
fleet1_sel_hist@retType <- "full"
fleet1_sel_hist@retMax <- 1
fleet1_sel_hist@Dmort <- 0

fleet1_sel_proj <- new("Fishery")
fleet1_sel_proj@title <- "Fleet 1 - Large Fish"
fleet1_sel_proj@vulType <- "logistic"
fleet1_sel_proj@vulParams <- c(10.0, 0.08)
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
fleet2_sel_proj@vulParams <- c(8.5, 0.12)
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
fleet3_sel_proj@vulParams <- c(7.0, 0.15)
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
# MANAGEMENT STRATEGIES - MULTIFLEET FLEET (3 fleets)
# ============================================================================

#global storage for Newton-Raphson diagnostics
nr_diagnostics_complex  <<- data.frame()

complexMP_TAC  <- function(phase, dataObject) {
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

  if(phase == 2) {
    #calculate the last year of the historical period
    yrHist <- TimeAreaObj@historicalYears + 1  # the end of historical period, e.g., 11
    n_years_avg <- 3
    # e.g., if j=12, start_year=9
    start_year <- max(2, j - n_years_avg) #ensures we never go before year 2 (prevent to going year 1 or earlier)

    cat(sprintf("\n=== TAC ALLOCATION: Year %d, Iter %d ===\n", j-1, k))

    #Step1: Calculate recent catch by fleet-area combination
    recent_catch_matrix <- matrix(0, nrow = areas, ncol = nfleets)

    for(m in 1:areas) {
      for(f in 1:nfleets) {
        catch_col_name <- paste0("fleet_", f, "_observed_catch_area_", m)

        #extract recent observed catches (from decision data) for this area and fleet
        recent_catches <- decisionData[[catch_col_name]][
          which(decisionData$k == k &
                  decisionData$j >= start_year &        #Years 9, 10, ...
                  decisionData$j < j) #excludes current year (cannot use data because data have not been collected yet!)
        ]
        #calculate mean catch for this fleet-area combo
        if(length(recent_catches) > 0 && !all(is.na(recent_catches))) {
          recent_catch_matrix[m, f] <- mean(recent_catches, na.rm = TRUE)
        } else {
          #fallback: use proportion of current retained biomass if no catch history
          recent_catch_matrix[m, f] <- 0.10 * RB_by_fleet[j, k, m, f]
        }
      }
    }

    cat("Recent 3-year average catches by fleet-area:\n")
    print(recent_catch_matrix)
    #example:
    #      Fleet1  Fleet2 fleet 3
    # Area1  45.2   30.1
    # Area2   8.5    5.7


    #Step2: Calculate total TAC based on sum of recent catches

    total_historical_catch <- sum(recent_catch_matrix)

    #apply overall management adjustment (e.g., 10% reduction for conservation)
    overall_TAC_multiplier <- 0.9
    total_TAC <- total_historical_catch * overall_TAC_multiplier

    cat(sprintf("\nTotal historical catch: %.2f\n", total_historical_catch))
    cat(sprintf("TAC multiplier: %.2f\n", overall_TAC_multiplier))
    cat(sprintf("Total TAC: %.2f\n", total_TAC))



    #Step3: allocate TAC by area based on historical proportions
    #this precserve the historical catch  distribution

    #two-stage allocation: Area first (spatial), then fleet within area (user groups)

    area_proportions <- rowSums(recent_catch_matrix) / total_historical_catch
    TAC_by_area <- total_TAC * area_proportions

    cat("\nTAC allocation by area:\n")
    for(m in 1:areas) {
      cat(sprintf("  Area %d: %.2f (%.1f%% of total)\n",
                  m, TAC_by_area[m], area_proportions[m] * 100))
    }

    #   example:
    #   Area 1: 67.8 (75.3% of total)
    #   Area 2: 12.8 (14.2% of total)

    # STEP 4: allocate area TAC among fleets based on historical proportions
    TAC_decisions <- data.frame() #initialize decision storage

    for(m in 1:areas) {
      #calculate fleet proportions within this area
      area_total_catch <- sum(recent_catch_matrix[m, ])

      if(area_total_catch > 0) {
        #use historical catch proportions
        fleet_props_in_area <- recent_catch_matrix[m, ] / area_total_catch
      } else {
        #fallback to biomass proportions if no catch history
        fleet_biomass <- sapply(1:nfleets, function(f) RB_by_fleet[j, k, m, f])
        fleet_props_in_area <- fleet_biomass / sum(fleet_biomass)
      }

      cat(sprintf("\nArea %d fleet allocation:\n", m))

      #allocate TAC to each fleet in this area
      for(f in 1:nfleets) {
        TAC_fleet_area <- TAC_by_area[m] * fleet_props_in_area[f]

        TAC_decisions <- rbind(TAC_decisions, data.frame(
          year = j,
          iteration = k,
          area = m,
          fleet = f,
          TAC = TAC_fleet_area,
          area_proportion = area_proportions[m],
          fleet_proportion_in_area = fleet_props_in_area[f],
          total_TAC = total_TAC,
          stringsAsFactors = FALSE
        ))

        cat(sprintf("  Fleet %d: %.2f (%.1f%% of area TAC)\n",
                    f, TAC_fleet_area, fleet_props_in_area[f] * 100))
      }
    }

    cat("=====================================\n\n")

    return(TAC_decisions)
  }


  if(phase == 3) {

    #extract TAC decisions from Phase 2
    TAC_decisions <- decisionAnnual[decisionAnnual$year == j &
                                      decisionAnnual$iteration == k, ]

    #initialize result storage
    F_results <- data.frame()

    #process each area independently because Newton-Raphson needs area-specific N and selectivity data
    for(m in 1:areas) {
      #get TACs for all fleets in this area
      fleet_TACs_this_area <- TAC_decisions$TAC[TAC_decisions$area == m]
      # Example for Area 1: [40.7, 27.1] (Fleet 1, Fleet 2)

      cat(sprintf("\n=== SOLVING Area %d ===\n", m))
      cat("Fleet TACs:", paste(round(fleet_TACs_this_area, 2), collapse = ", "), "\n")

      #prepare single-area, multi-fleet data structure
      #create temporary single-area N structure
      #solver expects N[[gtg]][age, year, area] where area dimension = 1

      N_temp <- lapply(1:lh$gtg, function(gtg_idx) {
        temp_array <- array(0, dim = c(dim(N[[gtg_idx]])[1], ## ages
                                       dim(N[[gtg_idx]])[2], # years
                                       1))  # 1 area only    #one area
        #copy data from area m into position 1
        temp_array[, , 1] <- N[[gtg_idx]][, , m]
        temp_array
      })

      #extract all fleet selectivities for this area
      #selGroup[[m]] contains all fleets for area m
      #wrap in list because solver expects [[area]][[fleet]] structure
      selGroup_temp <- list(selGroup[[m]])  # wrap in list for 1-area structure

      #verify selectivity structure
      cat("\n--- SELECTIVITY VERIFICATION ---\n")
      for(f in 1:nfleets) {
        cat(sprintf("Fleet %d GTG 1 Age 10 keep selectivity: %.4f\n",
                    f, selGroup_temp[[1]][[f]]$keep[[1]][10]))
      }
      cat("--------------------------------\n\n")

      # Call Newton-Raphson solver for this area
      F_vector <- solveTAC_to_F_fishSimGTG(
        j = j,
        k = k,
        TAC_targets = fleet_TACs_this_area,  # Vector: [fleet1_TAC, fleet2_TAC]
        N = N_temp,                          # abundance (area m only, in position 1)
        lh = lh,
        selGroup = selGroup_temp,            #selectivity (area m only, in position 1)
        M_rate = lh$LifeHistory@M,
        effort_F_by_fleet = numeric(0),  # Not used in TAC mode
        is_multifleet = TRUE,
        areas = 1,  # Solver sees 1 area (the current m area)
        nfleets = nfleets,
        TAC_type = "keep",
        control = list(maxiterF = 300, tolF = 1e-4)
      )

      cat("Solved F values:", paste(round(F_vector, 6), collapse = ", "), "\n")

      # Store results for each fleet in this area
      for(f in 1:nfleets) {
        # add to F_results data frame
        F_results <- rbind(F_results, data.frame(
          year = j,
          iteration = k,
          area = m,
          fleet = f,
          Flocal = F_vector[f]
        ))
        #store to track Newton-Raphson performance 9global varibale)
        nr_diagnostics_complex  <<- rbind(nr_diagnostics_complex , data.frame(
          year = j,
          iteration = k,
          area = m,
          fleet = f,
          TAC = fleet_TACs_this_area[f],
          initial_F_guess = ifelse(is.null(attr(F_vector, "initial_guess")),
                                   NA_real_, as.numeric(attr(F_vector, "initial_guess"))[f]),
          final_F = F_vector[f],
          converged = ifelse(is.null(attr(F_vector, "converged")),
                             FALSE, attr(F_vector, "converged")),
          nr_iterations = ifelse(is.null(attr(F_vector, "iterations")),
                                 NA_integer_, as.integer(attr(F_vector, "iterations"))),
          nr_error = ifelse(is.null(attr(F_vector, "final_error")),
                            NA_real_, as.numeric(attr(F_vector, "final_error"))[f]),
          predicted_catch = ifelse(is.null(attr(F_vector, "predicted_catch")),
                                   NA_real_, as.numeric(attr(F_vector, "predicted_catch"))[f]),
          target_catch = fleet_TACs_this_area[f],
          stringsAsFactors = FALSE
        ))
      }
    }

    return(F_results)
  }
}


strategy_complex <- new("Strategy")
strategy_complex@title <- "Complex Spatial-Temporal Strategy"
strategy_complex@projectionYears <- 7
strategy_complex@projectionName <- "complexMP_TAC"
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
    areas = c(1,2),
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
  fileName = "test3_complex_spatial_temporal_NR",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "complexMP"
)

cat("Complex simulation completed successfully\n")


# ============================================================================
# VALIDATION
# ============================================================================


result_complex <- readProjection(getwd(), "test3_complex_spatial_temporal_NR")


# Population outputs
result_complex$dynamics$multifleet$actual_catch_proportions

result_complex$dynamics$SB
result_complex$dynamics$recN
result_complex$dynamics$SPR

# obs model outputs (new structure data frame)
result_complex$HCR$decisionData
result_complex$HCR$decisionData$IDX_Survey_1
result_complex$HCR$decisionLocal

result_complex$HCR$decisionData$fleet_1_observed_catch_area_1
result_complex$HCR$decisionData$fleet_1_observed_catch_area_2
result_complex$HCR$decisionData$fleet_2_observed_catch_area_1
result_complex$HCR$decisionData$fleet_2_observed_catch_area_2


result_complex$HCR$decisionAnnual$TAC
result_complex$HCR$decisionAnnual

result_complex$HCR$decisionLocal
result_complex$HCR$decisionAnnual


#checkeando convergencia de Newton Raphson

nr_diag <- nr_diagnostics_complex
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
all(nr_diagnostics_complex$converged)

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




# Total catch across all areas
plot_catch_observations_both(result_complex,show_individual=TRUE)

plot_catch_observations_both(result_complex, areas = 1)
plot_catch_observations_both(result_complex, areas = 2)
plot_catch_observations_both(result_complex, areas = c(1, 2))



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



#===========================================================================================#
#========================== Exploring NR outputs and performance ===========================#
#===========================================================================================#

#check SB
plot_SB(result_complex,areas=1)
plot_SB(result_complex,areas=2)

plot_Ftotal_multi(result_complex,areas=c(1,2))
plot_Ftotal_multi(result_complex,areas=c(1))
plot_Ftotal_multi(result_complex,areas=c(2))

#cacth (add fleet)
plot_catchB_multi(result_complex,areas=1)
plot_catchB_multi(result_complex,areas=2)

#TAC (add fleet)
plot_TAC_by_area(result_complex, areas = 1)
plot_TAC_by_area(result_complex, areas = 2)
plot_TAC_by_area(result_complex, areas = 1)


plot_TAC(result_complex, areas = "all",show_fleets = TRUE)  # All areas in one plot
plot_TAC_by_fleet(result_complex)
plot_TAC(result_complex, areas = c(1,2), show_fleets = TRUE)


result_complex$HCR$decisionAnnual


#add fleet problem with these plots (problem: unused argument (areas_to_plot)
# Total catch across all areas (original behavior)
plot_catch_observations_both(result_complex)

# Area 1 only
plot_catch_observations_both(result_complex, areas = 1)
# Area 2 only
plot_catch_observations_both(result_complex, areas = 2)

# Both areas with faceting
plot_catch_observations_both(result_complex, areas = c(1, 2))

#combined catch obs:
plot_catch_observations_both(result_complex,show_individual = TRUE)


#calculate relative error (problem with area 2 -only)
nr_diag$relative_error <- abs(nr_diag$predicted_catch - nr_diag$target_catch) / nr_diag$target_catch
ggplot(nr_diag, aes(x = target_catch, y = predicted_catch, color = factor(area))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Target vs Predicted Catch - should be on 1:1 line")


#check for extreme values (extreme F values detected, convergencia failed)
extreme_F <- nr_diag[nr_diag$final_F > 2 | nr_diag$final_F < 0.0001, ]
if(nrow(extreme_F) > 0) {
  cat("WARNING: Extreme F values detected\n")
  print(extreme_F)
}

#realized catch vs target TAC
catchB <- result_complex$dynamics$catchB
TAC_decisions <- result_complex$HCR$decisionAnnual



#simple selectivity check - compare fleet selectivities directly
lh <- LHwrapper(result_complex$LifeHistoryObj, result_complex$TimeAreaObj)

#check multiple ages on the selectivity curve
test_ages <- c(1,2,3,4,5,6,7,8,9,10,11, 12,13,14,15,20,25,30,35,40)
ageClasses <- length(lh$L[[1]])

for(age in test_ages) {
  if(age <= ageClasses) {
    fleet1_sel <- selWrapper(lh, result_complex$TimeAreaObj,
                             FisheryObj = result_complex$MultifleetObj@fleet_selectivity_proj_list[[1]][[1]],
                             doPlot = FALSE)
    fleet2_sel <- selWrapper(lh, result_complex$TimeAreaObj,
                             FisheryObj = result_complex$MultifleetObj@fleet_selectivity_proj_list[[1]][[2]],
                             doPlot = FALSE)
    fleet3_sel <- selWrapper(lh, result_complex$TimeAreaObj,
                             FisheryObj = result_complex$MultifleetObj@fleet_selectivity_proj_list[[1]][[3]],
                             doPlot = FALSE)

    f1_keep <- fleet1_sel$keep[[1]][age]
    f2_keep <- fleet2_sel$keep[[1]][age]
    f3_keep <- fleet3_sel$keep[[1]][age]
    length_at_age <- lh$L[[1]][age]

    cat(sprintf("\nAge %d (Length %.1f cm):\n", age, length_at_age))
    cat(sprintf("  Fleet 1 (L50=9.0):  %.4f\n", f1_keep))
    cat(sprintf("  Fleet 2 (L50=8.8):  %.4f\n", f2_keep))
    cat(sprintf("  Fleet 3 (L50=7.8):  %.4f\n", f3_keep))
  }
}


#ploting sel for testing
ages <- 1:ageClasses
lengths <- lh$L[[1]]  # GTG 1 lengths

sel_data <- data.frame()
for(f in 1:3) {
  fleet_sel <- selWrapper(lh, result_complex$TimeAreaObj,
                          FisheryObj = result_complex$MultifleetObj@fleet_selectivity_proj_list[[1]][[f]],
                          doPlot = FALSE)

  sel_data <- rbind(sel_data, data.frame(
    Age = ages,
    Length = lengths,
    Selectivity = fleet_sel$keep[[1]],
    Fleet = paste0("Fleet ", f)
  ))
}


p <- ggplot(sel_data, aes(x = Length, y = Selectivity, color = Fleet)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = c(7.8, 8.8, 9.0), linetype = "dashed", alpha = 0.5) +
  labs(title = "Fleet Selectivity Comparison (Projection Period)",
       subtitle = paste("Linf =", lh$LifeHistory@Linf, "cm"),
       x = "Length (cm)",
       y = "Retention Probability") +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("Fleet 1" = "blue",
                                "Fleet 2" = "red",
                                "Fleet 3" = "green")) +
  annotate("text", x = 9.0, y = 0.55, label = "F1 L50", size = 3, angle = 90) +
  annotate("text", x = 8.8, y = 0.55, label = "F2 L50", size = 3, angle = 90) +
  annotate("text", x = 7.8, y = 0.55, label = "F3 L50", size = 3, angle = 90)

print(p)

#TAC vs catch comparison

proj_year <- 14
TAC_decisions <- result_complex$HCR$decisionAnnual

for(area in 1:2) {
  cat(sprintf("\nArea %d:\n", area))
  for(fleet in 1:3) {
    TAC_val <- TAC_decisions$TAC[TAC_decisions$year == proj_year &
                                   TAC_decisions$iteration == 1 &
                                   TAC_decisions$area == area &
                                   TAC_decisions$fleet == fleet]

    realized <- result_complex$dynamics$multifleet$catchB_by_fleet[proj_year, 1, area, fleet]
    error_pct <- abs(realized - TAC_val) / TAC_val * 100

    cat(sprintf("  Fleet %d: TAC=%.1f, Realized=%.1f, Error=%.1f%%\n",
                fleet, TAC_val, realized, error_pct))
  }
}
