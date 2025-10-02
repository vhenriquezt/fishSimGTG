#testing Newton - multi- fleet
#validation


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

# Shared random seed
test_seed <- 12345

# ============================================================================
# MANAGEMENT STRATEGIES - MULTIFLEET FLEET
# ============================================================================

# Multifleet strategy
# In phase 2, I am creating a 3 level hierarchical approach to
# split a TAC among areas and fleets

#global storage for Newton-Raphson diagnostics
nr_diagnostics_multifleet <<- data.frame()

multiCompMP_TAC <- function(phase, dataObject) {
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

    #Step1: Calculate Total TAC value from average of recent catches
    total_TAC <- 0

    for(m in 1:areas) {
      for(f in 1:nfleets) {
        catch_col_name <- paste0("fleet_", f, "_observed_catch_area_", m)

      #extract recent observed catches for this area and fleet
      recent_catches <- decisionData[[catch_col_name]][
        which(decisionData$k == k &
                decisionData$j >= start_year &        #Years 9, 10, ...
                decisionData$j < j) #excludes current year (cannot use data because data have not been collected yet!)
      ]

      if(length(recent_catches) > 0 && !all(is.na(recent_catches))) {
        avg_catch <- mean(recent_catches, na.rm = TRUE)
        total_TAC <- total_TAC + avg_catch
      } else {
        fallback <- 0.10 * RB_by_fleet[j, k, m, f]
        total_TAC <- total_TAC + fallback
      }
      }
    }
    #calculate total TAC (sum of averages across all fleet-area combinations)

    total_TAC <- 0

    for(m in 1:areas) {
      for(f in 1:nfleets) {
        catch_col_name <- paste0("fleet_", f, "_observed_catch_area_", m)

        recent_catches <- decisionData[[catch_col_name]][
          which(decisionData$k == k &
                  decisionData$j >= start_year &
                  decisionData$j < j)
        ]

        if(length(recent_catches) > 0 && !all(is.na(recent_catches))) {
          avg_catch <- mean(recent_catches, na.rm = TRUE)
          total_TAC <- total_TAC + avg_catch
        } else {
          #fallback for this fleet-area: use 10% of retained biomass
          fallback <- 0.10 * RB_by_fleet[j, k, m, f]
          total_TAC <- total_TAC + fallback
        }
      }
    }

    #we  can apply a overall TAC multiplier if we need toS
    overall_TAC_multiplier <- 1.0  # 1= no change, or 0.9= 10% reduction
    total_TAC <- total_TAC * overall_TAC_multiplier

    cat(sprintf("\n=== TAC ALLOCATION: Yr %d, Iter %d ===\n", j-1, k))
    cat(sprintf("Total TAC (sum of 3-yr averages): %.2f\n", total_TAC))


    #Step2: Allocate TAC by area

    area_proportions <- c(0.4, 0.6)  # Area 1: 40%, Area 2: 60%

    if(abs(sum(area_proportions) - 1.0) > 1e-6) {
      stop("Area proportions must sum to 1.0")
    }

    TAC_by_area <- total_TAC * area_proportions

    cat("area Allocations:\n")
    for(m in 1:areas) {
      cat(sprintf("  area %d: %.2f (%.1f%%)\n",
                  m, TAC_by_area[m], area_proportions[m] * 100))
    }

    #Step3: Allocate TAC by fleet

    fleet_proportions_by_area <- matrix(c(
      0.7, 0.3,  # Area 1: Fleet 1 = 70%, Fleet 2 = 30%
      0.5, 0.5   # Area 2: Fleet 1 = 50%, Fleet 2 = 50%
    ), nrow = areas, byrow = TRUE)

    for(m in 1:areas) {
      if(abs(sum(fleet_proportions_by_area[m, ]) - 1.0) > 1e-6) {
        stop(sprintf("fleet proportions for area %d must sum to 1.0", m))
      }
    }

    TAC_decisions <- data.frame()

    for(m in 1:areas) {
      cat(sprintf("\nArea %d breakdown (Total: %.2f):\n", m, TAC_by_area[m]))

      for(f in 1:nfleets) {
        fleet_proportion <- fleet_proportions_by_area[m, f]
        TAC_fleet_area <- TAC_by_area[m] * fleet_proportion

        TAC_decisions <- rbind(TAC_decisions, data.frame(
          year = j,
          iteration = k,
          area = m,
          fleet = f,
          TAC = TAC_fleet_area,
          area_proportion = area_proportions[m],
          fleet_proportion = fleet_proportion,
          total_TAC = total_TAC
        ))

        cat(sprintf("  Fleet %d: %.2f (%.1f%% of area TAC)\n",
                    f, TAC_fleet_area, fleet_proportion * 100))
      }
    }

    cat("=====================================\n\n")

    return(TAC_decisions)
  }


  #Vania edit's to match Bill's edits
  if(phase == 3) {

    # Extract TAC decisions for this year/iteration
    TAC_decisions <- decisionAnnual[decisionAnnual$year == j &
                                      decisionAnnual$iteration == k, ]

    # Initialize result storage
    F_results <- data.frame()

    # Process each area
    for(m in 1:areas) {
      # Get TACs for all fleets in this area
      fleet_TACs <- TAC_decisions$TAC[TAC_decisions$area == m]

      # Prepare area-specific data structures for the solver
      # Critical: Create 1-area, multi-fleet structure
      N_temp <- lapply(1:lh$gtg, function(gtg_idx) {
        temp_array <- array(0, dim = c(dim(N[[gtg_idx]])[1],
                                       dim(N[[gtg_idx]])[2],
                                       1))  # Only 1 area
        temp_array[, , 1] <- N[[gtg_idx]][, , m]  # Copy area m data
        temp_array
      })

      # Selectivity: Extract all fleets for this area
      selGroup_temp <- list(selGroup[[m]])  # List with 1 area containing all fleet selectivities


      #when solving for area 1:
      # selGroup_temp <- list(selGroup[[1]])
      # selGroup_temp[[1]][[1]] = Fleet 1 selectivity for area 1
      # selGroup_temp[[1]][[2]] = Fleet 2 selectivity for area 1



      # Check selectivity structure
      cat(sprintf("\n=== SOLVING Area %d ===\n", m))
      cat("Fleet TACs:", paste(round(fleet_TACs, 2), collapse = ", "), "\n")

      cat("\n--- SELECTIVITY VERIFICATION ---\n")
      for(f in 1:nfleets) {
        cat(sprintf("Fleet %d selectivity check:\n", f))
        cat(sprintf("  Original selGroup[[%d]][[%d]]$keep[[1]][10] = %.4f\n",
                    m, f, selGroup[[m]][[f]]$keep[[1]][10]))
        cat(sprintf("  Temp selGroup_temp[[1]][[%d]]$keep[[1]][10] = %.4f\n",
                    f, selGroup_temp[[1]][[f]]$keep[[1]][10]))
        cat(sprintf("  Match? %s\n",
                    ifelse(selGroup[[m]][[f]]$keep[[1]][10] == selGroup_temp[[1]][[f]]$keep[[1]][10],
                           "YES", "NO")))
      }
      cat("--------------------------------\n\n")



      # Call Newton-Raphson solver for this area
      F_vector <- solveTAC_to_F_fishSimGTG(
        j = j,
        k = k,
        TAC_targets = fleet_TACs,  # Vector of TACs, one per fleet
        N = N_temp,
        lh = lh,
        selGroup = selGroup_temp,
        M_rate = lh$LifeHistory@M,
        effort_F_by_fleet = numeric(0),  # Not used in TAC mode
        is_multifleet = TRUE,
        areas = 1,  # Solver sees 1 area
        nfleets = nfleets,
        TAC_type = "keep",
        control = list(maxiterF = 300, tolF = 1e-4)
      )

      cat("Solved F values:", paste(round(F_vector, 6), collapse = ", "), "\n")

      # Store results
      for(f in 1:nfleets) {
        F_results <- rbind(F_results, data.frame(
          year = j,
          iteration = k,
          area = m,
          fleet = f,
          Flocal = F_vector[f]
        ))

        nr_diagnostics_multifleet <<- rbind(nr_diagnostics_multifleet, data.frame(
          year = j,
          iteration = k,
          area = m,
          fleet = f,
          TAC = fleet_TACs[f],
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
          target_catch = fleet_TACs[f],
          stringsAsFactors = FALSE
        ))
      }
    }

    return(F_results)
  }
}


# Strategy objects
strategy_multi_comp <- new("Strategy")
strategy_multi_comp@title <- "Multifleet Comprehensive"
strategy_multi_comp@projectionYears <- 5
strategy_multi_comp@projectionName <- "multiCompMP_TAC"
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
  fileName = "test_multiCompMP_TAC",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "multiCompMP_TAC"
)
cat("Multifleet (2 fleeets) comprehensive simulation completed\n")

result_multi_comp  <- readProjection(getwd(), "test_multiCompMP_TAC")

# Population outputs

result_multi_comp$dynamics$SB
result_multi_comp$dynamics$recN
result_multi_comp$dynamics$SPR

# obs model outputs (new structure data frame)
result_multi_comp$HCR$decisionData
result_multi_comp$HCR$decisionData$IDX_Survey_1
result_multi_comp$HCR$decisionLocal

result_multi_comp$HCR$decisionData$fleet_1_observed_catch_area_1
result_multi_comp$HCR$decisionData$fleet_1_observed_catch_area_2
result_multi_comp$HCR$decisionData$fleet_2_observed_catch_area_1
result_multi_comp$HCR$decisionData$fleet_2_observed_catch_area_2


result_multi_comp$HCR$decisionAnnual$TAC
result_multi_comp$HCR$decisionAnnual

result_multi_comp$HCR$decisionLocal



#checkeando convergencia de Newton Raphson

nr_diag <- nr_diagnostics_multifleet
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

plot_discB(result_multi_comp)
plot_discB(result_multi_comp,areas=1)
plot_discB(result_multi_comp,areas=2)

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



#not working need to review
# Total catch across all areas
plot_catch_observations_both(result_multi_comp)

# Area 1 only
plot_catch_observations_both(result_multi_comp, areas = 1)
# Area 2 only
plot_catch_observations_both(result_multi_comp, areas = 2)

# Both areas with faceting
plot_catch_observations_both(result_multi_comp, areas = c(1, 2))

#combined catch obs:
plot_catch_observations_both(result_multi_comp,show_individual = TRUE)
#stop()


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



#===========================================================================================#
#========================== Exploring NR outputs and performance ===========================#
#===========================================================================================#

#check SB
#Area 1 (TAC reduced by xxx): SB should xxx
plot_SB(result_multi_comp,areas=1)
#Area 2 (TAC increased by xxx): SB should xxx
plot_SB(result_multi_comp,areas=2)

#problem with F area 2
plot_Ftotal_multi(result_multi_comp,areas=c(1,2))
plot_Ftotal_multi(result_multi_comp,areas=c(1))
plot_Ftotal_multi(result_multi_comp,areas=c(2))

#cacth (add fleet)
plot_catchB_multi(result_multi_comp,areas=1)
plot_catchB_multi(result_multi_comp,areas=2)

#TAC (add fleet)
plot_TAC(result_multi_comp, areas = 1, show_fleets = TRUE)
plot_TAC_by_area(result_multi_comp, areas = 2)
plot_TAC(result_multi_comp, areas = "all")  # All areas in one plot


#add fleet problem with these plots (problem: unused argument (areas_to_plot)
# Total catch across all areas (original behavior)
plot_catch_observations_both(result_multi_comp)

# Area 1 only
plot_catch_observations_both(result_multi_comp, areas = 1)
# Area 2 only
plot_catch_observations_both(result_multi_comp, areas = 2)

# Both areas with faceting
plot_catch_observations_both(result_multi_comp, areas = c(1, 2))

#combined catch obs:
plot_catch_observations_both(result_multi_comp,show_individual = TRUE)


#calculate relative error (problem with area 2 -only)
nr_diag$relative_error <- abs(nr_diag$predicted_catch - nr_diag$target_catch) / nr_diag$target_catch
ggplot(nr_diag, aes(x = target_catch, y = predicted_catch, color = factor(area))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Target vs Predicted Catch - should be on 1:1 line")


#check for extreme values (extreme F values detected, convergencia failed)
extreme_F <- nr_diag[nr_diag$final_F > 2 | nr_diag$final_F < 0.001, ]
if(nrow(extreme_F) > 0) {
  cat("WARNING: Extreme F values detected\n")
  print(extreme_F)
}

#realized catch vs target TAC
catchB <- result_multi_comp$dynamics$catchB
TAC_decisions <- result_multi_comp$HCR$decisionAnnual


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

