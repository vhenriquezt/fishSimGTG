#MP1: RFB (Reference Fishing Biommass) rule (ICES)


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
fleet1_sel_proj@vulParams <- c(8, 0.1)  # Same as base single fleet for 1-fleet test
fleet1_sel_proj@retType <- "full"
fleet1_sel_proj@retMax <- 1
fleet1_sel_proj@Dmort <- 0

fleet2_sel_proj <- new("Fishery")
fleet2_sel_proj@title <- "Fleet 2"
fleet2_sel_proj@vulType <- "logistic"
fleet2_sel_proj@vulParams <- c(12, 0.1)  # Different selectivity
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

# multifleet_rfb_MP: Multifleet strategy applying rfb rule


#global storage for Newton-Raphson diagnostics
nr_diagnostics_multifleet <<- data.frame()

multifleet_rfb_MP  <- function(phase, dataObject) {

  #Some configurations

  #survey index to use for biomass trend
  survey_index_name <- "IDX_Survey_1"   # it could be CPUE

  #fleet proportions: set to NULL for automatic (historical), or specify manually
  #example: c(0.65, 0.35) for 65% fleet 1, 35% fleet 2
  manual_fleet_proportions <- NULL

  #initial TAC: set to NULL to use last historical catch, or specify value
  initial_TAC_override <- NULL #initial_TAC_override <- 200  # uncomment to use fixed initial TAC

  #precautionary multiplier (maintain 95% probability above Blim)
  precautionary_m <- 0.95

  #number of years for historical catch averaging (for fleet proportions)
  n_years_for_proportions <- 3

  # Minimum capture length (Lc)
  # This should match fishery minimum size regulation
  Lc <- 8.4  # cm - used to calculate LF/M

  # LF_M calculation method:
  LF_M_method <- "beverton_holt"  # Options: "beverton_holt" -"traditional"



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
    } else {
      stop("rfb rule requires LengthCompObj for mean length calculation")
    }


    return(combined_data)
  }

  if(phase == 2) {

    # Purpose: calculate TAC using RFB formula and allocate to fleets/areas
    # Formula: TAC = TAC_prev × r × f × b × m
    # Where:
    #   r = biomass trend ratio (recent/historical index)
    #   f = length indicator (mean length / target length)
    #   b = biomass safeguard (protection when stock low)
    #   m = precautionary multiplier
    # Returns: Data frame with TAC for each fleet-area combination



    #calculate the last year of the historical period
    yrHist <- TimeAreaObj@historicalYears + 1  # the end of historical period, e.g., 11
    cat(sprintf("\n=== rfb RULE: Year %d, Iteration %d ===\n", j-1, k))

    #Step 0: get previous TAC or use the last historical catch


    # first projection year: use last historical catch or manual override
    # subsequent years: use previous year's TAC

    if(j == (yrHist + 1)) {
      #first projection year
      if(!is.null(initial_TAC_override)) {
        TAC_previous <- initial_TAC_override
        cat(sprintf("Using manual initial TAC: %.2f\n", TAC_previous))
      } else {
        #sum all fleet-area catches from last historical year
        TAC_previous <- 0
        for(m in 1:areas) {
          for(f in 1:nfleets) {
            catch_col <- paste0("fleet_", f, "_observed_catch_area_", m)
            last_hist_catch <- decisionData[[catch_col]][
              which(decisionData$k == k & decisionData$j == yrHist)
            ]
            if(length(last_hist_catch) > 0 && !is.na(last_hist_catch)) {
              TAC_previous <- TAC_previous + last_hist_catch
            }
          }
        }
        cat(sprintf("Initial TAC from last historical catch: %.2f\n", TAC_previous))
      }
    } else {
      #subsequent projection years: sum previous year's TACs
      prev_year_TACs <- decisionAnnual$TAC[decisionAnnual$year == (j-1) &
                                             decisionAnnual$iteration == k]
      TAC_previous <- sum(prev_year_TACs)
      cat(sprintf("TAC_previous from year %d: %.2f\n", j-2, TAC_previous))
    }


    #Step 1: calculate biomass index ratio (r)

    # Purpose: measure recent biomass trend from survey data
    # Index A = mean of last 2 survey values (recent)
    # Index B = mean of 3 preceding values (historical)
    # r = Index A / Index B
    # r > 1: increasing trend, r < 1: decreasing trend

    #extract index values for this iteration
    survey_data <- decisionData[decisionData$k == k, c("j", survey_index_name)]
    survey_data <- survey_data[!is.na(survey_data[[survey_index_name]]), ]

    #calculate Index A (mean of last 2 values)and  Index B (mean of 3 preceding values)
    Index_A <- mean(tail(survey_data[[survey_index_name]], 2), na.rm = TRUE)
    Index_B <- mean(tail(head(survey_data[[survey_index_name]], -2), 3), na.rm = TRUE)
    #biomass trend ratio
    r <- Index_A / Index_B
    cat(sprintf("Index A (last 2): %.2f | Index B (previous 3): %.2f | r = %.3f\n",
                Index_A, Index_B, r))

    #Step 2: calculate length indicator (f)

    # Purpose: measure if fish are being caught at optimal size
    # L_mean = current mean length in catch (from length composition)
    # LF_M = target length (length at F=M, MSY proxy)
    # f = L_mean / LF_M
    # f > 1: fishing larger fish, f < 1: fishing smaller fish

    #get previous year's length composition (most recent available)
    current_lc <- decisionData[decisionData$k == k & decisionData$j == (j-1), ]  # Changed from j to j-1
    bin_cols <- grep("_count_bin_", names(current_lc), value = TRUE)

      #get length bin info from LengthCompObj
    length_bin_width <- LengthCompObj@length_bin_width
    max_length <- max(unlist(lh$L))
    length_bins <- seq(0, max_length + length_bin_width, by = length_bin_width)
    bin_midpoints <- length_bins[-length(length_bins)] + length_bin_width/2

      #aggregate counts across all FD programs
    total_counts <- numeric(length(bin_midpoints))
    for(col in bin_cols) {
      if(!is.na(current_lc[[col]])) {
        bin_idx <- as.numeric(sub(".*_bin_", "", col))
        if(bin_idx <= length(total_counts)) {
          total_counts[bin_idx] <- total_counts[bin_idx] + current_lc[[col]]
        }
      }
    }

      #calculate mean length from length frecuency obs
    if(sum(total_counts) > 0) {
      L_mean <- sum(bin_midpoints * total_counts) / sum(total_counts)
    } else {
      stop("No length composition samples available for year ", j-1)
    }

      # MSY proxy length: LF=M - calculate target length (LF_M) using selected method
    if(LF_M_method == "traditional") {
      # LF_M = Linf × (1 - M/K)
      # the asymptotic length at F=M (Beverton & Holt 1957)
      LF_M <- lh$LifeHistory@Linf * (1 - lh$LifeHistory@M / lh$LifeHistory@K)
      cat(sprintf("Using LF_M = Linf × (1 - M/K)\n"))
    } else {
      # Beverton-Holt 1993 / Kell et al 2022: LF_M = 0.75×Lc + 0.25×Linf
      # accounts for minimum size regulations (Lc)
      # realistic target when fishery has size limits
      LF_M <- 0.75 * Lc + 0.25 * lh$LifeHistory@Linf
      cat(sprintf("Using Beverton-Holt (1993) LF_M = 0.75×Lc + 0.25×Linf\n"))
      cat(sprintf("  Lc (min capture length) = %.2f cm\n", Lc))
    }
      f_length <- L_mean / LF_M #(instead of LF_M could I use L50 mat? )

      cat(sprintf("L_mean: %.2f cm | LF_M: %.2f cm | f = %.3f\n",
                  L_mean, LF_M, f_length))

    #Step 3: calculate biomass Safeguard (b)
      # Purpose: Reduce TAC when stock approaches limit reference point
      # I_trigger = I_loss × 1.4 (where I_loss is lowest observed index)
      # b = min(I_last / I_trigger, 1.0)
      # When I_last < I_trigger: b < 1, TAC is reduced
      # When I_last ≥ I_trigger: b = 1, no reduction


      I_last <- tail(survey_data[[survey_index_name]], 1)
      I_loss <- min(survey_data[[survey_index_name]], na.rm = TRUE)
      I_trigger <- I_loss * 1.4
      b <- min(I_last / I_trigger, 1.0)

      cat(sprintf("I_last: %.2f | I_trigger: %.2f (I_loss × 1.4) | b = %.3f\n",
                  I_last, I_trigger, b))


    #Step 4: calculate precautionary multiplier (m)
      # Purpose: additional precautionary reduction (e.g., 0.95 = 5% reduction)
      # maintains higher probability of staying above limit reference points

      cat(sprintf("Precautionary multiplier m = %.2f\n", precautionary_m))

    #Step 5: Calculate Preliminary TAC

      # TAC_preliminary = TAC_previous × r × f × b × m
      # each component adjusts TAC based on stock status:
      # - r: biomass trend (up/down)
      # - f: size structure (optimal/suboptimal)
      # - b: safeguard (triggered/not triggered)
      # - m: precautionary reduction

      TAC_preliminary <- TAC_previous * r * f_length * b * precautionary_m

      cat(sprintf("TAC_preliminary = %.2f × %.3f × %.3f × %.3f × %.2f = %.2f\n",
                  TAC_previous, r, f_length, b, precautionary_m, TAC_preliminary))

    #Step 6: stability clause

      # Purpose: limit year-to-year TAC changes
      # When b ≥ 1 (stock above trigger):
      #   - Maximum increase: +20% of previous TAC
      #   - Maximum decrease: -30% of previous TAC
      # When b < 1 (stock below trigger):
      #   - No stability clause, allow full reduction from safeguard


      if(b >= 1.0) {
        TAC_max <- TAC_previous * 1.2   # +20% cap
        TAC_min <- TAC_previous * 0.7   # -30% floor
        TAC_final <- max(TAC_min, min(TAC_preliminary, TAC_max))
        stability_applied <- "YES (±20%/-30%)"
      } else {
        TAC_final <- TAC_preliminary
        stability_applied <- sprintf("NO (b=%.2f reduces TAC)", b)
      }

      cat(sprintf("TAC_final: %.2f (stability: %s)\n", TAC_final, stability_applied))



    #Step 7: Allocate TAC by fleet proportions

      # Purpose: allocate total TAC among fleets
      # Manual mode: use user-specified proportions
      # Automatic mode: calculate from recent historical catch patterns

      if(!is.null(manual_fleet_proportions)) {
        #user-specified proportions
        if(length(manual_fleet_proportions) != nfleets) {
          stop(sprintf("manual_fleet_proportions must have %d elements", nfleets))
        }
        if(abs(sum(manual_fleet_proportions) - 1.0) > 1e-6) {
          manual_fleet_proportions <- manual_fleet_proportions / sum(manual_fleet_proportions)
          warning("Fleet proportions normalized to sum to 1")
        }
        fleet_proportions <- manual_fleet_proportions
        cat(sprintf("\nUsing MANUAL fleet proportions: %s\n",
                    paste(round(fleet_proportions, 3), collapse = ", ")))
      } else {
        #calculate from historical catch
        start_year <- max(2, j - n_years_for_proportions)
        fleet_catches <- numeric(nfleets)

        for(f in 1:nfleets) {
          for(m in 1:areas) {
            catch_col <- paste0("fleet_", f, "_observed_catch_area_", m)
            recent_catches <- decisionData[[catch_col]][
              which(decisionData$k == k &
                      decisionData$j >= start_year &
                      decisionData$j < j)
            ]
            fleet_catches[f] <- fleet_catches[f] + sum(recent_catches, na.rm = TRUE)
          }
        }

        if(sum(fleet_catches) > 0) {
          fleet_proportions <- fleet_catches / sum(fleet_catches)
        } else {
          fleet_proportions <- rep(1/nfleets, nfleets)
        }

        cat(sprintf("\nUsing HISTORICAL fleet proportions (%d-yr avg): %s\n",
                    n_years_for_proportions,
                    paste(round(fleet_proportions, 3), collapse = ", ")))
      }

    #Step 8: Allocate TAC to fleet area combination

      # Purpose: distribute each fleet TAC among areas based on catch patterns
      # For each fleet:
      #   1. Calculate fleet's total TAC (total × fleet_proportion)
      #   2. Calculate area proportions within that fleet (from recent catch)
      #   3. Allocate fleet TAC to areas

      TAC_decisions <- data.frame()
      start_year <- max(2, j - n_years_for_proportions)

      for(f in 1:nfleets) {
        fleet_TAC_total <- TAC_final * fleet_proportions[f]

      #calculate area proportions for this fleet
        area_catches_for_fleet <- numeric(areas)
        for(m in 1:areas) {
          catch_col <- paste0("fleet_", f, "_observed_catch_area_", m)
          recent_catches <- decisionData[[catch_col]][
            which(decisionData$k == k &
                    decisionData$j >= start_year &
                    decisionData$j < j)
          ]
          area_catches_for_fleet[m] <- sum(recent_catches, na.rm = TRUE)
        }

      #calculate area proportions within this fleet
        if(sum(area_catches_for_fleet) > 0) {
          area_props <- area_catches_for_fleet / sum(area_catches_for_fleet)
        } else {
          area_props <- rep(1/areas, areas)
        }

        cat(sprintf("\n  Fleet %d (%.1f%% of total TAC):\n", f, fleet_proportions[f]*100))

        for(m in 1:areas) {
          TAC_fleet_area <- fleet_TAC_total * area_props[m]

          TAC_decisions <- rbind(TAC_decisions, data.frame(
            year = j,
            iteration = k,
            area = m,
            fleet = f,
            TAC = TAC_fleet_area,
            fleet_proportion = fleet_proportions[f],
            area_proportion_in_fleet = area_props[m],
            total_TAC = TAC_final,
            # RFB components for tracking
            r = r,
            f_length = f_length,
            b = b,
            m = precautionary_m,
            L_mean = L_mean,
            LF_M = LF_M,
            LF_M_method = LF_M_method,
            Lc_used = if(LF_M_method == "beverton_holt") Lc else NA,
            I_last = I_last,
            I_trigger = I_trigger,
            TAC_previous = TAC_previous,
            TAC_preliminary = TAC_preliminary,
            stability_applied = stability_applied,
            stringsAsFactors = FALSE
          ))

          cat(sprintf("    Area %d: TAC = %.2f (%.1f%% of fleet TAC)\n",
                      m, TAC_fleet_area, area_props[m]*100))
        }
      }

    cat("=====================================\n\n")

    return(TAC_decisions)
  }


  if(phase == 3) {

    # Purpose:use Newton-Raphson solver to find F values that achieve TACs
    # Process each area independently (different N and selectivity by area)
    # For each area:
    #   1. Extract TACs for all fleets in that area
    #   2. Prepare area-specific abundance (N) and selectivity data
    #   3. Call Newton-Raphson solver to find F vector [F1, F2, ..., Fn]
    #   4. Store results for population dynamics
    # Returns: Data frame with F for each year-iteration-area-fleet combination

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

      # Call Newton-Raphson Solver

      # Inputs:
      #   - TAC_targets: Vector of TACs for each fleet in this area
      #   - N_temp: Abundance arrays (area m only, in position 1)
      #   - selGroup_temp: Selectivity objects (area m only, in position 1)
      # Output:
      #   - F_vector: Vector of F values [F_fleet1, F_fleet2, ..., F_fleetn]
      # steps:
      #   1. Initialize F guess from TAC/vulnerable_biomass
      #   2. Iterate using Newton-Raphson:
      #      - Calculate predicted catch with current F
      #      - Calculate derivative (dCatch/dF)
      #      - Update F: F_new = F_old - error/(0.8 × derivative)
      #   3. Converge when |predicted - target| < tolerance


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

      # STEP 1: #extract NR checks before fleet loop
      #extract NR checks (before fleet loop)
      #since both fleets are TAC-managed, attributes have length = nfleets

      # The solver returns one F_vector with attributes attached.
      # These attributes contain data for all TAC-managed fleets.
      # We extract them once for clarity
      # Inside loop - just use simple variable names


      converged_flag <- attr(F_vector, "converged")
      nr_iters <- attr(F_vector, "iterations")
      initial_guess_vec <- attr(F_vector, "initial_guess")
      final_error_vec <- attr(F_vector, "final_error")
      predicted_catch_vec <- attr(F_vector, "predicted_catch")
      target_catch_vec <- attr(F_vector, "target_catch")

      #Step 2: identify TAC-managed fleets (both in this case) (determine which fleets are managed by TAC (versus effort))
      tac_managed_fleets <- which(!is.na(fleet_TACs_this_area))

      #validate attribute lengths match expectations
      expected_attr_length <- length(tac_managed_fleets)

      if(!is.null(initial_guess_vec) && length(initial_guess_vec) != expected_attr_length) {
        warning(sprintf("Area %d: initial_guess length mismatch. Expected %d, got %d",
                        m, expected_attr_length, length(initial_guess_vec)))
      }



      # STEP3: mapping fleet to attribute index
      #Store results for each fleet in this area
      for(f in 1:nfleets) {
        # add to F_results data frame
        F_results <- rbind(F_results, data.frame(
          year = j,
          iteration = k,
          area = m,
          fleet = f,
          Flocal = F_vector[f]
        ))

        # Diagnostics: mappin fleet number to attribute index

        if(f %in% tac_managed_fleets) {

          #find position in attribute vectors
          #for TAC-managed fleets: f=1: attr_index=1, f=2: attr_index=2
          attr_index <- which(tac_managed_fleets == f)


        #store to track Newton-Raphson performance 9global varibale)
        nr_diagnostics_multifleet <<- rbind(nr_diagnostics_multifleet, data.frame(
          year = j,
          iteration = k,
          area = m,
          fleet = f,
          management_type = "TAC",
          TAC = fleet_TACs_this_area[f],
          initial_F_guess = if(is.null(initial_guess_vec)) NA_real_
          else initial_guess_vec[attr_index],
          final_F = F_vector[f],
          converged = if(is.null(converged_flag)) FALSE else converged_flag,
          nr_iterations = if(is.null(nr_iters)) NA_integer_ else as.integer(nr_iters),
          nr_error = if(is.null(final_error_vec)) NA_real_
          else final_error_vec[attr_index],
          predicted_catch = if(is.null(predicted_catch_vec)) NA_real_
          else predicted_catch_vec[attr_index],
          target_catch = fleet_TACs_this_area[f],
          stringsAsFactors = FALSE
        ))

        } else {
          # effort-managed fleet (should not happen in RFB, but included anyway)
          nr_diagnostics_multifleet <<- rbind(nr_diagnostics_multifleet, data.frame(
            year = j,
            iteration = k,
            area = m,
            fleet = f,
            management_type = "effort",
            TAC = NA_real_,
            initial_F_guess = NA_real_,
            final_F = F_vector[f],
            converged = TRUE,
            nr_iterations = 0L,
            nr_error = NA_real_,
            predicted_catch = NA_real_,
            target_catch = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    cat("=====================================\n\n")

    return(F_results)
  }
}




# Strategy objects
strategy_rfb  <- new("Strategy")
strategy_rfb @title <- "rfb Rule with Multifleet"
strategy_rfb @projectionYears <- 5
strategy_rfb @projectionName <- "multifleet_rfb_MP"
strategy_rfb@projectionParams <- list()


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


result_multi_RFB <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StrategyObj = strategy_rfb,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_2fleet,
  IndexObj = multi_comprehensive_index,
  CatchObsObj = multi_comprehensive_catch,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test_multi_RFB_TAC",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "multifleet_rfb_MP"
)
cat("Multifleet (2 fleeets) comprehensive simulation completed\n")

result_multi_RFB  <- readProjection(getwd(), "test_multi_RFB_TAC")

# Population outputs
result_multi_RFB$dynamics$multifleet$actual_catch_proportions

result_multi_RFB$dynamics$SB
result_multi_RFB$dynamics$recN
result_multi_RFB$dynamics$SPR

# obs model outputs (new structure data frame)
result_multi_RFB$HCR$decisionData
result_multi_RFB$HCR$decisionData$IDX_Survey_1
result_multi_RFB$HCR$decisionLocal

result_multi_RFB$HCR$decisionData$fleet_1_observed_catch_area_1
result_multi_RFB$HCR$decisionData$fleet_1_observed_catch_area_2
result_multi_RFB$HCR$decisionData$fleet_2_observed_catch_area_1
result_multi_RFB$HCR$decisionData$fleet_2_observed_catch_area_2


result_multi_RFB$HCR$decisionAnnual$TAC
result_multi_RFB$HCR$decisionAnnual

result_multi_RFB$HCR$decisionLocal
result_multi_RFB$HCR$decisionAnnual


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
all(nr_diagnostics_multifleet$converged)

#obs: need to add units
plot_SB(result_multi_RFB)
plot_SB(result_multi_RFB, areas=1)
plot_SB(result_multi_RFB, areas=2)
plot_SB(result_multi_RFB, areas=c(1,2))


plot_catchB(result_multi_RFB)
plot_catchB(result_multi_RFB,areas=1)
plot_catchB(result_multi_RFB,areas=2)

plot_catchN(result_multi_RFB)
plot_catchN(result_multi_RFB, areas=1)
plot_catchN(result_multi_RFB, areas=2)

plot_discB(result_multi_RFB)
plot_discB(result_multi_RFB,areas=1)
plot_discB(result_multi_RFB,areas=2)

plot_discN(result_multi_RFB)
plot_discN(result_multi_comp,areas=1)
plot_discN(result_multi_comp,areas=2)

plot_catchB_multi(result_multi_RFB)
plot_catchB_multi(result_multi_RFB, areas=1)
plot_catchB_multi(result_multi_RFB, areas=2)

plot_catchN_multi(result_multi_RFB,show_individual = TRUE)
plot_catchN_multi(result_multi_RFB, areas=1)
plot_catchN_multi(result_multi_RFB, areas=2)

plot_SPR(result_multi_RFB)
plot_recN(result_multi_RFB)


#plot obs models (indices)
plot_survey_indices(result_multi_RFB)
plot_cpue_indices(result_multi_RFB)

plot_all_indices(result_multi_RFB)

#plot individual indices
plot_indices(result_multi_RFB,
             index_pattern = "IDX_CPUE.*Fleet_1",
             show_individual = TRUE,
             title = "Fleet 1 CPUE Only")

plot_indices(result_multi_RFB,
             index_pattern = "IDX_CPUE.*Fleet_2",
             show_individual = TRUE,
             title = "Fleet 2 CPUE Only")



#not working need to review
# Total catch across all areas
plot_catch_observations_both(result_multi_RFB,show_individual=TRUE)

plot_catch_observations_both(result_multi_RFB, areas = 1)
plot_catch_observations_both(result_multi_RFB, areas = 2)
plot_catch_observations_both(result_multi_RFB, areas = c(1, 2))



# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_multi_RFB,show_individual = TRUE)
plot_fishery_length_comp(result_multi_RFB, areas=1,show_individual = TRUE)
plot_fishery_length_comp(result_multi_RFB, areas=2,show_individual = TRUE)

plot_survey_length_comp(result_multi_RFB,show_individual = TRUE)
plot_survey_length_comp(result_multi_RFB, areas=1,show_individual = TRUE)
plot_survey_length_comp(result_multi_RFB, areas=2,show_individual = TRUE)



# NEW: Custom filtering for fleets and areas
plot_length_composition_by_area(result_multi_RFB,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets



plot_length_composition_by_area(result_multi_RFB,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets


plot_length_composition_by_area(result_multi_RFB,
                                program_pattern = "LC_Survey",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets



#===========================================================================================#
#========================== Exploring NR outputs and performance ===========================#
#===========================================================================================#

#check SB
plot_SB(result_multi_RFB,areas=1)
plot_SB(result_multi_RFB,areas=2)

plot_Ftotal_multi(result_multi_RFB,areas=c(1,2))
plot_Ftotal_multi(result_multi_RFB,areas=c(1))
plot_Ftotal_multi(result_multi_RFB,areas=c(2))

#cacth (add fleet)
plot_catchB_multi(result_multi_RFB,areas=1)
plot_catchB_multi(result_multi_RFB,areas=2)

#TAC (add fleet)
plot_TAC_by_area(result_multi_RFB, areas = 1)
plot_TAC_by_area(result_multi_RFB, areas = 2)
plot_TAC_by_area(result_multi_RFB, areas = 1)


plot_TAC(result_multi_RFB, areas = "all",show_fleets = TRUE)  # All areas in one plot
plot_TAC_by_fleet(result_multi_RFB)
plot_TAC(result_multi_RFB, areas = c(1,2), show_fleets = TRUE)


result_multi_RFB$HCR$decisionAnnual


#add fleet problem with these plots (problem: unused argument (areas_to_plot)
# Total catch across all areas (original behavior)
plot_catch_observations_both(result_multi_RFB)

# Area 1 only
plot_catch_observations_both(result_multi_RFB, areas = 1)
# Area 2 only
plot_catch_observations_both(result_multi_RFB, areas = 2)

# Both areas with faceting
plot_catch_observations_both(result_multi_RFB, areas = c(1, 2))

#combined catch obs:
plot_catch_observations_both(result_multi_RFB,show_individual = TRUE)


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
catchB <- result_multi_RFB$dynamics$catchB
TAC_decisions <- result_multi_RFB$HCR$decisionAnnual


#compare for each year/iteration/area
proj_start <- 12
nfleets <- 2
for(yr in proj_start:(proj_start+2)) {
  for(iter in 1:3) {
    for(area in 1:2) {
      cat(sprintf("\nYr %d, Iter %d, Area %d:\n", yr-1, iter, area))

      for(fleet in 1:nfleets) {
        # Fleet-specific TAC
        TAC_fleet <- TAC_decisions$TAC[TAC_decisions$year == yr &
                                         TAC_decisions$iteration == iter &
                                         TAC_decisions$area == area &
                                         TAC_decisions$fleet == fleet]

        # Fleet-specific realized catch
        realized_fleet <- result_multi_RFB$dynamics$multifleet$catchB_by_fleet[yr, iter, area, fleet]

        error <- abs(realized_fleet - TAC_fleet) / TAC_fleet * 100

        cat(sprintf("  Fleet %d: TAC=%.2f, Realized=%.2f, Error=%.1f%%\n",
                    fleet, TAC_fleet, realized_fleet, error))

        if(error > 5) {
          cat("    WARNING: Error >5%\n")
        }
      }
    }
  }
}

#simple selectivity check - compare fleet selectivities directly
lh <- LHwrapper(result_multi_RFB$LifeHistoryObj, result_multi_RFB$TimeAreaObj)

#check multiple ages on the selectivity curve
test_ages <- c(1,2,3,4,5,6,7,8,9,10,11, 12, 13,14, 15)

for(age in test_ages) {
  cat(sprintf("\nAge %d:\n", age))

  fleet1_sel <- selWrapper(lh, result_multi_RFB$TimeAreaObj,
                           FisheryObj = result_multi_RFB$MultifleetObj@fleet_selectivity_proj_list[[1]][[1]],
                           doPlot = FALSE)

  fleet2_sel <- selWrapper(lh, result_multi_RFB$TimeAreaObj,
                           FisheryObj = result_multi_RFB$MultifleetObj@fleet_selectivity_proj_list[[1]][[2]],
                           doPlot = FALSE)

  f1_keep <- fleet1_sel$keep[[1]][age]
  f2_keep <- fleet2_sel$keep[[1]][age]

  cat(sprintf("  Fleet 1 (L50=8):  %.4f\n", f1_keep))
  cat(sprintf("  Fleet 2 (L50=12): %.4f\n", f2_keep))
  cat(sprintf("  Different? %s\n", abs(f1_keep - f2_keep) > 0.001))



}
