#MP4: Index ratio (1 or multiple indices) + length indicators (none, one or multiple length indicators)
#Modified

# Here I am using 2 index ratio indicators (Survey + CPUE) and multiple length indicators

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
fleet1_sel_hist@vulParams <- c(11.2, 0.1)
fleet1_sel_hist@retType <- "full"
fleet1_sel_hist@retMax <- 1
fleet1_sel_hist@Dmort <- 0

fleet2_sel_hist <- new("Fishery")
fleet2_sel_hist@title <- "Fleet 2"
fleet2_sel_hist@vulType <- "logistic"
fleet2_sel_hist@vulParams <- c(8.2, 0.1)
fleet2_sel_hist@retType <- "full"
fleet2_sel_hist@retMax <- 1
fleet2_sel_hist@Dmort <- 0


# Fleet selectivities for multifleet tests (projection)
fleet1_sel_proj <- new("Fishery")
fleet1_sel_proj@title <- "Fleet 1"
fleet1_sel_proj@vulType <- "logistic"
fleet1_sel_proj@vulParams <- c(8, 0.1)
fleet1_sel_proj@retType <- "full"
fleet1_sel_proj@retMax <- 1
fleet1_sel_proj@Dmort <- 0

fleet2_sel_proj <- new("Fishery")
fleet2_sel_proj@title <- "Fleet 2"
fleet2_sel_proj@vulType <- "logistic"
fleet2_sel_proj@vulParams <- c(12, 0.1)
fleet2_sel_proj@retType <- "full"
fleet2_sel_proj@retMax <- 1
fleet2_sel_proj@Dmort <- 0


# Shared random seed
test_seed <- 12345

# ============================================================================
# MANAGEMENT STRATEGIES - MULTIFLEET FLEET
# ============================================================================

# multifleet MP: Index ratio

#global storage for Newton-Raphson diagnostics
nr_diagnostics_multifleet <<- data.frame()

multifleet_indexratio_length_MP_V2   <- function(phase, dataObject) {

  #New index configuration (as a list)

  index_config <- list(
    #   "IDX_Survey_1"                           # single FI survey
    #   c("IDX_Survey_1", "IDX_CPUE_1_Fleet_1")  # multiple indices
    #1. survey index to use for biomass trend
    survey_index_names = c("IDX_Survey_1", "IDX_CPUE_3_Fleet_1"),

    #2. observations for index ratio calculation
    n_recent_obs = 2,      #number of observations for recent period (Index A)
    n_historical_obs = 3,  #number of observations for historical reference (Index B)

    #3. how to combine multiple indices
    index_combination_method = "weighted_average",  # Options: "weighted_average", "multiplicative", "minimum"

    #4. Thr weigthing approach
    index_weights = c(1.0,0.8),  # match length of survey_index_names

    #5. Enable/disable index
    use_index_ratio = TRUE
  )


  #6. Length indicator configuration (this can be replaced by other L indicators - I just want to illustrate an example)

  length_config <- list(
    # enable/disable each indicator with TRUE or FALSE
    use_L_mean_L_mat = TRUE,    # L_mean (j-1) FD / L50 maturity
    use_L_mean_L_opt = FALSE,    # L_mean FD / L_optimal
    use_Pmature = TRUE,         # Proportion mature in catch (j-1)
    use_L_mean_LF_M = FALSE,     # L_mean FD / LF_M (Beverton-Holt)
    #Temporal trend indicator
    use_L_mean_ratio = TRUE,     #L_mean (j-1) / L_mean (average of previous 5 years)

    # weights for combining indicators (only for enabled indicators)
    weights = c(
      L_mean_L_mat = 0.25,
      #L_mean_L_opt = 0.25,
      Pmature = 0.25,
      #L_mean_LF_M = 0.25,
      L_mean_ratio = 0.30
    ),

    # target values (what we want each indicator to achieve)
    targets = list(
      L_mean_L_mat = 1.2,    # For example we may want mean length 20% above maturity
      #L_mean_L_opt = 1.0,    # or mean length at optimal
      Pmature = 0.50,        # or 50% of catch to be mature
      #L_mean_LF_M = 1.0      # or mean length at LF_M
      L_mean_ratio = 1.00    # new target no decline

    ),

    # penalty multipliers (TAC reduction when below target)
    # 1.0 = proportional penalty, >1.0 = more severe, <1.0 = less severe
    penalty_multipliers = c(
      L_mean_L_mat = 1.0,
      #L_mean_L_opt = 1.0,
      Pmature = 1.5,         #for example, more severe penalty for immature fish
      #L_mean_LF_M = 0.8,
      L_mean_ratio = 1.8     #strong penaly for decline
    ),

    # multipliers for TAC increase when above a target)
    # smaller than penalties to be more conservative
    reward_multipliers = c(
      L_mean_L_mat = 0.5,
      #L_mean_L_opt = 0.5,
      Pmature = 0.3,
      #L_mean_LF_M = 0.5
      L_mean_ratio = 0.4     #moderate reward if increasing

    ),

    # combination method for multiple length indicators
    combination_method = "weighted_average",  # Options: "weighted_average", "multiplicative", "minimum"

    # additional parameters
    Lc = 8.4,              # minimum length of capture (cm) - for LF_M calculation
    target_Pmature = 0.5, # target proportion mature (for Pmature, e.g.: target 50% mature)

    # Lookback period for L_mean_ratio calculation
    length_lookback_years = 5,  # compare current year to previous 5 years

    #enable/disable all indicators
    use_length_indicators = TRUE

  )


  #7. enable the use of biomass safeguard based on unfished spawning biomass
  #use_biomass_safeguard <- TRUE

  use_index_breaker <- TRUE
  breaker_index <- "IDX_CPUE_3_Fleet_1"  # MRIP index -choose monitoring index (CPUE or Survey)

  #8. trigger point as proportion of index
  breaker_lookback_years <- 5 # look back 5 years to find minimum index
  breaker_multiplier <- 0.50  # force 50% TAC when triggered

  #trigger_proportion <- 0.20  # Trigger at 20% of unfished biomass

  #9. precautionary multiplier
  precautionary_m <- 0.95 #5% reduction to maintain 95% probability above Blim

  #10. stability clauses
  max_increase <- 0.20  # +20% #Maximum annual TAC increase when stock is healthy
  max_decrease <- 0.30  # -30% # Maximum annual TAC decrease when stock is healthy
  stability_only_when_healthy <- TRUE #apply stability clause only when safeguard not triggered?

  # Maximum TAC increase based on length indicators
  max_tac_increase <- 1.2  # 20% maximum increase


  # It is important to mention here that:
  # Stability clause (7): Controls year-to-year TAC changes
  # max_increase = 0.20 = +20% cap
  # max_decrease = 0.30 = -30% floor
  # Only applies when b >= 1 (stock healthy)

  # The length indicator maximum (in combine_length_indicators)
  # max_tac_increase = 1.2 = caps individual length indicator adjustments at 120% (max increase 20%)
  # it is applied to each indicator adjustment factor before combining

  # Length indicators calculate f_length
  # - individual adjustments are capped at max_tac_increase (1.3)
  # - combined f_length could be anywhere from 0.01 to 1.3

  #TAC_preliminary <- TAC_previous × r × f_length × b × m

  #The stability clause comes later
  #TAC_final <- max(TAC_min, min(TAC_preliminary, TAC_max))

  # I have aligned the constrains!! both
  # - Length indicators can increase TAC by up to 20%
  # - Stability clause limits increases to only 20%

  # TAC increase is limited to 20% despite strong positive signals


  #11. How should TAC be split among fleets?
  #NULL = automatic (use historical catch proportions)
  #Vector = manual specification, e.g., c(0.65, 0.35) for 65%:35% split
  manual_fleet_proportions <- NULL #must sum to 1.0 if specified


  #12. Number of recent years for calculating fleet proportions (automatic mode)
  #only used when manual_fleet_proportions = NULL
  n_years_for_proportions <- 3 # 1= last year catch, 2= average two last years...etc)


  #13. Initial TAC:
  #NULL = use average of recent historical catches (see n_years_for_initial_TAC)
  #numeric value = override with specific TAC (e.g., 200 tonnes)
  initial_TAC_override <- NULL #initial_TAC_override <- 200  # uncomment to use fixed initial TAC

  #14. number of years to average for initial TAC calculation
  #only used when initial_TAC_override = NULL
  #examples:
  #   1 = use only last historical year catch
  #   3 = average last 3 years of historical catch
  #   5 = average last 5 years of historical catch
  n_years_for_initial_TAC <- 3

  #====================== End configuration section =================================#

  #unpack data object
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase == 1) {
    combined_data <- list()

    #needed for Index ratio (surveys and/or CPUE)
    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject) #returns 1-row tibble
      #add each column from the tibble to combined_data
      for(col_name in names(index_result)) {
        combined_data[[col_name]] <- index_result[[col_name]]
      }
    }

    #collect catch data (needed for fleet proportions and initial TAC)
    if(!is.null(CatchObsObj)) {
      catch_result <- calculate_single_CatchObs(dataObject)
      for(col_name in names(catch_result)) {
        combined_data[[col_name]] <- catch_result[[col_name]]
      }
    }

    #collect length composition data (if any length indicators enabled)
    any_length_enabled <- length_config$use_L_mean_L_mat ||
      length_config$use_L_mean_L_opt ||
      length_config$use_Pmature ||
      length_config$use_L_mean_LF_M||
      length_config$use_L_mean_ratio  #dded to check

    if(any_length_enabled) {
      if(is.null(LengthCompObj)) {
        stop("Length indicators are enabled but LengthCompObj is NULL. ",
             "This MP requires length composition data. ",
             "Provide LengthCompObj or disable all length indicators in configuration.")
      }

      lc_result <- calculate_single_LengthComp(dataObject)
      for(col_name in names(lc_result)) {
        combined_data[[col_name]] <- lc_result[[col_name]]
      }
    }

    return(combined_data)
  }

  if(phase == 2) {

    #Step 0: j index to real year

    #calculate the last year of the historical period
    yrHist <- TimeAreaObj@historicalYears + 1  # the end of historical period, e.g., 11
    # yrHist = end of historical period (e.g., year 11 in j-index)

    #real year for display (j-1 converts from j-index to simulation year)
    real_year <- j - 1

    # Example: j=11 = real_year = 10 (last historical year)
    # Example: j=12 = real_year = 11 (first projection year)

    # All observed data (FI and FD) uses previous year (j-1)

    obs_data_year  <- j - 1

    cat(sprintf("\n=== MULTI-INDEX-LENGTH MP: Year %d, Iteration %d ===\n", real_year, k))
    cat(sprintf("Configuration: %d indices, Length indicators: %s\n",
                length(index_config$survey_index_names),
                ifelse(length_config$use_length_indicators, "TRUE", "FALSE")))


    #Step 1: Determine previous TAC

    if(j == (yrHist + 1)) {
      # j = yrHist + 1 means we are at first projection year
      # e.g.,: if yrHist=11, then j=12 is first projection year (real year 11)

      if(!is.null(initial_TAC_override)) {
        #manual override: the user specified initial TAC
        TAC_previous <- initial_TAC_override
        cat(sprintf("Using manual initial TAC: %.2f\n", TAC_previous))

      } else {
        #Automatic calculation: average recent historical catches

        #calculate how many years back to start averaging
        #max(2, ...) ensures we never go before year 2 (year 1 is equilibrium)
        start_year <- max(2, yrHist - n_years_for_initial_TAC + 1)
        #e.g.: if yrHist=11 and n_years=3, start_year = max(2, 11-3+1) = 9
        # This means we average years j=9, j=10, and j=11 (real years 8, 9, 10)

        #initialize accumulator
        TAC_previous <- 0
        #Sum catches from all fleet-area combinations over the averaging window
        for(m in 1:areas) {
          for(f in 1:nfleets) {
            #column name for this fleet-area combination
            catch_col <- paste0("fleet_", f, "_observed_catch_area_", m)
            #e.g., "fleet_1_observed_catch_area_2"

            # Extract catches for this fleet-area in the averaging window
            avg_catches <- decisionData[[catch_col]][
              which(decisionData$k == k &
                      decisionData$j >= start_year &
                      decisionData$j <= yrHist)
            ]

            #e.g.,: if start_year=9, yrHist=11, this gets catches from j=9,10,11

            #add mean catch from this fleet-area to running total

            if(length(avg_catches) > 0 && !is.na(mean(avg_catches))) {
              TAC_previous <- TAC_previous + mean(avg_catches, na.rm = TRUE)
            }
          }
        }
        cat(sprintf("Initial TAC from averaging last %d years: %.2f\n",
                    n_years_for_initial_TAC, TAC_previous))
        cat(sprintf("  (Real years %d to %d, j-index %d to %d)\n",
                    start_year-1, yrHist-1, start_year, yrHist))
      }

      # Checking we successfully calculated a positive initial TAC
      # applies for both TAC:manual override or TAC:automatic calculation
      if(is.null(TAC_previous) || is.na(TAC_previous) || TAC_previous <= 0) {
        stop("cannot calculate initial TAC: TAC_previous = ", TAC_previous)
      }

    } else {
      # subsequent years projection
      # use previous year's TAC (sum across all fleet-area combinations)
      prev_year_TACs <- decisionAnnual$TAC[
        decisionAnnual$year == (j-1) &           #previous year in j-index
          decisionAnnual$iteration == k]           #current iteration

      #e.g.,: if j=13 (real year 12), we get TACs from j=12 (real year 11)

      #sum all fleet-area TACs from previous year
      TAC_previous <- sum(prev_year_TACs)
      cat(sprintf("TAC_previous from real year %d (j=%d): %.2f\n",
                  real_year-1, j-1, TAC_previous))
    }

    #Step 2: calculate biomass index ratio (r) for multiple indices

    # Purpose: measure recent biomass trend from survey data
    # Index A = mean of last 2 survey values (recent)
    # Index B = mean of 3 preceding values (historical)
    # r = Index A / Index B
    # r > 1: increasing trend, r < 1: decreasing trend

    # data: uses obs_data_year (j-1) - previous year data only

    r <- 1.0  # default if indices disabled
    index_ratios <- c()
    index_details <- list()

    if(index_config$use_index_ratio) {
      cat(sprintf("\n=== processing %d Index/Indices ===\n",
                  length(index_config$survey_index_names)))

      #calculate ratio for each index
      for(idx in 1:length(index_config$survey_index_names)) {
        index_name <- index_config$survey_index_names[idx]


        #validate index
        if(!index_name %in% names(decisionData)) {
          stop("Index '", index_name, "' not found in decisionData. ",
               "available indices: ", paste(grep("^IDX_", names(decisionData), value = TRUE), collapse = ", "))
        }



        # extract survey values up to determined maximum year
        survey_data <- decisionData[decisionData$k == k &
                                      decisionData$j <= obs_data_year,
                                    c("j", index_name)]
        # remove NA values (years without survey observations)
        survey_data <- survey_data[!is.na(survey_data[[index_name]]), ]

        #validating enough survey observations for trend calculation
        min_required <- index_config$n_recent_obs + index_config$n_historical_obs

        if(nrow(survey_data) < min_required) {
          stop("Index '", index_name, "': Insufficient data. Need ", min_required,
               " observations, found ", nrow(survey_data))
        }

        # calculate Index A (recent) - calculate Index B: mean of 3 values BEFORE the last 2
        Index_A <- mean(tail(survey_data[[index_name]], index_config$n_recent_obs), na.rm = TRUE)
        Index_B <- mean(tail(head(survey_data[[index_name]], -index_config$n_recent_obs),
                             index_config$n_historical_obs), na.rm = TRUE)

        # calculate ratio for this index
        r_index <- Index_A / Index_B
        index_ratios <- c(index_ratios, r_index)

        # store details
        index_details[[index_name]] <- list(
          Index_A = Index_A,
          Index_B = Index_B,
          ratio = r_index
        )

        cat(sprintf("Index %d (%s):\n", idx, index_name))
        cat(sprintf("  Index A (last %d): %.2f\n", index_config$n_recent_obs, Index_A))
        cat(sprintf("  Index B (previous %d): %.2f\n", index_config$n_historical_obs, Index_B))
        cat(sprintf("  Ratio: %.3f\n", r_index))
      }


      #combining multiple indices

      if(length(index_ratios) > 1) {
        cat(sprintf("\nCombining %d indices using '%s' method:\n",
                    length(index_ratios), index_config$index_combination_method))

        if(index_config$index_combination_method == "weighted_average") {
          # Validate weights
          if(length(index_config$index_weights) != length(index_ratios)) {
            stop("index_weights length must match number of indices")
          }
          weights <- index_config$index_weights / sum(index_config$index_weights)
          r <- sum(index_ratios * weights)
          cat(sprintf("  Weights: %s\n", paste(sprintf("%.2f", weights), collapse=", ")))
          cat(sprintf("  Combined r = %.3f\n", r))

        } else if(index_config$index_combination_method == "multiplicative") {
          r <- prod(index_ratios)
          cat(sprintf("  Product of ratios = %.3f\n", r))

        } else if(index_config$index_combination_method == "minimum") {
          r <- min(index_ratios)
          cat(sprintf("  Minimum ratio = %.3f\n", r))
        }
      } else {
        r <- index_ratios[1]
      }

      cat(sprintf("\nFinal biomass trend ratio (r): %.3f ", r))

      if(r > 1) {
        cat(sprintf("Stock INCREASING (r > 1)\n"))
      } else if(r < 1) {
        cat(sprintf("Stock DECLINING (r < 1)\n"))
      } else {
        cat(sprintf("Stock STABLE (r ~ 1)\n"))
      }
    } else {
      cat("\nIndex ratio: DISABLED (r = 1.0)\n")
    }



    #Step 3: length indicators

    f_length <- 1.0  # Default if disabled
    length_results <- list(L_mean = NA, indicators = list(), diagnostics = list())
    length_combination <- list(f_combined = 1.0, individual_adjustments = NULL,
                               method = "none", weights_used = NULL)

    if(length_config$use_length_indicators) {
      #check if any length indicators are enabled
      any_length_enabled <- length_config$use_L_mean_L_mat ||
        length_config$use_L_mean_L_opt ||
        length_config$use_Pmature ||
        length_config$use_L_mean_LF_M||
        length_config$use_L_mean_ratio

      if(any_length_enabled) {

        # get length bin setup first (needed for the last length indicator)
    length_bin_width <- LengthCompObj@length_bin_width
    max_length <- max(unlist(lh$L))
    length_bins <- seq(0, max_length + length_bin_width, by = length_bin_width)
    bin_midpoints <- length_bins[-length(length_bins)] + length_bin_width/2

    #get most recent length composition data
    current_lc <- decisionData[decisionData$k == k & decisionData$j == obs_data_year, ]

    cat(sprintf("\nLength indicator calculation:\n"))
    cat(sprintf("  Using length composition from real year %d (j=%d)\n",
                obs_data_year - 1, obs_data_year))

    #find all FD length bin columns
    bin_cols <- grep("LC_Fishery_.*_count_bin_", names(current_lc), value = TRUE)

    if(length(bin_cols) == 0) {
      stop("No FD length composition data available for year ", obs_data_year - 1, ". ",
           "This MP requires length composition data when length indicators are enabled. ",
           "Check LengthCompObj@survey_design configuration.")
    }

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



        # needed for proportion mature in catch (to get the probability of being mature at length of catch)
        calculate_maturity_ogive <- function(lengths, lh) {
          L50 <- lh$LifeHistory@L50
          L95 <- L50 + lh$LifeHistory@L95delta
          slope <- -(L95 - L50) / log(1/0.95 - 1)

          # apply logistic to get maturity probability at each length
          maturity <- 1 / (1 + exp(-(lengths - L50) / slope))
          return(maturity)
        }

        # calculate all length indicators
        calculate_length_indicators <- function(lh, length_config, bin_midpoints, total_counts) {

          indicators <- list() #  list of calculated indicator values
          diagnostics <- list() # reference used (L_mat, L_opt, etc.)

          # all indicators compare observed mean length to some reference length:
          # f > 1: fishing larger than reference
          # f < 1: fishing smaller than reference


          # calculate L_mean (needed for most indicators)
          # length indicators are based on this mean length from FD data
          # represents the average size of fish being caught
          if(sum(total_counts) > 0) {
            L_mean <- sum(bin_midpoints * total_counts) / sum(total_counts)

            # e.g., with 3 bins:
            # bins:    [8.5cm,  9.5cm,  10.5cm]
            # counts:  [10,     50,     30]
            # l_mean = (8.5×10 + 9.5×50 + 10.5×30) / (10+50+30)
            #          = (85 + 475 + 315) / 90
            #          = 875 / 90 = 9.72cm

          } else {
            stop("No length composition samples available for year ", obs_data_year - 1, ". ",
                 "This MP requires length composition data when length indicators are enabled.")
          }


          # 1. L_mean / L_mat (mean length relative to L50)
          # e.g.:
          #   f = 1.2:  fish 20% larger than 50% maturity
          #   f = 1.0:  fish right at 50% maturity
          #   f = 0.8:  fish 20% smaller than 50% maturity

          if(length_config$use_L_mean_L_mat) {
            L_mat <- lh$LifeHistory@L50
            indicators$L_mean_L_mat <- L_mean / L_mat
            diagnostics$L_mat <- L_mat
          } else {
            diagnostics$L_mat <- NULL
          }

          # 2. L_mean / L_opt (mean length relative to optimal length)
          # do we catch fish at the size that maximizes yield?
          #e.g.:
          #   f > 1.0: catching larger fish than optimal (less yield)
          #   f = 1.0: catching at optimal size (maximum yield)
          #   f < 1.0: catching smaller fish than optimal (growth overfishing)

          if(length_config$use_L_mean_L_opt) {
            MK_ratio <- lh$LifeHistory@M / lh$LifeHistory@K
            L_opt <- (3 * lh$LifeHistory@Linf) / (3 + MK_ratio)
            indicators$L_mean_L_opt <- L_mean / L_opt
            diagnostics$L_opt <- L_opt
          } else {
            diagnostics$L_opt <- NULL
          }

          # 3. Proportion mature in catch
          # what proportion of the catch is mature?- goal spawning Biomass proteccion
          # here we compare to a predetermined target

          #   f = 1.2: 20% above target
          #   f = 1.0:  at target
          #   f = 0.6: 40% below target

          #   e.g.,:
          #   bins:        [7.5cm,  8.5cm,  9.5cm,  10.5cm]
          #   counts:      [20,     60,     80,     40]
          #   maturity:    [0.15,   0.45,   0.75,   0.90]  (from ogive)
          #   Pmature = (20×0.15 + 60×0.45 + 80×0.75 + 40×0.90) / 200
          #           = (3 + 27 + 60 + 36) / 200 = 126/200 = 0.63
          #   If target = 0.5: f = 0.63/0.5 = 1.26 (26% above target)

          if(length_config$use_Pmature) {
            # Calculate maturity at each bin midpoint
            maturity_at_length <- calculate_maturity_ogive(bin_midpoints, lh)

            # proportion proportion mature in catch (weighted average)
            Pmature_value <- sum(total_counts * maturity_at_length) / sum(total_counts)

            # scale to target ("how close to target")
            # this makes it comparable to other indicators where f=1 is ideal
            indicators$Pmature <- Pmature_value / length_config$target_Pmature
            diagnostics$Pmature_value <- Pmature_value
          } else {
            diagnostics$Pmature_value <- NULL
          }

          # 4. L_mean / LF_M (Beverton-Holt as RFB)
          # LF_M (length when F=M, used in ICES RFB rule) (I know it is not optimal for BG, buthere we can replace by another L indicator)
          #   f > 1.0: catching larger than LF_M
          #   f = 1.0: catching at LF_M (MSY proxy)
          #   f < 1.0: catching smaller than LF_M


          if(length_config$use_L_mean_LF_M) {
            LF_M <- 0.75 * length_config$Lc + 0.25 * lh$LifeHistory@Linf
            indicators$L_mean_LF_M <- L_mean / LF_M
            diagnostics$LF_M <- LF_M
          } else {
            diagnostics$LF_M <- NULL
          }

          # 5. L_mean_ratio - Temporal trend in mean length
          # compares CURRENT YEAR (j-1) to PREVIOUS lookback YEARS
          # use to detect declining mean length in the most recent year
          # This is the ONLY indicator that uses a rolling window (5 years rolling window)
          # Current: year j-1 (obs_data_year)
          # Reference: average of years (j-1-lookback) through (j-2)

          # example with j=13, obs_data_year=12, lookback=5:
          #   current: year 12 L_mean
          #   reference: average L_mean from years 7, 8, 9, 10, 11
          #
          # Interpretation:
          #   f > 1.0: mean length INCREASING (good)
          #   f = 1.0: mean length STABLE
          #   f < 1.0: mean length DECLINING (bad)

          if(length_config$use_L_mean_ratio) {

            lookback <- length_config$length_lookback_years  # e.g., 5

            # CURRENT L_mean: already calculated above from obs_data_year (j-1)

            # REFERENCE L_mean: average of the PREVIOUS 'lookback' years
            # EXCLUDING the current year (obs_data_year) from the reference
            #
            # Example with obs_data_year = 12 and lookback = 5:
            #   Current: year 12
            #   Reference: years 7, 8, 9, 10, 11 (the 5 years BEFORE year 12)

            start_year_ref <- obs_data_year - lookback      # = 12 - 5 = 7
            end_year_ref <- obs_data_year - 1               # = 12 - 1 = 11


            # Verify getting exactly 'lookback' years
            n_ref_years <- end_year_ref - start_year_ref + 1
            if(n_ref_years != lookback) {
              warning(sprintf("Reference period has %d years, expected %d",
                              n_ref_years, lookback))
            }

            # extract length composition for reference period
            ref_lc <- decisionData[decisionData$k == k &
                                     decisionData$j >= start_year_ref &
                                     decisionData$j <= end_year_ref, ]

            # verify NOT including current year
            if(any(ref_lc$j == obs_data_year)) {
              stop("ERROR: Reference period incorrectly includes current observation year!")
            }

            # find FD length bins for reference period
            ref_bin_cols <- grep("LC_Fishery_.*_count_bin_", names(ref_lc), value = TRUE)

            # aggregate counts from reference period (5 years)
            ref_counts <- numeric(length(bin_midpoints))
            for(col in ref_bin_cols) {
              bin_idx <- as.numeric(sub(".*_bin_", "", col))
              if(bin_idx <= length(ref_counts)) {
                ref_counts[bin_idx] <- ref_counts[bin_idx] + sum(ref_lc[[col]], na.rm = TRUE)
              }
            }

            if(sum(ref_counts) > 0) {
              #we calculate reference L_mean (average over previous 'lookback' years)
              L_mean_reference <- sum(bin_midpoints * ref_counts) / sum(ref_counts)

              #maybe use here only the commercial lengths because they fish at deper waters
              # wider range of lengths

              indicators$L_mean_ratio <- L_mean / L_mean_reference
              diagnostics$L_mean_reference <- L_mean_reference
              diagnostics$L_mean_ref_years <- c(start_year_ref, end_year_ref)

              cat(sprintf("\n5. L_mean_ratio (Temporal Trend):\n"))
              cat(sprintf("   Current L_mean (year %d, j=%d): %.2f cm\n",
                          obs_data_year - 1, obs_data_year, L_mean))
              cat(sprintf("   Reference L_mean (years %d-%d, j=%d-%d): %.2f cm\n",
                          start_year_ref - 1, end_year_ref - 1,
                          start_year_ref, end_year_ref, L_mean_reference))
              cat(sprintf("   Number of reference years: %d\n", n_ref_years))
              cat(sprintf("   Ratio (current/reference): %.3f\n", indicators$L_mean_ratio))


              if(indicators$L_mean_ratio > 1.0) {
                cat(sprintf("    Status: Mean length INCREASING (good)\n"))
              } else if(indicators$L_mean_ratio < 1.0) {
                cat(sprintf("    Status: Mean length DECLINING (warning!)\n"))
              } else {
                cat(sprintf("    Status: Mean length STABLE\n"))
              }

            } else {
              stop(sprintf("No length composition data in reference period (years %d-%d)",
                           start_year_ref, end_year_ref))
            }
          } else {
            diagnostics$L_mean_reference <- NULL
            diagnostics$L_mean_ref_years <- NULL
          }


          return(list(
            indicators = indicators,       #list of f values for each enabled indicator
            L_mean = L_mean,               #observed mean length in catch
            diagnostics = diagnostics      #reference value used (L_mat, L_opt, etc.)
          ))
        }


        #combining length indicators
        #when using multiple length indicators, combine them into a single
        #multiplier (f_combined) for the TAC calculation
        #each indicator captures different aspects of size structure
        #L_mean/L_mat: Reproductive protection
        #L_mean/L_opt: Yield opti
        #Pmature: Spawning biomass
        #L_mean/LF_M: FMSY proxy
        #L_mean_ratio


        # There are 3 possible methods to combine the Length indicators:
        # 1. weighted_average: indicators based on importance (the user can allocate diff weigths)
        # 2. multiplicative:  all indicators must be good (more conservative)
        # 3. minimum: using the most pessimistic indicator (very conservative)


        combine_length_indicators <- function(indicators, length_config) {

          # incase where no indicators are enabled
          if(length(indicators) == 0) {
            return(list(
              f_combined = 1.0,   # no TAC adjustment
              individual_adjustments = NULL,
              method = "none",
              weights_used = NULL
            ))
          }

          indicator_names <- names(indicators)

          # calculate TAC adjustment for each indicator
          # for each indicator, we convert the observed/target ratio into a TAC multiplier
          # If indicator BELOW target (f < 1): apply penalty (reduce TAC)
          # If indicator ABOVE target (f > 1): apply rewrds  (increase TAC)

          # e.g.: penalty=1.0, reward=0.5 meaning:
          # 20% below target: 20% TAC reduction (more aggressive)
          # 20% above target: 10% TAC increase (more cautious)

          # a very conservative example for reduction:
          #   observed = 0.6, target = 1.0, penalty_mult = 1.5
          #   deviation = 1.0 - 0.6 = 0.4 (40% below)
          #   adjustment = 1 - 1.5 × (0.4/1.0) = 1 - 0.6 = 0.4
          #   TAC reduced by 60% (more than proportional)

          # a conservative example for increase:
          #   observed = 1.2, target = 1.0, reward_mult = 0.5
          #   deviation = 1.2 - 1.0 = 0.2 (20% above)
          #   adjustment = 1 + 0.5 × (0.2/1.0) = 1 + 0.1 = 1.1
          #   TAC increased by 10% (half of the improvement)


          adjustments <- sapply(indicator_names, function(name) {

            observed <- indicators[[name]]                            # e.g., f = 0.85 (15% below target)
            target <- length_config$targets[[name]]                   # e.g., target = 1.0
            penalty_mult <- length_config$penalty_multipliers[name]   # e.g., 1.0
            reward_mult <- length_config$reward_multipliers[name]     # e.g., 0.5

            if(observed < target) {
              # Below target: apply penalty
              deviation <- target - observed #1.0 - 0.85 = 0.15 (15% below)
              adjustment <- 1 - penalty_mult * (deviation / target) #1 - 1.0 × (0.15/1.0) = 1 - 0.15 = 0.85
              adjustment <- max(adjustment, 0.01)  # e.g., TAC reduced 15 - Floor at 1% of TAC
            } else {
              # Above target: apply reward
              deviation <- observed - target
              adjustment <- 1 + reward_mult * (deviation / target)
              adjustment <- min(adjustment, max_tac_increase)  # Cap at 130% of TAC (30% increase maximum)
            }

            return(adjustment)
          })

          # adjustments is now a named vector, e.g.:
          # c(L_mean_L_mat = 0.95, L_mean_L_opt = 1.10, Pmature = 0.85)

          # combine based on indicator method
          method <- length_config$combination_method

          #combine indicators based on their relative importance (weights)
          #all indicators contribute - higher weight = more influence on final TAC

          #For example (weighted_average)
          #L_mean/L_mat = 0.90 (10% below), weight = 0.25
          #L_mean/L_opt = 1.05 (5% above),  weight = 0.25
          #Pmature = 0.85 (15% below),      weight = 0.50 (higher concern)
          #
          #Normalized weights: [0.25, 0.25, 0.50] (already sum to 1)
          #f_combined = 0.25×0.90 + 0.25×1.05 + 0.50×0.85
          #              = 0.225 + 0.2625 + 0.425
          #              = 0.9125
          #   (1-0.9125 = 0.0875)TAC reduced by 8.75%


          if(method == "weighted_average") {
            # weighted average of adjustments
            weights <- length_config$weights[indicator_names]
            weights <- weights / sum(weights)     #weights are relative (normalized to sum to 1)
            f_combined <- sum(adjustments * weights)
            weights_used <- weights

            #For example (multiplicative)
            #f_combined = adjustment_1 × adjustment_2 × ... × adjustment_n
            #quite conservative: all indicators must be good for TAC increase
            #one poor indicator significantly reduces TAC
            #   L_mean/L_mat = 0.95 (5% reduction)
            #   L_mean/L_opt = 0.90 (10% reduction)
            #   Pmature = 0.85 (15% reduction)

            #   f_combined = 0.95 × 0.90 × 0.85 = 0.727
            #   TAC reduced by 27.3% (more reduction here compared when we allocated equal weights)

          } else if(method == "multiplicative") {
            # multiply all adjustments
            f_combined <- prod(adjustments)
            weights_used <- NULL


            #For example (minimum)- very pessimist (rebuilding or depleted stocks)
            # use the lowest adjustement
            #   L_mean/L_mat = 1.15 (15% increase possible)
            #   L_mean/L_opt = 1.08 (8% increase possible)
            #   Pmature = 0.75 (25% reduction needed)
            #   f_combined = min(1.15, 1.08, 0.75) = 0.75
            #   TAC reduced by 25%


          } else if(method == "minimum") {
            # most conservative (minimum adjustment)
            f_combined <- min(adjustments)
            weights_used <- NULL

          } else {
            stop("Unknown combination method: ", method)
          }

          return(list(
            f_combined = f_combined,                    #final combined multiplier for TAC
            individual_adjustments = adjustments,       #vector of individual indicator adjustments
            method = method,                            #method used for combination
            weights_used = weights_used                 #weights (NULL if not weighted_average)
          ))
        }


        # #get most recent length composition data
        # current_lc <- decisionData[decisionData$k == k & decisionData$j == obs_data_year, ]
        #
        # cat(sprintf("\nLength indicator calculation:\n"))
        # cat(sprintf("  Using length composition from real year %d (j=%d)\n",
        #             obs_data_year - 1, obs_data_year))
        #
        # #find all FD length bin columns
        # bin_cols <- grep("LC_Fishery_.*_count_bin_", names(current_lc), value = TRUE)
        #
        # if(length(bin_cols) == 0) {
        #   stop("No FD length composition data available for year ", obs_data_year - 1, ". ",
        #        "This MP requires length composition data when length indicators are enabled. ",
        #        "Check LengthCompObj@survey_design configuration.")
        # }
        #
        # # get length bin info
        # length_bin_width <- LengthCompObj@length_bin_width
        # max_length <- max(unlist(lh$L))
        # length_bins <- seq(0, max_length + length_bin_width, by = length_bin_width)
        # bin_midpoints <- length_bins[-length(length_bins)] + length_bin_width/2
        #
        # # aggregate counts across all FD programs
        # total_counts <- numeric(length(bin_midpoints))
        # for(col in bin_cols) {
        #   if(!is.na(current_lc[[col]])) {
        #     bin_idx <- as.numeric(sub(".*_bin_", "", col))
        #     if(bin_idx <= length(total_counts)) {
        #       total_counts[bin_idx] <- total_counts[bin_idx] + current_lc[[col]]
        #     }
        #   }
        # }
        # calculate all enabled length indicators
        length_results <- calculate_length_indicators(
          lh = lh,
          length_config = length_config,
          bin_midpoints = bin_midpoints,
          total_counts = total_counts
        )

        cat(sprintf("  L_mean = %.2f cm\n", length_results$L_mean))
        for(name in names(length_results$indicators)) {
          cat(sprintf("  %s = %.3f\n", name, length_results$indicators[[name]]))
        }

        # combine length indicators
        length_combination <- combine_length_indicators(
          indicators = length_results$indicators,
          length_config = length_config
        )

        f_length <- length_combination$f_combined

        cat(sprintf("\nLength indicator combination:\n"))
        cat(sprintf("  Method: %s\n", length_combination$method))
        if(!is.null(length_combination$weights_used)) {
          cat(sprintf("  Weights: %s\n",
                      paste(sprintf("%s=%.2f", names(length_combination$weights_used),
                                    length_combination$weights_used), collapse=", ")))
        }
        cat(sprintf("  Individual adjustments:\n"))
        for(name in names(length_combination$individual_adjustments)) {
          cat(sprintf("    %s: %.3f\n", name,
                      length_combination$individual_adjustments[name]))
        }
        cat(sprintf("  Combined f_length = %.3f\n", f_length))

      } else {
        # No length indicators enabled (any_length_enabled = FALSE)
        f_length <- 1.0
        length_results <- list(L_mean = NA, indicators = list(), diagnostics = list())
        length_combination <- list(f_combined = 1.0, individual_adjustments = NULL,
                                   method = "none", weights_used = NULL)
        cat("\nLength indicators: All disabled (f_length = 1.0)\n")
      }

    } else {

      # length_config$use_length_indicators = FALSE
      f_length <- 1.0
      length_results <- list(L_mean = NA, indicators = list(), diagnostics = list())
      length_combination <- list(f_combined = 1.0, individual_adjustments = NULL,
                                 method = "none", weights_used = NULL)
      cat("\nLength indicators: DISABLED (f_length = 1.0)\n")
    }


    #Step 4: calculate Safeguard (b)
    # Purpose: Reduce TAC when index approaches the index breaker
    if(use_index_breaker) {

      # get index history
      index_history <- decisionData[[breaker_index]][
        decisionData$k == k & decisionData$j <= obs_data_year
      ]
      index_history <- index_history[!is.na(index_history)]

      cat(sprintf("\n=== Index Breaker ===\n"))
      cat(sprintf("  Monitoring index: %s\n", breaker_index))
      cat(sprintf("  Index history: %d years available\n", length(index_history)))

      # need lookback + 1 years (e.g., 6 years for lookback=5)
      if(length(index_history) >= (breaker_lookback_years + 1)) {

        # current year index = obs_data_year (j-1) - last value
        I_current <- tail(index_history, 1)

        # previous 5 years (excluding current)
        I_previous <- tail(head(index_history, -1), breaker_lookback_years)
        I_min_reference <- min(I_previous)

        cat(sprintf("  Current index (year %d, j=%d): %.4f\n",
                    obs_data_year - 1, obs_data_year, I_current))
        cat(sprintf("  Min of previous %d years: %.4f\n",
                    breaker_lookback_years, I_min_reference))
        cat(sprintf("  Ratio (current/min): %.3f\n", I_current / I_min_reference))

        #breaker: if current < minimum --> trigger
        if(I_current < I_min_reference) {
          b <- breaker_multiplier
          cat(sprintf("\n BREAKER TRIGGERED!\n"))
          cat(sprintf("  Current (%.4f) < Min (%.4f)\n", I_current, I_min_reference))
          cat(sprintf("  Forcing b = %.2f (%.0f%% TAC reduction)\n",
                      b, (1 - b) * 100))
        } else {
          # not enough history yet
          b <- 1.0
          cat(sprintf("\n  Insufficient history (%d years available). Need %d years.\n",
                      length(index_history), breaker_lookback_years + 1))
          cat(sprintf("  No breaker applied (b = 1.0)\n"))
        }

      } else {
        b <- 1.0
        cat("\n breaker: DISABLED (b = 1.0)\n")
      }
    }



    #Step 5: apply precautionary multiplier (m)
    # Purpose: additional precautionary reduction (e.g., 0.95 = 5% reduction)
    # maintains higher probability of staying above limit reference points

    cat(sprintf("\nPrecautionary multiplier: m = %.2f\n", precautionary_m))
    if(precautionary_m < 1.0) {
      cat(sprintf(" Additional %.1f%% precautionary reduction\n",
                  (1 - precautionary_m) * 100))
    } else {
      cat("\n")
    }

    #Step 6: Calculate Preliminary TAC
    TAC_preliminary <- TAC_previous * r * f_length * b * precautionary_m

    cat(sprintf("\nPreliminary TAC calculation:\n"))
    cat(sprintf("  TAC_preliminary = TAC_previous × r × f_length × b × m\n"))
    cat(sprintf("  TAC_preliminary = %.2f × %.3f × %.3f × %.3f × %.2f\n",
                TAC_previous, r, f_length, b, precautionary_m))
    cat(sprintf("  TAC_preliminary = %.2f\n", TAC_preliminary))

    #Step 7: stability clause

    # Purpose: limit year-to-year TAC changes

    if(stability_only_when_healthy && b < 1.0) {

      #stock below trigger - allow full reduction
      TAC_final <- TAC_preliminary
      stability_applied <- sprintf("NO (safeguard triggered, b=%.2f)", b)
      cat(sprintf("\nStability clause: NO (safeguard takes priority)\n"))
      cat(sprintf("  TAC_final = %.2f\n", TAC_final))

    } else {
      #stock healthy or safeguard disabled, apply stability constraints
      TAC_max <- TAC_previous * (1 + max_increase)
      TAC_min <- TAC_previous * (1 - max_decrease)
      TAC_final <- max(TAC_min, min(TAC_preliminary, TAC_max))
      stability_applied <- sprintf("YES (+%.0f%%/-%.0f%%)", max_increase*100, max_decrease*100)


      cat(sprintf("\nStability clause: YES\n"))
      cat(sprintf("  TAC_min (-%.0f%%) = %.2f\n", max_decrease*100, TAC_min))
      cat(sprintf("  TAC_max (+%.0f%%) = %.2f\n", max_increase*100, TAC_max))
      cat(sprintf("  TAC_final = %.2f", TAC_final))

      if(TAC_final == TAC_preliminary) {
        cat(" (no constraint)\n")
      } else if(TAC_final == TAC_max) {
        cat(sprintf(" (CAPPED at +%.0f%%)\n", max_increase*100))
      } else {
        cat(sprintf(" (FLOORED at -%.0f%%)\n", max_decrease*100))
      }
    }

    #Step 8: Allocate TAC by fleet proportions

    # Purpose: allocate total TAC among fleets
    # Two modes:
    # - Manual: User-specified proportions
    # - Automatic: Calculate from recent historical catch patterns

    if(!is.null(manual_fleet_proportions)) {
      #MANUAL MODE: User specified fleet proportions

      if(length(manual_fleet_proportions) != nfleets) {
        stop(sprintf("manual_fleet_proportions must have %d elements", nfleets))
      }
      if(abs(sum(manual_fleet_proportions) - 1.0) > 1e-6) {
        manual_fleet_proportions <- manual_fleet_proportions / sum(manual_fleet_proportions)
        warning("Fleet proportions normalized to sum to 1")
      }
      fleet_proportions <- manual_fleet_proportions
      cat(sprintf("\nFleet allocation: MANUAL\n"))
      cat(sprintf("  Proportions: %s\n",
                  paste(sprintf("Fleet_%d=%.1f%%", 1:nfleets,
                                fleet_proportions*100), collapse = ", ")))

    } else {
      #AUTOMATIC MODE: Calculate from historical catch patterns

      # Data source: FD catches (uses recent years up to j-1)
      # cannot use current year catches (fishing has not occurred yet)


      #Determine averaging window
      start_year <- max(2, obs_data_year - n_years_for_proportions + 1)#ensures we do not go before year 2
      fleet_catches <- numeric(nfleets)

      for(f in 1:nfleets) {
        for(m in 1:areas) {
          #col name
          catch_col <- paste0("fleet_", f, "_observed_catch_area_", m)

          #extract catches in window (up to obs_data_year)
          recent_catches <- decisionData[[catch_col]][
            which(decisionData$k == k &
                    decisionData$j >= start_year &
                    decisionData$j <= obs_data_year)
          ]

          #add to fleet total
          fleet_catches[f] <- fleet_catches[f] + sum(recent_catches, na.rm = TRUE)
        }
      }

      #e.g., fleet_catches = [180, 120] (fleet 1 caught 180, fleet 2 caught 120)

      #calculate proportions

      if(sum(fleet_catches) > 0) {
        fleet_proportions <- fleet_catches / sum(fleet_catches)
      } else {
        #No historical catches, use equal proportions
        fleet_proportions <- rep(1/nfleets, nfleets)
      }

      cat(sprintf("\nFleet allocation: AUTOMATIC (%d-year average)\n",
                  n_years_for_proportions))
      cat(sprintf("  Years averaged: real years %d to %d (j=%d to %d)\n",
                  start_year-1, obs_data_year-1, start_year, obs_data_year))
      cat(sprintf("  Proportions: %s\n",
                  paste(sprintf("Fleet_%d=%.1f%%", 1:nfleets,
                                fleet_proportions*100), collapse = ", ")))
    }

    #Step 9: Allocate TAC to fleet area combination

    # Purpose: distribute each fleet TAC among areas based on catch patterns
    # For each fleet:
    #   1. Calculate fleet's total TAC (total × fleet_proportion)
    #   2. Calculate area proportions within that fleet (from recent catch)
    #   3. Allocate fleet TAC to areas

    #Initialize results storage
    TAC_decisions <- data.frame()
    start_year <- max(2, obs_data_year - n_years_for_proportions + 1)

    cat(sprintf("\nArea allocation within fleets:\n"))

    for(f in 1:nfleets) {
      #calculate this fleet's total TAC allocation
      fleet_TAC_total <- TAC_final * fleet_proportions[f]
      #e.g.: if TAC_final=100 and fleet_proportions[1]=0.6, fleet_TAC_total=60

      #calculate area proportions for this fleet
      #data source: FD catches (same window as fleet proportions)
      area_catches_for_fleet <- numeric(areas)
      for(m in 1:areas) {
        catch_col <- paste0("fleet_", f, "_observed_catch_area_", m)
        recent_catches <- decisionData[[catch_col]][
          which(decisionData$k == k &
                  decisionData$j >= start_year &
                  decisionData$j <= obs_data_year)
        ]
        area_catches_for_fleet[m] <- sum(recent_catches, na.rm = TRUE)
      }

      #e.g., fleet 1 caught [40, 20] in areas [1, 2]

      #calculate area proportions within this fleet
      if(sum(area_catches_for_fleet) > 0) {
        area_props <- area_catches_for_fleet / sum(area_catches_for_fleet)
        #eg.,: [40, 20] / 60 = [0.667, 0.333]
      } else {
        #no catches in this fleet, use equal proportions
        area_props <- rep(1/areas, areas)
      }

      cat(sprintf("  Fleet %d (%.1f%% of total):\n",
                  f, fleet_proportions[f]*100))

      #allocate fleet TAC to areas

      for(m in 1:areas) {
        TAC_fleet_area <- fleet_TAC_total * area_props[m]

        # e.g.,: fleet_TAC_total=60, area_props[1]=0.667
        #        TAC_fleet_area = 60 × 0.667 = 40

        # Store decision with full metadata

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
          n_indices_used = length(index_config$survey_index_names),
          index_combination_method = index_config$index_combination_method,

          #breaker components (replaces SB/B0)
          b = b,
          breaker_triggered = (b < 1.0),
          # Precautionary multiplier
          m = precautionary_m,

          #index details
          Index_A = Index_A,
          Index_B = Index_B,

          #length indicator components
          f_length = f_length,
          f_combination_method = length_combination$method,
          L_mean = length_results$L_mean,

          #individual length indicators (NA if not used)
          f_L_mean_L_mat = ifelse("L_mean_L_mat" %in% names(length_results$indicators),
                                  length_results$indicators$L_mean_L_mat, NA),
          f_L_mean_L_opt = ifelse("L_mean_L_opt" %in% names(length_results$indicators),
                                  length_results$indicators$L_mean_L_opt, NA),
          f_Pmature = ifelse("Pmature" %in% names(length_results$indicators),
                             length_results$indicators$Pmature, NA),
          f_L_mean_LF_M = ifelse("L_mean_LF_M" %in% names(length_results$indicators),
                                 length_results$indicators$L_mean_LF_M, NA),
          f_L_mean_ratio = ifelse("L_mean_ratio" %in% names(length_results$indicators),
                                  length_results$indicators$L_mean_ratio, NA),




          #reference
          L_mat = ifelse(!is.null(length_results$diagnostics$L_mat),
                         length_results$diagnostics$L_mat, NA),
          L_opt = ifelse(!is.null(length_results$diagnostics$L_opt),
                         length_results$diagnostics$L_opt, NA),
          LF_M = ifelse(!is.null(length_results$diagnostics$LF_M),
                        length_results$diagnostics$LF_M, NA),
          Pmature_value = ifelse(!is.null(length_results$diagnostics$Pmature_value),
                                 length_results$diagnostics$Pmature_value, NA),

          L_mean_reference = ifelse(!is.null(length_results$diagnostics$L_mean_reference),
                                    length_results$diagnostics$L_mean_reference, NA),


          # TAC calculation details

          # SB_ratio = if(use_biomass_safeguard) current_SB / B0 else NA,
          # trigger_proportion = if(use_biomass_safeguard) trigger_proportion else NA,
          TAC_previous = TAC_previous,
          TAC_preliminary = TAC_preliminary,
          stability_applied = stability_applied,

          # breaker info (replaces biomass safeguard)
          breaker_enabled = use_index_breaker,
          breaker_index_used = breaker_index,

          #safeguard_enabled = use_biomass_safeguard,
          obs_data_year_used = obs_data_year - 1,
          stringsAsFactors = FALSE
        ))

        cat(sprintf("    Area %d: TAC = %.2f (%.1f%% of fleet TAC)\n",
                    m, TAC_fleet_area, area_props[m]*100))
      }
    }

    cat("=====================================\n\n")

    # Return TAC decisions for all fleet-area combinations
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
      #these are the catch targets we need to achieve
      fleet_TACs_this_area <- TAC_decisions$TAC[TAC_decisions$area == m]
      # Example for Area 1 with 2 fleets: [40.7, 27.1]
      # fleet_TACs_this_area[1] = TAC for Fleet 1 in Area 1
      # fleet_TACs_this_area[2] = TAC for Fleet 2 in Area 1

      cat(sprintf("\n=== SOLVING Area %d ===\n", m))
      cat("Fleet TAC targets:",
          paste(sprintf("Fleet_%d=%.2f", 1:nfleets, fleet_TACs_this_area),
                collapse = ", "), "\n")

      # Prepare single-area, multi-fleet data structure for Solver
      # The Newton-Raphson solver is designed to work with one area at a time
      # We need to extract this area's data and present it as if it's the only area

      # Why this transformation?
      # - Original data: N[[gtg]][age, year, areas=1:n]
      # - Solver expects: N[[gtg]][age, year, areas=1]
      # - Solution: Extract area m's column and make it position 1


      #create temporary N structure containing only this area's abundance
      N_temp <- lapply(1:lh$gtg, function(gtg_idx) {
        temp_array <- array(0, dim = c(dim(N[[gtg_idx]])[1], ## ages
                                       dim(N[[gtg_idx]])[2], # years
                                       1))  # 1 area only    #one area
        #copy abundance from area m into position 1
        #Original: N[[gtg]][age, year, area=m]
        #New:  N_temp[[gtg]][age, year, area=1]
        temp_array[, , 1] <- N[[gtg_idx]][, , m]
        temp_array
      })

      #Result: N_temp contains area m's abundance in a 1-area structure

      #Extract selectivity for this area
      #Original structure: selGroup[[area=m]][[fleet=1:n]]
      #Solver expects:     selGroup[[area=1]][[fleet=1:n]]
      #Solution: Wrap area m's selectivity in a list to make it "area 1"
      selGroup_temp <- list(selGroup[[m]])  # wrap in list for 1-area structure

      #verify selectivity structure
      #ensure selectivity data was correctly extracted
      #and check sample values to confirm data structure is correct

      cat("\n--- SELECTIVITY VERIFICATION ---\n")
      for(f in 1:nfleets) {
        # Sample: GTG 1, Age 2, keep selectivity
        sample_sel <- selGroup_temp[[1]][[f]]$keep[[1]][5]
        cat(sprintf("Fleet %d | GTG 1 | Age 5 | Keep selectivity: %.4f\n",
                    f, sample_sel))
        # Example output: "Fleet 1 | GTG 1 | Age 5 | Keep selectivity: 0.8534"
        # This confirms:
        # - Fleet indexing is correct
        # - GTG indexing is correct
        # - Age indexing is correct
        # - Selectivity values are reasonable (0-1 range)
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

      # Result: F_vector = [F1, F2] - fishing mortality for each fleet
      # e.g.,: [0.142536, 0.095678] means Fleet 1 fishes at F=0.14, Fleet 2 at F=0.10

      cat("Solved F values:",
          paste(sprintf("Fleet_%d=%.6f", 1:nfleets, F_vector),
                collapse = ", "), "\n")

      # -------------------------------------------------------------------------
      # EXTRACT NEWTON-RAPHSON DIAGNOSTICS
      # -------------------------------------------------------------------------


      # STEP 1: #extract NR checks before fleet loop
      #extract NR checks (before fleet loop)
      #since both fleets are TAC-managed, attributes have length = nfleets

      # The solver returns one F_vector with attributes attached.
      # These attributes contain data for all TAC-managed fleets.
      # We extract them once for clarity
      # Inside loop - just use simple variable names


      converged_flag <- attr(F_vector, "converged") #TRUE if solver converged, FALSE if hit max iterations
      nr_iters <- attr(F_vector, "iterations")
      initial_guess_vec <- attr(F_vector, "initial_guess") #starting F values before iteration
      final_error_vec <- attr(F_vector, "final_error")     #relative error at convergence for each fleet
      predicted_catch_vec <- attr(F_vector, "predicted_catch") #actual catch achieved with solved F values
      target_catch_vec <- attr(F_vector, "target_catch")       #Original TAC targets

      #Step 2: identify TAC-managed fleets (both in this case) (determine which fleets are managed by TAC (versus effort))
      tac_managed_fleets <- which(!is.na(fleet_TACs_this_area))

      # Returns: [1, 2] if both fleets have TAC targets
      # Returns: [1] if only fleet 1 has TAC (fleet 2 effort-managed)

      #validate attribute lengths match expectations
      expected_attr_length <- length(tac_managed_fleets)

      if(!is.null(initial_guess_vec) && length(initial_guess_vec) != expected_attr_length) {
        warning(sprintf("Area %d: initial_guess length mismatch. Expected %d, got %d",
                        m, expected_attr_length, length(initial_guess_vec)))
      }


      # STEP3: mapping fleet to attribute index
      # 1. store F values for population dynamics (F_results)
      # 2. store convergence diagnostics for verification (nr_diagnostics_multifleet)


      for(f in 1:nfleets) {
        # add to F_results data frame
        # TAks 1: store F value for population dynamics
        # this data frame feeds back into evalMSE for next time step
        F_results <- rbind(F_results, data.frame(
          year = j,
          iteration = k,
          area = m,
          fleet = f,
          Flocal = F_vector[f] #solved fishing mortality
        ))


        #Task 2: store Newton-Raphson diagnostics (only for TAC-managed fleets)

        if(f %in% tac_managed_fleets) {

          #find position in attribute vectors
          #TAC-managed fleets may not be sequential if some are effort-managed
          #example: If fleet 1 and 3 are TAC-managed (fleet 2 is effort):
          #tac_managed_fleets = [1, 3]
          #For f=1: attr_index = 1
          #For f=3: attr_index = 2
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
            converged = if(is.null(converged_flag)) FALSE
            else converged_flag,
            nr_iterations = if(is.null(nr_iters)) NA_integer_
            else as.integer(nr_iters),
            nr_error = if(is.null(final_error_vec)) NA_real_
            else final_error_vec[attr_index],
            predicted_catch = if(is.null(predicted_catch_vec)) NA_real_
            else predicted_catch_vec[attr_index],
            target_catch = fleet_TACs_this_area[f],
            stringsAsFactors = FALSE
          ))
        }
      }

    }


    cat("=====================================\n\n")

    #return F values for all fleet-area combinations
    #this data frame is used by evalMSE to apply fishing mortality in population dynamics

    return(F_results)
  }
}


# Strategy objects
strategy_indexratio_length  <- new("Strategy")
strategy_indexratio_length @title <- "Multiple IRatio and Lindicators MP"
strategy_indexratio_length @projectionYears <- 5
strategy_indexratio_length @projectionName <- "multifleet_indexratio_length_MP_V2"
strategy_indexratio_length@projectionParams <- list()


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
    indexYears = seq(1, 15, 2),
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
    years = seq(1, 15, 1),
    sample_sizes = seq(100, 250, length.out = 15)
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

cat("\n=== INDEX RATIO - LENGTH MP TEST ===\n")


result_multi_indexratio_lengthV2<- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StrategyObj = strategy_indexratio_length,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_2fleet,
  IndexObj = multi_comprehensive_index,
  CatchObsObj = multi_comprehensive_catch,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test_multi_indexratio_length_MP",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "multifleet_indexratio_length_MP_V2"
)
cat(" multiple Index Ratio - lenght MP simulation completed\n")

# all debugging code reviewed!! and passed!!

result_multi_indexratio_lengthV2 <- readProjection(getwd(), "test_multi_indexratio_length_MP")


# Population outputs
result_multi_indexratio_lengthV2$dynamics$multifleet$actual_catch_proportions
result_multi_indexratio_lengthV2$dynamics$SB
result_multi_indexratio_lengthV2$dynamics$recN
result_multi_indexratio_lengthV2$dynamics$SPR

# obs model outputs (new structure data frame)
result_multi_indexratio_lengthV2$HCR$decisionData
result_multi_indexratio_lengthV2$HCR$decisionData$IDX_Survey_1
result_multi_indexratio_lengthV2$HCR$decisionLocal

result_multi_indexratio_lengthV2$HCR$decisionData$fleet_1_observed_catch_area_1
result_multi_indexratio_lengthV2$HCR$decisionData$fleet_1_observed_catch_area_2
result_multi_indexratio_lengthV2$HCR$decisionData$fleet_2_observed_catch_area_1
result_multi_indexratio_lengthV2$HCR$decisionData$fleet_2_observed_catch_area_2


result_multi_indexratio_lengthV2$HCR$decisionAnnual$TAC
result_multi_indexratio_lengthV2$HCR$decisionAnnual

result_multi_indexratio_lengthV2$HCR$decisionLocal
result_multi_indexratio_lengthV2$HCR$decisionAnnual


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
plot_SB(result_multi_indexratio_lengthV2)
plot_SB(result_multi_indexratio_lengthV2, areas=1)
plot_SB(result_multi_indexratio_lengthV2, areas=2)


plot_catchB(result_multi_indexratio_lengthV2)
plot_catchB(result_multi_indexratio_lengthV2,areas=1)
plot_catchB(result_multi_indexratio_lengthV2,areas=2)

plot_catchN(result_multi_indexratio_lengthV2)
plot_catchN(result_multi_indexratio_lengthV2, areas=1)
plot_catchN(result_multi_indexratio_lengthV2, areas=2)

plot_discB(result_multi_indexratio_lengthV2)
plot_discB(result_multi_indexratio_lengthV2,areas=1)
plot_discB(result_multi_indexratio_lengthV2,areas=2)

plot_discN(result_multi_indexratio_lengthV2)
plot_discN(result_multi_indexratio_lengthV2,areas=1)
plot_discN(result_multi_indexratio_lengthV2,areas=2)

plot_catchB_multi(result_multi_indexratio_lengthV2)
plot_catchB_multi(result_multi_indexratio_lengthV2, areas=1)
plot_catchB_multi(result_multi_indexratio_lengthV2, areas=2)

plot_catchN_multi(result_multi_indexratio_lengthV2,show_individual = TRUE)
plot_catchN_multi(result_multi_indexratio_lengthV2, areas=1)
plot_catchN_multi(result_multi_indexratio_lengthV2, areas=2)

plot_SPR(result_multi_indexratio_lengthV2)
plot_recN(result_multi_indexratio_lengthV2)


#plot obs models (indices)
plot_survey_indices(result_multi_indexratio_lengthV2)
plot_cpue_indices(result_multi_indexratio_lengthV2)

plot_all_indices(result_multi_indexratio_lengthV2)

#plot individual indices
plot_indices(result_multi_indexratio_lengthV2,
             index_pattern = "IDX_Survey_1",
             show_individual = TRUE,
             title = "FI survey")

# Total catch across all areas
plot_catch_observations_both(result_multi_indexratio_lengthV2,areas = c(1,2),show_individual=TRUE)
plot_catch_observations_both(result_multi_indexratio_lengthV2, areas = 1)
plot_catch_observations_both(result_multi_indexratio_lengthV2, areas = 2)


# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_multi_indexratio_lengthV2,show_individual = TRUE)
plot_fishery_length_comp(result_multi_indexratio_lengthV2, areas=1,show_individual = TRUE)
plot_fishery_length_comp(result_multi_indexratio_lengthV2, areas=2,show_individual = TRUE)

# plot_survey_length_comp(result_multi_indexratio_lengthV2,show_individual = TRUE)
# plot_survey_length_comp(result_multi_indexratio_lengthV2, areas=1,show_individual = TRUE)
# plot_survey_length_comp(result_multi_indexratio_lengthV2, areas=2,show_individual = TRUE)



# NEW: Custom filtering for fleets and areas
plot_length_composition_by_area(result_multi_indexratio_lengthV2,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets



plot_length_composition_by_area(result_multi_indexratio_lengthV2,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1,2),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets


# plot_length_composition_by_area(result_multi_RFB,
#                                 program_pattern = "LC_Survey",
#                                 area_filter = c(1,2),    # Specific areas
#                                 fleet_filter = c(2),
#                                 show_individual = TRUE)   # Specific fleets
#


#===========================================================================================#
#========================== Exploring NR outputs and performance ===========================#
#===========================================================================================#

#check SB
plot_SB(result_multi_indexratio_lengthV2,areas=1)
plot_SB(result_multi_indexratio_lengthV2,areas=2)

plot_Ftotal_multi(result_multi_indexratio_lengthV2,areas=c(1,2))
plot_Ftotal_multi(result_multi_indexratio_lengthV2,areas=c(1))
plot_Ftotal_multi(result_multi_indexratio_lengthV2,areas=c(2))

#cacth (add fleet)
plot_catchB_multi(result_multi_indexratio_lengthV2,areas=1)
plot_catchB_multi(result_multi_indexratio_lengthV2,areas=2)

#TAC (add fleet)
# plot_TAC_by_area(result_multi_RFB, areas = 1)
# plot_TAC_by_area(result_multi_RFB, areas = 2)
# plot_TAC_by_area(result_multi_RFB, areas = c(1,2))


plot_TAC(result_multi_indexratio_lengthV2, areas = "all",show_fleets = TRUE)  # All areas in one plot
#plot_TAC_by_fleet(result_indexratio)
plot_TAC(result_multi_indexratio_lengthV2, areas = c(1), show_fleets = TRUE)
plot_TAC(result_multi_indexratio_lengthV2, areas = c(2), show_fleets = TRUE)

result_multi_indexratio_lengthV2$HCR$decisionAnnual

# Area 1 only
plot_catch_observations_both(result_multi_indexratio_lengthV2, areas = 1)
# Area 2 only
plot_catch_observations_both(result_multi_indexratio_lengthV2, areas = 2)
# Both areas with faceting
plot_catch_observations_both(result_multi_indexratio_lengthV2, areas = c(1, 2))
#combined catch obs:
plot_catch_observations_both(result_multi_indexratio_lengthV2,areas = c(1, 2), show_individual = TRUE)


#calculate relative error (problem with area 2 -only)- fixed
nr_diag$relative_error <- abs(nr_diag$predicted_catch - nr_diag$target_catch) / nr_diag$target_catch
ggplot(nr_diag, aes(x = target_catch, y = predicted_catch, color = factor(area))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Target vs Predicted Catch - should be on 1:1 line")
#(given TAC target,find F that produces exactly that catch = passed)

#check for extreme values (extreme F values detected, convergencia failed)
extreme_F <- nr_diag[nr_diag$final_F > 2 | nr_diag$final_F < 0.001, ]
if(nrow(extreme_F) > 0) {
  cat("WARNING: Extreme F values detected\n")
  print(extreme_F)
}

#realized catch vs target TAC
catchB <- result_multi_indexratio_lengthV2$dynamics$catchB
TAC_decisions <- result_multi_indexratio_lengthV2$HCR$decisionAnnual


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
        realized_fleet <- result_multi_indexratio_lengthV2$dynamics$multifleet$catchB_by_fleet[yr, iter, area, fleet]

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
lh <- LHwrapper(result_multi_indexratio_lengthV2$LifeHistoryObj, result_multi_indexratio_lengthV2$TimeAreaObj)

#check multiple ages on the selectivity curve
test_ages <- c(1,2,3,4,5,6,7,8,9,10,11, 12, 13,14, 15)

for(age in test_ages) {
  cat(sprintf("\nAge %d:\n", age))

  fleet1_sel <- selWrapper(lh, result_multi_indexratio_lengthV2$TimeAreaObj,
                           FisheryObj = result_multi_indexratio_lengthV2$MultifleetObj@fleet_selectivity_proj_list[[1]][[1]],
                           doPlot = FALSE)

  fleet2_sel <- selWrapper(lh, result_multi_indexratio_lengthV2$TimeAreaObj,
                           FisheryObj = result_multi_indexratio_lengthV2$MultifleetObj@fleet_selectivity_proj_list[[1]][[2]],
                           doPlot = FALSE)

  f1_keep <- fleet1_sel$keep[[1]][age]
  f2_keep <- fleet2_sel$keep[[1]][age]

  cat(sprintf("  Fleet 1 (L50=8):  %.4f\n", f1_keep))
  cat(sprintf("  Fleet 2 (L50=12): %.4f\n", f2_keep))
  cat(sprintf("  Different? %s\n", abs(f1_keep - f2_keep) > 0.001))

}
