#MP1: RFB (Reference Fishing Biommass) rule (ICES)
# adding some reformulation to the Stragey - see below


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


# Shared random seed
test_seed <- 12345

# ============================================================================
# MANAGEMENT STRATEGIES - MULTIFLEET FLEET
# ============================================================================

# multifleet_rfb_MP: Multifleet strategy applying rfb rule
# adding some modifications and flexibility to account for survey timing before a
# a decision is made.

# Bill needs to review this to see if this makes sense.

# 1. It can happen that the survey (FI) from current real year is available (survey happened before TAC decision)
# so the index ratio is based taking into account the current observed year.

# 2. Survey (FI) from current real year is not available,
# for example it may after the TAC is defined or the data are not ready when the TAC need to be defined.

# 3. FD index: Always from previous year because the decision of TAC has not occurred yeat for current real year.
# The fishery has not occurred yet because no TAC has been defined.
# Catches and length comps MUST use previous year (fishery hasn't occurred yet)





#global storage for Newton-Raphson diagnostics
nr_diagnostics_multifleet <<- data.frame()

multifleet_rfb_MP  <- function(phase, dataObject) {

  # Some Indices Configurations
  # - Which survey index should be used for biomass trend calculation?
  # - Example: "IDX_Survey_1" uses the first fishery-independent survey

  #survey index to use for biomass trend
  survey_index_name <- "IDX_Survey_1"   # it could be CPUE (IDX_CPUE_1) or another survey

  #survey data availability
  # -current:     use current real year survey if available
  #               Use when: Survey occurs before TAC decision deadline
  # -previous:    always use previous real year survey
  #               Use when: Survey occurs late in year

  #IMPORTANT:     This ONLY affects FI (fishery-independent) survey indices
  #               FD (fishery-dependent) data ALWAYS uses previous year
  #               because fishing has not occurred yet in current year


  survey_data_timing <- "current"  # Options: "current" or "previous"


  #Other configurations
  #How should TAC be split among fleets?
  #NULL = automatic (use historical catch proportions)
  #Vector = manual specification, e.g., c(0.65, 0.35) for 65%:35% split
  #must sum to 1.0 if specified
  manual_fleet_proportions <- NULL


  #Fleet proportion calculation window
  #number of recent years to use for calculating fleet catch proportions
  #only used when manual_fleet_proportions = NULL
  n_years_for_proportions <- 3 # 1= last year catch, 2= average two last years...etc)


  #Initial TAC:
  #how should the first projection year's TAC be calculated?
  #NULL = use average of recent historical catches (see n_years_for_initial_TAC)
  #numeric value = override with specific TAC (e.g., 200 tonnes)
  initial_TAC_override <- NULL #initial_TAC_override <- 200  # uncomment to use fixed initial TAC

  #number of years to average for initial TAC calculation
  #only used when initial_TAC_override = NULL
  #examples:
  #   1 = use only last historical year catch
  #   3 = average last 3 years of historical catch
  #   5 = average last 5 years of historical catch
  n_years_for_initial_TAC <- 3

  #precautionary multiplier (maintain 95% probability above Blim)
  precautionary_m <- 0.95 #5% reduction to maintain 95% probability above Blim

  #Target length configuration (LF_M):
  #minimum capture length (cm) - used for Beverton-Holt LF_M calculation
  #should match fishery size regulations
  Lc <- 8.4  # cm - used to calculate LF/M (using this for now, we need to re-define)

  # We need to decide what approach to use here
  # LF_M calculation method has 2 options:
  # "beverton_holt" = 0.75×Lc + 0.25×Linf (accounts for size limits)
  # "direct_input"  = use value specified in LF_M_value
  LF_M_method <- "beverton_holt"

  #Direct LF_M value (only used if LF_M_method = "direct_input")
  #example: use L50 maturity as proxy
  #LF_M_value <- 8.4  # length at 50% maturity
  LF_M_value <- NULL

  #unpack data object
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

    #Step 0: j index to real year

    #calculate the last year of the historical period
    yrHist <- TimeAreaObj@historicalYears + 1  # the end of historical period, e.g., 11
    # yrHist = end of historical period (e.g., year 11 in j-index)

    #real year for display (j-1 converts from j-index to simulation year)
    real_year <- j - 1

    cat(sprintf("\n=== rfb RULE: Year %d, Iteration %d ===\n", real_year, k))

    # Example: j=11 = real_year = 10 (last historical year)
    # Example: j=12 = real_year = 11 (first projection year)


    #Step 1: Different data types have different availability

    # FI (Fishery-Independent) Survey Data:
    #   - "current":  can use current year (j) if survey already occurred
    #   - "previous": must use previous year (j-1) for complete annual data

    # FD (Fishery-Dependent) Data (catches, length compositions):
    #   - always use previous year (j-1)
    #   - because fishing hasn't occurred in current year yet (no TAC set)
    #   - cannot have current year fishery data before setting current year TAC


    # Determine maximum year index for FI survey data
    if(survey_data_timing == "current") {
      # assumes survey occurs before TAC decision in same real year
      survey_max_year <- j
      survey_timing_label <- "current (current year survey)"
      cat(sprintf("Survey data timing: %s\n", survey_timing_label))
      cat(sprintf("Using FI survey data through real year %d (j=%d)\n",
                  real_year, j))

    } else {  # "previous"
      # Waits for complete annual data before using in MP
      # Example: Use previous year data because real year survey is not available yet
      survey_max_year <- j - 1
      survey_timing_label <- "previous (previous year survey)"
      cat(sprintf("Survey data timing: %s\n", survey_timing_label))
      cat(sprintf("Using FI survey data through real year %d (j=%d)\n",
                  real_year - 1, j - 1))
    }

    # FD data ALWAYS uses previous year (j-1)
    # logical requirement
    fd_data_year <- j - 1
    cat(sprintf("Fishery-dependent data: ALWAYS previous year\n"))
    cat(sprintf("Using FD data (catches, length comps) from real year %d (j=%d)\n",
                real_year - 1, fd_data_year))
    cat(sprintf("Fishing has not occurred yet in real year %d\n\n", real_year))


    #Step 2: Determine previous TAC or use the last historical catch

    # first projection year: use last historical catch or manual override
    # subsequent years: use previous year's TAC

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

    #Step 3: calculate biomass index ratio (r)

    # Purpose: measure recent biomass trend from survey data
    # Index A = mean of last 2 survey values (recent)
    # Index B = mean of 3 preceding values (historical)
    # r = Index A / Index B
    # r > 1: increasing trend, r < 1: decreasing trend

    # data: FI survey (uses survey_max_year determined above)
    # - If "current": includes current year survey (j)
    # - If "previous": only previous year and earlier (j-1)

    # extract survey values up to determined maximum year
    survey_data <- decisionData[decisionData$k == k &
                                  decisionData$j <= survey_max_year,
                                c("j", survey_index_name)]
    # remove NA values (years without survey observations)
    survey_data <- survey_data[!is.na(survey_data[[survey_index_name]]), ]

    # calculate Index A: mean of LAST 2 survey observations
    Index_A <- mean(tail(survey_data[[survey_index_name]], 2), na.rm = TRUE)

    # calculate Index B: mean of 3 values BEFORE the last 2
    # first remove last 2 with head(..., -2), then take last 3
    Index_B <- mean(tail(head(survey_data[[survey_index_name]], -2), 3), na.rm = TRUE)

    #biomass trend ratio
    r <- Index_A / Index_B
    cat(sprintf("Index A (last 2): %.2f | Index B (previous 3): %.2f | r = %.3f\n",
                Index_A, Index_B, r))

    #for clarity, better to report which years were actually used
    survey_years_used <- survey_data$j - 1  # Convert to real years
    cat(sprintf("\nBiomass trend calculation (r):\n"))
    cat(sprintf("  Survey years available: %s (real years)\n",
                paste(survey_years_used, collapse = ", ")))
    cat(sprintf("Index A (last 2): %.2f\n", Index_A))
    cat(sprintf("Index B (previous 3): %.2f\n", Index_B))
    cat(sprintf("r = Index_A / Index_B = %.3f\n", r))
    if(r > 1) {
      cat(sprintf("Stock INCREASING (r > 1)\n"))
    } else if(r < 1) {
      cat(sprintf("Stock DECLINING (r < 1)\n"))
    } else {
      cat(sprintf("Stock STABLE (r ~ 1)\n"))
    }


    #Step 4: calculate length indicator (f)

    # Purpose: measure if fish are being caught at optimal size
    # L_mean = current mean length in catch (from length composition)
    # LF_M = target length (length at F=M, MSY proxy)
    # f = L_mean / LF_M
    # f > 1: fishing larger fish, f < 1: fishing smaller fish

    # Data source: FD length composition (ALWAYS uses fd_data_year = j-1)
    # Fishing has not occurred in current year yet

    # get most recent length composition data (from previous real year)
    current_lc <- decisionData[decisionData$k == k & decisionData$j == fd_data_year, ]

    cat(sprintf("\nLength indicator calculation (f):\n"))
    cat(sprintf("  Using length composition from real year %d (j=%d)\n",
                fd_data_year - 1, fd_data_year))

    #find all length bin columns in the data (FD data)
    bin_cols <- grep("LC_Fishery_.*_count_bin_", names(current_lc), value = TRUE)
    #e.g.:["LC_Fishery_1_Fleet_1_count_bin_1", "LC_Fishery_1_Fleet_1_count_bin_2", etc]

    if(length(bin_cols) == 0) {
      stop("No fishery-dependent length composition data available for year ",
           fd_data_year - 1)
    }


    #get length bin info from LengthCompObj
    length_bin_width <- LengthCompObj@length_bin_width  #e.g., 1 cm
    max_length <- max(unlist(lh$L))                     #maximum length across all GTGs
    length_bins <- seq(0, max_length + length_bin_width, by = length_bin_width)
    bin_midpoints <- length_bins[-length(length_bins)] + length_bin_width/2

    #aggregate counts across all FD programs
    total_counts <- numeric(length(bin_midpoints))
    for(col in bin_cols) {
      if(!is.na(current_lc[[col]])) {
        #extract bin index from column name
        bin_idx <- as.numeric(sub(".*_bin_", "", col))
        if(bin_idx <= length(total_counts)) {
          #add counts from this program to total
          total_counts[bin_idx] <- total_counts[bin_idx] + current_lc[[col]]
        }
      }
    }

      #calculate mean length from length frecuency obs
    if(sum(total_counts) > 0) {
      #weighted average: sum(midpoint × count) / sum(count)
      L_mean <- sum(bin_midpoints * total_counts) / sum(total_counts)
      #e.g.: bins [8.5, 9.5, 10.5], counts [10, 50, 30]
      #      L_mean = (8.5×10 + 9.5×50 + 10.5×30) / (10+50+30) = 9.72 cm
    } else {
      stop("No length composition samples available for year ", fd_data_year - 1)
    }

    # calculate target length (LF_M) using selected method

    if(LF_M_method == "direct_input") {
      # user provides specific target length (e.g., L50 maturity)
      if(is.null(LF_M_value)) {
        stop("LF_M_method = 'direct_input' requires LF_M_value to be specified")
      }
      LF_M <- LF_M_value
      cat(sprintf("  Target length method: DIRECT INPUT\n"))
      cat(sprintf("  LF_M = %.2f cm (user-specified)\n", LF_M))

    } else {
      # Beverton-Holt (1993) (default) - Kell et al 2022 formulation (LF_M = 0.75×Lc + 0.25×Linf)
      # accounts for minimum size regulations (Lc)
      LF_M <- 0.75 * Lc + 0.25 * lh$LifeHistory@Linf
      cat(sprintf("  Target length method: BEVERTON-HOLT\n"))
      cat(sprintf("  LF_M = 0.75×Lc + 0.25×Linf\n"))
      cat(sprintf("  Lc (min capture) = %.2f cm\n", Lc))
      cat(sprintf("  Linf = %.2f cm\n", lh$LifeHistory@Linf))
      cat(sprintf("  LF_M = %.2f cm\n", LF_M))
    }

      f_length <- L_mean / LF_M
      cat(sprintf("  L_mean (observed) = %.2f cm\n", L_mean))
      cat(sprintf("  f = L_mean / LF_M = %.3f\n", f_length))

    #Step 5: calculate biomass Safeguard (b)
      # Purpose: Reduce TAC when stock approaches limit reference point
      # I_trigger = I_loss × 1.4 (where I_loss is lowest observed index)
      # b = min(I_last / I_trigger, 1.0)
      # When I_last < I_trigger: b < 1, TAC is reduced
      # When I_last ≥ I_trigger: b = 1, no reduction


      #Data sourceE: Same FI survey data used for trend
      #get most recent survey observation
      I_last <- tail(survey_data[[survey_index_name]], 1)

      #Bill- we need to review this
      #find lowest ever observed survey value (proxy for Blim)
      I_loss <- min(survey_data[[survey_index_name]], na.rm = TRUE)

      #calculate trigger point (40% above lowest observed)
      I_trigger <- I_loss * 1.4

      #calculate safeguard multiplier (capped at 1.0)
      b <- min(I_last / I_trigger, 1.0)
      # If I_last > I_trigger: b = 1.0 (no reduction)
      # If I_last < I_trigger: b < 1.0 (proportional reduction)

      cat(sprintf("\nBiomass safeguard calculation (b):\n"))
      cat(sprintf("  I_last (most recent) = %.2f\n", I_last))
      cat(sprintf("  I_loss (lowest ever) = %.2f\n", I_loss))
      cat(sprintf("  I_trigger = I_loss × 1.4 = %.2f\n", I_trigger))
      cat(sprintf("  b = min(I_last/I_trigger, 1.0) = %.3f\n", b))
      if(b < 1.0) {
        cat(sprintf("SAFEGUARD TRIGGERED (b < 1): TAC will be reduced!\n"))
      } else {
        cat(sprintf("Stock above trigger (b = 1): No safeguard reduction\n"))
      }


    #Step 6: apply precautionary multiplier (m)
      # Purpose: additional precautionary reduction (e.g., 0.95 = 5% reduction)
      # maintains higher probability of staying above limit reference points

      cat(sprintf("\nPrecautionary multiplier: m = %.2f\n", precautionary_m))
      if(precautionary_m < 1.0) {
        cat(sprintf(" Additional %.1f%% precautionary reduction\n",
                    (1 - precautionary_m) * 100))
      }

    #Step 7: Calculate Preliminary TAC

      # Apply full RFB formula to previous TAC
      # TAC_preliminary = TAC_previous × r × f × b × m
      # Each component adjusts TAC based on stock status:
      # - r: biomass trend (up/down)
      # - f: size structure
      # - b: safeguard (triggered/not triggered)
      # - m: precautionary reduction

      TAC_preliminary <- TAC_previous * r * f_length * b * precautionary_m

      cat(sprintf("\nPreliminary TAC calculation:\n"))
      cat(sprintf("  TAC_preliminary = TAC_previous × r × f × b × m\n"))
      cat(sprintf("  TAC_preliminary = %.2f × %.3f × %.3f × %.3f × %.2f\n",
                  TAC_previous, r, f_length, b, precautionary_m))
      cat(sprintf("  TAC_preliminary = %.2f\n", TAC_preliminary))

    #Step 8: stability clause

      # Purpose: limit year-to-year TAC changes
      # When b ≥ 1 (stock above trigger):
      #   - Maximum increase: +20% of previous TAC
      #   - Maximum decrease: -30% of previous TAC
      # When b < 1 (stock below trigger):
      #   - No stability clause, allow full reduction from safeguard


      if(b >= 1.0) {
        #stock is healthy, apply stability constraints
        TAC_max <- TAC_previous * 1.2   # +20% cap
        TAC_min <- TAC_previous * 0.7   # -30% floor

        # Constrain TAC_preliminary within bounds
        TAC_final <- max(TAC_min, min(TAC_preliminary, TAC_max))
        # max(TAC_min, ...) ensures TAC ≥ 70% of previous
        # min(..., TAC_max) ensures TAC ≤ 120% of previous

        stability_applied <- "YES (±20%/-30%)"
        cat(sprintf("\nStability clause:\n"))
        cat(sprintf("  Stock above trigger (b ≥ 1) Apply stability limits\n"))
        cat(sprintf("  TAC_min (70%%) = %.2f\n", TAC_min))
        cat(sprintf("  TAC_max (120%%) = %.2f\n", TAC_max))
        cat(sprintf("  TAC_final = %.2f", TAC_final))

        if(TAC_final == TAC_preliminary) {
          cat(" (no constraint applied)\n")
        } else if(TAC_final == TAC_max) {
          cat(" (CAPPED at +20%%)\n")
        } else if(TAC_final == TAC_min) {
          cat(" (FLOORED at -30%%)\n")
        }


      } else {
        # Stock below trigger, allow full reduction
        TAC_final <- TAC_preliminary
        stability_applied <- sprintf("NO (b=%.2f)", b)

        cat(sprintf("\nStability clause:\n"))
        cat(sprintf("  Stock BELOW trigger (b < 1) NO stability limits\n"))
        cat(sprintf("  Safeguard reduction takes priority\n"))
        cat(sprintf("  TAC_final = %.2f\n", TAC_final))
      }


    #Step 9: Allocate TAC by fleet proportions

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
        start_year <- max(2, j - n_years_for_proportions)#ensures we do not go before year 2
        fleet_catches <- numeric(nfleets)

        for(f in 1:nfleets) {
          for(m in 1:areas) {
            #col name
            catch_col <- paste0("fleet_", f, "_observed_catch_area_", m)

            #extract catches in window (up to j-1, not including j)
            recent_catches <- decisionData[[catch_col]][
              which(decisionData$k == k &
                      decisionData$j >= start_year &
                      decisionData$j < j)  # < j, not <= j
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
                    start_year-1, j-2, start_year, j-1))
        cat(sprintf("  Proportions: %s\n",
                    paste(sprintf("Fleet_%d=%.1f%%", 1:nfleets,
                                  fleet_proportions*100), collapse = ", ")))
      }

    #Step 10: Allocate TAC to fleet area combination

      # Purpose: distribute each fleet TAC among areas based on catch patterns
      # For each fleet:
      #   1. Calculate fleet's total TAC (total × fleet_proportion)
      #   2. Calculate area proportions within that fleet (from recent catch)
      #   3. Allocate fleet TAC to areas

      #Initialize results storage
      TAC_decisions <- data.frame()
      start_year <- max(2, j - n_years_for_proportions)

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
                    decisionData$j < j) # < j, not <= j
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
            survey_data_timing = survey_data_timing,  # Record timing mode used
            survey_max_year_used = survey_max_year - 1,  # Real year
            fd_data_year_used = fd_data_year - 1,        # Real year
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
        # Sample: GTG 1, Age 10, keep selectivity
        sample_sel <- selGroup_temp[[1]][[f]]$keep[[1]][10]
        cat(sprintf("Fleet %d | GTG 1 | Age 10 | Keep selectivity: %.4f\n",
                    f, sample_sel))
        # Example output: "Fleet 1 | GTG 1 | Age 10 | Keep selectivity: 0.8534"
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

        } else {
          # effort-managed fleet (should not happen in RFB, but included anyway)
          # these fleets have F set directly from effort, no TAC solving needed
          nr_diagnostics_multifleet <<- rbind(nr_diagnostics_multifleet, data.frame(
            year = j,
            iteration = k,
            area = m,
            fleet = f,
            management_type = "effort",
            TAC = NA_real_,
            initial_F_guess = NA_real_,
            final_F = F_vector[f],
            converged = TRUE, #no solving needed, always "converged"
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

    #return F values for all fleet-area combinations
    #this data frame is used by evalMSE to apply fishing mortality in population dynamics

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



