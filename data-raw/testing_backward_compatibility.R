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

# Historical fishery - IDENTICAL for all tests
hist_fishery <- new("Fishery")
hist_fishery@title<-"Historical Fishery"
hist_fishery@vulType<-"logistic"
hist_fishery@vulParams<-c(9.0, 1.5)
hist_fishery@retType<-"full"
hist_fishery@retMax <- 1
hist_fishery@Dmort <- 0

# Projection fisheries - IDENTICAL for all tests
proj_fishery_area1 <- new("Fishery")
proj_fishery_area1@title<-"Proj Area 1"
proj_fishery_area1@vulType<-"logistic"
proj_fishery_area1@vulParams<-c(10.2, 0.1)
proj_fishery_area1@retType<-"logistic"
proj_fishery_area1@retParams <- c(10.2, 0.1)
proj_fishery_area1@retMax <- 1
proj_fishery_area1@Dmort <- 0

proj_fishery_area2 <- new("Fishery")
proj_fishery_area2@title<-"Proj Area 2"
proj_fishery_area2@vulType<-"logistic"
proj_fishery_area2@vulParams<-c(10.2, 0.1)
proj_fishery_area2@retType<-"logistic"
proj_fishery_area2@retParams <- c(10.2, 0.1)
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

  if(phase==1) return(list())
  if(phase==2) return(list())
  if(phase==3) {
    year = rep(j, areas)
    iteration = rep(k, areas)
    area = 1:areas
    Flocal = rep(0.15, areas)  # Constant F = 0.15

    return(list(year=year, iteration=iteration, area=area, Flocal=Flocal))
  }
}

# Strategy for multifleet (with fleet column)
simpleMP_multi <- function(phase, dataObject) {
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase==1) return(list())
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
# TEST SCENARIO 1: SINGLE FLEET (BASELINE)
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
  wd = getwd(),
  fileName = "validation_single_fleet",
  seed = validation_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
)
cat("Single fleet simulation completed\n")


# ============================================================================
# TEST SCENARIO 2: MULTIFLEET WITH 2 FLEETS (IDENTICAL SELECTIVITY)
# ============================================================================

cat("=========================================\n")
cat("RUNNING TEST 2: MULTIFLEET (2 IDENTICAL)\n")
cat("=========================================\n")


#Create IDENTICAL fleets for comparison
fleet1_identical <- new("Fishery")
fleet1_identical@vulType <- "logistic"
fleet1_identical@vulParams <- c(9.0, 1.5)  # SAME as hist_fishery
fleet1_identical@retType <- "full"
fleet1_identical@retMax <- 1
fleet1_identical@Dmort <- 0

fleet2_identical <- new("Fishery")
fleet2_identical@vulType <- "logistic"
fleet2_identical@vulParams <- c(9.0, 1.5)  # SAME as hist_fishery
fleet2_identical@retType <- "full"
fleet2_identical@retMax <- 1
fleet2_identical@Dmort <- 0

multifleet_identical <- new("Multifleet")
multifleet_identical@nfleets <- 2
multifleet_identical@fleet_proportions <- c(0.6, 0.4)
multifleet_identical@allocation_type <- "catch"  # Use effort allocation for comparison
multifleet_identical@fleet_selectivity_list <- list(fleet1_identical, fleet2_identical)

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

# ============================================================================
# TEST SCENARIO 3: MULTIFLEET WITH 1 FLEET (SHOULD MATCH SINGLE FLEET)
# ============================================================================

cat("======================================\n")
cat("RUNNING TEST 3: MULTIFLEET (1 FLEET)\n")
cat("======================================\n")

fleet1_single <- new("Fishery")
fleet1_single@vulType <- "logistic"
fleet1_single@vulParams <- c(9.0, 1.5)  # SAME as hist_fishery
fleet1_single@retType <- "full"
fleet1_single@retMax <- 1
fleet1_single@Dmort <- 0

multifleet_single <- new("Multifleet")
multifleet_single@nfleets <- 1
multifleet_single@fleet_proportions <- c(1.0)  # 100% to single fleet
multifleet_single@allocation_type <- "catch"
multifleet_single@fleet_selectivity_list <- list(fleet1_single)

strategy_multi1 <- new("Strategy")
strategy_multi1@title <- "Multifleet 1 Fleet Validation"
strategy_multi1@projectionYears <- 5
strategy_multi1@projectionName <- "simpleMP_multi"
strategy_multi1@projectionParams <- list()

result_multi1 <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  ProFisheryObj_list = proj_fishery_list,
  StrategyObj = strategy_multi1,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_single,
  wd = getwd(),
  fileName = "validation_multifleet_1fleet",
  seed = validation_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
)

cat("Multifleet (1 fleet) simulation completed\n")



# ============================================================================
# TEST SCENARIO 4: MULTIFLEET WITH DIFFERENT SELECTIVITIES
# ============================================================================


cat("============================================\n")
cat("RUNNING TEST 4: MULTIFLEET (DIFFERENT SELECT)\n")
cat("============================================\n")

fleet1_diff <- new("Fishery")
fleet1_diff@vulType <- "logistic"
fleet1_diff@vulParams <- c(8.0, 1.2)
fleet1_diff@retType <- "full"
fleet1_diff@retMax <- 1
fleet1_diff@Dmort <- 0

fleet2_diff <- new("Fishery")
fleet2_diff@vulType <- "logistic"
fleet2_diff@vulParams <- c(11.0, 2.0)
fleet2_diff@retType <- "full"
fleet2_diff@retMax <- 1
fleet2_diff@Dmort <- 0

multifleet_different <- new("Multifleet")
multifleet_different@nfleets <- 2
multifleet_different@fleet_proportions <- c(0.7, 0.3)
multifleet_different@allocation_type <- "catch"  # Use catch allocation
multifleet_different@fleet_selectivity_list <- list(fleet1_diff, fleet2_diff)

strategy_multi_diff <- new("Strategy")
strategy_multi_diff@title <- "Multifleet Different Selectivities"
strategy_multi_diff@projectionYears <- 5
strategy_multi_diff@projectionName <- "simpleMP_multi"
strategy_multi_diff@projectionParams <- list()

result_multi_diff <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  ProFisheryObj_list = proj_fishery_list,
  StrategyObj = strategy_multi_diff,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_different,
  wd = getwd(),
  fileName = "validation_multifleet_different",
  seed = validation_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE
)

cat("Multifleet (different selectivities) simulation completed\n")

# ============================================================================
# VALIDATION ANALYSIS
# ============================================================================

# load all
result_single <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_single_fleet")
result_multi2 <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_multifleet_2identical")
result_multi1 <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_multifleet_1fleet")
result_multi_diff <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_multifleet_different")


extract_complete_timeseries <- function(result, test_name) {

  years <- dim(result$dynamics$SB)[1]
  iterations <- dim(result$dynamics$SB)[2]
  areas <- dim(result$dynamics$SB)[3]
  year_seq <- 1:years

  cat(sprintf("Processing %s: %d years, %d iterations, %d areas\n",
              test_name, years, iterations, areas))

  complete_data <- data.frame()
  if(is.null(result$dynamics$multifleet)) {
    # Single fleet case
    cat("  Single fleet mode\n")
    for(area in 1:areas) {
      # average across iterations for each metric
      area_data <- data.frame(
        test = test_name,
        year = year_seq,
        area = area,
        fleet = 1,
        nfleets = 1,
        allocation_type = "single_fleet",
        fleet_proportion = 1.0,

        # population-level metrics (averaged across iterations)
        SB = rowMeans(result$dynamics$SB[, , area]),
        VB = rowMeans(result$dynamics$VB[, , area]),
        RB = rowMeans(result$dynamics$RB[, , area]),

        # total catch and mortality (averaged across iterations)
        total_catchB = rowMeans(result$dynamics$catchB[, , area]),
        total_catchN = rowMeans(result$dynamics$catchN[, , area]),
        total_discB = rowMeans(result$dynamics$discB[, , area]),
        total_discN = rowMeans(result$dynamics$discN[, , area]),
        total_F = rowMeans(result$dynamics$Ftotal[, , area]),

        # fleet-specific metrics (same as total for single fleet)
        fleet_catchB = rowMeans(result$dynamics$catchB[, , area]),
        fleet_catchN = rowMeans(result$dynamics$catchN[, , area]),
        fleet_discB = rowMeans(result$dynamics$discB[, , area]),
        fleet_discN = rowMeans(result$dynamics$discN[, , area]),
        fleet_F = rowMeans(result$dynamics$Ftotal[, , area]),

        # proportions (always 1.0 for single fleet)
        fleet_catchB_proportion = 1.0,
        fleet_F_proportion = 1.0
      )

      complete_data <- rbind(complete_data, area_data)
    }
  } else {
    # Multifleet case
    mf <- result$dynamics$multifleet
    nfleets <- mf$nfleets

    cat(sprintf("  Multifleet mode: %d fleets, allocation: %s\n",
                nfleets, mf$allocation_type))

    for(area in 1:areas) {
      # get total values for this area (averaged across iterations)
      total_catchB_ts <- rowMeans(result$dynamics$catchB[, , area])
      total_F_ts <- rowMeans(result$dynamics$Ftotal[, , area])

      for(fleet in 1:nfleets) {

        # get fleet-specific values (averaged across iterations)
        fleet_catchB_ts <- rowMeans(mf$catchB_by_fleet[, , area, fleet])
        fleet_F_ts <- rowMeans(mf$Ftotal_by_fleet[, , area, fleet])


        # CORRECT proportion calculation
        fleet_catchB_prop <- ifelse(total_catchB_ts > 0, fleet_catchB_ts / total_catchB_ts, 0)
        fleet_F_prop <- ifelse(total_F_ts > 0, fleet_F_ts / total_F_ts, 0)

        area_fleet_data <- data.frame(
          test = test_name,
          year = year_seq,
          area = area,
          fleet = fleet,
          nfleets = nfleets,
          allocation_type = mf$allocation_type,
          fleet_proportion = mf$fleet_proportions[fleet],

          # population-level metrics (same for all fleets in an area)
          SB = rowMeans(result$dynamics$SB[, , area]),
          VB = rowMeans(result$dynamics$VB[, , area]),
          RB = rowMeans(result$dynamics$RB[, , area]),

          # total catch and mortality (same for all fleets in an area)
          total_catchB = rowMeans(result$dynamics$catchB[, , area]),
          total_catchN = rowMeans(result$dynamics$catchN[, , area]),
          total_discB = rowMeans(result$dynamics$discB[, , area]),
          total_discN = rowMeans(result$dynamics$discN[, , area]),
          total_F = total_F_ts,

          # fleet-specific metrics (averaged across iterations)
          fleet_catchB = rowMeans(mf$catchB_by_fleet[, , area, fleet]),
          fleet_catchN = rowMeans(mf$catchN_by_fleet[, , area, fleet]),
          fleet_discB = rowMeans(mf$discB_by_fleet[, , area, fleet]),
          fleet_discN = rowMeans(mf$discN_by_fleet[, , area, fleet]),
          fleet_F = rowMeans(mf$Ftotal_by_fleet[, , area, fleet]),

          # fleet proportions over time
          fleet_catchB_proportion = fleet_catchB_prop,
          fleet_F_proportion = fleet_F_prop
        )

        complete_data <- rbind(complete_data, area_fleet_data)
      }
    }

    # add allocation details for multifleet
    if(!is.null(mf$target_catch_proportions)) {
      cat(sprintf("  Target catch proportions: %s\n",
                  paste(round(mf$target_catch_proportions, 3), collapse = ", ")))
    }
    if(!is.null(mf$actual_catch_proportions)) {
      cat(sprintf("  Actual catch proportions: %s\n",
                  paste(round(mf$actual_catch_proportions, 3), collapse = ", ")))
    }
    if(!is.null(mf$final_effort_proportions)) {
      cat(sprintf("  Final effort proportions: %s\n",
                  paste(round(mf$final_effort_proportions, 3), collapse = ", ")))
    }
  }

  return(complete_data)
}

# ============================================================================
# EXTRACT DATA FROM ALL SCENARIOS
# ============================================================================
# extract time series for all scenarios
complete_single <- extract_complete_timeseries(result_single, "Single_Fleet")
complete_multi2 <- extract_complete_timeseries(result_multi2, "Multifleet_2_Identical")
complete_multi1 <- extract_complete_timeseries(result_multi1, "Multifleet_1_Fleet")
complete_multi_diff <- extract_complete_timeseries(result_multi_diff, "Multifleet_Different")

#combine all scenarios
all_complete_data <- rbind(
  complete_single,
  complete_multi2,
  complete_multi1,
  complete_multi_diff
)

#add period classification
all_complete_data <- all_complete_data %>%
  mutate(
    period = case_when(
      year == 1 ~ "Equilibrium",
      year <= 11 ~ "Historical",
      year > 11 ~ "Projection"
    ),
    period_year = case_when(
      year == 1 ~ 0,
      year <= 11 ~ year - 1,
      year > 11 ~ year - 11
    )
  )


cat(sprintf("\nCombined dataset: %d rows, %d columns\n",
            nrow(all_complete_data), ncol(all_complete_data)))


# ============================================================================
# CREATE SUMMARY BY PERIOD
# ============================================================================

# summary by test, area, fleet, and period
period_summary <- all_complete_data %>%
  group_by(test, area, fleet, nfleets, allocation_type, fleet_proportion, period) %>%
  summarise(
    years_in_period = n(),

    # population metrics
    SB_mean = mean(SB, na.rm = TRUE),
    SB_final = last(SB),
    VB_mean = mean(VB, na.rm = TRUE),
    VB_final = last(VB),

    # total metrics
    total_catchB_mean = mean(total_catchB, na.rm = TRUE),
    total_catchB_final = last(total_catchB),
    total_F_mean = mean(total_F, na.rm = TRUE),
    total_F_final = last(total_F),

    # fleet-specific metrics
    fleet_catchB_mean = mean(fleet_catchB, na.rm = TRUE),
    fleet_catchB_final = last(fleet_catchB),
    fleet_catchN_mean = mean(fleet_catchN, na.rm = TRUE),
    fleet_F_mean = mean(fleet_F, na.rm = TRUE),
    fleet_F_final = last(fleet_F),

    # proportions
    fleet_catchB_prop_mean = mean(fleet_catchB_proportion, na.rm = TRUE),
    fleet_F_prop_mean = mean(fleet_F_proportion, na.rm = TRUE),

    .groups = "drop"
  )

# extract allocation results for multifleet scenarios
allocation_summary <- data.frame()

scenarios <- list(
  "Multifleet_2_Identical" = result_multi2,
  "Multifleet_1_Fleet" = result_multi1,
  "Multifleet_Different" = result_multi_diff
)

for(test_name in names(scenarios)) {
  result <- scenarios[[test_name]]

  if(!is.null(result$dynamics$multifleet)) {
    mf <- result$dynamics$multifleet

    for(fleet in 1:mf$nfleets) {
      allocation_summary <- rbind(allocation_summary, data.frame(
        test = test_name,
        fleet = fleet,
        nfleets = mf$nfleets,
        allocation_type = mf$allocation_type,
        original_fleet_proportion = mf$fleet_proportions[fleet],
        target_catch_proportion = if(!is.null(mf$target_catch_proportions)) mf$target_catch_proportions[fleet] else NA,
        actual_catch_proportion = if(!is.null(mf$actual_catch_proportions)) mf$actual_catch_proportions[fleet] else NA,
        final_effort_proportion = if(!is.null(mf$final_effort_proportions)) mf$final_effort_proportions[fleet] else NA,
        allocation_error = if(!is.null(mf$target_catch_proportions) && !is.null(mf$actual_catch_proportions)) {
          abs(mf$actual_catch_proportions[fleet] - mf$target_catch_proportions[fleet])
        } else NA
      ))
    }
  }
}

# ============================================================================
# EXPORT TO CSV FILES
# ============================================================================


#fix the calculation of proportions, I did that manually in excel

# Export complete time series
write.csv(all_complete_data, "multifleet_validation_complete_timeseries.csv", row.names = FALSE)
cat("Exported: multifleet_validation_complete_timeseries.csv\n")

# Export period summary
write.csv(period_summary, "multifleet_validation_period_summary.csv", row.names = FALSE)
cat("Exported: multifleet_validation_period_summary.csv\n")

# Export allocation summary
write.csv(allocation_summary, "multifleet_validation_allocation_summary.csv", row.names = FALSE)
cat("Exported: multifleet_validation_allocation_summary.csv\n")












# # extract key metrics for comparison
# extract_metrics <- function(result, test_name) {
#   list(
#     test = test_name,
#     SB_final = result$dynamics$SB[16, 1, 1],  # Final year, iter 1, area 1
#     VB_final = result$dynamics$VB[16, 1, 1],
#     catchB_final = result$dynamics$catchB[16, 1, 1],
#     Ftotal_final = result$dynamics$Ftotal[16, 1, 1],
#     SB_hist_mean = mean(result$dynamics$SB[2:11, 1, 1]),  # Historical period
#     catchB_hist_mean = mean(result$dynamics$catchB[2:11, 1, 1]),
#     Ftotal_hist_mean = mean(result$dynamics$Ftotal[2:11, 1, 1]),
#     has_multifleet = !is.null(result$dynamics$multifleet),
#     nfleets = if(!is.null(result$dynamics$multifleet)) result$dynamics$multifleet$nfleets else 1
#   )
# }
#
# # extract from all scenarios
# metrics_single <- extract_metrics(result_single, "Single Fleet")
# metrics_multi2 <- extract_metrics(result_multi2, "Multifleet (2 identical)")
# metrics_multi1 <- extract_metrics(result_multi1, "Multifleet (1 fleet)")
# metrics_multi_diff <- extract_metrics(result_multi_diff, "Multifleet (different)")
#
# # combine into data frame for analysis
# validation_df <- bind_rows(
#   data.frame(metrics_single),
#   data.frame(metrics_multi2),
#   data.frame(metrics_multi1),
#   data.frame(metrics_multi_diff)
# )
#
# print(validation_df)
#
#
#
# # ============================================================================
# # BACKWARD COMPATIBILITY TESTS
# # ============================================================================
#
# # test 1: Single Fleet vs Multifleet (1 fleet) - SHOULD BE IDENTICAL
# tolerance <- 1e-6
# cat("TEST 1: Single Fleet vs Multifleet (1 fleet)\n")
# cat("============================================\n")
#
# sb_match <- all(abs(result_single$dynamics$SB - result_multi1$dynamics$SB) < tolerance)
# vb_match <- all(abs(result_single$dynamics$VB - result_multi1$dynamics$VB) < tolerance)
# catch_match <- all(abs(result_single$dynamics$catchB - result_multi1$dynamics$catchB) < tolerance)
# f_match <- all(abs(result_single$dynamics$Ftotal - result_multi1$dynamics$Ftotal) < tolerance)
#
# cat(sprintf("  SB arrays match: %s\n", ifelse(sb_match, "PASS", "FAIL")))
# cat(sprintf("  VB arrays match: %s\n", ifelse(vb_match, "PASS", "FAIL")))
# cat(sprintf("  Catch arrays match: %s\n", ifelse(catch_match, "PASS", "FAIL")))
# cat(sprintf("  F arrays match: %s\n", ifelse(f_match, "PASS", "FAIL")))
#
# backward_compatible <- sb_match && vb_match && catch_match && f_match
# cat(sprintf("  OVERALL: %s\n", ifelse(backward_compatible, "BACKWARD COMPATIBLE", "COMPATIBILITY ISSUE")))
#
# # test 2: Multifleet (2 identical) should have similar total results to single fleet
# cat("\nTEST 2: Single Fleet vs Multifleet (2 identical)\n")
# cat("===============================================\n")
#
# # when fleets are identical with effort allocation, totals should be very similar
# sb_similar <- all(abs(result_single$dynamics$SB - result_multi2$dynamics$SB) < 0.01)
# catch_similar <- all(abs(result_single$dynamics$catchB - result_multi2$dynamics$catchB) < 0.01)
#
# cat(sprintf("  SB similarity (<1%% diff): %s\n", ifelse(sb_similar, "PASS", "FAIL")))
# cat(sprintf("  Catch similarity (<1%% diff): %s\n", ifelse(catch_similar, "PASS", "FAIL")))
#
# # ============================================================================
# # MULTIFLEET-SPECIFIC VALIDATION
# # ============================================================================
#
# # test multifleet outputs exist and are reasonable
# if(!is.null(result_multi2$dynamics$multifleet)) {
#   mf <- result_multi2$dynamics$multifleet
#
#   cat("Multifleet (2 identical) Structure:\n")
#   cat(sprintf("  Number of fleets: %d\n", mf$nfleets))
#   cat(sprintf("  Fleet proportions: [%.3f, %.3f]\n", mf$fleet_proportions[1], mf$fleet_proportions[2]))
#   cat(sprintf("  Allocation type: %s\n", mf$allocation_type))
#
#   # check that fleet catches sum to total
#   total_fleet_catch <- mf$catchB_by_fleet[16,1,1,1] + mf$catchB_by_fleet[16,1,1,2]
#   total_catch <- result_multi2$dynamics$catchB[16,1,1]
#   catch_sum_match <- abs(total_fleet_catch - total_catch) < tolerance
#
#   cat(sprintf("  Fleet catches sum to total: %s\n", ifelse(catch_sum_match, "PASS", "FAIL")))
#
#   # check fleet F values sum to total F
#   total_fleet_F <- mf$Ftotal_by_fleet[16,1,1,1] + mf$Ftotal_by_fleet[16,1,1,2]
#   total_F <- result_multi2$dynamics$Ftotal[16,1,1]
#   f_sum_match <- abs(total_fleet_F - total_F) < tolerance
#
#   cat(sprintf("  Fleet F values sum to total: %s\n", ifelse(f_sum_match, "PASS", "FAIL")))
# }
#
# # test catch allocation results for different selectivities
# if(!is.null(result_multi_diff$dynamics$multifleet)) {
#   mf_diff <- result_multi_diff$dynamics$multifleet
#
#   cat("\nMultifleet (different selectivities) Results:\n")
#   cat(sprintf("  Target catch proportions: [%.3f, %.3f]\n",
#               mf_diff$target_catch_proportions[1], mf_diff$target_catch_proportions[2]))
#   cat(sprintf("  Actual catch proportions: [%.3f, %.3f]\n",
#               mf_diff$actual_catch_proportions[1], mf_diff$actual_catch_proportions[2]))
#   cat(sprintf("  Final effort proportions: [%.3f, %.3f]\n",
#               mf_diff$final_effort_proportions[1], mf_diff$final_effort_proportions[2]))
#
#   # check if catch allocation achieved target (within tolerance)
#   catch_target_achieved <- all(abs(mf_diff$actual_catch_proportions - mf_diff$target_catch_proportions) < 0.05)
#   cat(sprintf("  Catch allocation within 5%%: %s\n", ifelse(catch_target_achieved, "PASS", "NEEDS TUNING")))
# }
#
# # ============================================================================
# # VISUAL VALIDATION PLOTS
# # ============================================================================
#
# # create comparison plots
# create_validation_plots <- function() {
#   # prepare data for plotting
#   years <- 1:16
#
#   # SB comparison
#   sb_data <- data.frame(
#     Year = rep(years, 4),
#     SB = c(result_single$dynamics$SB[,1,1],
#            result_multi1$dynamics$SB[,1,1],
#            result_multi2$dynamics$SB[,1,1],
#            result_multi_diff$dynamics$SB[,1,1]),
#     Test = rep(c("Single Fleet", "Multifleet (1)", "Multifleet (2 id)", "Multifleet (diff)"), each = 16)
#   )
#
#   p1 <- ggplot(sb_data, aes(x = Year, y = SB, color = Test, linetype = Test)) +
#     geom_line(size = 1) +
#     geom_vline(xintercept = 11.5, linetype = "dashed", alpha = 0.5) +
#     labs(title = "Spawning Biomass Comparison",
#          subtitle = "Vertical line = end of historical period",
#          y = "Spawning Biomass") +
#     theme_minimal()
#
#   # catch comparison
#   catch_data <- data.frame(
#     Year = rep(years, 4),
#     Catch = c(result_single$dynamics$catchB[,1,1],
#               result_multi1$dynamics$catchB[,1,1],
#               result_multi2$dynamics$catchB[,1,1],
#               result_multi_diff$dynamics$catchB[,1,1]),
#     Test = rep(c("Single Fleet", "Multifleet (1)", "Multifleet (2 id)", "Multifleet (diff)"), each = 16)
#   )
#
#   p2 <- ggplot(catch_data, aes(x = Year, y = Catch, color = Test, linetype = Test)) +
#     geom_line(size = 1) +
#     geom_vline(xintercept = 11.5, linetype = "dashed", alpha = 0.5) +
#     labs(title = "Catch Biomass Comparison",
#          subtitle = "Vertical line = end of historical period",
#          y = "Catch Biomass") +
#     theme_minimal()
#
#   # F comparison
#   f_data <- data.frame(
#     Year = rep(years, 4),
#     F = c(result_single$dynamics$Ftotal[,1,1],
#           result_multi1$dynamics$Ftotal[,1,1],
#           result_multi2$dynamics$Ftotal[,1,1],
#           result_multi_diff$dynamics$Ftotal[,1,1]),
#     Test = rep(c("Single Fleet", "Multifleet (1)", "Multifleet (2 id)", "Multifleet (diff)"), each = 16)
#   )
#
#   p3 <- ggplot(f_data, aes(x = Year, y = F, color = Test, linetype = Test)) +
#     geom_line(size = 1) +
#     geom_vline(xintercept = 11.5, linetype = "dashed", alpha = 0.5) +
#     labs(title = "Fishing Mortality Comparison",
#          subtitle = "Vertical line = end of historical period",
#          y = "Fishing Mortality (F)") +
#     theme_minimal()
#
#   return(list(sb = p1, catch = p2, f = p3))
# }
#
# plots <- create_validation_plots()
#
# # save plots
# ggsave("validation_SB_comparison.jpeg", plots$sb, width = 10, height = 6, dpi = 300)
# ggsave("validation_catch_comparison.jpeg", plots$catch, width = 10, height = 6, dpi = 300)
# ggsave("validation_F_comparison.jpeg", plots$f, width = 10, height = 6, dpi = 300)
#
# cat("Validation plots saved\n")
#











#debugging

# cat("================================================================\n")
# cat("TESTING MULTIFLEET DETECTION LOGIC\n")
# cat("================================================================\n\n")
#
# test_detection_with_real_objects <- function() {
#
#   # Test 1: NULL (single fleet)
#   cat("TEST 1: NULL MultifleetObj (Single Fleet)\n")
#   cat("-----------------------------------------\n")
#   MultifleetObj <- NULL
#   is_multifleet_current <- !is.null(MultifleetObj) && MultifleetObj@nfleets > 1
#   is_multifleet_fixed <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1
#
#   cat(sprintf("  MultifleetObj: NULL\n"))
#   cat(sprintf("  Current logic (>1):  %s\n", is_multifleet_current))
#   cat(sprintf("  Fixed logic (>=1):   %s\n", is_multifleet_fixed))
#   cat(sprintf("  Expected: FALSE (single fleet)\n"))
#   cat(sprintf("  Status: %s\n\n", ifelse(is_multifleet_current == FALSE, "CORRECT", "WRONG")))
#
#   # Test 2: 1 fleet multifleet
#   cat("TEST 2: 1 Fleet Multifleet\n")
#   cat("---------------------------\n")
#
#   # Create real Multifleet object with 1 fleet
#   fishery1 <- new("Fishery")
#   fishery1@vulType <- "logistic"
#   fishery1@vulParams <- c(9.0, 1.5)
#   fishery1@retType <- "full"
#   fishery1@retMax <- 1
#   fishery1@Dmort <- 0
#
#   MultifleetObj <- new("Multifleet")
#   MultifleetObj@nfleets <- 1
#   MultifleetObj@fleet_proportions <- c(1.0)
#   MultifleetObj@allocation_type <- "effort"
#   MultifleetObj@fleet_selectivity_list <- list(fishery1)
#
#   is_multifleet_current <- !is.null(MultifleetObj) && MultifleetObj@nfleets > 1
#   is_multifleet_fixed <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1
#
#   cat(sprintf("  nfleets: %d\n", MultifleetObj@nfleets))
#   cat(sprintf("  Current logic (>1):  %s\n", is_multifleet_current))
#   cat(sprintf("  Fixed logic (>=1):   %s\n", is_multifleet_fixed))
#   cat(sprintf("  Expected: TRUE (should use multifleet structures)\n"))
#   cat(sprintf("  Status: %s ← THIS IS THE ISSUE!\n\n",
#               ifelse(is_multifleet_current == TRUE, "CORRECT", "WRONG")))
#
#   # Test 3: 2 fleet multifleet
#   cat("TEST 3: 2 Fleet Multifleet\n")
#   cat("---------------------------\n")
#
#   fishery2 <- new("Fishery")
#   fishery2@vulType <- "logistic"
#   fishery2@vulParams <- c(11.0, 2.0)
#   fishery2@retType <- "full"
#   fishery2@retMax <- 1
#   fishery2@Dmort <- 0
#
#   MultifleetObj@nfleets <- 2
#   MultifleetObj@fleet_proportions <- c(0.6, 0.4)
#   MultifleetObj@fleet_selectivity_list <- list(fishery1, fishery2)
#
#   is_multifleet_current <- !is.null(MultifleetObj) && MultifleetObj@nfleets > 1
#   is_multifleet_fixed <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1
#
#   cat(sprintf("  nfleets: %d\n", MultifleetObj@nfleets))
#   cat(sprintf("  Current logic (>1):  %s\n", is_multifleet_current))
#   cat(sprintf("  Fixed logic (>=1):   %s\n", is_multifleet_fixed))
#   cat(sprintf("  Expected: TRUE (multifleet)\n"))
#   cat(sprintf("  Status: %s\n\n", ifelse(is_multifleet_current == TRUE, "CORRECT", "WRONG")))
#
#   return(list(
#     test1_correct = is_multifleet_current == FALSE,
#     test2_correct = is_multifleet_fixed == TRUE,  # Should be TRUE with fixed logic
#     test3_correct = is_multifleet_current == TRUE
#   ))
# }
#
# # Run the test
# results <- test_detection_with_real_objects()
#
# cat("================================================================\n")
# cat("SUMMARY\n")
# cat("================================================================\n\n")
#
# cat("PROBLEM IDENTIFIED:\n")
# cat("When nfleets = 1, current logic returns FALSE\n")
# cat("but multifleet data structures are already created!\n\n")
#
# cat("SOLUTION:\n")
# cat("Change detection logic from:\n")
# cat("  is_multifleet <- !is.null(MultifleetObj) && MultifleetObj@nfleets > 1\n")
# cat("To:\n")
# cat("  is_multifleet <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1\n\n")
#
# cat("WHAT HAPPENS WITH CURRENT LOGIC:\n")
# cat("1. MultifleetObj with 1 fleet created\n")
# cat("2. 4D arrays (histEffortDev) created for multifleet\n")
# cat("3. Detection logic says is_multifleet = FALSE\n")
# cat("4. Code tries to use single fleet logic on multifleet arrays\n")
# cat("5. rbind() fails due to column mismatch\n\n")
#
# cat("WHAT HAPPENS WITH FIXED LOGIC:\n")
# cat("1. MultifleetObj with 1 fleet created\n")
# cat("2. 4D arrays created for multifleet\n")
# cat("3. Detection logic says is_multifleet = TRUE\n")
# cat("4. Code uses multifleet logic consistently\n")
# cat("5. All functions return same column structure\n\n")
#
# cat("CONCLUSION: change > 1 to >= 1\n")
#
# cat("CHANGES:\n")
# cat("MSEWrappers.R (evalMSE function)\n")
# cat("managementStrategies.R (fixedStrategy function)\n")
# cat("MSEWrappers.R (runProjection function)\n\n")
#
#
#
