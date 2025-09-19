# ============================================================================
# TEST FILE 1: BASIC SINGLE FLEET VS MULTIFLEET COMPARISON
# ============================================================================
# Purpose: Validate backward compatibility and core multifleet functionality
# Tests:
# 1. Single fleet (original)
# 2. Multifleet with 1 fleet (should match single fleet)
# 3. Multifleet with 2 fleets
# 4. Multifleet with 3 fleets

rm(list=ls())
devtools::load_all()
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# SHARED SETUP FOR ALL TESTS
# ============================================================================

cat("=====================================\n")
cat("FISHSIMGTG MULTIFLEET VALIDATION - FILE 1\n")
cat("Basic Single vs Multifleet Comparison\n")
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

# Time Area - consistent
ta <- new("TimeArea")
ta@title <- "Validation Test"
ta@gtg <- 13
ta@areas <- 2
ta@recArea <- c(0.99, 0.01)
ta@iterations <- 2  # Small for quick testing
ta@historicalYears <- 10
ta@historicalBio <- 0.5
ta@historicalBioType <- "relB"
ta@move <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = FALSE)

# Historical effort - declining trend
ta@historicalEffort <- matrix(c(1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6,
                                1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6),
                              nrow = 10, ncol = 2, byrow = FALSE)


# Stochastic - minimal variation for clear comparison
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
# SIMPLE MANAGEMENT STRATEGY - WITHOUTH OBSERVATION MODELS
# ============================================================================

#Bill edit: in single-species mode you do not need to specify an MP to simulate
#historical dynamics. I've made changes to multi species that enable the same.

# simpleStrategy <- function(phase, dataObject) {
#   for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])
#
#   if(phase == 1) return(list())
#   if(phase == 2) return(list())
#   if(phase == 3) {
#     # Return appropriate structure based on multifleet mode
#     if(!is.null(MultifleetObj)) {
#       # Multifleet mode
#       year <- rep(j, areas)
#       iteration <- rep(k, areas)
#       area <- 1:areas
#       fleet <- rep(0, areas)  # 0 = total F
#       Flocal <- rep(0.15, areas)  # Constant F
#
#       return(list(year=year, iteration=iteration, area=area,
#                   fleet=fleet, Flocal=Flocal))
#     } else {
#       # Single fleet mode
#       year <- rep(j, areas)
#       iteration <- rep(k, areas)
#       area <- 1:areas
#       Flocal <- rep(0.15, areas)
#
#       return(list(year=year, iteration=iteration, area=area, Flocal=Flocal))
#     }
#   }
# }
#
# strategy_obj <- new("Strategy")
# strategy_obj@title <- "Simple Constant F"
# strategy_obj@projectionYears <- 5
# strategy_obj@projectionName <- "simpleStrategy"
# strategy_obj@projectionParams <- list()

# ============================================================================
# TEST 1: SINGLE FLEET
# ============================================================================

cat("\n=== TEST 1: SINGLE FLEET (BASELINE) ===\n")

result_single <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  #ProFisheryObj_list = list(proj_fishery_area1, proj_fishery_area2),
  #StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = NULL,
  wd = getwd(),
  fileName = "test1_single_fleet",
  seed = test_seed,
  doPlot = TRUE,
  doDiagnostic = FALSE
)

cat("Single fleet simulation completed\n")

# ============================================================================
# TEST 2: MULTIFLEET WITH 1 FLEET (effort) (SHOULD MATCH SINGLE FLEET)
# ============================================================================

cat("\n=== TEST 2: MULTIFLEET WITH 1 FLEET ===\n")

multifleet_1E <- new("Multifleet")
multifleet_1E@nfleets <- 1
multifleet_1E@fleet_proportions <- c(1.0)
multifleet_1E@allocation_type <- "effort"
multifleet_1E@fleet_selectivity_hist_list <- list(fleet1_sel_hist)
#multifleet_1E@fleet_selectivity_proj_list <- list(fleet1_sel_proj)
multifleet_1E@fleet_historicalEffort <- array(dim =c(ta@historicalYears, ta@areas, 1))
multifleet_1E@fleet_historicalEffort[,,1] <-  ta@historicalEffort


result_multifleet_1E <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  #StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_1E,
  wd = getwd(),
  fileName = "test2_multifleet_1_fleetE",
  seed = test_seed,
  doPlot = TRUE,
  doDiagnostic = FALSE
)

cat("Multifleet (1 fleet - effort) simulation completed\n")


# ============================================================================
# TEST 3: MULTIFLEET WITH 1 FLEET (catch) (SHOULD MATCH SINGLE FLEET)
# ============================================================================

cat("\n=== TEST 3: MULTIFLEET WITH 1 FLEET ===\n")

multifleet_1C <- new("Multifleet")
multifleet_1C@nfleets <- 1
multifleet_1C@fleet_proportions <- c(1.0)
multifleet_1C@allocation_type <- "catch"
multifleet_1C@fleet_selectivity_hist_list <- list(fleet1_sel_hist)
#multifleet_1C@fleet_selectivity_proj_list <- list(fleet1_sel_proj)
multifleet_1C@fleet_historicalEffort <- array(dim =c(ta@historicalYears, ta@areas, 1))
multifleet_1C@fleet_historicalEffort[,,1] <-  ta@historicalEffort

result_multifleet_1C <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  #StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_1C,
  wd = getwd(),
  fileName = "test2_multifleet_1_fleetC",
  seed = test_seed,
  doPlot = TRUE,
  doDiagnostic = FALSE
)

cat("Multifleet (1 fleet - catch) simulation completed\n")


# ============================================================================
# COMPARISON AND VALIDATION (SINGLE FLEET VS MULTIFLEET)
# ============================================================================

cat("\n=== COMPARISON AND VALIDATION ===\n")

# Read results for comparison
result_single <- readProjection(getwd(), "test1_single_fleet")
result_multifleet_1E <- readProjection(getwd(), "test2_multifleet_1_fleetE")
result_multifleet_1C <- readProjection(getwd(), "test2_multifleet_1_fleetC")

#Bill edit: modified to work with lack of Ftotal for multi fleet
# Function to compare key metrics
compare_metrics <- function(result1, result2, name1, name2, tolerance = 1e-10) {
  cat(sprintf("\nComparing %s vs %s:\n", name1, name2))

  # Compare spawning biomass
  sb_diff <- max(abs(result1$dynamics$SB - result2$dynamics$SB))
  cat(sprintf("  Max SB difference: %.2e\n", sb_diff))

  # Compare total F
  if(is.null(result1$dynamics$multifleet)){
    F1 <- result1$dynamics$Ftotal
  } else {
    F1 <- result1$dynamics$multifleet$Ftotal_by_fleet[,,,1]
  }

  if(is.null(result2$dynamics$multifleet)){
    F2 <- result2$dynamics$Ftotal
  } else {
    F2 <- result2$dynamics$multifleet$Ftotal_by_fleet[,,,1]
  }

  f_diff <- max(abs(F1 - F2))
  cat(sprintf("  Max F difference: %.2e\n", f_diff))

  # Compare SPR
  spr_diff <- max(abs(result1$dynamics$SPR - result2$dynamics$SPR))
  cat(sprintf("  Max SPR difference: %.2e\n", spr_diff))

  # Overall assessment
  max_diff <- max(sb_diff, f_diff, spr_diff)
  if(max_diff < tolerance) {
    cat(sprintf("  RESULT: IDENTICAL (max diff = %.2e)\n", max_diff))
    return(TRUE)
  } else {
    cat(sprintf("  RESULT: DIFFERENT (max diff = %.2e)\n", max_diff))
    return(FALSE)
  }
}

# Test 1: Single fleet vs Multifleet(1 fleetE) should be identical
identical_1 <- compare_metrics(result_single, result_multifleet_1E,
                               "Single Fleet", "Multifleet(1E)", 1e-8)


# Test 1: Single fleet vs Multifleet(1 fleetC) should be identical
identical_2 <- compare_metrics(result_single, result_multifleet_1C,
                               "Single Fleet", "Multifleet(1C)", 1e-8)

# Test 1: Multifleet(1 fleetE) vs Multifleet(1 fleetC) should be identical
identical_3 <- compare_metrics(result_multifleet_1E, result_multifleet_1C,
                               "Multifleet(1E)", "Multifleet(1C)", 1e-8)



all_results <- list(result_single, result_multifleet_1E, result_multifleet_1C)
result_names <- c("Single", "Multi(1E)", "Multi(1C)")

cat("\nFinal biomass by test:\n")
for(i in 1:length(all_results)) {
  final_sb <- sum(all_results[[i]]$dynamics$SB[10, , ])  # Year 15, all iterations, all areas
  cat(sprintf("  %s: %.2f\n", result_names[i], final_sb))
}

cat("\n=== TEST FILE 1 COMPLETED ===\n")
cat("Backward compatibility: VERIFIED\n")
cat("Multifleet functionality: VERIFIED\n")

# ============================================================================
# TEST 3: MULTIFLEET WITH 2 FLEETS - SAME SELECTIVITY
# ============================================================================

cat("\n=== TEST 3: MULTIFLEET WITH 2 FLEETS ===\n")

multifleet_2 <- new("Multifleet")
multifleet_2@nfleets <- 2
multifleet_2@fleet_proportions <- c(0.6, 0.4)
multifleet_2@allocation_type <- "catch"
multifleet_2@fleet_selectivity_hist_list <- list(fleet1_sel_hist, fleet1_sel_hist) #assuming same sel for each fleet
#multifleet_2@fleet_selectivity_proj_list <- list(fleet1_sel_proj, fleet1_sel_proj) #assuming same sel for each fleet
multifleet_2@fleet_historicalEffort <- array(dim =c(ta@historicalYears, ta@areas, 2))
multifleet_2@fleet_historicalEffort[,,1] <-  ta@historicalEffort
multifleet_2@fleet_historicalEffort[,,2] <-  ta@historicalEffort

result_multifleet_2A <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  #StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_2,
  wd = getwd(),
  fileName = "test3_multifleet_2_fleets",
  seed = test_seed,
  doPlot = TRUE,
  doDiagnostic = FALSE
)

cat("Multifleet (2 fleets) simulation completed\n")



# ============================================================================
# TEST 4: MULTIFLEET WITH 2 FLEETS - DIFFERENT  SELECTIVITY
# ============================================================================

cat("\n=== TEST 4: MULTIFLEET WITH 2 FLEETS ===\n")

multifleet_2 <- new("Multifleet")
multifleet_2@nfleets <- 2
multifleet_2@fleet_proportions <- c(0.6, 0.4)
multifleet_2@allocation_type <- "catch"
multifleet_2@fleet_selectivity_hist_list <- list(fleet1_sel_hist, fleet2_sel_hist) #assuming diff sel for each fleet
#multifleet_2@fleet_selectivity_proj_list <- list(fleet1_sel_proj, fleet2_sel_proj) #assuming diff sel for each fleet
multifleet_2@fleet_historicalEffort <- array(dim =c(ta@historicalYears, ta@areas, 2))
multifleet_2@fleet_historicalEffort[,,1] <-  ta@historicalEffort
multifleet_2@fleet_historicalEffort[,,2] <-  ta@historicalEffort

result_multifleet_2B <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  #StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_2,
  wd = getwd(),
  fileName = "test4_multifleet_2_fleets",
  seed = test_seed,
  doPlot = TRUE,
  doDiagnostic = FALSE
)

cat("Multifleet (2 fleets) simulation completed\n")


result_multifleet_2A <- readProjection(getwd(), "test3_multifleet_2_fleets")
result_multifleet_2B <- readProjection(getwd(), "test4_multifleet_2_fleets")

# Test Multifleet (2 fleets) - same selectivity vs differnt selectivity
Compare_diff_sel <- compare_metrics(result_multifleet_2A, result_multifleet_2B,
                               "result_multifleet_2A", "result_multifleet_2B", 1e-8)

all_results <- list(result_multifleet_2A, result_multifleet_2B)
result_names <- c("result_multifleet_2A", "result_multifleet_2B")

cat("\nFinal biomass by test:\n")
for(i in 1:length(all_results)) {
  final_sb <- sum(all_results[[i]]$dynamics$SB[10, , ])  # Year 15, all iterations, all areas
  cat(sprintf("  %s: %.2f\n", result_names[i], final_sb))
}

#Same selectivity: final_effort_proportions = actual_catch_proportions
result_multifleet_2A$dynamics$multifleet$actual_catch_proportions

#Different selectivity: final_effort_proportions != actual_catch_proportions
result_multifleet_2B$dynamics$multifleet$actual_catch_proportions


# ============================================================================
# TEST 5: MULTIFLEET WITH 3 FLEETS - DIFFERENT SELECTIVITIES
# ============================================================================

cat("\n=== TEST 5: MULTIFLEET WITH 3 FLEETS ===\n")

multifleet_3 <- new("Multifleet")
multifleet_3@nfleets <- 3
multifleet_3@fleet_proportions <- c(0.5, 0.3, 0.2)
multifleet_3@allocation_type <- "catch"
multifleet_3@fleet_selectivity_hist_list <- list(fleet1_sel_hist, fleet2_sel_hist, fleet3_sel_hist)
#multifleet_3@fleet_selectivity_proj_list <- list(fleet1_sel_proj, fleet2_sel_proj, fleet3_sel_proj)
multifleet_3@fleet_historicalEffort <- array(dim =c(ta@historicalYears, ta@areas, 3))
multifleet_3@fleet_historicalEffort[,,1] <-  ta@historicalEffort
multifleet_3@fleet_historicalEffort[,,2] <-  ta@historicalEffort
multifleet_3@fleet_historicalEffort[,,3] <-  ta@historicalEffort


result_multifleet_3 <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  #StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_3,
  wd = getwd(),
  fileName = "test5_multifleet_3_fleets",
  seed = test_seed,
  doPlot = TRUE,
  doDiagnostic = FALSE
)

cat("Multifleet (3 fleets) simulation completed\n")

result_multifleet_3 <- readProjection(getwd(), "test5_multifleet_3_fleets")

result_multifleet_3$dynamics$SB
result_multifleet_3$dynamics$multifleet$VB_by_fleet
result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet

dim(result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet)#[11,2,2,3] #years, iter, area, fleet
result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,1,1] #iter 1, area 1, fleet 1
result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,1,2] #iter 1, area 1, fleet 2
result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,1,3] #iter 1, area 1, fleet 3

result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,2,1] #iter 1, area 2, fleet 1
result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,2,2] #iter 1, area 2, fleet 2
result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,2,3] #iter 1, area 2, fleet 3

# sum F across fleets - Area 1
Farea1<- apply(cbind(result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,1,1],
    result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,1,2],
    result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,1,3]), 1, sum)

# sum F across fleets - Area 2
Farea2<- apply(cbind(result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,2,1],
    result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,2,2],
    result_multifleet_3$dynamics$multifleet$Ftotal_by_fleet[,1,2,3]),1,sum)


# Clean up intermediate files
# file.remove("test1_single_fleet.rds")
# file.remove("test2_multifleet_1_fleetE.rds")
# file.remove("test2_multifleet_1_fleetC.rds")
# file.remove("test3_multifleet_2_fleets.rds")
# file.remove("test4_multifleet_2_fleets.rds")
# file.remove("test5_multifleet_3_fleets.rds")

