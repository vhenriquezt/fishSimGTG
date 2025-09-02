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

# Strategy for multifleet (with fleet column) to maintain maintain data structure consistency with
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

result_single <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_single_fleet")
result_single$dynamics$SB # 1 set of vectors for each area for n iterations
result_single$dynamics$VB
result_single$dynamics$RB
result_single$dynamics$catchB
result_single$dynamics$catchN
result_single$dynamics$Ftotal # 1 set of vectors for each area for n iterations
result_single$dynamics$SPR # 1 set of vectors area-combined



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
result_multi2$dynamics$RB
result_multi2$dynamics$catchB
result_multi2$dynamics$catchN
result_multi2$dynamics$Ftotal
result_multi2$dynamics$SPR
result_multi2$dynamics$multifleet$Ftotal_by_fleet
result_multi2$dynamics$multifleet$catchB_by_fleet
result_multi2$dynamics$multifleet$final_effort_proportions
result_multi2$dynamics$multifleet$actual_catch_proportions
result_multi2$dynamics$multifleet$target_catch_proportions

# ============================================================================
# TEST SCENARIO 3: MULTIFLEET WITH 1 FLEET (SHOULD MATCH SINGLE FLEET)
# ============================================================================

cat("======================================\n")
cat("RUNNING TEST 3: MULTIFLEET (1 FLEET)\n")
cat("======================================\n")

fleet1_single_hist <- new("Fishery")
fleet1_single_hist@vulType <- "logistic"
fleet1_single_hist@vulParams <- c(10.2, 0.1)
fleet1_single_hist@retType <- "full"
fleet1_single_hist@retMax <- 1
fleet1_single_hist@Dmort <- 0

fleet1_single_proj <- new("Fishery")
fleet1_single_proj@vulType <- "logistic"
fleet1_single_proj@vulParams <- c(10.2, 0.1)
fleet1_single_proj@retType <- "full"
fleet1_single_proj@retMax <- 1
fleet1_single_proj@Dmort <- 0



multifleet_single <- new("Multifleet")
multifleet_single@nfleets <- 1
multifleet_single@fleet_proportions <- c(1.0)  # 100% to single fleet
multifleet_single@allocation_type <- "catch"
multifleet_single@fleet_selectivity_hist_list  <- list(fleet1_single_hist)
multifleet_single@fleet_selectivity_proj_list  <- list(fleet1_single_proj)



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

result_multi1 <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_multifleet_1fleet")

result_multi1$dynamics$SB
result_multi1$dynamics$VB
result_multi1$dynamics$RB
result_multi1$dynamics$catchB
result_multi1$dynamics$catchN
result_multi1$dynamics$Ftotal
result_multi1$dynamics$SPR
result_multi1$dynamics$multifleet$Ftotal_by_fleet
result_multi1$dynamics$multifleet$catchB_by_fleet


# ============================================================================
# COMPARING THE THREE PREVIOUS EXAMPLES
# ============================================================================

cat("=========================================================================\n")
cat("COMPARING: single fleet - multifleet (2 identical)- multifleet (1 fleet))\n")
cat("=========================================================================\n")

#Single Fleet vs Multifleet (1 fleet) - must be identical
identical(result_single$dynamics$SB, result_multi1$dynamics$SB)
identical(result_single$dynamics$VB, result_multi1$dynamics$VB)
identical(result_single$dynamics$catchB, result_multi1$dynamics$catchB)
identical(result_single$dynamics$Ftotal, result_multi1$dynamics$Ftotal)


# Check 1-fleet multifleet has correct fleet dimensions
dim(result_multi1$dynamics$multifleet$Ftotal_by_fleet)     # [years, iter, areas, 1]
dim(result_multi1$dynamics$multifleet$catchB_by_fleet)     # [years, iter, areas, 1]

# Fleet values should equal total values for 1-fleet case
identical(result_multi1$dynamics$Ftotal, result_multi1$dynamics$multifleet$Ftotal_by_fleet[,,,1])
identical(result_multi1$dynamics$catchB, result_multi1$dynamics$multifleet$catchB_by_fleet[,,,1])


#Multifleet (2 identical) fleet summation

# Fleet catches must sum to total
fleet_sum_catchB <- result_multi2$dynamics$multifleet$catchB_by_fleet[,,,1] +
  result_multi2$dynamics$multifleet$catchB_by_fleet[,,,2]


max(abs(fleet_sum_catchB - result_multi2$dynamics$catchB))


# Fleet F must sum to total
fleet_sum_F <- result_multi2$dynamics$multifleet$Ftotal_by_fleet[,,,1] +
  result_multi2$dynamics$multifleet$Ftotal_by_fleet[,,,2]

max(abs(fleet_sum_F - result_multi2$dynamics$Ftotal))  # Should be ~0

# Since selectivities are identical, population should be very similar
# (Small differences expected due to allocation algorithm)

# Check relative differences
rel_diff_SB <- abs(result_single$dynamics$SB - result_multi2$dynamics$SB) /
  result_single$dynamics$SB

max(rel_diff_SB, na.rm = TRUE)  # Should be < 5%

# Same for catches
rel_diff_catch <- abs(result_single$dynamics$catchB - result_multi2$dynamics$catchB) /
  result_single$dynamics$catchB

max(rel_diff_catch, na.rm = TRUE)  # Should be < 5%

# Check multifleet allocation worked
result_multi2$dynamics$multifleet$allocation_type
result_multi2$dynamics$multifleet$fleet_proportions          # [0.6, 0.4]
result_multi2$dynamics$multifleet$target_catch_proportions   # [0.6, 0.4]
result_multi2$dynamics$multifleet$actual_catch_proportions   # Should be close to [0.6, 0.4]
result_multi2$dynamics$multifleet$final_effort_proportions

# Check allocation accuracy
allocation_error <- abs(result_multi2$dynamics$multifleet$actual_catch_proportions -
                          result_multi2$dynamics$multifleet$target_catch_proportions)
max(allocation_error)  # Should be < 0.05 (5% tolerance)



  # Test 1: Single vs Multi1 (should be identical)
  cat("\n1. BACKWARD COMPATIBILITY (Single vs Multi1):\n")
  sb_identical <- identical(result_single$dynamics$SB, result_multi1$dynamics$SB)
  vb_identical <- identical(result_single$dynamics$VB, result_multi1$dynamics$VB)
  rb_identical <- identical(result_single$dynamics$RB, result_multi1$dynamics$RB)

  cat(sprintf("   SB identical: %s\n", sb_identical))
  cat(sprintf("   VB identical: %s\n", vb_identical))
  cat(sprintf("   RB identical: %s\n", rb_identical))

  if(!sb_identical) {
    max_sb_diff <- max(abs(result_single$dynamics$SB - result_multi1$dynamics$SB))
    cat(sprintf("   Max SB difference: %e\n", max_sb_diff))
  }

  # Test 2: Population similarity (Single vs Multi2)
  years <- dim(result_single$dynamics$SB)[1]
  areas <- dim(result_single$dynamics$SB)[3]
  final_year <- years

  cat("\n2. POPULATION SIMILARITY (Single vs Multi2):\n")
  for(area in 1:areas) {
    # Compare final year values (iteration 1)
    final_year <- years

    sb_single <- result_single$dynamics$SB[final_year, 1, area]
    sb_multi2 <- result_multi2$dynamics$SB[final_year, 1, area]
    sb_diff_pct <- abs(sb_single - sb_multi2) / sb_single * 100

    vb_single <- result_single$dynamics$VB[final_year, 1, area]
    vb_multi2 <- result_multi2$dynamics$VB[final_year, 1, area]
    vb_diff_pct <- abs(vb_single - vb_multi2) / vb_single * 100

    rb_single <- result_single$dynamics$RB[final_year, 1, area]
    rb_multi2 <- result_multi2$dynamics$RB[final_year, 1, area]
    rb_diff_pct <- abs(rb_single - rb_multi2) / rb_single * 100

    #print results for each area
    cat(sprintf("  Area %d:\n", area))
    cat(sprintf("    SB difference: %.2f%%\n", sb_diff_pct))
    cat(sprintf("    VB difference: %.2f%%\n", vb_diff_pct))
    cat(sprintf("    RB difference: %.2f%%\n", rb_diff_pct))
  }

  # Store results in a variable instead of return (since this isn't inside a function)
  comparison_results <- list(
    backward_compatible = sb_identical && vb_identical && rb_identical,
    max_differences = list(
      sb_single_multi2 = max(abs(result_single$dynamics$SB - result_multi2$dynamics$SB)),
      vb_single_multi2 = max(abs(result_single$dynamics$VB - result_multi2$dynamics$VB)),
      rb_single_multi2 = max(abs(result_single$dynamics$RB - result_multi2$dynamics$RB))
    )
  )

  print(comparison_results)


  #create some plots

create_biomass_plots <- function() {

  years <- dim(result_single$dynamics$SB)[1]
  areas <- dim(result_single$dynamics$SB)[3]
  year_seq <- 1:years

  #create data frame for plotting (using iteration 1, area 1)
  plot_data <- data.frame(
    Year = rep(year_seq, 9),  # 3 scenarios × 3 metrics
    Value = c(
      # SB values
      result_single$dynamics$SB[, 1, 1],
      result_multi1$dynamics$SB[, 1, 1],
      result_multi2$dynamics$SB[, 1, 1],
      # VB values
      result_single$dynamics$VB[, 1, 1],
      result_multi1$dynamics$VB[, 1, 1],
      result_multi2$dynamics$VB[, 1, 1],
      # RB values
      result_single$dynamics$RB[, 1, 1],
      result_multi1$dynamics$RB[, 1, 1],
      result_multi2$dynamics$RB[, 1, 1]
    ),
    Scenario = rep(rep(c("Single_Fleet", "Multi_1Fleet", "Multi_2Identical"), each = years), 3),
    Metric = rep(c("Spawning_Biomass", "Vulnerable_Biomass", "Retained_Biomass"), each = years * 3)
  )

  #add period classification
  plot_data <- plot_data %>%
    mutate(
      Period = case_when(
        Year == 1 ~ "Equilibrium",
        Year <= 11 ~ "Historical",
        Year > 11 ~ "Projection"
      )
    )

  #create the plot
  p1 <- ggplot(plot_data, aes(x = Year, y = Value, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    geom_vline(xintercept = 1.5, linetype = "dotted", alpha = 0.7, color = "gray") +
    geom_vline(xintercept = 11.5, linetype = "dashed", alpha = 0.7, color = "red") +
    facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    labs(
      title = "Biomass Metrics Comparison - Area 1, Iteration 1",
      subtitle = "Dotted line = End Equilibrium | Dashed line = End Historical",
      x = "Year",
      y = "Biomass",
      color = "Scenario",
      linetype = "Scenario"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 10, face = "bold")
    )

  #create difference plot (Multi2 - Single)
  diff_data <- data.frame(
    Year = rep(year_seq, 3),
    Difference = c(
      result_multi2$dynamics$SB[, 1, 1] - result_single$dynamics$SB[, 1, 1],
      result_multi2$dynamics$VB[, 1, 1] - result_single$dynamics$VB[, 1, 1],
      result_multi2$dynamics$RB[, 1, 1] - result_single$dynamics$RB[, 1, 1]
    ),
    Metric = rep(c("SB_Difference", "VB_Difference", "RB_Difference"), each = years)
  )

  p2 <- ggplot(diff_data, aes(x = Year, y = Difference)) +
    geom_line(color = "blue", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 11.5, linetype = "dashed", alpha = 0.7, color = "gray") +
    facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    labs(
      title = "Difference: Multi2 - Single Fleet",
      subtitle = "Should be small differences (identical selectivities)",
      x = "Year",
      y = "Difference"
    ) +
    theme_minimal()

  return(list(comparison = p1, differences = p2))
}

#run the analysis and create plots
plots <- create_biomass_plots()

#display plots
print(plots$comparison)
print(plots$differences)



# ============================================================================
# TEST SCENARIO 4: MULTIFLEET WITH DIFFERENT SELECTIVITIES
# ============================================================================

cat("============================================\n")
cat("RUNNING TEST 4: MULTIFLEET (DIFFERENT SELECT)\n")
cat("============================================\n")

fleet1_diff_hist <- new("Fishery")
fleet1_diff_hist@vulType <- "logistic"
fleet1_diff_hist@vulParams <- c(8.0, 1.2)
fleet1_diff_hist@retType <- "full"
fleet1_diff_hist@retMax <- 1
fleet1_diff_hist@Dmort <- 0

fleet2_diff_hist <- new("Fishery")
fleet2_diff_hist@vulType <- "logistic"
fleet2_diff_hist@vulParams <- c(11.0, 2.0)
fleet2_diff_hist@retType <- "full"
fleet2_diff_hist@retMax <- 1
fleet2_diff_hist@Dmort <- 0

fleet1_diff_proj <- new("Fishery")
fleet1_diff_proj@vulType <- "logistic"
fleet1_diff_proj@vulParams <- c(9.0, 1.2)
fleet1_diff_proj@retType <- "full"
fleet1_diff_proj@retMax <- 1
fleet1_diff_proj@Dmort <- 0

fleet2_diff_proj <- new("Fishery")
fleet2_diff_proj@vulType <- "logistic"
fleet2_diff_proj@vulParams <- c(12.0, 2.0)
fleet2_diff_proj@retType <- "full"
fleet2_diff_proj@retMax <- 1
fleet2_diff_proj@Dmort <- 0




multifleet_different <- new("Multifleet")
multifleet_different@nfleets <- 2
multifleet_different@fleet_proportions <- c(0.7, 0.3)
multifleet_different@allocation_type <- "catch"  # Use catch allocation
multifleet_different@fleet_selectivity_hist_list  <- list(fleet1_diff_hist, fleet2_diff_hist)
multifleet_different@fleet_selectivity_proj_list  <- list(fleet1_diff_proj, fleet2_diff_proj)




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

result_multi_diff <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_multifleet_different")
result_multi_diff$dynamics$multifleet$final_effort_proportions
result_multi_diff$dynamics$multifleet$actual_catch_proportions

# ============================================================================
# VALIDATION ANALYSIS
# ============================================================================

# load all
result_single <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_single_fleet")
result_multi2 <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_multifleet_2identical")
result_multi1 <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_multifleet_1fleet")
result_multi_diff <- readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "validation_multifleet_different")


result_multi_diff$dynamics$multifleet




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











