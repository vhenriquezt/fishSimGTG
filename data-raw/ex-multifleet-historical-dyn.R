
# ============================================================================
# example multifleet historical dyn
# ============================================================================
# initial testings
rm(list=ls())
#options(error = traceback)
#devtools::install()
devtools::load_all()
# #devtools::document()
library(ggplot2)
#library(fishSimGTG)
library(dplyr)

# Create simple examples of each class to understand their structure
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
lh_obj@recSD<-0 #Run with no rec var'n to see deterministic trends
lh_obj@recRho<-0
lh_obj@isHermaph<-FALSE
lh_obj@R0<-10000



#Time area set up
# The TimeArea class controls the size and structure of the simulation
ta <- new("TimeArea")
ta@title = "Test"
ta@gtg = 13
ta@areas = 2
ta@recArea = c(0.99, 0.01)
ta@iterations = 2
ta@historicalYears = 10
ta@historicalBio = 0.5
ta@historicalBioType = "relB" # or SPR
ta@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
#Matrix of fishing effort multipliers
#Dimensions: [historicalYears × areas]
#Values are multipliers of equilibrium fishing rate
#e.g., 1.5 = 50% more fishing than equilibrium
# ta@historicalEffort <- matrix(c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6),
#                              nrow = 10, ncol = 2, byrow = FALSE)


ta@historicalEffort <- matrix(c(1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6,
                                1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6),
                              nrow = 10, ncol = 2, byrow = FALSE)


dim(ta@historicalEffort)



# #Fishery ste up
hist_fishery <- new("Fishery")
hist_fishery@title<-"Test"
hist_fishery@vulType<-"logistic"
hist_fishery@vulParams<-c(10.2,2)
hist_fishery@retType<-"full"
hist_fishery@retMax <- 1
hist_fishery@Dmort <- 0

proj_fishery_area1 <- new("Fishery")
proj_fishery_area1@title<-"Proj Area 1"
proj_fishery_area1@vulType<-"logistic"
proj_fishery_area1@vulParams<-c(10.2,0.1)
proj_fishery_area1@retType<-"logistic"
proj_fishery_area1@retParams <- c(10.2, 0.1)
proj_fishery_area1@retMax <- 1
proj_fishery_area1@Dmort <- 0

proj_fishery_area2 <- new("Fishery")
proj_fishery_area2@title<-"Proj Area 2"
proj_fishery_area2@vulType<-"logistic"
proj_fishery_area2@vulParams<-c(10.2,0.1)
proj_fishery_area2@retType<-"logistic"
proj_fishery_area2@retParams <- c(10.2, 0.1)
proj_fishery_area2@retMax <- 1
proj_fishery_area2@Dmort <- 0

proj_fishery_list <- list(proj_fishery_area1, proj_fishery_area2)  # FIXED: Create the list


#stochastic object set up
stochastic_obj <- new("Stochastic")
stochastic_obj@historicalBio = c(0.3, 0.6)
stochastic_obj@Steep= c(0.45, 0.75)


lh <- LHwrapper(LifeHistoryObj = lh_obj, TimeAreaObj = ta)

# Create multifleet object
fishery1 <- new("Fishery")
fishery1@vulType <- "logistic"
fishery1@vulParams <- c(9.0, 1.5)
fishery1@retType <- "full"
fishery1@retMax <- 1
fishery1@Dmort <- 0

fishery2 <- new("Fishery")
fishery2@vulType <- "logistic"
fishery2@vulParams <- c(11.0, 2.0)
fishery2@retType <- "full"
fishery2@retMax <- 1
fishery2@Dmort <- 0

multifleet_obj <- new("Multifleet")
multifleet_obj@nfleets <- 2
multifleet_obj@fleet_proportions <- c(0.6, 0.4)
multifleet_obj@allocation_type <- "catch"
multifleet_obj@fleet_selectivity_list <- list(fishery1, fishery2)



#adding a simple strategy (based on the template)
simpleMP <- function(phase, dataObject) {
  # Unpack dataObject
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase==1){
    # Phase 1: No observations needed for this simple MP
    return(list())
  }

  if(phase==2){
    # Phase 2: No complex analysis needed for this simple MP
    return(list())
  }

  if(phase==3){
    # Phase 3: Return constant F for all areas
    # Use a simple constant F = 0.15 for projection years

    year = rep(j, areas)
    iteration = rep(k, areas)
    area = 1:areas
    Flocal = rep(0.15, areas)  # Constant F = 0.15

    return(list(year=year, iteration=iteration, area=area, Flocal=Flocal))
  }
}


#simple strategy
strategy_obj <- new("Strategy")
strategy_obj@title <- "Simple Fixed F Strategy"
strategy_obj@projectionYears <- 5
strategy_obj@projectionName <- "simpleMP"  # use the management procedure for projections
strategy_obj@projectionParams <- list()


lh <- LHwrapper(LifeHistoryObj = lh_obj, TimeAreaObj = ta)
cat("Testing basic functionality...\n")


#test basic functionality
sel_test <- selWrapper(lh, ta, hist_fishery, doPlot = FALSE)

#test equilibrium
eq_test <- solveD(lh, sel_test, doFit = TRUE, D_type = "relB", D_in = 0.4)

# Test multifleet equilibrium using selectivity list
hist_sel_list <- lapply(1:multifleet_obj@nfleets, function(f) {
  selWrapper(lh, ta, FisheryObj = multifleet_obj@fleet_selectivity_list[[f]], doPlot = FALSE)
})

eq_multifleet <- solveD_multifleet2(
  lh = lh,
  sel_list = hist_sel_list,
  doFit = TRUE,
  D_type = "relB",
  D_in = 0.4,
  fleet_proportions = multifleet_obj@fleet_proportions,
  allocation_type = multifleet_obj@allocation_type
)
cat("Multifleet equilibrium works\n")
cat(sprintf("  - Total F: %.4f\n", eq_multifleet$Feq))
cat(sprintf("  - Fleet F: [%.4f, %.4f]\n", eq_multifleet$F_by_fleet[1], eq_multifleet$F_by_fleet[2]))
cat(sprintf("  - Achieved depletion: %.4f\n", eq_multifleet$D))


wd <- getwd()
cat("Working directory:", wd, "\n")

#testing single fleet
single_result <- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  HistFisheryObj = hist_fishery,
  ProFisheryObj_list = proj_fishery_list,
  StrategyObj = strategy_obj,
  StochasticObj = stochastic_obj,
  MultifleetObj = NULL,
  wd = wd,
  fileName = "single_fleet_test",
  seed = 123,
  doPlot = FALSE,
  doDiagnostic = FALSE
)

single_result<-readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG", "single_fleet_test")
single_result$dynamics$SB
single_result$dynamics$catchB
single_result$dynamics$catchN
single_result$HCR$decisionLocal
single_result$HCR$decisionData
single_result$dynamics$Ftotal





stop()
#
# #selectivity list
# hist_sel_list <- lapply(1:MultifleetObj@nfleets, function(f) {
#   selWrapper(lh, ta, FisheryObj = MultifleetObj@fleet_selectivity_list[[f]], doPlot = FALSE)
# })
#
#
# # ============================================================================
# # CREATE STRATEGY OBJECT
# # ============================================================================
#
# strategy_obj <- new("Strategy")
# strategy_obj@title <- "Fixed Strategy Test"
# strategy_obj@projectionYears <- 5
# strategy_obj@projectionName <- "fixedStrategy"
# strategy_obj@projectionParams <- list()  # No special parameters for fixed strategy
#
#
#
#
#
# # ============================================================================
# # RUN FULL SIMULATION - SINGLE FLEET
# # ============================================================================
#
# # Create working directory for outputs
# wd <- tempdir()
# cat("Working directory:", wd, "\n")
#
# # Run single fleet simulation
#
#   single_result <- runProjection(
#     LifeHistoryObj = lh_obj,
#     TimeAreaObj = ta_obj,
#     HistFisheryObj = hist_fishery,
#     ProFisheryObj_list = proj_fishery_list,
#     StrategyObj = strategy_obj,
#     StochasticObj = stochastic_obj,
#     MultifleetObj = NULL,  # Single fleet mode
#     wd = wd,
#     fileName = "single_fleet_test",
#     seed = 123,
#     doPlot = FALSE,
#     doDiagnostic = FALSE
#   )
#
#
# # ============================================================================
# # RUN FULL SIMULATION - MULTIFLEET
# # ============================================================================
#
#   multifleet_result <- runProjection(
#     LifeHistoryObj = lh_obj,
#     TimeAreaObj = ta_obj,
#     HistFisheryObj = hist_fishery,
#     ProFisheryObj_list = proj_fishery_list,
#     StrategyObj = strategy_obj,
#     StochasticObj = stochastic_obj,
#     MultifleetObj = multifleet_obj,  # Multifleet mode
#     wd = wd,
#     fileName = "multifleet_test",
#     seed = 123,
#     doPlot = FALSE,
#     doDiagnostic = FALSE
#   )
#
#
#
#
# # ============================================================================
# # TEST BACKWARD COMPATIBILITY
# # ============================================================================
#
# # Test that multifleet with 1 fleet gives same results as single fleet
# single_fleet_multi <- new("Multifleet")
# single_fleet_multi@nfleets <- 1
# single_fleet_multi@fleet_proportions <- c(1.0)
# single_fleet_multi@allocation_type <- "effort"
# single_fleet_multi@fleet_selectivity_list <- list(hist_fishery)
#
#
#   compat_result <- runProjection(
#     LifeHistoryObj = lh_obj,
#     TimeAreaObj = ta_obj,
#     HistFisheryObj = hist_fishery,
#     ProFisheryObj_list = proj_fishery_list,
#     StrategyObj = strategy_obj,
#     StochasticObj = stochastic_obj,
#     MultifleetObj = single_fleet_multi,  # Single fleet via multifleet
#     wd = wd,
#     fileName = "compatibility_test",
#     seed = 123,
#     doPlot = FALSE,
#     doDiagnostic = FALSE
#   )
#
#   compat_final_ssb <- compat_result$dynamics$SB[nrow(compat_result$dynamics$SB), 1, 1]
#
#
