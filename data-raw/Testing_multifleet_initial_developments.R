#testing_multifleet_ init_conditions
rm(list=ls())
devtools::load_all()
library(ggplot2)
library(dplyr)

# Create life history
lh_obj <- new("LifeHistory")
lh_obj@title <- "Kole"
lh_obj@speciesName <- "Ctenochaetus strigosus"
lh_obj@Linf <- 17.7
lh_obj@K <- 0.423
lh_obj@t0 <- -0.51
lh_obj@L50 <- 8.4
lh_obj@L95delta <- 1.26
lh_obj@M <- 0.08
lh_obj@L_type <- "FL"
lh_obj@L_units <- "cm"
lh_obj@LW_A <- 0.046
lh_obj@LW_B <- 2.85
lh_obj@Steep <- 0.54
lh_obj@recSD <- 0
lh_obj@recRho <- 0
lh_obj@isHermaph <- FALSE
lh_obj@R0 <- 10000

# Create TimeArea
ta <- new("TimeArea")
ta@title = "Test"
ta@gtg = 13

# Create fishery
fishery1 <- new("Fishery")
fishery1@title <- "Test Fishery1"
fishery1@vulType <- "logistic"
fishery1@vulParams <- c(10.2, 2)
fishery1@retType <- "full"
fishery1@retMax <- 1
fishery1@Dmort <- 0

fishery2 <- new("Fishery")
fishery2@title <- "Test Fishery2"
fishery2@vulType <- "logistic"
fishery2@vulParams <- c(11.2, 2)
fishery2@retType <- "full"
fishery2@retMax <- 1
fishery2@Dmort <- 0

fishery3 <- new("Fishery")
fishery3@title <- "Test Fishery3"
fishery3@vulType <- "logistic"
fishery3@vulParams <- c(9.2, 2)
fishery3@retType <- "full"
fishery3@retMax <- 1
fishery3@Dmort <- 0

# Create life history wrapper
lh <- LHwrapper(lh_obj, ta)

# Create selectivity wrapper
sel1 <- selWrapper(lh, ta, fishery1, doPlot = FALSE)
sel2 <- selWrapper(lh, ta, fishery2, doPlot = FALSE)
sel3 <- selWrapper(lh, ta, fishery3, doPlot = FALSE)

#=================================================================================================#
# Example 1: Run solveD (single fleet) doFit = FALSE- Fixed F scenarios (F = 0.2)
#=================================================================================================#
base_single1 <- solveD(lh, sel1, doFit = FALSE, F_in = 0.2)
base_single1$Feq
base_single1$D
base_single1$B0
base_single1$SB
base_single1$VB


# Run solveD_multifleet (two identical fleets, 50:50), doFit = FALSE - Fixed F scenarios (F = 0.2)
result_multi1 <- solveD_multifleet(lh, list(sel1, sel1), doFit = FALSE, F_in = 0.2,
                                  fleet_proportions = c(0.5, 0.5))
result_multi1$Feq
result_multi1$fleet_proportions
result_multi1$F_by_fleet
result_multi1$D
result_multi1$B0
result_multi1$SB
result_multi1$VB

# # Compare the results
# tol <- 1e-10
# feq_match1 <- abs(result_single1$Feq - result_multi1$Feq) < tol
# d_match1 <- abs(result_single1$D - result_multi1$D) < tol
# catch_match1 <- abs(result_single1$catchB - result_multi1$catchB) < tol
# vb_match1 <- abs(result_single1$VB - result_multi1$VB) < tol
#
# #other checks
# fleet_f_correct1 <- abs(result_multi1$F_by_fleet[1] - result_single1$Feq/2) < tol &&
#   abs(result_multi1$F_by_fleet[2] - result_single1$Feq/2) < tol
#
# fleet_catch_sum1 <- abs(sum(result_multi1$catchB_by_fleet) - result_multi1$catchB) < tol
#
# test1_pass <- feq_match1 && d_match1 && catch_match1 && vb_match1 && fleet_f_correct1 && fleet_catch_sum1


#============================================================================================================#
# Example 2: Estimating F eq (finding the F that produce X depletion - doFit = TRUE)
#============================================================================================================#
base_single2 <- solveD(lh, sel1, doFit = TRUE, D_type = "relB", D_in = 0.4)
base_single2$Feq
base_single2$D
base_single2$SPR
base_single2$B0
base_single2$Req
base_single2$VB
base_single2$SB

# Run solveD_multifleet (two identical fleets, 50:50), doFit = TRUE
result_multi2 <- solveD_multifleet(lh, list(sel1, sel1), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                   fleet_proportions = c(0.5, 0.5))
result_multi2$Feq
result_multi2$D
result_multi2$SPR
result_multi2$B0
result_multi2$Req

result_multi2$nfleets
result_multi2$F_by_fleet
result_multi2$YPR_by_fleet
result_multi2$VB
result_multi2$SB

# # Simple comparison
# feq_match2 <- abs(result_single2$Feq - result_multi2$Feq) < tol
# d_match2 <- abs(result_single2$D - result_multi2$D) < tol
# test2_pass <- feq_match2 && d_match2
#===================================================================================================#

# Example 3: testing with different F proportions (same selectivity)
# Comparing the output with result_single1
# Not estimating Feq

base_single3 <- solveD(lh, sel1, doFit = FALSE, F_in = 0.2)
base_single3$Feq
base_single3$D
base_single3$B0
base_single3$SB
base_single3$VB

result_multi3 <- solveD_multifleet(lh, list(sel1, sel1), doFit = FALSE, F_in = 0.2,
                                  fleet_proportions = c(0.7, 0.3))
result_multi3$Feq
result_multi3$fleet_proportions
result_multi3$F_by_fleet
result_multi3$D
result_multi3$B0
result_multi3$SB
result_multi3$VB

# # Check if results match
# tol <- 1e-10
# feq_match3 <- abs(result_single1$Feq - result_multi3$Feq) < tol
# d_match3 <- abs(result_single1$D - result_multi3$D) < tol
# catch_match3 <- abs(result_single1$catchB - result_multi3$catchB) < tol
# vb_match3 <- abs(result_single1$VB - result_multi3$VB) < tol

# #other checks
# fleet_catch_sum3 <- abs(sum(result_multi3$catchB_by_fleet) - result_multi3$catchB) < tol
#
# test3_pass <- feq_match3 && d_match3 && catch_match3 && vb_match3 && fleet_catch_sum3
#============================================================================================#

# Example 4: testing with different F proportions (same selectivity)
# NOW testing fit to depletion
# Not estimating Feq

base_single4 <- solveD(lh, sel1, doFit = TRUE, D_type = "relB", D_in = 0.4)
base_single4$Feq
base_single4$D
base_single4$SPR
base_single4$B0
base_single4$Req
base_single4$VB
base_single4$SB

result_multi4 <- solveD_multifleet(lh, list(sel1, sel1), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                   fleet_proportions = c(0.7, 0.3))
result_multi4$Feq
result_multi4$D
result_multi4$SPR
result_multi4$B0
result_multi4$Req

result_multi4$nfleets
result_multi4$F_by_fleet
result_multi4$YPR_by_fleet
#result_multi4$sel_list
result_multi4$VB
result_multi4$SB

#===================================================================================================#

# Example 5: testing with 50:50 F proportions and different selectivity)
# Create fishery

result_multi5 <- solveD_multifleet(lh, list(sel1, sel2), doFit = FALSE, F_in = 0.2,
                                   fleet_proportions = c(0.5, 0.5))
result_multi5$Feq
result_multi5$fleet_proportions
result_multi5$F_by_fleet
result_multi5$D
result_multi5$B0
result_multi5$SB
result_multi5$VB
#====================================================================================#
# Example 6: testing with 50:50 F proportions and different selectivity)
# fitting

result_multi6 <- solveD_multifleet(lh, list(sel1, sel2), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                   fleet_proportions = c(0.5, 0.5))
result_multi6$Feq
result_multi6$fleet_proportions
result_multi6$F_by_fleet
result_multi6$D
result_multi6$B0
result_multi6$SB
result_multi6$VB

#====================================================================================#
# Example 7: testing with 30:70 F proportions and different selectivity)
# fitting

result_multi7 <- solveD_multifleet(lh, list(sel1, sel2), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                   fleet_proportions = c(0.3, 0.7))
result_multi7$Feq
result_multi7$fleet_proportions
result_multi7$F_by_fleet
result_multi7$D
result_multi7$B0
result_multi7$SB
result_multi7$VB

#====================================================================================#
# Example 8: testing with 30:70 F proportions and different selectivity)
# fitting
# adding more fleets

result_multi8 <- solveD_multifleet(lh, list(sel1, sel2, sel3), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                   fleet_proportions = c(0.3, 0.5,0.2))
result_multi8$Feq
result_multi8$fleet_proportions
result_multi8$F_by_fleet
result_multi8$D
result_multi8$B0
result_multi8$SB
result_multi8$VB



# Comparing the results
scenarios <- list(
  base_single1 = base_single1,
  result_multi1 = result_multi1,
  base_single2 = base_single2,
  result_multi2 = result_multi2,
  base_single3 = base_single3,
  result_multi3 = result_multi3,
  base_single4 = base_single4,
  result_multi4 = result_multi4,
  result_multi5 = result_multi5,
  result_multi6 = result_multi6,
  result_multi7 = result_multi7,
  result_multi8 = result_multi8
)

param_names <- c("Feq", "D", "B0", "SB", "VB", "SPR", "Req", "fleet_proportions", "F_by_fleet")

# Function to extract parameter, returns NA if not present
extract_param <- function(obj, pname, digits = 3) {
  out <- tryCatch(obj[[pname]], error = function(e) NA)
  if(is.null(out)) return(NA)
  # If it's numeric and a vector, round it before collapsing
  if(is.numeric(out) && length(out) > 1) {
    out <- round(out, digits = digits)
    paste(out, collapse = ", ")
  } else if(is.numeric(out) && length(out) == 1) {
    round(out, digits = digits)
  } else {
    out  # e.g. NA, or character
  }
}

summary_table <- sapply(scenarios, function(obj) {
  sapply(param_names, extract_param, obj=obj, digits=3)
})

summary_df <- as.data.frame(summary_table)
rownames(summary_df) <- param_names

print(summary_df)
write.csv(summary_df, "multifleet_solveD_examaples.csv")


##### MAIN CHANGES (SUMMARY)

# Explanation main modifications to solveD() to create solveD_multifleet

#CHANGE 1:
# BEFORE (solveD)
solveD(lh, sel, ...)                    # single selectivity object

# AFTER (solveD_multifleet)
solveD_multifleet(lh, sel_list, ...)    # list of selectivity objects


#CHANGE 2:
# adding a new fleet setup and proportions
nfleets <- length(sel_list)
fleet_proportions <- rep(1/nfleets, nfleets)  # how to distribute total F
F_by_fleet <- Ftotal * fleet_proportions      # individual fleet F values


#CHANGE 3:  Total Mortality Calculation
# BEFORE (single fleet):
Z = M + F * selectivity

# AFTER (multifleet - shared mortality):
total_removal_sel <- sapply(1:totalSteps, function(age) {
  sum(sapply(1:nfleets, function(f) {
    F_by_fleet[f] * sel_list[[f]]$removal[[x]][age]  # Sum all fleet contributions
  }))
})
Z = M + total_removal_sel

# CHANGE 4: Fleet-Specific Calculations
# Calculate catches, biomass, etc. for each fleet individually

#CHANGE 5: F-weighted combined vulnerability - used only to provide a metrics of total VB
combined_vul_by_gtg <- lapply(1:lh$gtg, function(x) {
  combined_vul <- rep(0, length(sel_list[[1]]$vul[[x]]))
  for(f in 1:nfleets) {
    combined_vul <- combined_vul + (F_by_fleet[f] / Feq_total) * sel_list[[f]]$vul[[x]]
  }
  combined_vul
})

