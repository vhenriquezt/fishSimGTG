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
base11 <- solveD(lh, sel1, doFit = FALSE, F_in = 0.2)

multi1 <- solveD_multifleet(lh, list(sel1, sel1), doFit = FALSE, F_in = 0.2,
                                 fleet_proportions = c(0.5, 0.5))

iter_eff_1 <- solveD_multifleet2(lh, list(sel1, sel1), doFit = FALSE, F_in = 0.2,
                                fleet_proportions = c(0.5, 0.5),
                                allocation_type="effort")

iter_ct_1 <- solveD_multifleet2(lh, list(sel1, sel1), doFit = FALSE, F_in = 0.2,
                                 fleet_proportions = c(0.5, 0.5),
                                 allocation_type="catch")

comparison_df1 <- data.frame(
  Single_Fleet1 = round(c(base11$Feq, base11$D, base11$SPR, base11$SB, base11$catchB, base11$YPR), 3),
  Multi_Original1 = round(c(multi1$Feq, multi1$D, multi1$SPR, multi1$SB, multi1$catchB, multi1$YPR), 3),
  Multi2_Effort1 = round(c(iter_eff_1$Feq, iter_eff_1$D, iter_eff_1$SPR, iter_eff_1$SB, iter_eff_1$catchB, iter_eff_1$YPR), 3),
  Multi2_Catch1 = round(c(iter_ct_1$Feq, iter_ct_1$D, iter_ct_1$SPR, iter_ct_1$SB, iter_ct_1$catchB, iter_ct_1$YPR), 3)
)

rownames(comparison_df1) <- c("Feq", "D", "SPR", "SB", "catchB", "YPR")

print(comparison_df1)

#============================================================================================================#
# Example 2: Estimating F eq (finding the F that produce X depletion - doFit = TRUE)
#============================================================================================================#
base2 <- solveD(lh, sel1, doFit = TRUE, D_type = "relB", D_in = 0.4)

multi2 <- solveD_multifleet(lh, list(sel1, sel1), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                   fleet_proportions = c(0.5, 0.5))


iter_eff_2 <- solveD_multifleet2(lh, list(sel1, sel1), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                   fleet_proportions = c(0.5, 0.5),
                                   allocation_type="effort")

iter_ct_2 <- solveD_multifleet2(lh, list(sel1, sel1), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                   fleet_proportions = c(0.5, 0.5),
                                   allocation_type="catch")

comparison_df2 <- data.frame(
  Single_Fleet2 = round(c(base2$Feq, base2$D, base2$SPR, base2$SB, base2$catchB, base2$YPR), 3),
  Multi_Original2 = round(c(multi2$Feq, multi2$D, multi2$SPR, multi2$SB, multi2$catchB, multi2$YPR), 3),
  Multi2_Effort2 = round(c(iter_eff_2$Feq, iter_eff_2$D, iter_eff_2$SPR, iter_eff_2$SB, iter_eff_2$catchB, iter_eff_2$YPR), 3),
  Multi2_Catch2 = round(c(iter_ct_2$Feq, iter_ct_2$D, iter_ct_2$SPR, iter_ct_2$SB, iter_ct_2$catchB, iter_ct_2$YPR), 3)
)

rownames(comparison_df2) <- c("Feq", "D", "SPR", "SB", "catchB", "YPR")

print(comparison_df2)


#============================================================================================================#
# Example 3: Estimating F eq (finding the F that produce X depletion - doFit = TRUE)- different proportions
#============================================================================================================#
base3 <- solveD(lh, sel1, doFit = TRUE, D_type = "relB", D_in = 0.4)

multi3 <- solveD_multifleet(lh, list(sel1, sel1), doFit = TRUE, D_type = "relB", D_in = 0.4,
                            fleet_proportions = c(0.2, 0.8))


iter_eff_3 <- solveD_multifleet2(lh, list(sel1, sel1), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                 fleet_proportions = c(0.2, 0.8),
                                 allocation_type="effort")

iter_ct_3 <- solveD_multifleet2(lh, list(sel1, sel1), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                fleet_proportions = c(0.2, 0.8),
                                allocation_type="catch")

comparison_df3 <- data.frame(
  Single_Fleet3 = round(c(base3$Feq, base3$D, base3$SPR, base3$SB, base3$catchB, base3$YPR), 3),
  Multi_Original3 = round(c(multi3$Feq, multi3$D, multi3$SPR, multi3$SB, multi3$catchB, multi3$YPR), 3),
  Multi2_Effort3 = round(c(iter_eff_3$Feq, iter_eff_3$D, iter_eff_3$SPR, iter_eff_3$SB, iter_eff_3$catchB, iter_eff_3$YPR), 3),
  Multi2_Catch3 = round(c(iter_ct_3$Feq, iter_ct_3$D, iter_ct_3$SPR, iter_ct_3$SB, iter_ct_3$catchB, iter_ct_3$YPR), 3)
)

rownames(comparison_df3) <- c("Feq", "D", "SPR", "SB", "catchB", "YPR")

print(comparison_df3)


#============================================================================================================#
# Example 4: Estimating F eq (finding the F that produce X depletion - doFit = TRUE)- different selectivities
#============================================================================================================#
base4 <- solveD(lh, sel1, doFit = TRUE, D_type = "relB", D_in = 0.4)

multi4 <- solveD_multifleet(lh, list(sel1, sel2), doFit = TRUE, D_type = "relB", D_in = 0.4,
                            fleet_proportions = c(0.5, 0.5))


iter_eff_4 <- solveD_multifleet2(lh, list(sel1, sel2), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                 fleet_proportions = c(0.5, 0.5),
                                 allocation_type="effort")

iter_ct_4 <- solveD_multifleet2(lh, list(sel1, sel2), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                fleet_proportions = c(0.5, 0.5),
                                allocation_type="catch")

comparison_df4 <- data.frame(
  Single_Fleet4 = round(c(base4$Feq, base4$D, base4$SPR, base4$SB, base4$catchB, base4$YPR), 3),
  Multi_Original4 = round(c(multi4$Feq, multi4$D, multi4$SPR, multi4$SB, multi4$catchB, multi4$YPR), 3),
  Multi2_Effort4 = round(c(iter_eff_4$Feq, iter_eff_4$D, iter_eff_4$SPR, iter_eff_4$SB, iter_eff_4$catchB, iter_eff_4$YPR), 3),
  Multi2_Catch4 = round(c(iter_ct_4$Feq, iter_ct_4$D, iter_ct_4$SPR, iter_ct_4$SB, iter_ct_4$catchB, iter_ct_4$YPR), 3)
)

rownames(comparison_df4) <- c("Feq", "D", "SPR", "SB", "catchB", "YPR")

print(comparison_df4)


#============================================================================================================#
# Example 5: Estimating F eq (finding the F that produce X depletion - doFit = TRUE)- different props and different selectivities
#============================================================================================================#
base5 <- solveD(lh, sel1, doFit = TRUE, D_type = "relB", D_in = 0.4)

multi5 <- solveD_multifleet(lh, list(sel1, sel2), doFit = TRUE, D_type = "relB", D_in = 0.4,
                            fleet_proportions = c(0.7, 0.3))


iter_eff_5 <- solveD_multifleet2(lh, list(sel1, sel2), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                 fleet_proportions = c(0.7, 0.3),
                                 allocation_type="effort")

iter_ct_5 <- solveD_multifleet2(lh, list(sel1, sel2), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                fleet_proportions = c(0.7, 0.3),
                                allocation_type="catch")

comparison_df5 <- data.frame(
  Single_Fleet5 = round(c(base5$Feq, base5$D, base5$SPR, base5$SB, base5$catchB, base5$YPR), 3),
  Multi_Original5 = round(c(multi5$Feq, multi5$D, multi5$SPR, multi5$SB, multi5$catchB, multi5$YPR), 3),
  Multi2_Effort5 = round(c(iter_eff_5$Feq, iter_eff_5$D, iter_eff_5$SPR, iter_eff_5$SB, iter_eff_5$catchB, iter_eff_5$YPR), 3),
  Multi2_Catch5 = round(c(iter_ct_5$Feq, iter_ct_5$D, iter_ct_5$SPR, iter_ct_5$SB, iter_ct_5$catchB, iter_ct_5$YPR), 3)
)

rownames(comparison_df5) <- c("Feq", "D", "SPR", "SB", "catchB", "YPR")

print(comparison_df5)



#============================================================================================================#
# Example 6: Estimating F eq (finding the F that produce X depletion - doFit = TRUE)- different props and sel
#            Testing the fucntion with 3 fleets
#============================================================================================================#
base6 <- solveD(lh, sel1, doFit = TRUE, D_type = "relB", D_in = 0.4)

multi6 <- solveD_multifleet(lh, list(sel1, sel2, sel3), doFit = TRUE, D_type = "relB", D_in = 0.4,
                            fleet_proportions = c(0.35, 0.45, 0.2))


iter_eff_6 <- solveD_multifleet2(lh, list(sel1, sel2, sel3), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                 fleet_proportions = c(0.35, 0.45, 0.2),
                                 allocation_type="effort")

iter_ct_6 <- solveD_multifleet2(lh, list(sel1, sel2, sel3), doFit = TRUE, D_type = "relB", D_in = 0.4,
                                fleet_proportions = c(0.35, 0.45, 0.2),
                                allocation_type="catch")

comparison_df6 <- data.frame(
  Single_Fleet6 = round(c(base6$Feq, base6$D, base6$SPR, base6$SB, base6$catchB, base6$YPR), 3),
  Multi_Original6 = round(c(multi6$Feq, multi6$D, multi6$SPR, multi6$SB, multi6$catchB, multi6$YPR), 3),
  Multi2_Effort6 = round(c(iter_eff_6$Feq, iter_eff_6$D, iter_eff_6$SPR, iter_eff_6$SB, iter_eff_6$catchB, iter_eff_6$YPR), 3),
  Multi2_Catch6 = round(c(iter_ct_6$Feq, iter_ct_6$D, iter_ct_6$SPR, iter_ct_6$SB, iter_ct_6$catchB, iter_ct_6$YPR), 3)
)

rownames(comparison_df6) <- c("Feq", "D", "SPR", "SB", "catchB", "YPR")

print(comparison_df6)



#cerate an csv

examples_list <- list(
  "Example1_IdenticalFleets_Fixed" = list(
    Single_Fleet = base11,
    Multi_Original = multi1,
    Multi2_Effort = iter_eff_1,
    Multi2_Catch = iter_ct_1
  ),

  "Example2_IdenticalFleets_FitDepletion" = list(
    Single_Fleet = base2,
    Multi_Original = multi2,
    Multi2_Effort = iter_eff_2,
    Multi2_Catch = iter_ct_2
  ),

  "Example3_IdenticalFleets_DiffProps" = list(
    Single_Fleet = base3,
    Multi_Original = multi3,
    Multi2_Effort = iter_eff_3,
    Multi2_Catch = iter_ct_3
  ),

  "Example4_DiffSelectivity_50-50" = list(
    Single_Fleet = base4,
    Multi_Original = multi4,
    Multi2_Effort = iter_eff_4,
    Multi2_Catch = iter_ct_4
  ),

  "Example5_DiffSelectivity_70-30" = list(
    Single_Fleet = base5,
    Multi_Original = multi5,
    Multi2_Effort = iter_eff_5,
    Multi2_Catch = iter_ct_5
  ),

  "Example6_ThreeFleets_DiffSel" = list(
    Single_Fleet = base6,
    Multi_Original = multi6,
    Multi2_Effort = iter_eff_6,
    Multi2_Catch = iter_ct_6
  )
)
create_multifleet_comparison_csv <- function(examples_list, filename = "multifleet_comparison_results.csv") {

  #extract fleet-specific data with separate columns
  extract_fleet_data <- function(result) {
    fleet_data <- list()

    # N fleets
    if(!is.null(result$nfleets)) {
      fleet_data[["nfleets"]] <- result$nfleets
    } else {
      fleet_data[["nfleets"]] <- 1
    }

    # allocation type
    if(!is.null(result$allocation_type)) {
      fleet_data[["allocation_type"]] <- result$allocation_type
    } else if(!is.null(result$nfleets)) {
      fleet_data[["allocation_type"]] <- "effort"  # Multi_Original is effort-based
    } else {
      fleet_data[["allocation_type"]] <- "single_fleet"
    }

    # F by fleet (create 3 columns, fill with NA if needed)
    if(!is.null(result$F_by_fleet)) {
      for(i in 1:3) {
        if(i <= length(result$F_by_fleet)) {
          fleet_data[[paste0("F_fleet_", i)]] <- round(result$F_by_fleet[i], 4)
        } else {
          fleet_data[[paste0("F_fleet_", i)]] <- NA
        }
      }
    } else {
      fleet_data[["F_fleet_1"]] <- round(result$Feq, 4)
      fleet_data[["F_fleet_2"]] <- NA
      fleet_data[["F_fleet_3"]] <- NA
    }


    if(!is.null(result$fleet_proportions)) {
      for(i in 1:3) {
        if(i <= length(result$fleet_proportions)) {
          fleet_data[[paste0("input_prop_fleet_", i)]] <- round(result$fleet_proportions[i], 3)
        } else {
          fleet_data[[paste0("input_prop_fleet_", i)]] <- NA
        }
      }
    } else {
      fleet_data[["input_prop_fleet_1"]] <- 1.0
      fleet_data[["input_prop_fleet_2"]] <- NA
      fleet_data[["input_prop_fleet_3"]] <- NA
    }

    if(!is.null(result$target_catch_proportions)) {
      for(i in 1:3) {
        if(i <= length(result$target_catch_proportions)) {
          fleet_data[[paste0("target_catch_prop_fleet_", i)]] <- round(result$target_catch_proportions[i], 3)
        } else {
          fleet_data[[paste0("target_catch_prop_fleet_", i)]] <- NA
        }
      }
    } else {
      fleet_data[["target_catch_prop_fleet_1"]] <- NA
      fleet_data[["target_catch_prop_fleet_2"]] <- NA
      fleet_data[["target_catch_prop_fleet_3"]] <- NA
    }


    # catch by fleet
    if(!is.null(result$catchB_by_fleet)) {
      for(i in 1:3) {
        if(i <= length(result$catchB_by_fleet)) {
          fleet_data[[paste0("catchB_fleet_", i)]] <- round(result$catchB_by_fleet[i], 1)
        } else {
          fleet_data[[paste0("catchB_fleet_", i)]] <- NA
        }
      }
    } else {
      fleet_data[["catchB_fleet_1"]] <- round(result$catchB, 1)
      fleet_data[["catchB_fleet_2"]] <- NA
      fleet_data[["catchB_fleet_3"]] <- NA
    }

    # final effort proportions (calculated)
    if(!is.null(result$final_effort_proportions)) {
      for(i in 1:3) {
        if(i <= length(result$final_effort_proportions)) {
          fleet_data[[paste0("final_effort_prop_fleet_", i)]] <- round(result$final_effort_proportions[i], 3)
        } else {
          fleet_data[[paste0("final_effort_prop_fleet_", i)]] <- NA
        }
      }
    } else {
      fleet_data[["final_effort_prop_fleet_1"]] <- NA
      fleet_data[["final_effort_prop_fleet_2"]] <- NA
      fleet_data[["final_effort_prop_fleet_3"]] <- NA
    }

    # actual catch proportions
    if(!is.null(result$actual_catch_proportions)) {
      for(i in 1:3) {
        if(i <= length(result$actual_catch_proportions)) {
          fleet_data[[paste0("actual_catch_prop_fleet_", i)]] <- round(result$actual_catch_proportions[i], 3)
        } else {
          fleet_data[[paste0("actual_catch_prop_fleet_", i)]] <- NA
        }
      }
    } else {
      fleet_data[["actual_catch_prop_fleet_1"]] <- NA
      fleet_data[["actual_catch_prop_fleet_2"]] <- NA
      fleet_data[["actual_catch_prop_fleet_3"]] <- NA
    }

    return(fleet_data)
  }

  # Create data frame
  all_results <- data.frame()

  for(example_name in names(examples_list)) {
    example <- examples_list[[example_name]]

    for(method_name in names(example)) {
      result <- example[[method_name]]


      row_data <- data.frame(
        Example = example_name,
        Method = method_name,
        Feq = round(result$Feq, 4),
        D = round(result$D, 3),
        SPR = round(result$SPR, 3),
        SB = round(result$SB, 1),
        catchB = round(result$catchB, 1),
        YPR = round(result$YPR, 3),
        stringsAsFactors = FALSE
      )

      #fleet-specific data
      fleet_data <- extract_fleet_data(result)
      fleet_df <- data.frame(fleet_data, stringsAsFactors = FALSE)

      # combine
      complete_row <- cbind(row_data, fleet_df)

      # bind
      all_results <- rbind(all_results, complete_row)
    }
  }

  # Write CSV
  write.csv(all_results, filename, row.names = FALSE)
  cat("Results saved to:", filename, "\n")
  cat("Number of rows:", nrow(all_results), "\n")
  cat("Number of columns:", ncol(all_results), "\n")

  return(all_results)
}

comprehensive_results <- create_multifleet_comparison_csv(examples_list, "multifleet_all_examples.csv")

print("Column names:")
print(colnames(comprehensive_results))
print("\nFirst few rows:")
print(head(comprehensive_results[, c("Example", "Method", "allocation_type", "F_fleet_1", "F_fleet_2")]))




iter_ct_1
iter_ct_1$fleet_proportions
iter_ct_1$target_catch_proportions
iter_ct_1$final_effort_proportions
iter_ct_1$actual_catch_proportions

iter_ct_2
iter_ct_2$fleet_proportions
iter_ct_2$target_catch_proportions
iter_ct_2$final_effort_proportions
iter_ct_2$actual_catch_proportions

iter_ct_3
iter_ct_3$fleet_proportions
iter_ct_3$target_catch_proportions
iter_ct_3$final_effort_proportions
iter_ct_3$actual_catch_proportions

iter_ct_4
iter_ct_4$fleet_proportions
iter_ct_4$target_catch_proportions
iter_ct_4$final_effort_proportions
iter_ct_4$actual_catch_proportions

iter_ct_5
iter_ct_5$fleet_proportions
iter_ct_5$target_catch_proportions
iter_ct_5$final_effort_proportions
iter_ct_5$actual_catch_proportions

iter_ct_6
iter_ct_6$fleet_proportions
iter_ct_6$target_catch_proportions
iter_ct_6$final_effort_proportions
iter_ct_6$actual_catch_proportions

