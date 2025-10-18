#Black grouper OM1 - test


rm(list=ls())
devtools::load_all()
library(ggplot2)
library(dplyr)
library(tidyr)


# Life History
lh_obj <- new("LifeHistory")
lh_obj@title <- "black grouper"
lh_obj@speciesName <- "Mycteroperca bonaci"
lh_obj@Linf <- 136.66
lh_obj@K <- 0.1338
lh_obj@t0 <- -0.4568
lh_obj@L50 <- 85.56
lh_obj@L95delta <- 11.41
lh_obj@M <- 0.195 # (0.09, 0.30)
lh_obj@L_type<-"TL"
lh_obj@L_units<-"cm"
lh_obj@LW_A <- 0.00001023602 #cm,kg #0.00001023602
lh_obj@LW_B <- 3.0722
lh_obj@Steep <- 0.72 #(0.58, 0.86)
lh_obj@recSD <- 0.4  # (0.2, 0.6)
lh_obj@recRho <- 0
lh_obj@isHermaph<-TRUE
lh_obj@H50<- 120.8
lh_obj@H95delta<-18.40
lh_obj@R0<-10000

# Time Area
ta <- new("TimeArea")
ta@title <- "Validation Test"
ta@gtg <- 13
ta@areas <- 2
ta@recArea <- c(0.999, 0.001)
ta@iterations <- 6
ta@historicalYears <- 25
ta@historicalBio <- 0.45 #c(0.3, 0.6)
ta@historicalBioType <- "relB"
ta@move <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = FALSE)

# Historical effort - increasing trend
colA1 <- seq(from = 1.5, to = 2.5, length.out = 25)
colA2 <- seq(from = 1, to = 1.5, length.out = 25)

ta@historicalEffort <- matrix(c(colA1, colA2), nrow = 25, ncol = 2, byrow = FALSE)


# Stochastic
stochastic_obj <- new("Stochastic")
stochastic_obj@historicalBio <- c(0.3, 0.6)
stochastic_obj@M<-c(0.09, 0.30)
stochastic_obj@Steep <- c(0.58, 0.86)
stochastic_obj@recSD<- c(0.2, 0.6)




#Selectivity: Fishery independent survey (RVC)
survey1_sel_hist <- new("Fishery")
survey1_sel_hist@title <- "Research Survey Historical"
survey1_sel_hist@vulType <- "explog" # dome-shaped
survey1_sel_hist@vulParams <- c(0.4,24.2,0.08) #dome-shaped - highest peak 35–40 cm -values above 0.2 produce a strongly dome-shaped (max value 0.5)
survey1_sel_hist@retType <- "full"
survey1_sel_hist@retMax <- 1
survey1_sel_hist@Dmort <- 0


survey1_sel_proj <- new("Fishery")
survey1_sel_proj@title <- "Research Survey Projection"
survey1_sel_proj@vulType <- "explog"
survey1_sel_proj@vulParams <- c(0.4,24.2,0.08)
survey1_sel_proj@retType <- "full"
survey1_sel_proj@retMax <- 1
survey1_sel_proj@Dmort <- 0



#Selectivity: fisheries (recreational - commercial)
#Recreational MRIP
fleet1_sel_hist <- new("Fishery")
fleet1_sel_hist@title <- "Fleet 1"
fleet1_sel_hist@vulType <- "explog"
fleet1_sel_hist@vulParams <- c(0.28,38.3,0.14)
fleet1_sel_hist@retType <- "logistic"
fleet1_sel_hist@retParams<-c(60.96,1)
fleet1_sel_hist@retMax <- 1
fleet1_sel_hist@Dmort <- 0.175 #9 -26%

fleet1_sel_proj <- new("Fishery")
fleet1_sel_proj@title <- "Fleet 1"
fleet1_sel_proj@vulType <- "explog"
fleet1_sel_proj@vulParams <- c(0.28,38.3,0.14)
fleet1_sel_proj@retType <- "logistic"
fleet1_sel_proj@retParams<-c(60.96,1)
fleet1_sel_proj@retMax <- 1
fleet1_sel_proj@Dmort <- 0.175 #9 -26%

#Commercial longline
fleet2_sel_hist <- new("Fishery")
fleet2_sel_hist@title <- "Fleet 2"
fleet2_sel_hist@vulType <- "logistic"
fleet2_sel_hist@vulParams <- c(86.27, 24.57) #95%= 98.9 #mini size 60. 96 cm
fleet2_sel_hist@retType <- "logistic"
fleet2_sel_hist@retParams<-c(60.96,1)
fleet2_sel_hist@retMax <- 1
fleet2_sel_hist@Dmort <- 0.375

fleet2_sel_proj <- new("Fishery")
fleet2_sel_proj@title <- "Fleet 2"
fleet2_sel_proj@vulType <- "logistic"
fleet2_sel_proj@vulParams <- c(86.27, 24.57) #95%= 98.9 #mini size 60. 96 cm
fleet2_sel_proj@retType <- "logistic"
fleet2_sel_proj@retParams<-c(60.96,1)
fleet2_sel_proj@retMax <- 1
fleet2_sel_proj@Dmort <- 0.375 #Sedar 48 recommended 25-50% -  midpoint assumed of 37.5%


#---Visualize LH
#To simply display to the console
lhOut<-LHwrapper(lh_obj, ta, doPlot = TRUE)
lhOut

#survey
selWrapper(lh = lhOut, ta, FisheryObj = survey1_sel_hist, doPlot = TRUE)

#recreational
selWrapper(lh = lhOut, ta, FisheryObj = fleet1_sel_hist, doPlot = TRUE)

#commercial
selWrapper(lh = lhOut, ta, FisheryObj = fleet2_sel_hist, doPlot = TRUE)

# #save?
# selWrapper(lh = lhOut, ta, FisheryObj = fleet2_sel_hist, doPlot = TRUE, wd = here(), imageName = "Vulnerability", dpi = 300)


# Shared random seed
test_seed <- 123456



#simple strategy for testing
# Multifleet strategy
multifleet_BG <- function(phase, dataObject) {
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase == 1) {
    combined_data <- list()

    #if indexObj exist
    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject) #returns 1-row tibble
      #add each column from the tibble to combined_data

      cat(sprintf("=== Year %d (j=%d), Iteration %d ===\n", j-1, j, k))
      cat("Columns:", ncol(index_result), "\n")
      cat("Names:", paste(sort(names(index_result)), collapse=", "), "\n\n")


      for(col_name in names(index_result)) {
        combined_data[[col_name]] <- index_result[[col_name]]
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

  if(phase == 2) return(list())


  #Vania edit's to match Bill's edits
  if(phase == 3) {

    Flocal <- data.frame()

    #creating one row per area-fleet combination
    for(m in 1:areas) {
      for(f in 1:nfleets) {

        #set fleet-specific F
        if(f == 1) {
          F_value <- 0.3  # Fleet 1 (Recreational) - higher F
        } else if(f == 2) {
          F_value <- 0.05  # Fleet 2 (Commercial) - lower F
        }
        #row contain: [year, iteration, area, fleet, F_value]
        Flocal <- rbind(Flocal, c(j, k, m, f, F_value))
      }
    }

    #results in 4 rows: (1,1), (1,2), (2,1), (2,2) for 2 areas × 2 fleets
    #                   (A1,F1) (A1,F2)

    return(list(
      year = Flocal[,1],           # [j, j, j, j]
      iteration = Flocal[,2],      # [k, k, k, k]
      area = Flocal[,3],           # [1, 1, 2, 2]
      fleet = Flocal[,4],          # [1, 2, 1, 2]
      Flocal = Flocal[,5]          # [0.05, 0.05, 0.05, 0.05]
    ))
  }
}



# Strategy objects
strategy_BG  <- new("Strategy")
strategy_BG @title <- "multifleet-test_BG"
strategy_BG @projectionYears <- 10
strategy_BG @projectionName <- "multifleet_BG"
strategy_BG@projectionParams <- list()


# ============================================================================
# MULTIFLEET (2 FLEETS) OBJECT
# ============================================================================

multifleet_2fleet <- new("Multifleet")
multifleet_2fleet@nfleets <- 2
multifleet_2fleet@fleet_proportions <- c(0.6312, 0.3688)
multifleet_2fleet@allocation_type <- "catch"
multifleet_2fleet@fleet_selectivity_hist_list <- list(fleet1_sel_hist, fleet2_sel_hist)

#proj
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
# OBSERVATION MODELS FOR BLACK GROUPER
# ============================================================================

# INDICES
multi_comprehensive_index <- new("Index")
multi_comprehensive_index@indexID <- "Indices_black_grouper"
multi_comprehensive_index@title <- "Indices_black_grouper"
multi_comprehensive_index@useWeight <- FALSE

# Survey selectivities (same as single fleet)
multi_comprehensive_index@selectivity_hist_list <- list(survey1_sel_hist)
multi_comprehensive_index@selectivity_proj_list <- list(survey1_sel_proj)

multi_comprehensive_index@survey_design <- list(
  # 1. FI Survey (RVC)
  list(
    indextype = "FI",
    areas = c(1),
    indexYears = seq(1, 35, 1),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.5,     # assuming midyear for now
    q_hist_bounds = c(1.56e-06, 2.41e-06),  #assuming low catchability range - does not observe the entire stock
    q_proj_bounds = c(1.56e-06, 2.41e-06),
    hyperstability_hist_bounds = c(0.8, 1.2),
    hyperstability_proj_bounds = c(0.8, 1.2),
    obsError_CV_hist_bounds = c(0.086, 0.208),
    obsError_CV_proj_bounds = c(0.086, 0.208)
  ),

  # 1. FD index (MRIP)
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    indexYears = seq(1, 35, 1),
    q_hist_bounds = c(2.62e-06, 4.07e-06),
    q_proj_bounds = c(2.62e-06, 4.07e-06),
    hyperstability_hist_bounds = c(0.8, 1.2),
    hyperstability_proj_bounds = c(0.8, 1.2),
    obsError_CV_hist_bounds = c(0.12, 0.30),
    obsError_CV_proj_bounds = c(0.12, 0.30)
  )

)

#LENGTH COMPOSITION
multi_comprehensive_lcomp <- new("LCompObs")
multi_comprehensive_lcomp@indexID <- "MultifleetLComp_BG"
multi_comprehensive_lcomp@title <- "Fleet-Specific-Length-Comp-BG"
multi_comprehensive_lcomp@length_bin_width <- 1
multi_comprehensive_lcomp@selectivity_hist_list <- list(survey1_sel_hist)
multi_comprehensive_lcomp@selectivity_proj_list <- list(survey1_sel_proj)

multi_comprehensive_lcomp@survey_design <- list(
  # 1. Fleet 1 fishery length composition - MRIP
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    years = seq(1, 35, 1),
    sample_sizes = rep(100, 35)
  ),

  # 2. Fleet 2 fishery length composition - Longline
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(1),
    years = seq(1, 35, 1),
    sample_sizes = rep(100, 35)
  ),

  # 3. Survey 1 length composition - RVC
  list(
    indextype = "FI",
    areas = c(1),
    years = seq(1,35,1),
    sample_sizes = rep(150, 35),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.5
  )

)

cat("\n=== Testing with obs models ===\n")



result_OM1_BG<- runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StrategyObj = strategy_BG,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_2fleet,
  IndexObj = multi_comprehensive_index,
  CatchObsObj = NULL,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = getwd(),
  fileName = "test_OM1_BG",
  seed = test_seed,
  doPlot = FALSE,
  doDiagnostic = FALSE,
  customToCluster = "multifleet_BG"
)
cat(" multi-fleet BG OM1 - simulation completed\n")

# all debugging code reviewed!! and passed!!



result_OM1_BG <- readProjection(getwd(), "test_OM1_BG")
# my setup create relSb last historical ~ average 0.2
result_OM1_BG$dynamics$relSB
result_OM1_BG$dynamics$relSB[26, ]  # Year 26 = end of historical
mean(result_OM1_BG$dynamics$relSB[26, ])  # Average across iterations

dim(result_OM1_BG$dynamics$multifleet$Ftotal_by_fleet)#[years, iter, area, fleet]
result_OM1_BG$dynamics$multifleet$Ftotal_by_fleet

# Population outputs
result_OM1_BG$dynamics$multifleet$actual_catch_proportions
result_OM1_BG$dynamics$SB
result_OM1_BG$dynamics$recN
result_OM1_BG$dynamics$SPR

# obs model outputs (new structure data frame)
result_OM1_BG$HCR$decisionData
result_OM1_BG$HCR$decisionData$IDX_Survey_1
result_OM1_BG$HCR$decisionLocal

result_OM1_BG$HCR$decisionData$fleet_1_observed_catch_area_1
result_OM1_BG$HCR$decisionData$fleet_1_observed_catch_area_2
result_OM1_BG$HCR$decisionData$fleet_2_observed_catch_area_1
result_OM1_BG$HCR$decisionData$fleet_2_observed_catch_area_2


result_OM1_BG$HCR$decisionAnnual$TAC
result_OM1_BG$HCR$decisionAnnual

result_OM1_BG$HCR$decisionLocal
result_OM1_BG$HCR$decisionAnnual


#checkeando convergencia de Newton Raphson

# nr_diag <- nr_diagnostics_multifleet
# #write.csv(nr_diag, "nr_diag.csv")
#
# #check for failures
# failures <- nr_diag[!nr_diag$converged, ]
# if(nrow(failures) > 0) {
#   cat("WARNING: Found", nrow(failures), "convergence failures:\n")
#   print(failures)
# } else {
#   cat("All solver calls converged successfully!\n")
# }
# View(nr_diag)
# all(nr_diagnostics_multifleet$converged)

#obs: need to add units
plot_SB(result_OM1_BG)
plot_SB(result_OM1_BG, areas=1)
plot_SB(result_OM1_BG, areas=2)


plot_catchB(result_OM1_BG)
plot_catchB(result_OM1_BG,areas=1)
plot_catchB(result_OM1_BG,areas=2)

plot_catchN(result_OM1_BG)
plot_catchN(result_OM1_BG, areas=1)
plot_catchN(result_OM1_BG, areas=2)

# plot_discB(result_OM1_BG)
# plot_discB(result_OM1_BG,areas=1)
# plot_discB(result_OM1_BG,areas=2)

plot_discN(result_OM1_BG)
plot_discN(result_OM1_BG,areas=1)
plot_discN(result_OM1_BG,areas=2)

plot_catchB_multi(result_OM1_BG)
plot_catchB_multi(result_OM1_BG, areas=1)
plot_catchB_multi(result_OM1_BG, areas=2)

plot_catchN_multi(result_OM1_BG,show_individual = TRUE)
plot_catchN_multi(result_OM1_BG, areas=1)
plot_catchN_multi(result_OM1_BG, areas=2)

plot_SPR(result_OM1_BG)
plot_recN(result_OM1_BG)


#plot obs models (indices)
plot_survey_indices(result_OM1_BG)
plot_cpue_indices(result_OM1_BG)

plot_all_indices(result_OM1_BG)

#plot individual indices
plot_indices(result_OM1_BG,
             index_pattern = "IDX_Survey_1",
             show_individual = FALSE,
             title = "FI survey")

plot_indices(result_OM1_BG,
             index_pattern = "IDX_CPUE_2",
             show_individual = FALSE,
             title = "FD index")





# Total catch across all areas
plot_catch_observations_both(result_OM1_BG,areas = c(1),show_individual=TRUE)
plot_catch_observations_both(result_OM1_BG, areas = 1)
plot_catch_observations_both(result_OM1_BG, areas = 2)


# plot LC obs models
# NEW: Area-specific functions (median across iterations are dispayed)
plot_fishery_length_comp(result_OM1_BG,show_individual = TRUE)
# plot_fishery_length_comp(result_OM1_BG, areas=1,show_individual = TRUE)
# plot_fishery_length_comp(result_OM1_BG, areas=1,show_individual = TRUE)

# plot_survey_length_comp(result_OM1_BG,show_individual = TRUE)
# plot_survey_length_comp(result_OM1_BG, areas=1,show_individual = TRUE)
# plot_survey_length_comp(result_OM1_BG, areas=1,show_individual = TRUE)



# NEW: Custom filtering for fleets and areas
plot_length_composition_by_area(result_OM1_BG,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1),    # Specific areas
                                fleet_filter = c(1),
                                show_individual = TRUE)   # Specific fleets



plot_length_composition_by_area(result_OM1_BG,
                                program_pattern = "LC_Fishery",
                                area_filter = c(1),    # Specific areas
                                fleet_filter = c(2),
                                show_individual = TRUE)   # Specific fleets


plot_length_composition_by_area(result_OM1_BG,
                                program_pattern = "LC_Survey",
                                area_filter = c(1),    # Specific areas
                                show_individual = TRUE)   # Specific fleets



# #===========================================================================================#
# #========================== Exploring NR outputs and performance ===========================#
# #===========================================================================================#
#
# #check SB
plot_SB(result_OM1_BG,areas=1)
plot_SB(result_OM1_BG,areas=2)

plot_Ftotal_multi(result_OM1_BG,areas=c(1,2))
plot_Ftotal_multi(result_OM1_BG,areas=c(1))
plot_Ftotal_multi(result_OM1_BG,areas=c(2))
#
# #cacth (add fleet)
# plot_catchB_multi(result_OM1_BG,areas=1)
# plot_catchB_multi(result_OM1_BG,areas=2)
#
# #TAC (add fleet)
# # plot_TAC_by_area(result_OM1_BG, areas = 1)
# # plot_TAC_by_area(result_OM1_BG, areas = 2)
# # plot_TAC_by_area(result_OM1_BG, areas = c(1,2))
#
#
# plot_TAC(result_OM1_BG, areas = "all",show_fleets = TRUE)  # All areas in one plot
# #plot_TAC_by_fleet(result_OM1_BG)
# plot_TAC(result_OM1_BG, areas = c(1), show_fleets = TRUE)
# plot_TAC(result_OM1_BG, areas = c(2), show_fleets = TRUE)
#
# result_OM1_BG$HCR$decisionAnnual
#
# # Area 1 only
# plot_catch_observations_both(result_multi_indexratio_lengthV2, areas = 1)
# # Area 2 only
# plot_catch_observations_both(result_multi_indexratio_lengthV2, areas = 2)
# # Both areas with faceting
# plot_catch_observations_both(result_multi_indexratio_lengthV2, areas = c(1, 2))
# #combined catch obs:
# plot_catch_observations_both(result_multi_indexratio_lengthV2,areas = c(1, 2), show_individual = TRUE)
#
#
# #calculate relative error (problem with area 2 -only)- fixed
# nr_diag$relative_error <- abs(nr_diag$predicted_catch - nr_diag$target_catch) / nr_diag$target_catch
# ggplot(nr_diag, aes(x = target_catch, y = predicted_catch, color = factor(area))) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, color = "red") +
#   labs(title = "Target vs Predicted Catch - should be on 1:1 line")
# #(given TAC target,find F that produces exactly that catch = passed)
#
# #check for extreme values (extreme F values detected, convergencia failed)
# extreme_F <- nr_diag[nr_diag$final_F > 2 | nr_diag$final_F < 0.001, ]
# if(nrow(extreme_F) > 0) {
#   cat("WARNING: Extreme F values detected\n")
#   print(extreme_F)
# }
#
# #realized catch vs target TAC
# catchB <- result_multi_indexratio_lengthV2$dynamics$catchB
# TAC_decisions <- result_multi_indexratio_lengthV2$HCR$decisionAnnual
#
#
# #compare for each year/iteration/area
# proj_start <- 12
# nfleets <- 2
# for(yr in proj_start:(proj_start+2)) {
#   for(iter in 1:3) {
#     for(area in 1:2) {
#       cat(sprintf("\nYr %d, Iter %d, Area %d:\n", yr-1, iter, area))
#
#       for(fleet in 1:nfleets) {
#         # Fleet-specific TAC
#         TAC_fleet <- TAC_decisions$TAC[TAC_decisions$year == yr &
#                                          TAC_decisions$iteration == iter &
#                                          TAC_decisions$area == area &
#                                          TAC_decisions$fleet == fleet]
#
#         # Fleet-specific realized catch
#         realized_fleet <- result_multi_indexratio_lengthV2$dynamics$multifleet$catchB_by_fleet[yr, iter, area, fleet]
#
#         error <- abs(realized_fleet - TAC_fleet) / TAC_fleet * 100
#
#         cat(sprintf("  Fleet %d: TAC=%.2f, Realized=%.2f, Error=%.1f%%\n",
#                     fleet, TAC_fleet, realized_fleet, error))
#
#         if(error > 5) {
#           cat("    WARNING: Error >5%\n")
#         }
#       }
#     }
#   }
# }
#
# #simple selectivity check - compare fleet selectivities directly
# lh <- LHwrapper(result_multi_indexratio_lengthV2$LifeHistoryObj, result_multi_indexratio_lengthV2$TimeAreaObj)
#
# #check multiple ages on the selectivity curve
# test_ages <- c(1,2,3,4,5,6,7,8,9,10,11, 12, 13,14, 15)
#
# for(age in test_ages) {
#   cat(sprintf("\nAge %d:\n", age))
#
#   fleet1_sel <- selWrapper(lh, result_multi_indexratio_lengthV2$TimeAreaObj,
#                            FisheryObj = result_multi_indexratio_lengthV2$MultifleetObj@fleet_selectivity_proj_list[[1]][[1]],
#                            doPlot = FALSE)
#
#   fleet2_sel <- selWrapper(lh, result_multi_indexratio_lengthV2$TimeAreaObj,
#                            FisheryObj = result_multi_indexratio_lengthV2$MultifleetObj@fleet_selectivity_proj_list[[1]][[2]],
#                            doPlot = FALSE)
#
#   f1_keep <- fleet1_sel$keep[[1]][age]
#   f2_keep <- fleet2_sel$keep[[1]][age]
#
#   cat(sprintf("  Fleet 1 (L50=8):  %.4f\n", f1_keep))
#   cat(sprintf("  Fleet 2 (L50=12): %.4f\n", f2_keep))
#   cat(sprintf("  Different? %s\n", abs(f1_keep - f2_keep) > 0.001))
#
# }


# # exploring the obs models
# #plot to explore if indices follow stock trends
#
# # Survey plot - abundance of fish < 10 years old
# plot_survey_young <- function(result, iter = 1, max_age = 10) {
#
#   #survey index
#   index_data <- result$HCR$decisionData
#   survey <- index_data$IDX_Survey_1[index_data$k == iter]
#
#   #calculate abundance of young fish only (ages 1 to max_age)
#   #years 2-36 to match survey years 1-35
#   true_young_abundance <- sapply(2:36, function(year) {
#     sum(sapply(1:length(result$dynamics$N), function(gtg) {
#       sum(result$dynamics$N[[gtg]][1:max_age, year, 1])
#     }))
#   })
#
#   #remove NAs
#   valid_idx <- !is.na(survey)
#   survey_clean <- survey[valid_idx]
#   young_abund_clean <- true_young_abundance[valid_idx]
#   years <- (1:35)[valid_idx]
#
#   #scale to 0-1
#   survey_scaled <- (survey_clean - min(survey_clean)) / (max(survey_clean) - min(survey_clean))
#   young_scaled <- (young_abund_clean - min(young_abund_clean)) / (max(young_abund_clean) - min(young_abund_clean))
#
#   plot(years, young_scaled, type="l", lwd=3, col="steelblue",
#        ylim=c(0,1), xlab="Year", ylab="Scaled (0-1)",
#        main=paste("RVC Survey vs Young Fish Abundance (ages 1-", max_age, ")", sep=""),
#        las=1)
#
#   points(years, survey_scaled, pch=19, col="coral", cex=1.5)
#   lines(years, survey_scaled, lty=2, col="coral", lwd=2)
#   abline(v=25, lty=3, col="gray50", lwd=2)
#
#   legend("topleft",
#          legend=c(paste("True Abundance (ages 1-", max_age, ")", sep=""), "RVC Survey"),
#          col=c("steelblue", "coral"), lwd=3, pch=c(NA, 19), bty="n", cex=1.1)
#
#   cor_val <- cor(survey_scaled, young_scaled)
#   text(5, 0.1, paste("Cor =", round(cor_val, 2)), cex=1.3, font=2)
#
#   return(cor_val)
# }
#
# # CPUE plot - abundance of fish < 10 years old
# plot_cpue_young <- function(result, iter = 1, max_age = 10) {
#
#   #CPUE index
#   index_data <- result$HCR$decisionData
#   cpue <- index_data$IDX_CPUE_2_Fleet_1[index_data$k == iter]
#
#   #calculate abundance of young fish only (ages 1 to max_age)
#   true_young_abundance <- sapply(2:36, function(year) {
#     sum(sapply(1:length(result$dynamics$N), function(gtg) {
#       sum(result$dynamics$N[[gtg]][1:max_age, year, 1])
#     }))
#   })
#
#   #remove NAs
#   valid_idx <- !is.na(cpue)
#   cpue_clean <- cpue[valid_idx]
#   young_abund_clean <- true_young_abundance[valid_idx]
#   years <- (1:35)[valid_idx]
#
#   #scale to 0-1
#   cpue_scaled <- (cpue_clean - min(cpue_clean)) / (max(cpue_clean) - min(cpue_clean))
#   young_scaled <- (young_abund_clean - min(young_abund_clean)) / (max(young_abund_clean) - min(young_abund_clean))
#
#   plot(years, young_scaled, type="l", lwd=3, col="steelblue",
#        ylim=c(0,1), xlab="Year", ylab="Scaled (0-1)",
#        main=paste("MRIP CPUE vs Young Fish Abundance (ages 1-", max_age, ")", sep=""),
#        las=1)
#
#   points(years, cpue_scaled, pch=19, col="darkorange", cex=1.5)
#   lines(years, cpue_scaled, lty=2, col="darkorange", lwd=2)
#   abline(v=25, lty=3, col="gray50", lwd=2)
#
#   legend("topleft",
#          legend=c(paste("True Abundance (ages 1-", max_age, ")", sep=""), "MRIP CPUE"),
#          col=c("steelblue", "darkorange"), lwd=3, pch=c(NA, 19), bty="n", cex=1.1)
#
#   cor_val <- cor(cpue_scaled, young_scaled)
#   text(5, 0.1, paste("Cor =", round(cor_val, 2)), cex=1.3, font=2)
#
#   return(cor_val)
# }
#
#
# par(mfrow=c(1,2), mar=c(4,4,3,1))
# plot_survey_young(result_OM1_BG, iter = 1, max_age = 7)
# plot_cpue_young(result_OM1_BG, iter = 1, max_age = 7)
# par(mfrow=c(1,1))
