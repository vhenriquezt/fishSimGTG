#Modifying this example with multiple blocks of selectivity to
#mimic walleye history
devtools::load_all()

# Life History parameters
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
lh_obj@LW_A <- 0.00001023602 #cm,kg
lh_obj@LW_B <- 3.0722
lh_obj@Steep <- 0.72 #(0.58, 0.86)
lh_obj@recSD <- 0  # (0.2, 0.6)
lh_obj@recRho <- 0
lh_obj@isHermaph<-TRUE
lh_obj@H50<- 120.8
lh_obj@H95delta<-18.40
lh_obj@R0<-10000
lh_obj@Walpha_units<-"kg"

# Time Area
ta <- new("TimeArea")
ta@title <- "Validation Test"
ta@gtg <- 13
ta@areas <- 2
ta@recArea <- c(0.99, 0.01)
ta@iterations <- 1
ta@historicalYears <- 56
ta@historicalBio <- 0.45 #c(0.3, 0.6)
ta@historicalBioType <- "relB"
ta@move <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = FALSE)

# Stochastic
stochastic_obj <- new("Stochastic")
stochastic_obj@historicalBio <- c(0.2, 0.6) #Dustin's suggestion
stochastic_obj@M<-c(0.09, 0.30)
stochastic_obj@Steep <- c(2.6, 1.42)#parameters from beta distribution



#---------------------------------------------------------
# SELECTIVITY/RETENTION BY REGULATION PERIOD FOR WALLEYE #
#---------------------------------------------------------

# Block 1: 1970-1998 - NO SIZE LIMIT
# Vulnerability same, retention = full retention
sel_block1_noLimit <- new("Fishery")
sel_block1_noLimit@title <- "1970-1998: No size limit"
sel_block1_noLimit@vulType <- "logistic"
sel_block1_noLimit@vulParams <- c(40, 2)  # Gear selectivity unchanged
sel_block1_noLimit@retType <- "full"      # Keep everything caught
sel_block1_noLimit@retMax <- 1
sel_block1_noLimit@Dmort <- 0.175

# Block 2: 1999-2013 - SLOT LIMIT 40-60 cm
# Fish < 40 cm are released, fish 40-60 cm are kept, fish > 60 cm are released
sel_block2_slot40_60 <- new("Fishery")
sel_block2_slot40_60@title <- "1999-2013: Slot 40-60 cm"
sel_block2_slot40_60@vulType <- "logistic"
sel_block2_slot40_60@vulParams <- c(40, 2)
sel_block2_slot40_60@retType <- "slotLimit"  # Slot limit retention
sel_block2_slot40_60@retParams <- c(40, 60)  # minimum and maximum length
sel_block2_slot40_60@retMax <- 1
sel_block2_slot40_60@Dmort <- 0.175

# Block 3: 2014-2021 - MINIMUM SIZE 46 cm
# Fish < 46 cm are released, fish >= 46 cm are kept
sel_block3_min46 <- new("Fishery")
sel_block3_min46@title <- "2014-2021: Min size 46 cm"
sel_block3_min46@vulType <- "logistic"
sel_block3_min46@vulParams <- c(40, 2)
sel_block3_min46@retType <- "logistic"
sel_block3_min46@retParams <- c(46, 2)  # L50 = 46 cm
sel_block3_min46@retMax <- 1
sel_block3_min46@Dmort <- 0.175

# Block 4: 2022-2025 - SLOT LIMIT 40-45 cm (narrow slot)
sel_block4_slot40_45 <- new("Fishery")
sel_block4_slot40_45@title <- "2022-2025: Slot 40-45 cm"
sel_block4_slot40_45@vulType <- "logistic"
sel_block4_slot40_45@vulParams <- c(40, 2)
sel_block4_slot40_45@retType <- "slotLimit"
sel_block4_slot40_45@retParams <- c(40, 45)  #narrow slot 40-45
sel_block4_slot40_45@retMax <- 1
sel_block4_slot40_45@Dmort <- 0.175


multifleet_2fleet <- new("Multifleet")
multifleet_2fleet@nfleets <- 2
multifleet_2fleet@fleet_proportions <- c(0.6312, 0.3688)
multifleet_2fleet@allocation_type <- "catch"

# Defining year blocks for each fleet
# Fleet 1 (e.g., Recreational): All 4 regulation periods apply
# Fleet 2 (e.g., Commercial):   All 4 regulation periods apply

multifleet_2fleet@fleet_block_hist_list <- list(
  # Fleet 1 blocks
  list(
    1:29,    # Block 1: 1970-1998 (years 1-29)
    30:44,   # Block 2: 1999-2013 (years 30-44)
    45:52,   # Block 3: 2014-2021 (years 45-52)
    53:56    # Block 4: 2022-2025 (years 53-56)
  ),
  # Fleet 2 blocks (same periods, but could have different selectivity)
  list(
    1:29,    # Block 1: 1970-1998
    30:44,   # Block 2: 1999-2013
    45:52,   # Block 3: 2014-2021
    53:56    # Block 4: 2022-2025
  )
)


multifleet_2fleet@fleet_selectivity_hist_list <- list(
  # Fleet 1: Recreational - same regulations for all periods
  list(
    sel_block1_noLimit,    # 1970-1998: No size limit
    sel_block2_slot40_60,  # 1999-2013: Slot 40-60
    sel_block3_min46,      # 2014-2021: Min 46 cm
    sel_block4_slot40_45   # 2022-2025: Slot 40-45
  ),
  # Fleet 2: Commercial - same regulations (for now, it could be different)
  list(
    sel_block1_noLimit,    # 1970-1998: No size limit
    sel_block2_slot40_60,  # 1999-2013: Slot 40-60
    sel_block3_min46,      # 2014-2021: Min 46 cm
    sel_block4_slot40_45   # 2022-2025: Slot 40-45
  )
)

#--------------------------------------------------------#
# ADDING PROJECTION TO CREATE A SIMPLE STARTEGY OBJ      #
# SO WE CAN SIMULATE OBS MODELS                          #
#--------------------------------------------------------#

# Projection selectivity: [[area]][[fleet]] (using current regulation)
multifleet_2fleet@fleet_selectivity_proj_list <- list(
  # Area 1
  list(sel_block4_slot40_45, sel_block4_slot40_45),
  # Area 2
  list(sel_block4_slot40_45, sel_block4_slot40_45)
)

multifleet_2fleet@fleet_historicalEffort <- array(dim = c(ta@historicalYears, ta@areas, 2))
multifleet_2fleet@fleet_historicalEffort[,,1] <- matrix(1:1, nrow = 56, ncol = 2, byrow = FALSE)
multifleet_2fleet@fleet_historicalEffort[,,2] <- matrix(1:1, nrow = 56, ncol = 2, byrow = FALSE)


#----------
# INDICES
#Selectivity: Fishery independent survey (RVC)
survey1_sel_hist <- new("Fishery")
survey1_sel_hist@title <- "Research Survey Historical"
survey1_sel_hist@vulType <- "explogFlex" # dome-shaped
survey1_sel_hist@vulParams <- c(0.058,24.2,0.9,24)
survey1_sel_hist@retType <- "full"
survey1_sel_hist@retMax <- 1
survey1_sel_hist@Dmort <- 0

multi_comprehensive_index <- new("Index")
multi_comprehensive_index@indexID <- "Indices_BG_OM1"
multi_comprehensive_index@title <- "Indices_BG_OM1"
multi_comprehensive_index@useWeight <- FALSE
multi_comprehensive_index@selectivity_hist_list <- list(survey1_sel_hist)
multi_comprehensive_index@selectivity_proj_list <- list(survey1_sel_hist)

multi_comprehensive_index@survey_design <- list(
  # 1. FI Survey (RVC)
  list(
    indextype = "FI",
    areas = c(1),
    indexYears = seq(1, 61, 1), # 56 hist + 5 proj
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.5,     # assuming midyear for now
    q_hist_bounds = c(1.985e-06, 1.985e-06),  # I took the average (1.56e-06, 2.41e-06)assuming low catchability range - does not observe the entire stock
    q_proj_bounds = c(1.985e-06, 1.985e-06),
    hyperstability_hist_bounds = c(1, 1),
    hyperstability_proj_bounds = c(1, 1),
    obsError_CV_hist_bounds = c(0.086, 0.208),
    obsError_CV_proj_bounds = c(0.086, 0.208)
  ),

  # 1. FD index (MRIP)
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    indexYears = seq(1, 61, 1),
    q_hist_bounds = c(3.345e-06, 3.345e-06), #i took the avergae (2.62e-06, 4.07e-06)
    q_proj_bounds = c(3.345e-06, 3.345e-06),
    hyperstability_hist_bounds = c(1, 1),
    hyperstability_proj_bounds = c(1, 1),
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
multi_comprehensive_lcomp@selectivity_proj_list <- list(survey1_sel_hist)

multi_comprehensive_lcomp@survey_design <- list(
  # 1. Fleet 1 fishery length composition - MRIP
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    years = seq(1, 61, 1),
    sample_sizes = rep(100, 61)
  ),

  # 2. Fleet 2 fishery length composition - Longline
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(1),
    years = seq(1, 61, 1),
    sample_sizes = rep(100, 61)
  ),

  # 3. Survey 1 length composition - RVC
  list(
    indextype = "FI",
    areas = c(1),
    years = seq(1,61,1),
    sample_sizes = rep(150, 61),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.5
  )

)

#--------------------------------------------------------#
# ADDING A SIMPLE STARTEGY                               #
#--------------------------------------------------------#


# add a simple strategy to explore obs models after adding sel blocks
simple_selblock_test <- function(phase, dataObject) {
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

    # if(!is.null(CatchObsObj)) {
    #   catch_result <- calculate_single_CatchObs(dataObject)
    #   for(col_name in names(catch_result)) {
    #     combined_data[[col_name]] <- catch_result[[col_name]]
    #   }
    # }

    if(!is.null(LengthCompObj)) {
      lc_result <- calculate_single_LengthComp(dataObject)
      for(col_name in names(lc_result)) {
        combined_data[[col_name]] <- lc_result[[col_name]]
      }
    }


    return(combined_data)
  }
  # Phase 2: No TAC decisions needed for testing
  if(phase == 2) return(list())


  if(phase == 3) {

    Flocal <- data.frame()

    #creating one row per area-fleet combination
    for(m in 1:areas) {
      for(f in 1:nfleets) {
        #row contain: [year, iteration, area, fleet, F_value]
        Flocal <- rbind(Flocal, c(j, k, m, f, 0.05))  # Conservative F = 0.05
      }
    }

    # Results in 4 rows: (A1,F1), (A1,F2), (A2,F1), (A2,F2) for 2 areas × 2 fleets

    return(list(
      year = Flocal[,1],           # [j, j, j, j]
      iteration = Flocal[,2],      # [k, k, k, k]
      area = Flocal[,3],           # [1, 1, 2, 2]
      fleet = Flocal[,4],          # [1, 2, 1, 2]
      Flocal = Flocal[,5]          # [0.05, 0.05, 0.05, 0.05]
    ))
  }
}

strategy_multifleet <- new("Strategy")
strategy_multifleet@title <- "Multifleet Block Test"
strategy_multifleet@projectionYears <- 5
strategy_multifleet@projectionName <- "simple_selblock_test"
strategy_multifleet@projectionParams <- list()




#dir.create("P:/Fork_fish_Sim_GTG/fishSimGTG/Test")
#--------
  runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_2fleet,
  StrategyObj = strategy_multifleet,
  IndexObj = multi_comprehensive_index,
  CatchObsObj = NULL,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = paste0(getwd(), "/Test"),
  fileName = "Test_Blocks",
  seed = 1,
  doPlot = TRUE,
  doDiagnostic = FALSE
)


test1_block <- readProjection("Test","Test_Blocks")


lh <- LHwrapper(LifeHistoryObj = lh_obj, TimeAreaObj = ta)
selOut_fleet1 <- selWrapper(lh = lh, TimeAreaObj = ta, FisheryObj = fleet1_sel_hist, doPlot = TRUE)
selOut_fleet2 <- selWrapper(lh = lh, TimeAreaObj = ta, FisheryObj = fleet2_sel_hist, doPlot = TRUE)


#Plots
plot_SB_total_modified(test1_block, save_plot=FALSE,show_individual = FALSE)
plot_recN_modified(test1_block, save_plot=TRUE,show_individual = FALSE)
plot_SPR(test1_block,save_plot=TRUE,show_individual = FALSE)
plot_catchB_total_modified(test1_block, save_plot=TRUE,show_individual = FALSE)
plot_catchN_total_modified(test1_block, save_plot=TRUE,show_individual = FALSE)
plot_discN_total_modified(test1_block, save_plot=TRUE,show_individual = FALSE)
#Fishing Mortality - Area 1 (standardized)
#plot_Ftotal_modified(test1_block, save_plot=TRUE,show_individual = FALSE)
#TAC by Fleet - standardized by CATCH (large fonts)
#plot_TAC_total_modified(test1_block, save_plot=TRUE,show_individual = FALSE)
plot_indices_modified(test1_block,
                      index_pattern = "IDX_Survey_1",
                      title = "FI survey",
                      save_plot = TRUE,show_individual = TRUE)

plot_indices_modified(test1_block,
                      index_pattern = "IDX_CPUE_2",
                      title = "FD index",
                      save_plot = TRUE,show_individual = TRUE)



plot_length_composition_by_area_modified(test1_block,
                                         program_pattern = "LC_Fishery",
                                         area_filter = c(1),    # Specific areas
                                         fleet_filter = c(1),
                                         show_individual = TRUE,
                                         save_plot=TRUE,
                                         filename="LC_Fishery_Fleet1.png")   # Specific fleets

plot_length_composition_by_area_modified(test1_block,
                                         program_pattern = "LC_Fishery",
                                         area_filter = c(1),    # Specific areas
                                         fleet_filter = c(2),
                                         show_individual = TRUE,
                                         save_plot=TRUE,
                                         filename="LC_Fishery_Fleet2.png")   # Specific fleets





