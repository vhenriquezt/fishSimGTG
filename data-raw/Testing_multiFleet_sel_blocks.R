
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
ta@historicalYears <- 25
ta@historicalBio <- 0.45 #c(0.3, 0.6)
ta@historicalBioType <- "relB"
ta@move <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = FALSE)

# Stochastic
stochastic_obj <- new("Stochastic")
stochastic_obj@historicalBio <- c(0.2, 0.6) #Dustin's suggestion
stochastic_obj@M<-c(0.09, 0.30)
stochastic_obj@Steep <- c(2.6, 1.42)#parameters from beta distribution


#Selectivity: fisheries (recreational - commercial)
#Recreational
fleet1_sel_hist <- new("Fishery")
fleet1_sel_hist@title <- "Fleet 1"
fleet1_sel_hist@vulType <- "logistic"
fleet1_sel_hist@vulParams <- c(40,2)
fleet1_sel_hist@retType <- "logistic"
fleet1_sel_hist@retParams<-c(60.96,1)
fleet1_sel_hist@retMax <- 1
fleet1_sel_hist@Dmort <- 0.175 #9 -26%

#Commercial longline
fleet2_sel_hist <- new("Fishery")
fleet2_sel_hist@title <- "Fleet 2"
fleet2_sel_hist@vulType <- "logistic"
fleet2_sel_hist@vulParams <- c(86.27, 24.57) #95%= 98.9 #mini size 60. 96 cm
fleet2_sel_hist@retType <- "logistic"
fleet2_sel_hist@retParams<-c(60.96,1)
fleet2_sel_hist@retMax <- 1
fleet2_sel_hist@Dmort <- 0.375

multifleet_2fleet <- new("Multifleet")
multifleet_2fleet@nfleets <- 2
multifleet_2fleet@fleet_proportions <- c(0.6312, 0.3688)
multifleet_2fleet@allocation_type <- "catch"


multifleet_2fleet@fleet_block_hist_list <- list(  #length nfleets
  list(13:25, 1:12), #length n breaks
  list(1:20, 21:25) #length n breaks
)

multifleet_2fleet@fleet_selectivity_hist_list <- list(  #length nfleets
  list(fleet1_sel_hist, fleet2_sel_hist), #length n breaks
  list(fleet2_sel_hist, fleet1_sel_hist) #length n breaks
)

multifleet_2fleet@fleet_historicalEffort <- array(dim = c(ta@historicalYears, ta@areas, 2))
multifleet_2fleet@fleet_historicalEffort[,,1] <- matrix(1:1, nrow = 25, ncol = 2, byrow = FALSE)
multifleet_2fleet@fleet_historicalEffort[,,2] <- matrix(1:1, nrow = 25, ncol = 2, byrow = FALSE)


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

# Survey selectivities (same as single fleet)
multi_comprehensive_index@selectivity_hist_list <- list(survey1_sel_hist)
#multi_comprehensive_index@selectivity_proj_list <- list(survey1_sel_proj)

multi_comprehensive_index@survey_design <- list(
  # 1. FI Survey (RVC)
  list(
    indextype = "FI",
    areas = c(1),
    indexYears = seq(1, 75, 1),
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
    indexYears = seq(1, 75, 1),
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
#multi_comprehensive_lcomp@selectivity_proj_list <- list(survey1_sel_proj)

multi_comprehensive_lcomp@survey_design <- list(
  # 1. Fleet 1 fishery length composition - MRIP
  list(
    indextype = "FD",
    fleet_id = 1,
    areas = c(1),
    years = seq(1, 75, 1),
    sample_sizes = rep(100, 75)
  ),

  # 2. Fleet 2 fishery length composition - Longline
  list(
    indextype = "FD",
    fleet_id = 2,
    areas = c(1),
    years = seq(1, 75, 1),
    sample_sizes = rep(100, 75)
  ),

  # 3. Survey 1 length composition - RVC
  list(
    indextype = "FI",
    areas = c(1),
    years = seq(1,75,1),
    sample_sizes = rep(150, 75),
    selectivity_hist_idx = 1,
    selectivity_proj_idx = 1,
    survey_timing = 0.5
  )

)

#--------
runProjection(
  LifeHistoryObj = lh_obj,
  TimeAreaObj = ta,
  StochasticObj = stochastic_obj,
  MultifleetObj = multifleet_2fleet,
  IndexObj = multi_comprehensive_index,
  CatchObsObj = NULL,
  LengthCompObj = multi_comprehensive_lcomp,
  wd = paste0(getwd(), "/Test"),
  fileName = "Test",
  seed = 1,
  doPlot = TRUE,
  doDiagnostic = FALSE
)



