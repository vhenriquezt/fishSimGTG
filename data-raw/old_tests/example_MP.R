# This first sections create the objects required by fishSimGTG and
# run the simulation using runProjection()

#clean the memory
rm(list=ls())
#load the modified fucntion
#source("obs_models_fishSimGTG_updated.R")
#load fishSimGTG
#library(fishSimGTG)

# used to test the fucntions to see if everything is working

devtools::load_all()
library(ggplot2)

# Create simple examples of each class to understand their structure
lh <- new("LifeHistory")
ta <- new("TimeArea")
fishery <- new("Fishery")
strategy <- new("Strategy")
stochastic <- new("Stochastic")

lh@title<-"Kole"
lh@speciesName<-"Ctenochaetus strigosus"
lh@Linf<-17.7
lh@K<-0.423
lh@t0<- -0.51
lh@L50<-8.4
lh@L95delta<-1.26
lh@M<-0.08
lh@L_type<-"FL"
lh@L_units<-"cm"
lh@LW_A<-0.046
lh@LW_B<-2.85
lh@Steep<-0.54
lh@recSD<-0 #Run with no rec var'n to see deterministic trends
lh@recRho<-0
lh@isHermaph<-FALSE
lh@R0<-10000

#Fishery ste up
fishery@title<-"Test"
fishery@vulType<-"logistic"
fishery@vulParams<-c(10.2,2) #Approx. knife edge
fishery@retType<-"full"
fishery@retMax <- 1
fishery@Dmort <- 0

#Time area set up
# The TimeArea class controls the size and structure of the simulation
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
ta@historicalEffort <- matrix(c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6),
                              nrow = 10, ncol = 2, byrow = FALSE)
dim(ta@historicalEffort)

#stochastic object set up
stochastic@historicalBio = c(0.3, 0.6)
stochastic@Steep= c(0.45, 0.75)

#projection fishery object
ProFisheryObj<-new("Fishery")
slotNames(ProFisheryObj)
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(10.2,0.1)
ProFisheryObj@retType<-"logistic"
ProFisheryObj@retParams <- c(10.2, 0.1)
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0



testMP <- function(phase, dataObject) {

  #------------------------------------
  #Unpacking: required - do not delete
  #------------------------------------
  #Unpack dataObject - objects passed from the OM to the MP
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])


  #------------------------------------------
  #Book keeping - be aware of year j indexing!
  #------------------------------------------
  #OM years are indexed 1:(1 + TimeAreaObj@historicalYears + StrategyObj@projectionYears)

  #Indexes specified within Class Strategy are typically indexed 1:projectionYears
  #To obtain the current projection year:
  #yr <- j - TimeAreaObj@historicalYears - 1

  #To obtain the terminal year of the historical period
  #yrHist <- TimeAreaObj@historicalYears + 1


  #-----------------------------------------------------------------
  #Summary of all available objects unpacked and passed from the OM
  #-----------------------------------------------------------------
  #j, current time step in OM
  #k, current iteration in OM
  #is, initial stock (equilibrium) for iteration k, year 1, contains values returned by solveD
  #lh, LHwrapper object for iteration k, contains values returned by LHwrapper
  #areas, contains TimeAreaObj@areas
  #ageClasses, contains the number of age classes. Obtained from lh$ageClasses. Caution: age 0 is always the first age class.
  #N, abundance-at-age for iteration k. A list of length lh$gtg that contains elements array(dim=c(ageClasses, years, areas))
  #selGroup, selectivity for iteration k and step j. A list of length areas with elements contains values returned by selWrapper
  #selHist, selectivity for iteration k for historical time period. A list of length areas with elements contains values returned by selWrapper
  #selPro, selectivity for iteration k for projection time period. A list of length areas with elements contains values returned by selWrapper
  #SB, spawing biomass at the beginnning of the year, array(dim=c(years, iterations, areas))
  #VB, vulnerable biomass at the beginning of the year, array(dim=c(years, iterations, areas))
  #RB, retained biomass at the beginning of the year, array(dim=c(years, iterations, areas))
  #catchN, annual catch in numbers in aggregate across gtg, array(dim=c(years, iterations, areas))
  #catchB, annual catch in biomass in aggregate across gtg, array(dim=c(years, iterations, areas))
  #discN, dead discards in numbers in aggregate across gtg, array(dim=c(years, iterations, areas))
  #discB, dead discards in biomass in aggregate across gtg, array(dim=c(years, iterations, areas))
  #Ftotal, fishing mortality rate, array(dim=c(years, iterations, areas))
  #SPR, spawning potential ratio, for entire population, array(dim=c(years, iterations))
  #relB, relative spawning biomass (depletion) for entire population, array(dim=c(years, iterations))
  #recN, annual age-0 recruitment in numbers, array(dim=c(years, iterations))
  #decisionData, optional data frame containing annual sampling details. (see phase 1 below)
  #decisionAnnual, optional data frame containing analysis and HCR info (see phase 2 below)
  #decisionLocal, mandatory data frame containing fishing mortality by area for year j and iteration k, (see phase 3 below)
  #RdevMatrix, recruitment deviations, contains values returned by recDev
  #Ddev, initial relative biomass deviations, contains values returned by bioDev
  #TimeAreaObj
  #StrategyObj


  #--------------------------------------
  #Phases of an MP
  #--------------------------------------
  #An MP can have up to 3 phases. Phase 3 is mandatory, phases 1 & 2 are optional
  #The OM will call this MP 3 times in each time step:
  #Phase 1: Observation model or simulated sampling/observations of the OM.
  #Phase 2: Analysis and harvest control rule
  #Phase 3: Calculate fishing mortality to be passed back to the OM.

  ########
  #Phase 1 - OPTIONAL. Observation process in year j
  #Phase 1 is called by the OM once for each iteration, k, and year, j.
  #Used to define imperfect sampling in year j
  #Returned list is appended to data frame: decisionData
  if(phase==1){
    #User defines computations needed.



    #User defines variable names for returned list, as this info will be used in Phase 2 (HCR)
    #return(list()) #e.g., return(list(year = j, iteration = k, CPUE = 1.4))
    return(calculate_single_Index(dataObject))
  }

  ########
  #Phase 2 - OPTIONAL. Analysis, decision-rule and/or HCR
  #Phase 2 is called by the OM once for each iteration, k, and year, j.
  #Used to define how catches will be regulated
  #Returned list is appended to data frame: decisionAnnual
  if(phase==2){
    #User defines computations needed (likely will rely on info saved to decisionData)

    #User defines variable names for returned list, as this info will be used in Phase 3
    return(list()) #e.g., return(list(year = j, iteration = k, TAC = 1400))

  }

  ########
  #Phase 3 - MANDATORY. Calculation of fishing mortality to be used by the OM
  #Phase 3 is called by the OM once for each iteration, k, and year, j.
  #Used to compute fishing mortality for the called k and j.
  #Returned list is appended to data frame: decisionLocal and retrieved by OM to calcualte catchN, catchB, etc.
  if(phase==3){

    #Must contain these elements, each a vector of length areas:
    year = rep(j, areas) #Time step to which fishing mortality applies
    iteration = rep(k, areas) #Iteration to which fishing mortality applies
    area = 1:areas #Area to which fishing mortality applies
    Flocal = rep(0.2, areas) #Calculation of fishing mortality (typically a complex computation is needed, e.g., to covert TAC to F)

    #Must return this list structure:
    return(list(year=year, iteration=iteration, area=area,  Flocal=Flocal))
  }
}


strategy@title <- "testing Obs"           # add a descriptive name
strategy@projectionYears <- 10                  #future projection
strategy@projectionName <- "testMP" # Which function to use - next: the parameters of projectionStrategy



#---------------------------------------------------------------------------------------#
# Observation model 1: Simulation of one CPUE (biomass) index covering both areas       #
#---------------------------------------------------------------------------------------#

# Create a new Index class
cpueB1 <- new("Index")
cpueB1@indexID <- "One cpueB1 - all_areas"
cpueB1@title <- "One cpueB1 - all_areas"
cpueB1@useWeight <- TRUE                # CPUE in biomass

# Define the survey design with cpue specific parameters
cpueB1@survey_design <- list(
  list(
    indextype = "FD",                          # new: add indextype to EACH index design
    areas = c(1, 2),                           # this CPUE covers both areas (1, 2)
    indexYears = 2:21,                         # all years have data (collect data every year)
    q_hist_bounds = c(0.00008, 0.00012),       # historical period - qbounds: 0.00008-0.00012
    q_proj_bounds = c(0.0001, 0.0003),         # projection period - qbounds: 0.0001-0.0003
    hyperstability_hist_bounds = c(0.85, 1.05),# historical - hypersbounds: 0.85-1.05
    hyperstability_proj_bounds = c(0.85, 1.05),# projection - hypersbounds: 0.85-1.05
    obsError_CV_hist_bounds = c(0.25, 0.35),   # historical CV bounds: 0.25-0.35
    obsError_CV_proj_bounds = c(0.2, 0.3)      # projection CV bounds: 0.2-0.3 (improved precision)
  )
)

# FD indices (CPUE) - empty list for FD data
cpueB1@selectivity_hist_list <- list()
cpueB1@selectivity_proj_list <- list()


# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB1,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_Bcpue1",
  doPlot = TRUE
)

X1<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_Bcpue1")
X1$HCR$decisionData
X1$HCR$decisionData$CPUE_1

#--------------------------------------#
# Plot function for observation models #
#--------------------------------------#

#testing plot fucntion
#tibble_data=X1$HCR$decisionData


plotIndex_tibble <- function(tibble_data, save_plot = FALSE,
                             filename = "index_plot.jpeg") {

  # extract data from the tibble
  title <- unique(tibble_data$title)   #unique gets unique values (removes duplicates)
  historical_end <- unique(tibble_data$historical_end)

  # identify CPUE/Survey columns (the actual index values)
  # find the columns names
  # find the pattern name with grep (must start (^) with CPUE_|Survey_    and have one or more digists(\\d+) $= ends here)
  index_columns <- grep("^(CPUE_|Survey_)\\d+$", names(tibble_data), value = TRUE)

  # create a empty data frame to prepare data for plotting
  all_data <- data.frame()


  # index_col (a column) = eg. "CPUE_1"
  for(index_col in index_columns) {
    # get the corresponding columns
    indextype_col <- paste0(index_col, "_indextype")
    areas_col <- paste0(index_col, "_areas")
    indexyears_col <- paste0(index_col, "_indexYears")

    # extract indexYears for this index
    index_years_str <- unique(tibble_data[[indexyears_col]])[1]
    #  "2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21"
    index_years <- as.numeric(unlist(strsplit(index_years_str, "_")))  # Splits string by "_"
    #2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
    # get index type and areas for panel naming
    index_type <- unique(tibble_data[[indextype_col]])[1]
    areas_str <- unique(tibble_data[[areas_col]])[1]

    # create panel name
    panel_name <- paste(index_col, "-", "Area(s)", areas_str)

    # prepare data for this index
    index_data <- tibble_data %>%
      select(k, j, all_of(index_col)) %>%
      rename(
        iteration = k,
        year = j,
        value = !!sym(index_col) #!!sym=creates a "symbol" representing the column name CPUE_1- !! =unquote = CPUE_1
      ) %>%
      #filter(!is.na(value)) %>%
      mutate(
        panel = panel_name,
        type = "Iterations",
        iteration_label = paste0("Iter_", iteration)
      )

    # calculate median across iterations for each year
    # drops remove the grouping after summarising
    median_data <- index_data %>%
      group_by(year, panel) %>%
      summarise(
        value = median(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      #text labels created using later for plotting
      mutate(
        type = "Median",
        iteration = "Median",
        iteration_label = "Median"
      )

    # set NA for years not in indexYears (to match original behavior)
    all_years <- min(tibble_data$j):max(tibble_data$j)
    missing_years <- setdiff(all_years, index_years)

    if(length(missing_years) > 0) {
      missing_data <- data.frame(
        year = missing_years,
        value = NA,
        panel = panel_name,
        type = "Median",
        iteration = "Median",
        iteration_label = "Median"
      )
      median_data <- rbind(median_data, missing_data)
    }

    # combine iteration and median data
    combined_data <- rbind(
      index_data %>% select(year, value, panel, type, iteration_label),
      median_data %>% select(year, value, panel, type, iteration_label)
    )

    all_data <- rbind(all_data, combined_data)
  }
  #Y-axis label based on data type
  use_weight <- unique(tibble_data$useWeight)[1]
  final_indextype <- unique(tibble_data$final_indextype)[1]

  if(final_indextype == "FD") {
    # Fishery Dependent (CPUE)
    y_label <- if(use_weight) "CPUE (biomass)" else "CPUE (numbers)"
  } else if(final_indextype == "FI") {
    # Fishery Independent (Survey)
    y_label <- if(use_weight) "Survey (biomass)" else "Survey (numbers)"
  } else {
    # Mixed types
    y_label <- if(use_weight) "Index (biomass)" else "Index (numbers)"
  }



  # create the plot (matching original style)
  p <- ggplot(all_data, aes(x = year, y = value)) +
    geom_line(data = subset(all_data, type == "Iterations"),
              aes(group = iteration_label),
              color = "steelblue", alpha = 0.6, size = 0.5) +
    geom_point(data = subset(all_data, type == "Iterations" & !is.na(value)),
               color = "steelblue", alpha = 0.7, size = 1) +
    geom_line(data = subset(all_data, type == "Median"),
              color = "black", size = 1.2) +
    geom_point(data = subset(all_data, type == "Median" & !is.na(value)),
               color = "black", size = 1.5) +
    facet_wrap(~ panel, scales = "free_y") +
    geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
    labs(title = title, x = "Year", y = y_label) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )


  if(save_plot) {
    ggsave(filename, plot = p, width = 12, height = 8, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}



p1 <- plotIndex_tibble(X1$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/index_plot_CPUEB1.jpeg")

p1
#---------------------------------------------------------------------------------------#
# Observation model 2: Simulation of two CPUE (biomass) index covering both areas       #
#---------------------------------------------------------------------------------------#
cpueB2 <- new("Index")
cpueB2@indexID <- "Two_CPUE_programs_all_areas"
cpueB2@title <- "Two CPUE Programs - Both Cover All Areas"
cpueB2@useWeight <- TRUE                # CPUE in biomass

cpueB2@survey_design <- list(

  # CPUE 1
  list(
    indextype = "FD",
    areas = c(1, 2),                              # covers all areas
    indexYears = 2:21,                            # annual data collection
    q_hist_bounds = c(0.0001, 0.00015),           # historical: 0.0001-0.00015
    q_proj_bounds = c(0.00015, 0.0003),           # projection: 0.00015-0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),     # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),   # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),        # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)       # projection CV: 0.15-0.25 (better precision)
  ),

  # CPUE 2
  list(
    indextype = "FD",
    areas = c(1, 2),                                     # also covers all areas
    indexYears = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),  # data collected every other year
    q_hist_bounds = c(0.0001, 0.00015),                  # historical: 0.0001- 0.00015
    q_proj_bounds = c(0.00015, 0.0003),                  # projection: 0.00015- 0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),            # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),          # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),               # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)              # projection CV: 0.15-0.25 (some improvement)
  )
)

# FD indices (CPUE) - empty selectivity lists for FD data
cpueB2@selectivity_hist_list <- list()
cpueB2@selectivity_proj_list <- list()

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB2,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_Bcpue2",
  doPlot = TRUE
)

X2<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_Bcpue2")
X2$HCR$decisionData
X2$HCR$decisionData$areas_total
X2$HCR$decisionData$CPUE_1_areas
X2$HCR$decisionData$CPUE_1
X2$HCR$decisionData$CPUE_2_areas
X2$HCR$decisionData$CPUE_2


p2 <- plotIndex_tibble(X2$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/index_plot_CPUEB2.jpeg")

p2
#-----------------------------------------------------------------------
# Observation model 3: Two CPUE Biomass - Only in Area 1
#-----------------------------------------------------------------------
cpueB3 <- new("Index")
cpueB3@indexID <- "Two_CPUE_programs_only_area1"
cpueB3@title <- "Two CPUE Programs - Cover only Area 1"
cpueB3@useWeight <- TRUE                # CPUE in biomass
cpueB3@survey_design <- list(
  # CPUE 1
  list(
    indextype = "FD",
    areas = c(1),                              # covers  area 1
    indexYears = 2:21,                      # annual data collection
    q_hist_bounds = c(0.0001, 0.00015),        # historical: 0.0001-0.00015
    q_proj_bounds = c(0.00015, 0.0003),        # projection: 0.00015-0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),  # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),# projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),     # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)    # projection CV: 0.15-0.25 (better precision)
  ),

  # CPUE 2
  list(
    indextype = "FD",
    areas = c(1),                                        # also covers Area 1
    indexYears = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),  # data collected every other year
    q_hist_bounds = c(0.0001, 0.00015),                  # historical: 0.0001- 0.00015
    q_proj_bounds = c(0.00015, 0.0003),                  # projection: 0.00015- 0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),            # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),          # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),               # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)              # projection CV: 0.15-0.25 (some improvement)
  )
)
cpueB3@selectivity_hist_list <- list()
cpueB3@selectivity_proj_list <- list()

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB3,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_Bcpue3",
  doPlot = TRUE
)

X3<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_Bcpue3")
X3$HCR$decisionData$CPUE_1
X3$HCR$decisionData$CPUE_2

p3 <- plotIndex_tibble(X3$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/index_plot_CPUEB3.jpeg")

p3


#---------------------------------------------------------------------------------------#
# Observation model 4: Simulation of two CPUE (numbers) index covering both areas       #
#---------------------------------------------------------------------------------------#
cpueN4 <- new("Index")
cpueN4@indexID <- "Two_CPUE_programs_all_areas"
cpueN4@title <- "Two CPUE Programs - Both Cover All Areas"
cpueN4@useWeight <- FALSE                # CPUE in numbers

cpueN4@survey_design <- list(

  # CPUE 1
  list(
    indextype = "FD",
    areas = c(1, 2),                              # covers all areas
    indexYears = 2:21,                            # annual data collection
    q_hist_bounds = c(0.0001, 0.00015),           # historical: 0.0001-0.00015
    q_proj_bounds = c(0.00015, 0.0003),           # projection: 0.00015-0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),     # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),   # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),        # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)       # projection CV: 0.15-0.25 (better precision)
  ),

  # CPUE 2
  list(
    indextype = "FD",
    areas = c(1, 2),                                     # also covers all areas
    indexYears = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),  # data collected every other year
    q_hist_bounds = c(0.0001, 0.00015),                  # historical: 0.0001- 0.00015
    q_proj_bounds = c(0.00015, 0.0003),                  # projection: 0.00015- 0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),            # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),          # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),               # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)              # projection CV: 0.15-0.25 (some improvement)
  )
)

# FD indices (CPUE) - empty selectivity lists for FD data
cpueN4@selectivity_hist_list <- list()
cpueN4@selectivity_proj_list <- list()

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueN4,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_Ncpue4",
  doPlot = TRUE
)

X4<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_Ncpue4")
X4$HCR$decisionData
X4$HCR$decisionData$areas_total
X4$HCR$decisionData$CPUE_1_areas
X4$HCR$decisionData$CPUE_1
X4$HCR$decisionData$CPUE_2_areas
X4$HCR$decisionData$CPUE_2

p4 <- plotIndex_tibble(X4$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/index_plot_CPUEN4.jpeg")

p4

#---------------------------------------------------------------------------------------#
# Observation model 5: Simulation of two Survey (biomass) index covering both areas     #
#---------------------------------------------------------------------------------------#

SurveyB5 <- new("Index")
SurveyB5@indexID <- "Two_survey_programs_all_areas"
SurveyB5@title <- "Two survey Programs - Both Cover All Areas"
SurveyB5@useWeight <- TRUE                # survey in biomass
SurveyB5@survey_design <- list(

  # survey 1
  list(
    indextype = "FI",
    areas = c(1, 2),                    # covers all areas
    indexYears = 2:21,               # annual data collection
    survey_timing = 1,
    selectivity_hist_idx = 1,           # uses survey_sel1_hist
    selectivity_proj_idx = 1,           # uses survey_sel1_proj
    q_hist_bounds = c(0.0001, 0.00015),       # historical: 0.0001-0.00015
    q_proj_bounds = c(0.00015, 0.0003),       # projection: 0.00015-0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),     # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),   # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),        # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)       # projection CV: 0.15-0.25 (better precision)
  ),

  # survey 2
  list(
    indextype = "FI",
    areas = c(1, 2),                    # also covers all areas
    indexYears = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),  # data collected every other year
    survey_timing = 1,
    selectivity_hist_idx = 2,           # uses survey_sel2_hist
    selectivity_proj_idx = 2,           # uses survey_sel2_proj
    q_hist_bounds = c(0.0001, 0.00015),       # historical: 0.0001- 0.00015
    q_proj_bounds = c(0.00015, 0.0003),      # projection: 0.00015- 0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),     # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),     # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),       # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)        # projection CV: 0.15-0.25 (some improvement)
  )
)

# Survey Selectivity 1 (hist)
survey_sel1_hist <- new("Fishery")
survey_sel1_hist@title <- "Survey 1 - Historical"
survey_sel1_hist@vulType <- "logistic"
survey_sel1_hist@vulParams <- c(6, 1.0)
survey_sel1_hist@retType <- "full"
survey_sel1_hist@retMax <- 1.0
survey_sel1_hist@Dmort <- 0.0

# Survey Selectivity 2 (hist)
survey_sel2_hist <- new("Fishery")
survey_sel2_hist@title <- "Survey 2 - Historical"
survey_sel2_hist@vulType <- "logistic"
survey_sel2_hist@vulParams <- c(12, 2.0)
survey_sel2_hist@retType <- "full"
survey_sel2_hist@retMax <- 1.0
survey_sel2_hist@Dmort <- 0.0

# Survey Selectivity 1 (proj)
survey_sel1_proj <- new("Fishery")
survey_sel1_proj@title <- "Survey 1 - Projection"
survey_sel1_proj@vulType <- "logistic"
survey_sel1_proj@vulParams <- c(5, 0.8)
survey_sel1_proj@retType <- "full"
survey_sel1_proj@retMax <- 1.0
survey_sel1_proj@Dmort <- 0.0

# Survey Selectivity 2 (proj)
survey_sel2_proj <- new("Fishery")
survey_sel2_proj@title <- "Survey 2 - Projection"
survey_sel2_proj@vulType <- "logistic"
survey_sel2_proj@vulParams <- c(11, 1.8)   # slightly different
survey_sel2_proj@retType <- "full"
survey_sel2_proj@retMax <- 1.0
survey_sel2_proj@Dmort <- 0.0


# Step 3: FI indices (survey) -  selectivity lists
SurveyB5@selectivity_hist_list <- list(survey_sel1_hist, survey_sel2_hist)
SurveyB5@selectivity_proj_list <- list(survey_sel1_proj, survey_sel2_proj)

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = SurveyB5,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_Bsurvey5",
  doPlot = TRUE
)

X5<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_Bsurvey5")
X5$HCR$decisionData

p5 <- plotIndex_tibble(X5$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/index_plot_SurveyB5.jpeg")

p5

#---------------------------------------------------------------------------------------#
# Observation model 6: Simulation of two Survey (numbers) index covering both areas     #
#---------------------------------------------------------------------------------------#

SurveyN6 <- new("Index")
SurveyN6@indexID <- "Two_survey_programs_all_areas"
SurveyN6@title <- "Two survey Programs - Both Cover All Areas"
SurveyN6@useWeight <- FALSE                # survey in numbers
SurveyN6@survey_design <- list(

  # survey 1
  list(
    indextype = "FI",
    areas = c(1, 2),                    # covers all areas
    indexYears = 2:21,               # annual data collection
    survey_timing = 1,
    selectivity_hist_idx = 1,           # uses survey_sel1_hist
    selectivity_proj_idx = 1,           # uses survey_sel1_proj
    q_hist_bounds = c(0.0001, 0.00015),       # historical: 0.0001-0.00015
    q_proj_bounds = c(0.00015, 0.0003),       # projection: 0.00015-0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),     # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),   # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),        # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)       # projection CV: 0.15-0.25 (better precision)
  ),

  # survey 2
  list(
    indextype = "FI",
    areas = c(1, 2),                    # also covers all areas
    indexYears = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),  # data collected every other year
    survey_timing = 1,
    selectivity_hist_idx = 2,           # uses survey_sel2_hist
    selectivity_proj_idx = 2,           # uses survey_sel2_proj
    q_hist_bounds = c(0.0001, 0.00015),       # historical: 0.0001- 0.00015
    q_proj_bounds = c(0.00015, 0.0003),      # projection: 0.00015- 0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),     # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),     # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),       # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)        # projection CV: 0.15-0.25 (some improvement)
  )
)

# Survey Selectivity 1 (hist)
survey_sel1_hist <- new("Fishery")
survey_sel1_hist@title <- "Survey 1 - Historical"
survey_sel1_hist@vulType <- "logistic"
survey_sel1_hist@vulParams <- c(6, 1.0)
survey_sel1_hist@retType <- "full"
survey_sel1_hist@retMax <- 1.0
survey_sel1_hist@Dmort <- 0.0

# Survey Selectivity 2 (hist)
survey_sel2_hist <- new("Fishery")
survey_sel2_hist@title <- "Survey 2 - Historical"
survey_sel2_hist@vulType <- "logistic"
survey_sel2_hist@vulParams <- c(12, 2.0)
survey_sel2_hist@retType <- "full"
survey_sel2_hist@retMax <- 1.0
survey_sel2_hist@Dmort <- 0.0

# Survey Selectivity 1 (proj)
survey_sel1_proj <- new("Fishery")
survey_sel1_proj@title <- "Survey 1 - Projection"
survey_sel1_proj@vulType <- "logistic"
survey_sel1_proj@vulParams <- c(5, 0.8)
survey_sel1_proj@retType <- "full"
survey_sel1_proj@retMax <- 1.0
survey_sel1_proj@Dmort <- 0.0

# Survey Selectivity 2 (proj)
survey_sel2_proj <- new("Fishery")
survey_sel2_proj@title <- "Survey 2 - Projection"
survey_sel2_proj@vulType <- "logistic"
survey_sel2_proj@vulParams <- c(11, 1.8)   # slightly different
survey_sel2_proj@retType <- "full"
survey_sel2_proj@retMax <- 1.0
survey_sel2_proj@Dmort <- 0.0


# Step 3: FI indices (survey) -  selectivity lists
SurveyN6@selectivity_hist_list <- list(survey_sel1_hist, survey_sel2_hist)
SurveyN6@selectivity_proj_list <- list(survey_sel1_proj, survey_sel2_proj)

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = SurveyN6,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_Nsurvey6",
  doPlot = TRUE
)

X6<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_Nsurvey6")
X6$HCR$decisionData

p6 <- plotIndex_tibble(X6$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/index_plot_SurveyN6.jpeg")

p6


#---------------------------------------------------------------------------------------#
# Observation model 7: Simulation of CPUE (FD) + Survey (FI)  index covering both areas (biomass)    #
#---------------------------------------------------------------------------------------#

mixed_biomass <- new("Index")
mixed_biomass@indexID <- "Mixed_CPUE_and_Survey_Biomass"
mixed_biomass@title <- "Mixed CPUE and  Survey Biomass"
mixed_biomass@useWeight <- TRUE                # biomass

# Define MIXED survey designs
mixed_biomass@survey_design <- list(

  # CPUE (FD)
  list(
    indextype = "FD",
    areas = c(1, 2),                    # covers both areas
    indexYears = 2:21,               # annual data collection
    q_hist_bounds = c(0.0001, 0.00015),       # historical period
    q_proj_bounds = c(0.00012, 0.0002),       # projection period
    hyperstability_hist_bounds = c(0.9, 1.1),     # historical
    hyperstability_proj_bounds = c(0.85, 1.05),   # projection
    obsError_CV_hist_bounds = c(0.2, 0.3),        # Historical CV
    obsError_CV_proj_bounds = c(0.15, 0.25)       # Projection CV (improved precision)
    # Note: No selectivity parameters needed for FD
  ),

  # Survey (FI)
  list(
    indextype = "FI",
    areas = c(1, 2),                    # covers both areas
    indexYears = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),
    survey_timing = 1,
    selectivity_hist_idx = 1,           # uses survey_sel1_hist
    selectivity_proj_idx = 1,           # uses survey_sel1_proj
    q_hist_bounds = c(0.00008, 0.00012),      # historical period
    q_proj_bounds = c(0.0001, 0.00016),       # projection period
    hyperstability_hist_bounds = c(0.95, 1.05),   # historical
    hyperstability_proj_bounds = c(0.9, 1.1),     # projection
    obsError_CV_hist_bounds = c(0.15, 0.25),      # Historical CV
    obsError_CV_proj_bounds = c(0.1, 0.2)         # Projection CV (better precision)
  )
)
# Survey Selectivity 1 (hist)
survey_sel1_hist <- new("Fishery")
survey_sel1_hist@title <- "Survey 1 - Historical"
survey_sel1_hist@vulType <- "logistic"
survey_sel1_hist@vulParams <- c(6, 1.0)
survey_sel1_hist@retType <- "full"
survey_sel1_hist@retMax <- 1.0
survey_sel1_hist@Dmort <- 0.0

# Survey Selectivity 1 (proj)
survey_sel1_proj <- new("Fishery")
survey_sel1_proj@title <- "Survey 1 - Projection"
survey_sel1_proj@vulType <- "logistic"
survey_sel1_proj@vulParams <- c(5, 0.8)
survey_sel1_proj@retType <- "full"
survey_sel1_proj@retMax <- 1.0
survey_sel1_proj@Dmort <- 0.0

# survey selectivity
mixed_biomass@selectivity_hist_list <- list(survey_sel1_hist)
mixed_biomass@selectivity_proj_list <- list(survey_sel1_proj)

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = mixed_biomass,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_BMixed7",
  doPlot = TRUE
)

X7<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_BMixed7")
X7$HCR$decisionData$historical_end

p7 <- plotIndex_tibble(X7$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/index_plot_BMixed7.jpeg")

p7




# Function to create year / iteration table

create_year_iteration_table <- function(data, index_col) {
  # get subset with non-NA values
  subset_data <- data[!is.na(data[[index_col]]), c("k", "j", index_col)]
  # get unique years and iterations
  years <- sort(unique(subset_data$j))
  iterations <- sort(unique(subset_data$k))
  # create empty matrix
  result_matrix <- matrix(NA,
                          nrow = length(years),
                          ncol = length(iterations),
                          dimnames = list(paste0("Year_", years),
                                          paste0("Iter_", iterations)))

  # fill the matrix
  for(i in 1:nrow(subset_data)) {
    year_pos <- which(years == subset_data$j[i])
    iter_pos <- which(iterations == subset_data$k[i])
    result_matrix[year_pos, iter_pos] <- subset_data[[index_col]][i]
  }

  return(result_matrix)
}

X7$HCR$decisionData$CPUE_1
X7$HCR$decisionData$Survey_2

cpue_data <- create_year_iteration_table(X7$HCR$decisionData, "CPUE_1")
survey_data <- create_year_iteration_table(X7$HCR$decisionData, "Survey_2")



