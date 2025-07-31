rm(list=ls())
#devtools::install()
devtools::load_all()
library(ggplot2)
#library(fishSimGTG)
library(dplyr)

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

    #Need to combine obs (e.g., index and catch) in a list before returning
    combined_data <- list()

    # get index obs data and add to combined_data list
    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject)
      # add all index columns
      for(col_name in names(index_result)) {
        combined_data[[col_name]] <- index_result[[col_name]]
      }
    }

    # get catch obs data and add to combined_data list
    if(!is.null(CatchObsObj)) {
      catch_result <- calculate_single_CatchObs(dataObject)
      # add all catch columns
      for(col_name in names(catch_result)) {
        combined_data[[col_name]] <- catch_result[[col_name]]
      }
    }
    # adding the new obs model
    if(!is.null(LengthCompObj)) {
      lc_result <- calculate_single_LengthComp(dataObject)
      # add all length composition columns
      for(col_name in names(lc_result)) {
        combined_data[[col_name]] <- lc_result[[col_name]]
      }
    }

    #User defines variable names for returned list, as this info will be used in Phase 2 (HCR)
    #return(list()) #e.g., return(list(year = j, iteration = k, CPUE = 1.4))
    return(combined_data)
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
    indexYears = 1:20,                         # all years have data (collect data every year)
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

#-----------------------------------------------#
#------------NEW--------------------------------#
#-----------------------------------------------#
# Create a new Catchobs class
catch_obs1 <- new("CatchObs",
                      catchID = "Fishery Catch",
                      title = "Fishery catch observations",
                      areas = c(1, 2),  # Both areas like CPUE
                      catchYears = 1:20,  # Same years as CPUE
                      reporting_rates = c(seq(0.6, 1, length=10), rep(1, 10)),  # Improving historical, non systematic bias in projection
                      obs_CVs = matrix(c(rep(0.2, 20), rep(0.4, 20)), ncol=2))  # CV bounds: 0.2-0.4 for all years


# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB1,
  CatchObsObj = catch_obs1,
  customToCluster = "testMP",
  wd = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs",
  fileName = "test_run_MP_with_catch_obs",
  doPlot = TRUE
)

out1<-readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG/outputs", "test_run_MP_with_catch_obs")
out1$HCR$decisionData
out1$HCR$decisionData$CPUE_1
out1$HCR$decisionData$observed_catch
out1$dynamics$Ftotal


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
    # to match (j-1) changes
    user_years <- (all_years - 1)  # convert j to user year format

    missing_user_years  <- setdiff(user_years, index_years)
    missing_years <- missing_user_years + 1  # convert back to simulation years (j format)

    #Remove year 1 from missing_years since it's initial conditions
    missing_years <- missing_years[missing_years > 1]

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

  # fixing to match the new (j-1) index
  sim_years <- sort(unique(all_data$year))  # get actual years in the data (e.g., 2, 3, 4, ..., 21)
  user_years <- sim_years - 1  # convert to user-friendly years (e.g., 1, 2, 3, ..., 20)

  # convert historical_end to user-friendly scale for plotting
  user_historical_end <- historical_end - 1




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

    #new
    scale_x_continuous(
      breaks = sim_years,    # use actual simulation years for positioning
      labels = user_years    # but show user-friendly labels
    ) +

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



#--------------------------------------------#
# Plot function for catch observation models #
#--------------------------------------------#

plotCatch_tibble <- function(tibble_data, save_plot = FALSE,
                             filename = "catch_plot.jpeg") {

  # extract data from the tibble
  title <- unique(tibble_data$catchID)[1]
  historical_end <- unique(tibble_data$end_hist)

  # only plot observed_catch
  index_columns <- "observed_catch"

  # create a empty data frame
  all_data <- data.frame()

  # process observed_catch only
  for(index_col in index_columns) {

    # define panel name
    areas_str <- unique(tibble_data$areas_included)[1]
    panel_name <- paste("Observed catch -", "Area(s)", areas_str)

    # years with catch data
    index_years <- sort(unique(tibble_data$j[!is.na(tibble_data[[index_col]])]))

    # prepare data same structure as index plot
    index_data <- tibble_data %>%
      select(k, j, all_of(index_col)) %>%
      rename(
        iteration = k,
        year = j,
        value = !!sym(index_col)
      ) %>%
      mutate(
        panel = panel_name,
        type = "Iterations",
        iteration_label = paste0("Iter_", iteration)
      )

    # calculate median across iterations for each year as index plot
    median_data <- index_data %>%
      group_by(year, panel) %>%
      summarise(
        value = median(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        type = "Median",
        iteration = "Median",
        iteration_label = "Median"
      )

    # set NA for years not in catchYears (as index plot)
    all_years <- min(tibble_data$j):max(tibble_data$j)

    #NEW to match (j-1)
    #convert simulation years to user years for comparison
    user_years <- (all_years - 1)  # convert j to user year format
    index_user_years <- index_years - 1  # convert index_years to user years for comparison

    #find missing user years then convert back to simulation years
    missing_user_years <- setdiff(user_years, index_user_years)
    missing_years <- missing_user_years + 1  # Convert back to simulation years (j format)


    # Remove year 1 from missing_years since it is initial conditions
    missing_years <- missing_years[missing_years > 1]

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

    # combine iteration and median data (as index plot)
    combined_data <- rbind(
      index_data %>% select(year, value, panel, type, iteration_label),
      median_data %>% select(year, value, panel, type, iteration_label)
    )

    all_data <- rbind(all_data, combined_data)
  }

  # Y-axis label for catch
  y_label <- "Observed Catch (biomass)"

  #NEW
  # Calculate dynamic x-axis labels (same as index plot)
  sim_years <- sort(unique(all_data$year))  # gets actual years in the data (e.g., 2, 3, 4, ..., 21)
  user_years <- sim_years - 1  # convert to user-friendly years (e.g., 1, 2, 3, ..., 20)

  # convert historical_end to user-friendly scale for plotting
  user_historical_end <- historical_end - 1



  # create the plot (as index plot)
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
    # dynamic x-axis with user-friendly labels
    scale_x_continuous(
      breaks = sim_years,    # use actual simulation years for positioning
      labels = user_years    # but show user-friendly labels
    ) +

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




p1 <- plotIndex_tibble(out1$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/index_plot_CPUEB1.jpeg")

p1


p2 <- plotCatch_tibble (out1$HCR$decisionData,
                           save_plot = TRUE,
                             filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/plot_catchobs1.jpeg")
p2

#higer F in the projection period exaplain the increases in catchs
out1$dynamics$Ftotal


#---------------------------------------------------------------------------------------#
# EXAMPLE 2                                                                             #
# Observation model 2: Simulation of two CPUE (biomass) index covering both areas       #
# adding catch observtion model                                                         #
#---------------------------------------------------------------------------------------#

#Define the new Index class
cpueB2 <- new("Index")
cpueB2@indexID <- "Two_CPUE_programs_all_areas"
cpueB2@title <- "Two CPUE Programs - Both Cover All Areas"
cpueB2@useWeight <- TRUE                # CPUE in biomass

cpueB2@survey_design <- list(

  # CPUE 1
  list(
    indextype = "FD",
    areas = c(1, 2),                              # covers all areas
    indexYears = 1:20,                            # annual data collection
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


# Define CVs for CatchObs class
total_years <- strategy@projectionYears+ta@historicalYears
end_hist <- ta@historicalYears
cv_periods <- matrix(nrow = 20, ncol = 2)  #nrow= years with data, ncol lower, upper bound

for(pos in 1:20) {
  actual_year <- pos

  if(actual_year <= end_hist) {
    # Historical period (years 1-10)
    cv_periods[pos, ] <- c(0.2, 0.4)
  } else {
    # Projection period (years 11-20)
    cv_periods[pos, ] <- c(0.1, 0.2)
  }
}


#Define the new CatchObs class
catch_obs2 <- new("CatchObs",
                  catchID = "period_specific_CV",
                  title = "Different CV bounds by period",
                  areas = c(1, 2),
                  catchYears = 1:total_years,
                  #improving reporting rates over time
                  reporting_rates = seq(0.5, 0.9, length.out = total_years),
                  #more precise catch record for the projection period
                  obs_CVs = cv_periods)


# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB2,
  CatchObsObj = catch_obs2,
  customToCluster = "testMP",
  wd = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs",
  fileName = "test_run_MP_with_catch_obs2",
  doPlot = TRUE
)

out2<-readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG/outputs", "test_run_MP_with_catch_obs2")
out2$HCR$decisionData
out2$HCR$decisionData$CPUE_1
out2$HCR$decisionData$CPUE_2
out2$HCR$decisionData$observed_catch

out2$dynamics$Ftotal


p1 <- plotIndex_tibble(out2$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/index_plot_CPUEB2.jpeg")

p1

p2 <- plotCatch_tibble (out1$HCR$decisionData,
                        save_plot = TRUE,
                        filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/plot_catchobs2.jpeg")
p2



#-----------------------------------------------------------------------
# Observation model 3: Two CPUE Biomass - Only in Area 1
# adding catch observtion model
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
    indexYears = 1:20,                      # annual data collection
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


#Populating CVs for catchobs object
cv_yearly <- matrix(nrow = 20, ncol = 2)
# define a CV for each year individualy
cv_yearly[1, ] <- c(0.2, 0.4)   #  user year 1: CV 0.2-0.4
cv_yearly[2, ] <- c(0.1, 0.3)   #  user year 2
cv_yearly[3, ] <- c(0.15, 0.25) #  user year 3
# user years 4-10: CV 0.15-0.25
cv_yearly[4:10, 1] <- 0.15      # min CV 0.15
cv_yearly[4:10, 2] <- 0.25      # max CV 0.25
# user years 11-20 (projection period): CV 0.18-0.3
cv_yearly[11:20, 1] <- 0.18
cv_yearly[11:20, 2] <- 0.3


catch_obs3 <- new("CatchObs",
                  catchID = "period_specific_CV",
                  title = "Different CV bounds by period",  # only in area 1 catch data are collected
                  areas = c(1),
                  catchYears = 1:total_years,
                  #improving reporting rates over time
                  # historical: overreporting, projection: perfect reporting
                  reporting_rates = c(seq(1.2, 2.5, length=10), rep(1,10)), #length= length catchYears
                  #more precise catch record for the projection period
                  obs_CVs = cv_yearly)


# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB3,
  CatchObsObj = catch_obs3,
  customToCluster = "testMP",
  wd = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs",
  fileName = "test_run_MP_with_catch_obs3",
  doPlot = TRUE
)

out3<-readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG/outputs", "test_run_MP_with_catch_obs3")
out3$HCR$decisionData

out3$HCR$decisionData$CPUE_1
out3$HCR$decisionData$CPUE_2
out3$HCR$decisionData$observed_catch
out3$HCR$decisionData$observed_catch_area_1



p1 <- plotIndex_tibble(out3$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/index_plot_CPUEB3.jpeg")

p1

p2 <- plotCatch_tibble (out3$HCR$decisionData,
                        save_plot = TRUE,
                        filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/plot_catchobs3.jpeg")
p2



#---------------------------------------------------------------------------------------#
# Observation model 4: Simulation of two CPUE (numbers) index covering both areas       #
# adding catch observtion model
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
    indexYears = 1:total_years,                            # annual data collection
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


#Populating CVs for catchobs object
cv_block <- matrix(nrow = 10, ncol = 2)

cv_block[1:5, 1] <- 0.3
cv_block[1:5, 2] <- 0.5

cv_block[6:10, 1] <- 0.2
cv_block[6:10, 2] <- 0.3


catch_obs4 <- new("CatchObs",
                  catchID = "period_specific_CV",
                  title = "Different CV bounds by period",  # only in area 2 catch data are collected
                  areas = c(2),
                  catchYears = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),
                  #improving reporting rates over time
                  # systematic overreporting
                  reporting_rates = c(rep(1.1,10)), #length= length catchYears
                  obs_CVs = cv_block)

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueN4,
  CatchObsObj = catch_obs4,
  customToCluster = "testMP",
  wd = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs",
  fileName = "test_run_MP_with_catch_obs4",
  doPlot = TRUE
)

out4<-readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG/outputs", "test_run_MP_with_catch_obs4")
out4$HCR$decisionData
out4$HCR$decisionData$CPUE_1
out4$HCR$decisionData$CPUE_2
out4$HCR$decisionData$observed_catch
out4$HCR$decisionData$observed_catch_area_2


p1 <- plotIndex_tibble(out4$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/index_plot_CPUEN4.jpeg")
p1

p2 <- plotCatch_tibble (out4$HCR$decisionData,
                        save_plot = TRUE,
                        filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/plot_catchobs4.jpeg")
p2

#---------------------------------------------------------------------------------------#
# Observation model 5: Simulation of two Survey (biomass) index covering both areas     #
# and catch observation model
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
    indexYears = 1:total_years,               # annual data collection
    survey_timing = 0.1,                # near the beginning of the year
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
    survey_timing = 0.5,                #half year
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


#Populating CVs for catchobs object
cv_block <- matrix(nrow = 20, ncol = 2)
# define a CV for each year individualy
cv_block[1:20, 1] <- 0.1
cv_block[1:20, 2] <- 0.2


catch_obs5 <- new("CatchObs",
                  catchID = "One_CV",
                  title = "One CV",
                  areas = c(1,2),
                  catchYears = c(1:total_years),
                  #improving reporting rates over time
                  # 100% reporting
                  reporting_rates = c(rep(1,20)), #length= length catchYears
                  obs_CVs = cv_block)

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = SurveyB5,
  CatchObsObj = catch_obs5,
  customToCluster = "testMP",
  wd = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs",
  fileName = "test_run_MP_with_catch_obs5",
  doPlot = TRUE
)

out5<-readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG/outputs", "test_run_MP_with_catch_obs5")
out5$HCR$decisionData
out5$HCR$decisionData$Survey_1
out5$HCR$decisionData$Survey_2
out5$HCR$decisionData$observed_catch   # sum of both areas
out5$HCR$decisionData$observed_catch_area_1
out5$HCR$decisionData$observed_catch_area_2

p1 <- plotIndex_tibble(out5$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/index_plot_SurveyB5.jpeg")
p1

p2 <- plotCatch_tibble (out5$HCR$decisionData,
                        save_plot = TRUE,
                        filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/plot_catchobs5.jpeg")
p2


#---------------------------------------------------------------------------------------#
# Observation model 6: Simulation of CPUE (FD) + Survey (FI)  index covering both areas (biomass)    #
# catch observation model
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
    indexYears = 1:total_years,               # annual data collection
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
    survey_timing = 0.5,                # half year
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


#Populating CVs for catchobs object
cv_blocks <- matrix(nrow = 20, ncol = 2)
# define a CV for each year individualy
cv_blocks[1:10, 1] <- 0.25
cv_blocks[1:10, 2] <- 0.4
cv_blocks[11:15, 1] <- 0.2
cv_blocks[11:15, 2] <- 0.3
cv_blocks[16:20, 1] <- 0.1
cv_blocks[16:20, 2] <- 0.2


catch_obs6 <- new("CatchObs",
                  catchID = "Blocks_CV",
                  title = "Blocks CV",
                  areas = c(1,2),
                  catchYears = c(1:total_years),
                  #improving reporting rates over time
                  # overeporting
                  reporting_rates = c(rep(1.8,20)), #length= length catchYears
                  obs_CVs = cv_blocks)


# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = mixed_biomass,
  CatchObsObj = catch_obs6,
  customToCluster = "testMP",
  wd = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs",
  fileName = "test_run_MP_with_mixed_index_and_catch_obs6",
  doPlot = TRUE
)

out6<-readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG/outputs", "test_run_MP_with_mixed_index_and_catch_obs6")
out6$HCR$decisionData$CPUE_1
out6$HCR$decisionData$Survey_2
out6$HCR$decisionData$observed_catch



p1 <- plotIndex_tibble(out6$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/index_plot_mixedB6.jpeg")
p1

p2 <- plotCatch_tibble (out6$HCR$decisionData,
                        save_plot = TRUE,
                        filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/plot_catchobs6.jpeg")
p2




#simple function to organize data -  to create year / iteration table (if needed)

create_year_iteration_table <- function(data, index_col) {
  # get subset with non-NA values
  subset_data <- data[!is.na(data[[index_col]]), c("k", "j", index_col)]

  # get unique years and iterations
  sim_years <- sort(unique(subset_data$j))
  iterations <- sort(unique(subset_data$k))

  # convert to user years for labeling
  user_years <- sim_years - 1

  # create empty matrix
  result_matrix <- matrix(NA,
                          nrow = length(sim_years),
                          ncol = length(iterations),
                          dimnames = list(paste0("Year_", user_years),
                                          paste0("Iter_", iterations)))

  # fill the matrix
  for(i in 1:nrow(subset_data)) {
    year_pos <- which(sim_years  == subset_data$j[i])
    iter_pos <- which(iterations == subset_data$k[i])
    result_matrix[year_pos, iter_pos] <- subset_data[[index_col]][i]
  }

  return(result_matrix)
}

out6$HCR$decisionData$CPUE_1
out6$HCR$decisionData$Survey_2
out6$HCR$decisionData$observed_catch

cpue_data <- create_year_iteration_table(out6$HCR$decisionData, "CPUE_1")
survey_data <- create_year_iteration_table(out6$HCR$decisionData, "Survey_2")
catch_data <- create_year_iteration_table(out6$HCR$decisionData, "observed_catch")



####################################################################################
#####################New example including LC ######################################
####################################################################################

#---------------------------------------------------------------------------------------#
# Observation model 7: Simulation of CPUE (FD) + Survey (FI)  index covering both areas (biomass)    #
# catch observation model
# length composition obs model
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
    indexYears = 1:total_years,               # annual data collection
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
    survey_timing = 0.5,                #half year
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


#Populating CVs for catchobs object
cv_blocks <- matrix(nrow = 20, ncol = 2)
# define a CV for each year individualy
cv_blocks[1:10, 1] <- 0.25   # position 1-10: User Years 1-10 (Historical): CV 0.25-0.4
cv_blocks[1:10, 2] <- 0.4
cv_blocks[11:15, 1] <- 0.2
cv_blocks[11:15, 2] <- 0.3
cv_blocks[16:20, 1] <- 0.1
cv_blocks[16:20, 2] <- 0.2


catch_obs7 <- new("CatchObs",
                  catchID = "Blocks_CV",
                  title = "Blocks CV",
                  areas = c(1,2),
                  catchYears = c(1:total_years),
                  #improving reporting rates over time
                  # overeporting
                  reporting_rates = c(rep(1.8,20)), #length= length catchYears
                  obs_CVs = cv_blocks)

#NEW addition

# create length composition object
length_comp7 <- new("LCompObs")
length_comp7@indexID <- "Mixed_LC_FD_FI"
length_comp7@title <- "Mixed FD and FI Length Composition"
#length_comp7@length_cv <- 0.10         # 10% CV for individual length variability
length_comp7@length_bin_width <- 1     # 1cm bins

# define survey design for length composition
length_comp7@survey_design <- list(

  # Fishery-dependent length composition (uses catch data)
  list(
    indextype = "FD",
    areas = c(1, 2),                     # Both areas (same as CPUE)
    years = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),  # every other year
    sample_sizes = c(50, 60, 70, 80, 90, 100, 110, 120, 130, 140)  # increasing sample sizes
  ),

  # Fishery-independent survey length composition - REUSING previous selectivity
  list(
    indextype = "FI",
    areas = c(1, 2),                     # both areas (same as Survey)
    years = c(6, 8, 10, 12, 14, 16, 18, 20),            # every 2 years
    sample_sizes = c(200, 250, 300, 350, 300, 350, 300, 350),  # large sample sizes
    selectivity_hist_idx = 1,            # REUSES survey_sel1_hist (same as  Index)
    selectivity_proj_idx = 1,             # REUSES survey_sel1_proj (same as  Index)
    survey_timing = 0.5      # half of year
  )
)

# REUSE the same selectivity objects defined for the Index
length_comp7@selectivity_hist_list <- list(survey_sel1_hist)  # SAME as mixed_biomass
length_comp7@selectivity_proj_list <- list(survey_sel1_proj)  # SAME as mixed_biomass


# Now run the simulation with the three observation models
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = mixed_biomass,
  CatchObsObj = catch_obs7,
  LengthCompObj = length_comp7,      # NEW: Length composition object
  customToCluster = "testMP",
  wd = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs",
  fileName = "test_run_MP_with_mixed_index_obs_catch_obs_length_obs",
  doPlot = TRUE
)

out7<-readProjection("P:/Fork_fish_Sim_GTG/fishSimGTG/outputs", "test_run_MP_with_mixed_index_obs_catch_obs_length_obs")
out7$HCR$decisionData$CPUE_1
out7$HCR$decisionData$Survey_2
out7$HCR$decisionData$observed_catch
out7$HCR$decisionData$Fishery_1_count_bin_10
out7$HCR$decisionData$Fishery_1_count_bin_13
out7$HCR$decisionData$Survey_2_count_bin_12
out7$HCR$decisionData$Survey_2_count_bin_18


p1 <- plotIndex_tibble(out7$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/index_plot_mixedB7.jpeg")
p1

p2 <- plotCatch_tibble (out7$HCR$decisionData,
                        save_plot = TRUE,
                        filename = "P:/Fork_fish_Sim_GTG/fishSimGTG/outputs/plot_catchobs7.jpeg")
p2


#plotting length composition
plotLengthComp_distributions <- function(tibble_data, program_name,
                                         years_to_plot = NULL,
                                         show_iterations = TRUE,
                                         save_plot = FALSE,
                                         filename = "length_comp_distributions.jpeg") {

  # get number of length bins and basic info
  n_bins <- unique(tibble_data$n_length_bins)[1]
  length_bin_width <- unique(tibble_data$length_bin_width)[1]
  title <- unique(tibble_data$title)[1]

  # get years with data
  sample_size_col <- paste0(program_name, "_sample_size")
  years_with_data <- sort(unique(tibble_data$j[!is.na(tibble_data[[sample_size_col]])]))

  #convert to user years
  user_years_with_data <- years_with_data - 1

  #which years to plot
  if(is.null(years_to_plot)) {
    years_to_plot <- user_years_with_data  # all years with data
  }

  # convert user years back to simulation years for data extraction
  sim_years_to_plot <- years_to_plot + 1

  #create data frame for plotting
  plot_data <- data.frame()

  for(sim_year in sim_years_to_plot) {
    if(!sim_year %in% years_with_data) {
      cat("No data for year", sim_year-1, "(user year)\n")
      next
    }

    #get data for this year
    year_data <- tibble_data[tibble_data$j == sim_year, ]

    if(nrow(year_data) == 0) next

    #extract number for each iteration
    for(iter in unique(year_data$k)) {
      iter_data <- year_data[year_data$k == iter, ]

      if(nrow(iter_data) == 0) next

      #get length composition numbers
      numbers <- numeric(n_bins)
      for(bin in 1:n_bins) {
        col_name <- paste0(program_name, "_count_bin_", bin)
        if(col_name %in% names(iter_data) && !is.na(iter_data[[col_name]][1])) {
          numbers[bin] <- iter_data[[col_name]][1]
        }
      }

      # include if there is actual data (not all zeros)
      if(sum(numbers) > 0) {
        #create length bins (center of each bin)
        length_bins <- seq(length_bin_width/2, n_bins * length_bin_width - length_bin_width/2, by = length_bin_width)

        iteration_df <- data.frame(
          length_bin = length_bins,
          numbers = numbers,
          user_year = sim_year - 1,  #convert to user year
          iteration = iter,
          year_label = paste("Year", sim_year - 1)
        )

        plot_data <- rbind(plot_data, iteration_df)
      }
    }
  }

  #calculate median across iterations for each year and length bin
  median_data <- plot_data %>%
    group_by(user_year, length_bin, year_label) %>%
    summarise(
      median_numbers = median(numbers, na.rm = TRUE),
      .groups = "drop"
    )

  #convert year label to ordered factor (so the years are shown organizd)
  plot_data$year_label <- factor(plot_data$year_label,
                                 levels = paste("Year", sort(unique(plot_data$user_year))))
  median_data$year_label <- factor(median_data$year_label,
                                   levels = paste("Year", sort(unique(median_data$user_year))))



  #create the plot
  p <- ggplot()

  if(show_iterations) {
    #ndividual iterations as thin lines
    p <- p + geom_line(data = plot_data,
                       aes(x = length_bin, y = numbers, group = interaction(user_year, iteration)),
                       color = "lightblue", alpha = 0.3, size = 0.3)
  }

  #median line
  p <- p + geom_line(data = median_data,
                     aes(x = length_bin, y = median_numbers),
                     color = "darkblue", size = 1.2) +
    geom_point(data = median_data,
               aes(x = length_bin, y = median_numbers),
               color = "darkblue", size = 0.8) +
    facet_wrap(~ year_label, scales = "free_y") +
    labs(
      title = paste("Length Composition Distributions -", program_name),
      subtitle = paste("Data from:", title),
      x = "Length (cm)",
      y = "Numbers"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  if(save_plot) {
    ggsave(filename, plot = p, width = 12, height = 8, dpi = 300)
    cat("Length composition plot saved as:", filename, "\n")
  }

  return(p)
}

p1 <- plotLengthComp_distributions(
  out7$HCR$decisionData,
  program_name = "Fishery_1",
  years_to_plot = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),  # User years
  show_iterations = TRUE,
  save_plot = TRUE,
  filename = "FD_length_distributions.jpeg"
)
p1

p2 <- plotLengthComp_distributions(
  out7$HCR$decisionData,
  program_name = "Survey_2",
  years_to_plot = c(6, 8, 10, 12, 14, 16, 18, 20),  # User years
  show_iterations = TRUE,
  save_plot = TRUE,
  filename = "FI_length_distributions.jpeg"
)
p2








