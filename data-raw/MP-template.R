
templateMP <- function(phase, dataObject) {

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
  yr <- j - TimeAreaObj@historicalYears - 1

  #To obtain the terminal year of the historical period
  yrHist <- TimeAreaObj@historicalYears + 1


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
    return(list()) #e.g., return(list(year = j, iteration = k, CPUE = 1.4))
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
