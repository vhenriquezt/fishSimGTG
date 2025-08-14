


#----------------------------------------------------------------------
#Fixed time period - for fixed harvest time period (not for MSE time period)
#----------------------------------------------------------------------

#Roxygen header
#'Historical fishing pressure
#'
#' @param phase Management procedures are coded in three phases: 1 - data collection, 2 - a decision making process, 3 - conversion of that process into annual F
#' @param dataObject The needed inputs to the management procedure
#' @export

# new addition: updated fixedStrategy function for multifleet support
fixedStrategy<-function(phase, dataObject){

  #Unpack dataObject (adding MultifleetObj)
  j <- areas <- k <- TimeAreaObj <- is <- histEffortDev <- MultifleetObj <- NULL
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  #new addition: detecting multifleet mode
  is_multifleet <- !is.null(MultifleetObj) && MultifleetObj@nfleets > 1

  if(is_multifleet) {
    nfleets <- MultifleetObj@nfleets
    fleet_proportions <- MultifleetObj@fleet_proportions
    F_eq_by_fleet <- is$F_by_fleet
  } else {
    nfleets <- 1
    fleet_proportions <- c(1.0)
    F_eq_by_fleet <- c(is$Feq)
  }


  #Booking keeping for year for items in TimeAreaObj
  # e.g., simulation year j=3 needs historical effort from year yr=2
  yr <- j - 1  #the simulation year (starting from 2 in time dynamics) and yr is the index for historical effort arrays (starting from 1)

  if(phase==3){

    if(is_multifleet) {

      Flocal<-data.frame()

      for (m in 1:areas) {
        #total F for this area (will be distributed among fleets)- this approach for now, for testing
        total_F_area <- 0
        for(f in 1:nfleets) {
          # fleet-specific F calculation
          F_fleet <- F_eq_by_fleet[f] * TimeAreaObj@historicalEffort[yr,m] * histEffortDev[j,k,m,f]
          total_F_area <- total_F_area + F_fleet
          #add fleet-specific row (for tracking purposes)
          Flocal<-rbind(Flocal, c(j, k, m, f, F_fleet))
        }
        # add total F row (for backward compatibility)
        Flocal<-rbind(Flocal, c(j, k, m, 0, total_F_area))  # fleet=0 indicates total
      }

      colnames(Flocal) <- c("year", "iteration", "area", "fleet", "Flocal")
      return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3],
                  fleet=Flocal[,4], Flocal=Flocal[,5]))

  } else {

    #otherwise continue with original implementation
    #Create a temp data frame of fishing mortalities by area
    Flocal<-data.frame()
    for (m in 1:areas){
    Flocal<-rbind(Flocal, c(j, k, m, TimeAreaObj@historicalEffort[yr,m]*is$Feq*histEffortDev[j,k,m]))
    }
    return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3],  Flocal=Flocal[,4]))
  }
}
}

#-------------------------------------------------------------------
#Projection modeling - no harvest control rule, simple projections
#-------------------------------------------------------------------

#Roxygen header
#' Static projections of combinations of effort changes, bag limit, and/or spatial closures. Also used in combination with ProFisheryObj (e.g. size limit change) to project temporal dynamics of size limits.
#'
#' The Strategy object should be specified as follows. (1) Set Strategy@projectionYears to the number of forward projection years you wish to simulate. (2) Strategy@projectionName = "projectionStrategy".
#' (3)  Strategy@projectionParams should be a list with two items. First items is a vector of length areas containing bag limit. For no bag limit use -99.
#' The bag limit should be thought of as take per unit time (e.g. day) and basically acts like a CPUE threshold.
#' The effect of bag limit is calculated against the historical CPUE (e.g. in same units of take per unit time) in the Stochastic object historicalCPUE. Make sure that the bag limit and historicalCPUE are consistent with historicalCPUEType (e.g., biomass or abundance (numbers)) and this parameter is used in determing the effect of fish biomass or abundance on CPUE.
#' The second item in the list is a matrix of nrows = projectionYears and ncols = areas that contains value multipiers of initial equilibrium fishing effort. This allows projection of effort reduction and of marine reserves via setting effort to 0.
#' @param phase Management procedures are coded in three phases: 1 - data collection, 2 - a decision making process, 3 - conversion of that process into annual F
#' @param dataObject The needed inputs to the management procedure
#' @importFrom stats dpois ppois
#' @export

projectionStrategy<-function(phase, dataObject){

  #Unpack dataObject
  j <- TimeAreaObj <- areas <- StrategyObj <- is <- k <- StochasticObj <- lh <- N <- selHist <- Cdev <- Edev <- selGroup <- Ftotal <- NULL
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  #Book keeping year for items in StrategyObj
  yr <- j - TimeAreaObj@historicalYears - 1

  #Book keeping year for terminal year of historical period, if available
  yrHist <- TimeAreaObj@historicalYears + 1

  if(phase==3){

    Flocal<-data.frame()
    for (m in 1:areas){

      bag <- StrategyObj@projectionParams[['bag']][m]

      if(bag == -99){

        #Apply Flocal
        Ftmp<-StrategyObj@projectionParams[['effort']][yr,m]*Ftotal[yrHist,k,m]*Edev[k]
        Flocal<-rbind(Flocal, c(j, k, m, Ftmp))

      } else {

        if(StrategyObj@projectionParams[['CPUEtype']] == "retN") {

          #Initial equilibrium vulnerable N
          Nvul<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist[[m]]$keep[[x]])))

          #Specify assumed initial lambda
          lambdaInitial <- Cdev[k]

          #Solove for q
          q<-lambdaInitial/Nvul

          #Get current F multiplier
          lambda<-q*sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup[[m]]$keep[[x]])))
          probs<-c(dpois(0:(bag-1), lambda), 1-ppois(bag-1,lambda))
          probs<-probs/sum(probs)
          nm<-sum(0:bag*probs)
          Fmult<-min(nm/lambda, 1.0)

          #Apply Flocal
          Ftmp<-StrategyObj@projectionParams[['effort']][yr,m]*Ftotal[yrHist,k,m]*Fmult*Edev[k]
          Flocal<-rbind(Flocal, c(j, k, m, Ftmp))

        }

        if(StrategyObj@projectionParams[['CPUEtype']] == "retB") {

          #Initial equilibrium vulnerable N
          Nvul<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist[[m]]$keep[[x]]*lh$W[[x]])))

          #Specify assumed initial lambda
          lambdaInitial <- Cdev[k]

          #Solove for q
          q<-lambdaInitial/Nvul

          #Get current F multiplier
          lambda<-q*sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup[[m]]$keep[[x]]*lh$W[[x]])))
          probs<-c(dpois(0:(bag-1), lambda), 1-ppois(bag-1,lambda))
          probs<-probs/sum(probs)
          nm<-sum(0:bag*probs)
          Fmult<-min(nm/lambda, 1.0)

          #Apply Flocal
          Ftmp<-StrategyObj@projectionParams[['effort']][yr,m]*Ftotal[yrHist,k,m]*Fmult*Edev[k]
          Flocal<-rbind(Flocal, c(j, k, m, Ftmp))

        }
      }
    }
    return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3],  Flocal=Flocal[,4]))
  }
}

#--------------------------------------------------------------------------
#Constant catch - TAC based on average catch of final x years of historical
#--------------------------------------------------------------------------

#Roxygen header
#' MP for setting a TAC based on average catch of final x years of historical.
#'
#' The Strategy object should be specified as follows. (1) Set Strategy@projectionYears to the number of forward projection years you wish to simulate. (2) Strategy@projectionName = "constantCatchStrategy".
#' (3)  Strategy@projectionParams should be a named list with one item called aveYrs, where aveYrs is the number of years to average catch across, including and proceeding the final historical year.
#' @param phase Management procedures are coded in three phases: 1 - data collection, 2 - a decision making process, 3 - conversion of that process into annual F
#' @param dataObject The needed inputs to the management procedure
#' @export

constantCatchStrategy<-function(phase, dataObject){

  #Unpack dataObject
  j <- TimeAreaObj <- areas <- StrategyObj <- is <- k <- StochasticObj <- lh <- N <- selHist <- Cdev <- Edev <- selGroup <- Ftotal <- NULL
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  ########
  if(phase==1){
    #User defines computations needed.
    combined_data <- list()

    # get catch obs data and add to combined_data list
    if(!is.null(CatchObsObj)) {
      catch_result <- calculate_single_CatchObs(dataObject)
      # add all catch columns
      for(col_name in names(catch_result)) {
        combined_data[[col_name]] <- catch_result[[col_name]]
      }
    }
    return(combined_data)
  }

  ########
  #Phase 2 - OPTIONAL. Analysis, decision-rule and/or HCR
  if(phase==2){

    #Book keeping year for terminal year of historical period, if available
    yrHist <- TimeAreaObj@historicalYears + 1
    frYear <- TimeAreaObj@historicalYears + 2 - StrategyObj@projectionParams[['aveYrs']]

    TACvec<-vector()
    for(m in 1:areas){
      TACvec[m]<-mean(decisionData[[paste0("observed_catch_area_", m)]][which(decisionData$k == k & decisionData$j >= frYear & decisionData$j <= yrHist)])
    }

    #User defines variable names for returned list, as this info will be used in Phase 3
    return(list(year=rep(j, areas), iteration=rep(k,areas), area=1:areas, TAC=TACvec, frYear = rep(frYear, areas), yrHist = rep(yrHist, areas)))
  }

  ########
  #Phase 3
  if(phase==3){
    return(solveFfromTAC(dataObject))
  }
}


#----------------------------------------
#Calculate F from TAC - utility function
#----------------------------------------
#Roxygen header
#' Calculate F from TAC - utility function.
#' @param dataObject The needed inputs to the management procedure
#' @export

# need an example here, removed to create documentation
solveFfromTAC <- function(dataObject){
  #Unpack dataObject
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  Flocal<-vector()
  for(m in 1:areas){
    TACtmp<-decisionAnnual$TAC[decisionAnnual$year==j & decisionAnnual$iteration==k & decisionAnnual$area==m]
    data<-list(N=N, selGroup=selGroup, lh = lh, TACtmp=TACtmp)
    min.RSS<-function(logFmort, data){
      Fmort<-exp(logFmort)
      catchBtmp <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Fmort*selGroup[[m]]$keep[[x]]/(Fmort*selGroup[[m]]$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Fmort*selGroup[[m]]$removal[[x]] - lh$LifeHistory@M))*N[[x]][,j,m])))
      (catchBtmp-data$TACtmp)^2
    }
    Flocal[m]<-ifelse(data$TACtmp==0, 0.0, exp(optimize(min.RSS, lower = -14, upper = 1.1, maximum=FALSE, data=data)$minimum))
  }
  return(list(year=rep(j,areas), iteration=rep(k,areas), area=1:areas, Flocal=Flocal))
}

