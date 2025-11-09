


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

  #Bill edit: re-wrote most of this function

  #Unpack dataObject (adding MultifleetObj)
  j <- areas <- k <- TimeAreaObj <- is <- histEffortDev <- MultifleetObj <- NULL
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  #Booking keeping for year for items in TimeAreaObj
  yr <- j - 1


  if(phase==3){

    #Single fleet
    if(!is_multifleet){
      #Create a temp data frame of fishing mortalities by area
      Flocal<-data.frame()
      for (m in 1:areas) Flocal<-rbind(Flocal, c(j, k, m, TimeAreaObj@historicalEffort[yr,m]*is$Feq*histEffortDev[j,k,m]))
      return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3],  Flocal=Flocal[,4]))
    }

    #Multi-fleet
    if(is_multifleet){
      Flocal<-data.frame()
      for (m in 1:areas){
        for(f in 1:nfleets){
          Flocal<-rbind(Flocal, c(j, k, m, f, MultifleetObj@fleet_historicalEffort[yr,m,f]*is$F_by_fleet[f]*histEffortDev[j,k,m,f]))
        }
      }
      return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3], fleet = Flocal[,4],  Flocal=Flocal[,5]))
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

#' Function to convert Total Allowable Catch (TAC) targets into F rates for multiple fleets using Newton-Raphson iteration, adapted from OpenMSE and https://api.admb-project.org/baranov_8cpp_source.html code for fishSimGTG
#'
#' @param j Current simulation year
#' @param k Current iteration
#' @param TAC_targets Vector of TAC targets [fleet1_TAC, fleet2_TAC, ...]. NA = effort-managed
#' @param N List of abundance arrays by GTG: N[[gtg]][age, year, area]
#' @param lh Life history object from LHwrapper
#' @param selGroup Selectivity object: selGroup[[area]][[fleet]] (multifleet) or selGroup[[area]] (single fleet)
#' @param M_rate Natural mortality rate (scalar)
#' @param effort_F_by_fleet Vector of F values for effort-managed fleets (length = nfleets)
#' @param is_multifleet Logical indicating multifleet mode
#' @param areas Number of areas
#' @param nfleets Number of fleets (1 for single fleet)
#' @param TAC_type Character: "keep" (retained catch) or "removals" (total removals)
#' @param control List with algorithm settings (maxiterF, tolF)
#' @return Vector of F values [fleet1_F, fleet2_F, ...] or single value for single fleet
#' @export

solveTAC_to_F_fishSimGTG <- function(j, k, TAC_targets, N, lh, selGroup, M_rate,
                                     effort_F_by_fleet, is_multifleet, areas, nfleets,
                                     TAC_type = "keep", control = NULL) {

  #control parameters
  maxiterF <- 300 #same as openMSE
  tolF <- 1e-4    #same as openMSE
  tiny <- 1e-15   #same as openMSE (if catch is super tiny, do not need to iterate (same as openMSE))

  #some validation
  # if control NULL - all the rest skipped
  if(!is.null(control)) {
    if (!is.null(control$maxiterF) && is.numeric(control$maxiterF)) {
      maxiterF <- as.integer(control$maxiterF)
    }
    if (!is.null(control$tolF) && is.numeric(control$tolF)) {
      tolF <- control$tolF
    }
  }

  if (length(TAC_targets) != nfleets) {
    stop("TAC_targets length must equal nfleets")
  }

  if (!TAC_type %in% c("keep", "removals")) {
    stop("TAC_type must be 'keep' or 'removals'")
  }

  cat("SOLVER DEBUG: areas =", areas, "nfleets =", nfleets, "\n")
  cat("SOLVER DEBUG: TAC_targets =", TAC_targets, "\n")
  cat("SOLVER DEBUG: Sample N (GTG 1, Age 5, Area 1) =", N[[1]][5, j, 1], "\n")
  cat("SOLVER DEBUG: Sample selectivity (GTG 1, Age 5) =",
      if(is_multifleet) selGroup[[1]][[1]]$keep[[1]][5] else selGroup[[1]]$keep[[1]][5], "\n")


  ct <- TAC_targets  # TAC targets for current iteration


  #calculating initial guesses - matching GTG structure
  ft <- sapply(1:nfleets, function(f) {
    if(is.na(ct[f])) return(effort_F_by_fleet[f])  #changed - use predetermined F to work when effort based

    #calculate total vulnerable biomass
    total_vuln_biomass <- 0

    for (gtg in 1:lh$gtg) {
      for (age in 1:lh$ageClasses) {
        for (area in 1:areas) {
          #access selectivity correctly based on fleet structure
          if (is_multifleet) {
            selectivity <- selGroup[[area]][[f]]$vul[[gtg]][age]
          } else {
            selectivity <- selGroup[[area]]$vul[[gtg]][age]
          }

          #add this GTG-age-area combination to total
          total_vuln_biomass <- total_vuln_biomass +
            N[[gtg]][age, j, area] * lh$W[[gtg]][age] * selectivity
        }
      }
    }

    if (total_vuln_biomass < tiny) return(tiny)


    # Initial guess: F ~ TAC / VulnerableBiomass
    # CHANGED:
    # Removed min(1, ...) cap that prevented TAC adjustment from triggering (SEE BELOW WHY - following Martell example)
    # BEFORE: guess <- min(1, ct[f] / total_vuln_biomass) #Calculate intitial F, but if it's greater than 1.0, cap it at 1.0
    # NOW: Allow F to be calculated without artificial ceiling

    guess <-  ct[f] / total_vuln_biomass

    # Diagnostic
    cat(sprintf("SOLVER: Fleet %d, TAC=%.2f, VulnBiomass=%.2f, Initial_F=%.6f\n",
                f, ct[f], total_vuln_biomass, guess))

    return(guess)
  })

  #store initial guess before Newton-Raphson iterations for diagnostics
  #BEFORE ANY ADJUSTMENT
  ft_initial <- ft

  #Show initial guesses
  cat("\n=== show intitial guesses ===\n")
  for(f in 1:nfleets) {
    if(!is.na(ct[f])) {
      cat(sprintf("Fleet %d: TAC=%.2f, Initial F=%.4f\n",
                  f, ct[f], ft[f]))
    }
  }

  # ============================================================================
  # TAC ADJUSTMENT FOR BIOLOGICAL REALISM (FOLLOWING MARTELL CODE'S IDEA)
  # ============================================================================
  # If initial F exceeds biological limit, reduce TAC to what's achievable at the limit

  # PROBLEM: If TAC is too high for available biomass, the solver tries to
  #          achieve it by increasing F indefinitely, hitting the biological cap
  #          (F=3.0) and staying there, causing population collapse.

  # SOLUTION: Following Martell's ADMB approach: if initial F exceeds the biological limit, we:
  #           1. Cap F at the biological limit (F = 3.0)
  #           2. Calculate what catch (new target) is ACHIEVABLE at this F
  #           3. ADJUST THE TAC DOWNWARD to this achievable level
  #           4. Solve for the adjusted (realistic) TAC

  # RESULT: F converges at< = 3.0) instead of staying stuck
  #         at the cap trying to achieve an impossible target.

  # This is equivalent to Martell code (lines 117-129 in ADMB version):
  #   if(1.-ft(i)<minsurv) {            // If survival < 5%
  #     ft(i)=1.-minsurv;               // Cap F (exploitation rate)
  #     ctmp(i)=ft(i)*ba*V(i)*...;      // Recalculate achievable catch
  #   }
  #   ct=ctmp;                          // Use adjusted TAC as target

  max_F_bio <- 3.0  # this is my F cap (bio limit) - (roughly 95% exploitation, 5% survival)

  for(f in 1:nfleets) {
    if(!is.na(ct[f]) && ft[f] > max_F_bio) {  # Is 1.0 > 3.0?  NO! so this block does not work wit the cap=1 in initial guess
     #if this block does not get executed:
     #I will do the following example:
      # initial_F = 800 / 200 = 4.0
      # ft[f] = min(1, 4.0) = 1.0  # F is capped at 1.0
      # if (ft[f] > 3.0) { # Is 1.0 > 3.0?  NO!
      #   # This block NEVER EXECUTES
      #   # TAC stays at 800 kg
      #   }
      # Then NR will start at F=1 with a target TAC (800) that is imposible to achive
      # ANd NR will continue forever trying to achieve the 800 TAC
      cat(sprintf("\n Fleet %d: Initial F=%.4f exceeds F biological limit (%.2f)\n",
                  f, ft[f], max_F_bio))
      cat(sprintf("    Original TAC: %.2f kg\n", ct[f]))

      #Cap F at bio limit (like Martell's: ft(i)=1.-minsurv)
      ft[f] <- max_F_bio

      #Calculate achievable catch at F = max_F_bio using Baranov equation (like Martell's: ctmp(i)=ft(i)*ba*V(i)*...)
      achievable_catch <- 0

      for (gtg in 1:lh$gtg) {
        for (age in 1:lh$ageClasses) {
          for (area in 1:areas) {

            # Get selectivity based on TAC type
            if (is_multifleet) {
              sel <- if(TAC_type == "keep") selGroup[[area]][[f]]$keep[[gtg]][age]
              else selGroup[[area]][[f]]$removal[[gtg]][age]
            } else {
              sel <- if(TAC_type == "keep") selGroup[[area]]$keep[[gtg]][age]
              else selGroup[[area]]$removal[[gtg]][age]
            }

            # Calculate biomass for this GTG-age-area combination
            biomass <- N[[gtg]][age, j, area] * lh$W[[gtg]][age]

            # Calculate Z with capped F
            Z_approx <- M_rate + ft[f] * sel
            Z_approx <- max(Z_approx, tiny)

            # Baranov equation: calculate catch at F=max_F_bio
            achievable_catch <- achievable_catch +
              ft[f] * sel / Z_approx * (1 - exp(-Z_approx)) * biomass
          }
        }
      }

      cat(sprintf("    Achievable catch at F=%.2f: %.2f kg\n", max_F_bio, achievable_catch))
      cat(sprintf("    Reducing TAC to %.1f%% of original\n", 100*achievable_catch/ct[f]))

      # ADJUST TAC DOWNWARD to achievable level (like Martell's: ct=ctmp)
      # This is the KEY modification - here we change the target instead of forcing high F
      ct[f] <- achievable_catch

      warning(sprintf("Fleet %d: TAC reduced to %.2f kg due to biological constraints",
                      f, achievable_catch))
    }
  }

  # Print adjusted TAC targets that will be used in Newton-Raphson

  cat("\n=== Adjusted TAC targets (used in Newton-Raphson) ===\n")
  for(f in 1:nfleets) {
    if(!is.na(ct[f])) {
      cat(sprintf("Fleet %d: TAC=%.2f, Starting F=%.4f\n",
                  f, ct[f], ft[f]))
    }
  }

  # ============================================================================
  # END MARTELL AMBD CODE APPROACH)
  # ============================================================================



  # Early termination check: # If TACs or F values are negligible,
  # no fishing occurs - skip iteration

  #if NA ct or ct too tiny or ft too tiny - early termination check
  if (all(ct[!is.na(ct)] <= tiny) || all(ft <= tiny)) {
    return(if (!is_multifleet && nfleets == 1) tiny else rep(tiny, nfleets))
  }

  #Initialize variables
  pct <- numeric(nfleets)  # Predicted catches (calculated each iteration)
  dct <- numeric(nfleets)  # Derivatives: dCatch/dF (calculated each iteration)

  #Newton-Rapson iteration
  converged <- FALSE        # Convergence - Non convergence flag
  iteration_count <- 0      # Track number of iterations used



  # #initializing arrays: [age, area, fleet]
  # Fmat <- Fmat_ret <- predC <- array(NA, c(lh$ageClasses, areas, nfleets))
  # dct <- pct <- numeric(nfleets)  # Derivatives and predicted catches

  #====================here we start NR calculations=======================#
  #Method: F_new = F_guess - (error / derivative)
  #error = predicted_catch - target_catch
  #derivative = dCatch/dF

  for(iter in 1:maxiterF) {
    iteration_count <- iter

    # diagnostic output (first 3 iterations and every 50th)
    if(iter <= 3 || iter %% 50 == 0) {
      cat(sprintf("\n--- Iteration %d ---\n", iter))
      for(f in 1:nfleets) {
        if(!is.na(ct[f])) {
          cat(sprintf("  Fleet %d: F=%.4f\n", f, ft[f]))
        }
      }
    }

    #reset predicted catches and derivatives for this iter
    pct[] <- 0
    dct[] <- 0

    #calculate predicted catches and derivatives using fishSimGTG approach
    # loop through all GTGs, ages, and areas to calculate:
    # Total predicted catch for each fleet using current F values
    # Derivative of catch with respect to F for each fleet

    for (gtg in 1:lh$gtg) {
      for (age in 1:lh$ageClasses) {
        for (area in 1:areas) {

          #calculate GTG-specific total Z
          # Z = M + sum_across_fleets(F_fleet * selectivity_fleet)
          total_F_mortality <- 0
          for (ff in 1:nfleets) {
            if (is.na(ct[ff])) {
              #effort-managed fleet: use predetermined F
              fleet_F <- effort_F_by_fleet[ff]
            } else {
              #TAC-managed fleet: use current F guess
              fleet_F <- ft[ff]
            }

            #access removal selectivity
            if (is_multifleet) {
              removal_sel <- selGroup[[area]][[ff]]$removal[[gtg]][age]
            } else {
              removal_sel <- selGroup[[area]]$removal[[gtg]][age]
            }

            total_F_mortality <- total_F_mortality + fleet_F * removal_sel
          }

          Z_gtg <- M_rate + total_F_mortality
          Z_gtg <- max(Z_gtg, tiny)  #avoid division by zero

          #calculate biomass for this GTG-age-area combination
          biomass_gtg <- N[[gtg]][age, j, area] * lh$W[[gtg]][age]

          #skip if no biomass
          if (biomass_gtg < tiny) next  #skip if no biomass

          #calculate fleet-specific catches and derivatives
          for (f in 1:nfleets) {

            #skip effort-managed fleets (F is predetermined, not solved)
            if (is.na(ct[f])) next

            #get appropriate selectivity based on TAC type
            #TAC Types:
            #"keep" = landed catch
            #"removals" = total mortality (including dead discards)

            #keep Selectivity: matches selGroup[[m]]$keep[[l]] for retained catch
            if (TAC_type == "keep") {
              if (is_multifleet) {
                sel_f <- selGroup[[area]][[f]]$keep[[gtg]][age]
              } else {
                sel_f <- selGroup[[area]]$keep[[gtg]][age]
              }
            } else {  # removals
              #removal Selectivity: matches selGroup[[m]]$removal[[l]] for total removals
              if (is_multifleet) {
                sel_f <- selGroup[[area]][[f]]$removal[[gtg]][age]
              } else {
                sel_f <- selGroup[[area]]$removal[[gtg]][age]
              }
            }

            # Baranov catch equation for this GTG (following fishSimGTG structure)
            # C = F*sel/Z * (1-exp(-Z)) * N * W
            catch_component <- ft[f] * sel_f / Z_gtg * (1 - exp(-Z_gtg)) * biomass_gtg

            #accumulate if finite (valid calculation)
            if (is.finite(catch_component)) {
              #add to fleet totals
              pct[f] <- pct[f] + catch_component

              # Calculate derivative component
              # d(C)/d(F) = sel/Z * (1-exp(-Z)) * biomass - F*sel/Z^2 * (1-exp(-Z)) * biomass + F*sel/Z * exp(-Z) * sel * biomass
              exp_neg_Z <- exp(-Z_gtg)
              one_minus_exp <- 1 - exp_neg_Z

              #the three derivative terms - see my excel
              term1 <- sel_f / Z_gtg * one_minus_exp * biomass_gtg
              term2 <- ft[f] * sel_f / (Z_gtg^2) * one_minus_exp * biomass_gtg * sel_f
              term3 <- ft[f] * sel_f / Z_gtg * exp_neg_Z * sel_f * biomass_gtg

              derivative_component <- term1 - term2 + term3

              # accumulate if finite
              if (is.finite(derivative_component)) {
                #add to fleet totals
                dct[f] <- dct[f] + derivative_component
              }
            }
          }#fleet
        }#areas
      }#age classess
    }#gtg

    #identify TAC-managed fleets
    #added: to work when effort based - remember: effort-managed fleets (TAC = NA) use predetermined F values
    tac_managed <- !is.na(ct)


    #check convergence and update and numerical issues
    #before using derivatives, verify they're valid
    #invalid derivatives indicate numerical instability

    if (any(!is.finite(dct[tac_managed]))) {
      warning("Derivative became infinite or NaN. Stopping.")
      converged <- FALSE
      break
    }

    #check if derivatives are too small (indicates flat catch-F relationship)
    #very small derivatives can cause enormous Newton-Raphson steps
    #Rememeber: step = (error/derivative)

    if (any(abs(dct[tac_managed]) < 1e-6)) {
      warning("Derivative became too small (NR unstable) Stop iteration.")
      converged <- FALSE  # Mark as FAILED now, not converged
      break               #exist the loop inmediately- this prevent F explosion
    }


    #Calculate error and check steps
    error <- pct - ct      # Error: predicted catch - target catch

    #check if NR step would be unreasonably large
    #large steps indicate numerical issues (instability)

    potential_steps <- abs(error[tac_managed] / dct[tac_managed])
    if (any(potential_steps > 5, na.rm = TRUE)) {
      warning(sprintf("Newton-Raphson step would be too large (%.2f). Stopping.",
                      max(potential_steps, na.rm = TRUE)))
      converged <- FALSE
      break
    }



    #Newton-Raphson update with damping factor

    # The 0.8 factor is a damping parameter that makes each step smaller
    # This helps prevent overshooting and improves stability
    # Formula: F_new = F_old - 0.8 × (error / derivative)
    # The 0.8 multiplies the step size, making updates more conservative
    # Basically it avoids huhe jumps in F
    #OBS: before that was incorrecly applied only to the derivative (as openMSE)
    #but that was making steps larger instead of smaller (more overshooting and probab;y less stable)


    #changed: only update TAC-managed fleets
    for(f in 1:nfleets) {
      if(tac_managed[f]) {  # Only if TAC-managed
        #ft[f] <- ft[f] - error[f] / (0.8 * dct[f]) #incorrect damping effect
        #ft[f] <- ft[f] - error[f] / dct[f]         #removing damping effect
        ft[f] <- ft[f] - 0.8*(error[f] / dct[f])    #placing damping correctly to effectivelity control the step
      }
      #effort-managed fleets: ft[f] stays unchanged
    }

    #ensure F values remain positive
    ft <- pmax(ft, tiny)

    #cap F at bio limit (F = 3) to prevent unrealistic fishing mortality
    ft[tac_managed] <- pmin(ft[tac_managed], 3)

    # Warning if F approaches the cap during iterations
    # this suggests TAC may still be too high despite adjustment
    if (any(ft[tac_managed] > 2.5)) {
      cat(sprintf("\n  WARNING at iter %d: F approaching biological limit\n", iter))
      for(f in which(tac_managed & ft > 2.5)) {
        cat(sprintf("  Fleet %d: F=%.4f (Target TAC=%.2f, Predicted=%.2f)\n",
                    f, ft[f], ct[f], pct[f]))
      }
      cat(" this suggests TAC may still be too high for available biomass\n")
    }

    #Check convergence - now after capping
    relative_error <- abs(error / pmax(ct, tiny))
    relative_error[!tac_managed] <- 0  # Zero out for effort fleets

    if (all(relative_error[tac_managed] < tolF)) {
      converged <- TRUE

      #debug
      cat(sprintf("\nCONVERGED at iteration %d\n", iter))
      for(f in 1:nfleets) {
        if(tac_managed[f]) {
          cat(sprintf("  Fleet %d: F=%.4f, Target=%.2f, Predicted=%.2f (%.1f%% error)\n",
                      f, ft[f], ct[f], pct[f], 100*relative_error[f]))
        }
      }
      #end debug

      break
    }

  } #close iter loop: End Newton-Raphson iteration loop


  #convergence diagnostics and warnings

  if (!converged) {
    warning(paste("Newton-Raphson did not converge after", maxiterF, "iterations."))

    cat("\n CONVERGENCE FAILED\n")
    for(f in 1:nfleets) {
      if(tac_managed[f]) {
        cat(sprintf("  Fleet %d: Final F=%.4f, Target=%.2f, Predicted=%.2f (%.1f%% error)\n",
                    f, ft[f], ct[f], pct[f], 100*abs(pct[f]-ct[f])/ct[f]))
      }
    }
  }

  # Check for large final errors (even if converged)
  final_error <- abs(pct - ct) / pmax(ct, tiny)
  if (any(final_error[!is.na(ct)] > 0.1)) {
    warning(paste("Large differences between predicted and target catches.",
                  "Final errors (%):", paste(round(final_error * 100, 1), collapse = ", ")))
  }


  #add convergence information as attributes
  attr(ft, "converged") <- converged
  attr(ft, "iterations") <- iteration_count
  attr(ft, "final_error") <- final_error[!is.na(ct)]
  attr(ft, "predicted_catch") <- pct[!is.na(ct)]
  attr(ft, "target_catch") <- ct[!is.na(ct)]
  attr(ft, "initial_guess") <- ft_initial[!is.na(ct)]

  if (!is_multifleet && nfleets == 1) {
    # For single fleet, return the scalar F but preserve attributes
    F_value <- ft[1]
    attr(F_value, "converged") <- converged
    attr(F_value, "iterations") <- iteration_count
    attr(F_value, "final_error") <- final_error[!is.na(ct)]
    attr(F_value, "predicted_catch") <- pct[!is.na(ct)]
    attr(F_value, "target_catch") <- ct[!is.na(ct)]
    attr(F_value, "initial_guess") <- ft_initial[1]
    return(F_value)
  } else {
    return(ft)
  }

}
