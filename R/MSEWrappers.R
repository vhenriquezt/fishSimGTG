# all multifleet modifications try to preserve the  existing functionality

#---------------------------------------
#Evaluate MSE
#---------------------------------------

#Roxygen header
#'Population dynamics wrapper called by runProjection
#'
#'Contains population dynamics equations. Should not be run directly, instead called by runProjection
#'
#' @param inputObject  A list of objects passed from runProjection
#' @export

evalMSE<-function(inputObject){

  #evalMSE: maintain the same parameters in the multifleet version

  #------------------
  #Unpack dataObject
  #------------------
  #new addition: adding "MultifleetObj" to unpacking so MSEeval can access to the object
  TimeAreaObj <- StrategyObj <- LifeHistoryObj <- HistFisheryObj <- ProFisheryObj_list <- iterations <- iter <- Ddev <- Edev <- LHdev <- Sdev <- Cdev <- Edev <- histEffortDev <- RdevMatrix <- doDiagnostic <- MultifleetObj <- NULL
  for(r in 1:NROW(inputObject)) assign(names(inputObject)[r], inputObject[[r]])

  #new addition: detect multifleet mode (maintain backward compatibility)
  #the code has now two paths: single and multifleet
  is_multifleet <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1

  if(is_multifleet) {
    nfleets <- MultifleetObj@nfleets
    fleet_proportions <- MultifleetObj@fleet_proportions
    cat("Multifleet mode detected with", nfleets, "fleets\n")
  } else {
    nfleets <- 1
    fleet_proportions <- c(1.0)
    cat("Single fleet mode detected\n")
  }

  #same as before
  #controlRuleYear: Defines which years use management strategy
  #years: Total simulation years (equilibrium + historical + projection)
  #areas: Number of spatial areas
  controlRuleYear<-c(FALSE, rep(FALSE,(TimeAreaObj@historicalYears)), rep(TRUE, ifelse(is(StrategyObj, "Strategy")  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)))
  years <- 1 + TimeAreaObj@historicalYears + ifelse(is(StrategyObj, "Strategy")  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)
  areas <- TimeAreaObj@areas

  #--------------
  #Arrays setup (same as before)
  #--------------
  SB<-array(dim=c(years, iterations, areas))
  VB<-array(dim=c(years, iterations, areas))
  RB<-array(dim=c(years, iterations, areas))
  catchN<-array(dim=c(years, iterations, areas))
  catchB<-array(dim=c(years, iterations, areas))
  discN<-array(dim=c(years, iterations, areas))
  discB<-array(dim=c(years, iterations, areas))
  Ftotal<-array(dim=c(years, iterations, areas))
  SPR<-array(dim=c(years, iterations))
  relSB<-array(dim=c(years, iterations))
  recN<-array(dim=c(years, iterations))

  #new addition: for multifleet 4D arrays [years, iterations, areas, fleets]
  if(is_multifleet) {
    Ftotal_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    catchB_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    catchN_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    discB_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    discN_by_fleet <- array(dim=c(years, iterations, areas, nfleets))

  #new addition: initialize with NA
    Ftotal_by_fleet[] <- NA
    catchB_by_fleet[] <- NA
    catchN_by_fleet[] <- NA
    discB_by_fleet[] <- NA
    discN_by_fleet[] <- NA
  }


  #Optional exports for diagnostic mode (remain unchaged)
  Nexport<-NULL
  catchNageExport<-NULL
  Zexport<-NULL

  #-----------------------------------------------
  #Setup recording of management strategy details
  #-----------------------------------------------
  decisionData<-data.frame()
  decisionAnnual<-data.frame()
  decisionLocal<-data.frame()

  #-------------------------------------------
  #Setup capture of benchmarks
  #-------------------------------------------
  ref<-array(dim = c(iterations, 10))

  #-------------------------------------------
  #Deteministic LH and Sel, if present : detect which parameters have stochastic components
  #-------------------------------------------
  LHList<-names(LHdev[!unlist(lapply(LHdev, is.null))])
  selListHist<-names(Sdev$hist[!unlist(lapply(Sdev$hist, is.null))])
  selListPro<-lapply(1:TimeAreaObj@areas, function(x){
    names(Sdev$pro[[x]][!unlist(lapply(Sdev$pro[[x]], is.null))])
  })
  #Is deterministic, so save time by make calculations only once.
  if(NROW(LHList) == 0 & NROW(selListHist) == 0 & NROW(unlist(selListPro)) == 0){
    lh<-LHwrapper(LifeHistoryObj, TimeAreaObj)
    ageClasses <- lh$ageClasses
    if(!is.null(lh) & lh$LifeHistory@Steep < 0.21) lh$LifeHistory@Steep <- 0.21
    if(!is.null(lh) & lh$LifeHistory@Steep > 1) lh$LifeHistory@Steep <- 1

    #new addition: setup selectivity structure for multifleet (deterministic)
    #change from selHist[[area]] to selHist[[area]][[fleet]] structure
    #each fleet has it own sel
    #Hist. sel: uses MultifleetObj@fleet_selectivity_list[[f]] for each fleet
    if(is_multifleet) {
    #   selHist<-lapply(1:TimeAreaObj@areas, function(area){
    #     lapply(1:nfleets, function(f) {
    #       selWrapper(lh, TimeAreaObj, FisheryObj = MultifleetObj@fleet_selectivity_list[[f]], doPlot = FALSE)
    #     })
    #   })

      #changed:
      selHist <- lapply(1:TimeAreaObj@areas, function(area){
        lapply(1:nfleets, function(f) {
          selWrapper(lh, TimeAreaObj,
                     FisheryObj = MultifleetObj@fleet_selectivity_hist_list[[f]],
                     doPlot = FALSE)
        })
      })




      # modif: projection selectivity handling "THIS NEED TO BE IMPROVED"

      # the structure of selPro[[area]][[fleet]]
      # selPro<-lapply(1:TimeAreaObj@areas, function(area){
      #   lapply(1:nfleets, function(f) {
      #     #If ProFisheryObj_list exists, use area-specific projection fishery
      #     if(!is.null(ProFisheryObj_list) && length(ProFisheryObj_list) >= area) {
      #       #each fleet has it own sel per area
      #       selWrapper(lh, TimeAreaObj, FisheryObj = ProFisheryObj_list[[area]], doPlot = FALSE)
      #     } else {
      #       #otherwise fall back to hist fleet sel.
      #       selWrapper(lh, TimeAreaObj, FisheryObj = MultifleetObj@fleet_selectivity_list[[f]], doPlot = FALSE)
      #     }
      #   })
      # })

      #changed- fleet 1-sel1 - area 1/ fleet 1-sel1 - area 2
      #         fleet 2-sel2 - area 1/ fleet 2-sel2 - area 2
      selPro <- lapply(1:TimeAreaObj@areas, function(area){
        lapply(1:nfleets, function(f) {
          # Use fleet-specific projection selectivity
          selWrapper(lh, TimeAreaObj,
                     FisheryObj = MultifleetObj@fleet_selectivity_proj_list[[f]],
                     doPlot = FALSE)
        })
      })

      #Note: I am using area 1 and fleet 1 for this calculation (it could be changed)
      refCalc<-gtgYPRWrapper_Fonly(lh=lh, sel=selHist[[1]][[1]])

    } else {

    #continue with original single fleet strucuture (remain unchanged)
    selHist<-lapply(1:TimeAreaObj@areas, function(x){
      selWrapper(lh, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = FALSE)
    })
    selPro<-lapply(1:TimeAreaObj@areas, function(x){
      selWrapper(lh, TimeAreaObj, FisheryObj = ProFisheryObj_list[[x]], doPlot = FALSE)
    })
    refCalc<-gtgYPRWrapper_Fonly(lh=lh, sel=selHist[[1]]) #unchanged (sel area 1)
    }
    #both paths will create refCalc, so we can use it once here:
    for(k in iter[1]:iter[2]) ref[k, ]<-as.matrix(refCalc$sim)[1,]
    colnames(ref)<-names(refCalc$sim)
  }

  #------------------------------
  #Run simulator of k iterations
  #------------------------------
  if(!is.null(hostName) & !is.null(waitName)){
    waitName$show()
  }
  #step through iterations k
  for(k in iter[1]:iter[2]){
    #print(k)

    #-----------------------------------------------------------------
    #Setup iteration-specific life history & selectivity (if present)
    #-----------------------------------------------------------------
    #check if stochastic parameters exist
    if(NROW(LHList) > 0 | NROW(selListHist) > 0 | NROW(unlist(selListPro)) > 0){
      #Stochastic LH and Fishery objects
      #applies iteration-specific stochastic values to life history and fishery objects
      #creates _TMP objects with stochastic parameters for this iteration
      #LH
      LifeHistoryObj_TMP<-LifeHistoryObj
      if(NROW(LHList) > 0){
        for(x in 1:NROW(LHList)) slot(LifeHistoryObj_TMP, LHList[x]) <- LHdev[[LHList[x]]][k]
      }
      #Hist sel
      HistFisheryObj_TMP<-HistFisheryObj
      if(NROW(selListHist) > 0){
        for(x in 1:NROW(selListHist)) slot(HistFisheryObj_TMP, selListHist[x]) <- Sdev$hist[[selListHist[x]]][k,]
      }
      #Pro sel
      ProFisheryObj_TMP<-lapply(1:TimeAreaObj@areas, function(x){
        TMP<-ProFisheryObj_list[[x]]
        if(NROW(selListPro[[x]]) > 0){
          for(y in 1:NROW(selListPro[[x]])) slot(TMP, selListPro[[x]][y]) <- Sdev$pro[[x]][[selListPro[[x]][y]]][k,]
        }
        TMP
      })
      #setup
      lh<-LHwrapper(LifeHistoryObj_TMP, TimeAreaObj)
      ageClasses <- lh$ageClasses
      if(!is.null(lh) & lh$LifeHistory@Steep < 0.21) lh$LifeHistory@Steep <- 0.21
      if(!is.null(lh) & lh$LifeHistory@Steep > 1) lh$LifeHistory@Steep <- 1

      #new addition: handle selectivity for both single and multifleet
      #multifleet stochastic
      if(is_multifleet) {

        #multifleet selectivity with stochastic parameters
        #maintains [[area]][[fleet]] structure even with stochasticity
        #simplified approach for now: all fleets use same stochastic fishery object
        #fleet differences come only from the base selectivity curves (MultifleetObj@fleet_selectivity_list)
        #stochastic variation affects all fleets equally
        selHist<-lapply(1:TimeAreaObj@areas, function(area){
          lapply(1:nfleets, function(f) {
            selWrapper(lh, TimeAreaObj, FisheryObj = HistFisheryObj_TMP, doPlot = FALSE)
          })
        })
        selPro<-lapply(1:TimeAreaObj@areas, function(area){
          lapply(1:nfleets, function(f) {
            if(!is.null(ProFisheryObj_list) && length(ProFisheryObj_list) >= area) {
              selWrapper(lh, TimeAreaObj, FisheryObj = ProFisheryObj_TMP[[area]], doPlot = FALSE)
            } else {
              selWrapper(lh, TimeAreaObj, FisheryObj = HistFisheryObj_TMP, doPlot = FALSE)
            }
          })
        })

      refCalc<-gtgYPRWrapper_Fonly(lh=lh, sel=selHist[[1]][[1]]) # area 1 and fleet 1 for benchmarks (for now)
      } else {
        #same as before (single fleet stochastic path)
        selHist<-lapply(1:TimeAreaObj@areas, function(x){
          selWrapper(lh, TimeAreaObj, FisheryObj = HistFisheryObj_TMP, doPlot = FALSE)
        })
        selPro<-lapply(1:TimeAreaObj@areas, function(x){
          selWrapper(lh, TimeAreaObj, FisheryObj = ProFisheryObj_TMP[[x]], doPlot = FALSE)
        })
        refCalc<-gtgYPRWrapper_Fonly(lh=lh, sel=selHist[[1]])
      }
      ref[k, ]<-as.matrix(refCalc$sim)[1,]
      colnames(ref)<-names(refCalc$sim)
    }

    #-----------------------------------------
    #Initial equilibrium - year 1 (modification for multifleet approach)
    #-----------------------------------------
    #new addition
    if(is_multifleet) {
      # create selectivity list for multifleet equilibrium
      hist_sel_list <- lapply(1:nfleets, function(f) {
        selHist[[1]][[f]]  # Use area 1 selectivity for equilibrium
      })

    #new addition: multifleet
    is <- solveD_multifleet2(lh = lh, sel_list = hist_sel_list,doFit = TRUE,D_type = TimeAreaObj@historicalBioType,
                               D_in = Ddev[k], fleet_proportions = fleet_proportions,
                               allocation_type = MultifleetObj@allocation_type)

    #extract both total F and fleet-specific F

    # adding new:  store the final proportions to report after runProjection()
    final_effort_proportions <- is$final_effort_proportions
    target_catch_proportions <- is$target_catch_proportions
    actual_catch_proportions <- is$actual_catch_proportions
    allocation_type <- is$allocation_type

    #debugging
    cat("=== evalMSE EQUILIBRIUM DEBUG (iteration", k, ") ===\n")
    cat("is$allocation_type:", is$allocation_type, "\n")
    cat("is$final_effort_proportions:", is$final_effort_proportions, "\n")
    cat("is$actual_catch_proportions:", is$actual_catch_proportions, "\n")
    cat("is$target_catch_proportions:", is$target_catch_proportions, "\n")
    cat("stored final_effort_proportions:", final_effort_proportions, "\n")
    cat("stored actual_catch_proportions:", actual_catch_proportions, "\n")
    cat("Are stored values equal?", identical(final_effort_proportions, actual_catch_proportions), "\n")
    cat("===============================================\n")



    total_Feq <- is$Feq                # total F for population
    F_eq_by_fleet <- is$F_by_fleet     # fleet-specific F values

    } else {
      # single fleet equilibrium (original code unchanged)
      is<-solveD(lh, sel = selHist[[1]], doFit = TRUE, D_type = TimeAreaObj@historicalBioType, D_in = Ddev[k])
      total_Feq <- is$Feq
      F_eq_by_fleet <- c(is$Feq)

      # set NULL values for single fleet
      final_effort_proportions <- NULL
      target_catch_proportions <- NULL
      actual_catch_proportions <- NULL
      allocation_type <- NULL
    }

    #Burn-in to calibrate N by area, noting effect of movement (this is for area distribution)
    #remain the same
    Ntmp <- list()
    yrsTmp <- (ageClasses*4)
    for (l in 1:lh$gtg){
      Ntmp[[l]]<-array(dim=c(ageClasses, yrsTmp, areas))
      for (m in 1:areas){
        Ntmp[[l]][,1,m]<-is$N[[l]]*TimeAreaObj@recArea[m]
      }
    }
    for (j in 1: yrsTmp){
      for (l in 1:lh$gtg){
        #Cohort equations + recruitment
        if(j< yrsTmp){
          P<-matrix(nrow=ageClasses*areas, ncol=ageClasses*areas)


          for(m in 1:areas){
            #new addition:ASSUMPTION FLEET 1 is a representative fleet (THIS IS IMPORTANT)
            # the selectivity here is important  because it affects which ages survive during burn in
            # affecting the age structure in each area (ASK Bill , what we should do here)
            if(is_multifleet) {
              # Fleet 1 is the "main" fleet, equivalent to single fleet
              sel_removal <- selHist[[m]][[1]]$removal[[l]]

              #debugging
              # cat("DEBUG: m=", m, ", l=", l, ", length of sel_removal=", length(sel_removal), "\n")
              # if(length(sel_removal) == 0) {
              #   cat("ERROR: sel_removal has length 0!\n")
              #   print(selHist[[m]][[1]])
              # }

            } else {
              # Single fleet
              sel_removal <- selHist[[m]]$removal[[l]]
            }

            #S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=is$Feq, S_in=selHist[[m]]$removal[[l]] )
            #new addition:
            S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=is$Feq, S_in=sel_removal)

            # #debugging:
            # cat("About to call SurvMat: j=", j, ", k=", k, ", m=", m, ", l=", l, "\n")
            # cat("F_in =", Ftotal[j,k,m], "\n")
            # cat("selGroup structure for area", m, ":\n")
            # if(is_multifleet) {
            #   cat("Multifleet mode - checking selGroup[[", m, "]]$removal[[", l, "]]\n")
            #   if(is.null(selGroup[[m]]$removal[[l]])) {
            #     cat("ERROR: selGroup[[", m, "]]$removal[[", l, "]] is NULL!\n")
            #   } else {
            #     cat("Length:", length(selGroup[[m]]$removal[[l]]), "\n")
            #   }
            # } else {
            #   cat("Single fleet mode - checking selGroup[[", m, "]]$removal[[", l, "]]\n")
            #   if(is.null(selGroup[[m]]$removal[[l]])) {
            #     cat("ERROR: selGroup[[", m, "]]$removal[[", l, "]] is NULL!\n")
            #   } else {
            #     cat("Length:", length(selGroup[[m]]$removal[[l]]), "\n")
            #   }
            # }



            #remain unchanged
            rows<-c((m-1)*dim(Ntmp[[l]])[1]+1,m*dim(Ntmp[[l]])[1])
            cols<-c(1,(dim(Ntmp[[l]])[1]*areas))
            P[rows[1]:rows[2],cols[1]:cols[2]]<- MoveMat(Surv_in=S, Move_in=TimeAreaObj@move, area_in=m)
          }
          tmp<-matrix(as.vector(Ntmp[[l]][,j,]), nrow=(dim(Ntmp[[l]])[1]*areas), ncol=1)
          Ntmp[[l]][,(j+1),]<-matrix(P%*%tmp, nrow=dim(Ntmp[[l]])[1], ncol=areas, byrow=FALSE)
          Ntmp[[l]][1,(j+1),]<- is$Req*lh$recProb[l]*TimeAreaObj@recArea
        }
      }
    }

    #Specify iniitial conditions to start sims in iteration k
    N<-list()
    for(l in 1:lh$gtg){
      N[[l]]<-array(dim=c(ageClasses, years, areas))
    }
    for (l in 1:lh$gtg){
      for (m in 1:areas){
        N[[l]][,1,m]<-Ntmp[[l]][,yrsTmp,m]
        N[[l]][,2,m]<-Ntmp[[l]][,yrsTmp,m]
      }
    }

    #Initialize catch-at-age arrays and Z array
    catchNage<-list()
    Z<-list()
    for(l in 1:lh$gtg){
      catchNage[[l]]<-array(dim=c(ageClasses, years, areas))
      Z[[l]]<-array(dim=c(ageClasses, years, areas))
    }

    #Arrays (remain unchanged)
    for(m in 1:areas) SB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]][,1,m]*lh$mat[[x]]*lh$W[[x]])[2:ageClasses])))
    SPR[1,k]<-(sum(SB[1,k,])/is$Req)/(is$B0/lh$LifeHistory@R0)
    relSB[1,k]<-sum(SB[1,k,])/is$B0
    recN[1,k]<-is$Req

    #new addition:
    if(is_multifleet) {
      #initialize fleet-specific catch arrays for year 1
      #structure: catchNage_by_fleet[[fleet]][[gtg]][age, year, area]
      catchNage_by_fleet <- list()
      for(f in 1:nfleets) {
        catchNage_by_fleet[[f]] <- list()
        for(l in 1:lh$gtg) {
          catchNage_by_fleet[[f]][[l]] <- array(dim=c(ageClasses, years, areas))
        }
      }
    }
    #main area loop: calculate Year 1 catches and biomass by area
    for(m in 1:areas){
      if(is_multifleet) {

        #new addition: MULTIFLEET YEAR 1 CALCULATIONS
        #fleet-specific calculations for year 1
        for(f in 1:nfleets) {
          #store fleet-specific F in 4D array [year, iteration, area, fleet]
          Ftotal_by_fleet[1,k,m,f] <- F_eq_by_fleet[f]

          # fleet-specific catches for year 1
          for(l in 1:lh$gtg){
            #debugging:
            # cat("=== Year 1 Multifleet Debug ===\n")
            # cat("m=", m, ", f=", f, ", l=", l, "\n")
            # cat("F_eq_by_fleet:", F_eq_by_fleet, "\n")
            # cat("nfleets:", nfleets, "\n")
            #
            # # Checking if selHist structure is OK
            # if(is.null(selHist[[m]])) {
            #   cat("ERROR: selHist[[", m, "]] is NULL\n")
            # } else if(is.null(selHist[[m]][[1]])) {
            #   cat("ERROR: selHist[[", m, "]][[1]] is NULL\n")
            # } else {
            #   cat("selHist[[", m, "]][[1]]$removal[[", l, "]] length:", length(selHist[[m]][[1]]$removal[[l]]), "\n")
            # }






            #calculate total Z from all fleet contributions
            #Z = M + sum_across_fleets(F_fleet * selectivity_fleet)
            # NEVER combine selectivities - each fleet contributes independently
            total_fishing_mortality <- sapply(1:ageClasses, function(age) {
              sum(sapply(1:nfleets, function(ff) {
                F_eq_by_fleet[ff] * selHist[[m]][[ff]]$removal[[l]][age]
              }))
            })

            #debugging:
            # cat("total_fishing_mortality range:", range(total_fishing_mortality, na.rm = TRUE), "\n")
            # if(any(is.na(total_fishing_mortality))) {
            #   cat("WARNING: total_fishing_mortality contains NA values\n")
            # }

            #total mortality shared by all fleets: Z = M + total_fishing_mortality
            Z[[l]][,1,m] <- total_fishing_mortality + lh$LifeHistory@M

            #debugging:
            # cat("Z range:", range(Z[[l]][,1,m], na.rm = TRUE), "\n")
            # cat("===============================\n")

            #fleet-specific catch (using using Baranov equation with shared Z)
            catchNage_by_fleet[[f]][[l]][,1,m] <- F_eq_by_fleet[f] * selHist[[m]][[f]]$keep[[l]] /
              Z[[l]][,1,m] * (1-exp(-Z[[l]][,1,m])) * N[[l]][,1,m]
          }

          # fleet totals: sum across GTGs and ages for each fleet
          catchN_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(catchNage_by_fleet[[f]][[x]][,1,m])))
          catchB_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*catchNage_by_fleet[[f]][[x]][,1,m])))
          discN_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(F_eq_by_fleet[f]*selHist[[m]][[f]]$discard[[x]]/(Z[[l]][,1,m])*(1-exp(-Z[[l]][,1,m]))*N[[x]][,1,m])))
          discB_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*F_eq_by_fleet[f]*selHist[[m]][[f]]$discard[[x]]/(Z[[l]][,1,m])*(1-exp(-Z[[l]][,1,m]))*N[[x]][,1,m])))
        }

        #area totals (sum across fleets)
        #VB now sums across all fleets (each fleet contributes to vulnerable biomass)
        # VB[1,k,m] <- sum(sapply(1:nfleets, function(f) {
        #   sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist[[m]][[f]]$vul[[x]]*lh$W[[x]])))
        # }))

        #changed:
        #VB uses maximum vulnerability across fleets (otherwise I would double-counting VB)
        VB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) {
          max_vuln_by_age <- sapply(1:lh$ageClasses, function(age) {
            max(sapply(1:nfleets, function(f) selHist[[m]][[f]]$vul[[x]][age]))
          })
          sum(N[[x]][,1,m] * max_vuln_by_age * lh$W[[x]])
        }))



        #RB now sums fleet-specific catches
        RB[1,k,m] <- sum(catchB_by_fleet[1,k,m,1:nfleets], na.rm = TRUE)
        #Ftotal now sums all fleet F values
        Ftotal[1,k,m] <- sum(F_eq_by_fleet)  # Total F

        #Total catchNage across fleets: Sum fleet-specific catches to create total catch-at-age
        #create total catchNage across fleets for backward compatibility
        for(l in 1:lh$gtg){
          catchNage[[l]][,1,m] <- rowSums(sapply(1:nfleets, function(f) catchNage_by_fleet[[f]][[l]][,1,m]), na.rm = TRUE)
        }

        # single fleet (remain unchanged for single fleet)
      } else {
      VB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist[[m]]$vul[[x]]*lh$W[[x]])))
      RB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist[[m]]$keep[[x]]*lh$W[[x]])))
      Ftotal[1,k,m] <- is$Feq

      #Calculate catch matrices by gtg and age
      for(l in 1:lh$gtg){
        Z[[l]][,1,m] <- Ftotal[1,k,m]*selHist[[m]]$removal[[l]] + lh$LifeHistory@M
        catchNage[[l]][,1,m] <- Ftotal[1,k,m]*selHist[[m]]$keep[[l]]/(Z[[l]][,1,m])*(1-exp(-Z[[l]][,1,m]))*N[[l]][,1,m]
      }
      }

      #total catch and discards (multiffleet)
      if(is_multifleet) {
        #sum across multifleets
        catchN[1,k,m] <- sum(catchN_by_fleet[1,k,m,1:nfleets], na.rm = TRUE)
        catchB[1,k,m] <- sum(catchB_by_fleet[1,k,m,1:nfleets], na.rm = TRUE)
        discN[1,k,m] <- sum(discN_by_fleet[1,k,m,1:nfleets], na.rm = TRUE)
        discB[1,k,m] <- sum(discB_by_fleet[1,k,m,1:nfleets], na.rm = TRUE)


        } else {
          #remain unchanged for single fleet

      catchN[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(catchNage[[x]][,1,m])))
      catchB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*catchNage[[x]][,1,m])))
      discN[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal[1,k,m]*selHist[[m]]$discard[[x]]/(Ftotal[1,k,m]*selHist[[m]]$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[1,k,m]*selHist[[m]]$removal[[x]]-lh$LifeHistory@M))*N[[x]][,1,m])))
      discB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal[1,k,m]*selHist[[m]]$discard[[x]]/(Ftotal[1,k,m]*selHist[[m]]$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[1,k,m]*selHist[[m]]$removal[[x]]-lh$LifeHistory@M))*N[[x]][,1,m])))
        }
    }

    #--------------------
    #Time dynamics
    #--------------------
    # No modifications so far
    for (j in 2:years){

      #Selgroup
      if(controlRuleYear[j]) selGroup <- selPro
      if(!controlRuleYear[j]) selGroup <- selHist

      #Annual regulation decisions - phase 2
      dataObject<-c(list(j=j,
                         k=k,
                         is=is,
                         lh = lh,
                         areas = areas,
                         ageClasses = ageClasses,
                         N=N,
                         Z=Z,
                         catchNage=catchNage,
                         selGroup = selGroup,
                         selHist = selHist,
                         selPro = selPro,
                         SB=SB,
                         VB=VB,
                         RB=RB,
                         catchN=catchN,
                         catchB=catchB,
                         discN=discN,
                         discB=discB,
                         Ftotal=Ftotal,
                         SPR=SPR,
                         relSB=relSB,
                         recN=recN,
                         decisionData=decisionData,
                         decisionAnnual=decisionAnnual,
                         decisionLocal=decisionLocal
      ),
      inputObject
      )
      if(controlRuleYear[j]) decisionAnnual<-rbind(decisionAnnual, do.call(get(StrategyObj@projectionName), list(phase=2, dataObject)))

      #Localized F at each location - phase 3
      dataObject<-c(list(j=j,
                         k=k,
                         is=is,
                         lh = lh,
                         areas = areas,
                         ageClasses = ageClasses,
                         N=N,
                         Z=Z,
                         catchNage=catchNage,
                         selGroup = selGroup,
                         selHist = selHist,
                         selPro = selPro,
                         SB=SB,
                         VB=VB,
                         RB=RB,
                         catchN=catchN,
                         catchB=catchB,
                         discN=discN,
                         discB=discB,
                         Ftotal=Ftotal,
                         SPR=SPR,
                         relSB=relSB,
                         recN=recN,
                         decisionData=decisionData,
                         decisionAnnual=decisionAnnual,
                         decisionLocal=decisionLocal
      ),
      inputObject
      )
      if(controlRuleYear[j]) { decisionLocal<-rbind(decisionLocal, do.call(get(StrategyObj@projectionName), list(phase=3, dataObject)))
      } else { decisionLocal<-rbind(decisionLocal, do.call(fixedStrategy, list(phase=3, dataObject)))}  #fixedStrategy() need modification for multifleet

      #SB and recruits
      for(m in 1:areas) SB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]][,j,m]*lh$mat[[x]]*lh$W[[x]])[2:ageClasses])))
      Rtmp<-recruit(LifeHistoryObj = lh$LifeHistory, B0=is$B0, stock=sum(SB[j,k,]))
      SPR[j,k]<-(sum(SB[j,k,])/Rtmp)/(is$B0/lh$LifeHistory@R0)
      relSB[j,k]<-sum(SB[j,k,])/is$B0
      recN[j,k]<-Rtmp*RdevMatrix[j,k]
      for (l in 1:lh$gtg) N[[l]][1,j,]<- Rtmp*lh$recProb[l]*TimeAreaObj@recArea*RdevMatrix[j,k]

      #Arrays
      # Modifications need to be added here , this section is moved to line 626-643

      # for(m in 1:areas){
      #   VB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup[[m]]$vul[[x]]*lh$W[[x]])))
      #   RB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup[[m]]$keep[[x]]*lh$W[[x]])))
      #   xRow<-which(decisionLocal$year==j & decisionLocal$iteration==k & decisionLocal$area==m)
      #   Ftotal[j,k,m] <- decisionLocal$Flocal[xRow]
      #
      #   for(l in 1:lh$gtg){
      #     Z[[l]][,j,m] <- Ftotal[j,k,m]*selGroup[[m]]$removal[[l]] + lh$LifeHistory@M
      #     catchNage[[l]][,j,m] <- Ftotal[j,k,m]*selGroup[[m]]$keep[[l]]/(Z[[l]][,j,m])*(1-exp(-Z[[l]][,j,m]))*N[[l]][,j,m]
      #   }
      #
      #   catchN[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(catchNage[[x]][,j,m])))
      #   catchB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*catchNage[[x]][,j,m])))
      #   discN[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal[j,k,m]*selGroup[[m]]$discard[[x]]/(Ftotal[j,k,m]*selGroup[[m]]$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[j,k,m]*selGroup[[m]]$removal[[x]]-lh$LifeHistory@M))*N[[x]][,j,m])))
      #   discB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal[j,k,m]*selGroup[[m]]$discard[[x]]/(Ftotal[j,k,m]*selGroup[[m]]$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[j,k,m]*selGroup[[m]]$removal[[x]]-lh$LifeHistory@M))*N[[x]][,j,m])))
      # }

      #new addition:check if we have the multifleet mode activated
      if(is_multifleet) {
        for(m in 1:areas) { # main area loop for multfleet
      #initialize fleet-specific F values for this area
          F_by_fleet_current <- numeric(nfleets)

      # we are in the historical period (before management)
          if(!controlRuleYear[j]) {
      #historical period - fleet specific F scaling (multiplier)
      #each fleet scale independently from equilibrium F
            yr <- j - 1  # year index for historical effort (the time dyn loop starts at 2) this is to index correctly the TimeAreaObj@historicalEffort

              for(f in 1:nfleets) {
                #hist effort for each fleet * multiplier eff (scaling)
                #fleet f's equilibrium F from Year 1 times area specific effor multiplier* deviation
              F_by_fleet_current[f] <- F_eq_by_fleet[f] *
                TimeAreaObj@historicalEffort[yr,m] *
                histEffortDev[j,k,m,f]  # this is a 4D array now
      #store in 4D fleet array
              Ftotal_by_fleet[j,k,m,f] <- F_by_fleet_current[f]
            }
          } else {

          # New projection period - Use management strategy results

            # decisionLocal: df stores the F decisions (Flocal) made by management startegies for each year, iteration, and area
            xRow <- which(decisionLocal$year==j & decisionLocal$iteration==k &
                            decisionLocal$area==m & decisionLocal$fleet==0) # search in the df to find the row that match the year, iter, and m, xRow return the row that provide the combined value

            if(length(xRow) == 0) {
              # fallback: sum fleet-specific F values if total not available
              fleet_rows <- which(decisionLocal$year==j & decisionLocal$iteration==k &
                                    decisionLocal$area==m & decisionLocal$fleet > 0)
              total_F_from_strategy <- sum(decisionLocal$Flocal[fleet_rows])
            } else {
              total_F_from_strategy <- decisionLocal$Flocal[xRow]
            }

            # scale each fleet proportionally
            # for now, assuming proportional scaling (will improve later)
            for(f in 1:nfleets) {
              F_by_fleet_current[f] <- total_F_from_strategy * fleet_proportions[f]
              Ftotal_by_fleet[j,k,m,f] <- F_by_fleet_current[f]
            }
            }

              # # SINGLE FLEET - Original logic
              # xRow <- which(decisionLocal$year==j & decisionLocal$iteration==k & decisionLocal$area==m)
              # total_F_from_strategy <- decisionLocal$Flocal[xRow]

            #   # For single fleet, just set the total F
            #   Ftotal[j,k,m] <- F_by_fleet_current
            # #}



          # calculate total F for this area (sum across fleets)
          # this total F is used for population dynamics (survival, movement)
          # also maintain backward compatibi;ity
          Ftotal[j,k,m] <- sum(F_by_fleet_current) #sum all fleet-specific F values to get the total fishing pressure

          # calculate Z and catches for each GTG
          for(l in 1:lh$gtg) {
            # CRITICAL - Calculate total fishing mortality across all fleets
            # Z = M + sum_fleets(F_fleet * selectivity_fleet)
            # DO NOT COMBINE SELECTIVITIES - each fleet contributes separately
            # regardless of historical vs projection: the total F is calculated in this way
            total_fishing_mortality <- sapply(1:ageClasses, function(age) {
              fleet_mortality_sum <- 0   # sum contribution from all fleets
              for(f in 1:nfleets) {
                # Each fleet contributes: F_fleet * selectivity_fleet
                # for example:
                # F_by_fleet_current = [0.24, 0.16]* selGroup[[area]][[1]]$removal[[1]][5] = 0.8 # Age 5, GTG 1
                fleet_mortality_sum <- fleet_mortality_sum +
                  F_by_fleet_current[f] * selGroup[[m]][[f]]$removal[[l]][age]
              }
              return(fleet_mortality_sum)
            })

            # Z= natural + all fleet fishing mortalities
            # Z[gtg][age,year,area]
            Z[[l]][,j,m] <- total_fishing_mortality + lh$LifeHistory@M

            # fleet-specific catches using shared Z
            for(f in 1:nfleets) {
              # fleet-specific catch using Baranov equation with shared Z
              catchNage_by_fleet[[f]][[l]][,j,m] <-
                F_by_fleet_current[f] * selGroup[[m]][[f]]$keep[[l]] /
                Z[[l]][,j,m] * (1-exp(-Z[[l]][,j,m])) * N[[l]][,j,m]
            }

            # catchNage is sum across all fleets (for backward compativbility)
            # maintains compatibility with existing single-fleet code
            catchNage[[l]][,j,m] <- rowSums(sapply(1:nfleets, function(f)
              catchNage_by_fleet[[f]][[l]][,j,m]), na.rm = TRUE)
          }

          # fleet-specific totals
          for(f in 1:nfleets) {
            catchN_by_fleet[j,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x)
              sum(catchNage_by_fleet[[f]][[x]][,j,m])))

            catchB_by_fleet[j,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x)
              sum(lh$W[[x]] * catchNage_by_fleet[[f]][[x]][,j,m])))

            discN_by_fleet[j,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x)
              sum(F_by_fleet_current[f] * selGroup[[m]][[f]]$discard[[x]] /
                    Z[[l]][,j,m] * (1-exp(-Z[[l]][,j,m])) * N[[x]][,j,m])))

            discB_by_fleet[j,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x)
              sum(lh$W[[x]] * F_by_fleet_current[f] * selGroup[[m]][[f]]$discard[[x]] /
                    Z[[l]][,j,m] * (1-exp(-Z[[l]][,j,m])) * N[[x]][,j,m])))
          }

          #calculate area totals for existing arrays (backward compatibility)
          # each fleet sees different vulnerable biomass
          # VB[j,k,m] <- sum(sapply(1:nfleets, function(f) {
          #   sum(sapply(1:lh$gtg, FUN=function(x)
          #     sum(N[[x]][,j,m] * selGroup[[m]][[f]]$vul[[x]] * lh$W[[x]])))
          # }))

          #changed:
          # VB uses maximum vulnerability across fleets (no double-counting)
          VB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) {
            max_vuln_by_age <- sapply(1:lh$ageClasses, function(age) {
              max(sapply(1:nfleets, function(f) selGroup[[m]][[f]]$vul[[x]][age]))
            })
            sum(N[[x]][,j,m] * max_vuln_by_age * lh$W[[x]])
          }))



          RB[j,k,m] <- sum(catchB_by_fleet[j,k,m,1:nfleets], na.rm = TRUE)

          catchN[j,k,m] <- sum(catchN_by_fleet[j,k,m,1:nfleets], na.rm = TRUE)
          catchB[j,k,m] <- sum(catchB_by_fleet[j,k,m,1:nfleets], na.rm = TRUE)
          discN[j,k,m] <- sum(discN_by_fleet[j,k,m,1:nfleets], na.rm = TRUE)
          discB[j,k,m] <- sum(discB_by_fleet[j,k,m,1:nfleets], na.rm = TRUE)
        }

      } else {

        #Arrays
        for(m in 1:areas){
          VB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup[[m]]$vul[[x]]*lh$W[[x]])))
          RB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup[[m]]$keep[[x]]*lh$W[[x]])))
          xRow<-which(decisionLocal$year==j & decisionLocal$iteration==k & decisionLocal$area==m)
          Ftotal[j,k,m] <- decisionLocal$Flocal[xRow]

          for(l in 1:lh$gtg){
            Z[[l]][,j,m] <- Ftotal[j,k,m]*selGroup[[m]]$removal[[l]] + lh$LifeHistory@M
            catchNage[[l]][,j,m] <- Ftotal[j,k,m]*selGroup[[m]]$keep[[l]]/(Z[[l]][,j,m])*(1-exp(-Z[[l]][,j,m]))*N[[l]][,j,m]
          }

          catchN[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(catchNage[[x]][,j,m])))
          catchB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*catchNage[[x]][,j,m])))
          discN[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal[j,k,m]*selGroup[[m]]$discard[[x]]/(Ftotal[j,k,m]*selGroup[[m]]$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[j,k,m]*selGroup[[m]]$removal[[x]]-lh$LifeHistory@M))*N[[x]][,j,m])))
          discB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal[j,k,m]*selGroup[[m]]$discard[[x]]/(Ftotal[j,k,m]*selGroup[[m]]$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[j,k,m]*selGroup[[m]]$removal[[x]]-lh$LifeHistory@M))*N[[x]][,j,m])))
        }
      }



      #Next year abundance, move through each gtg (no changes needed)
      for (l in 1:lh$gtg){
        if(j<years){
          P<-matrix(nrow=ageClasses*areas, ncol=ageClasses*areas)
          for(m in 1:areas){
            # no changes needed: it uses total Ftotal[j,k,m] which is correct for both single and multifleet
            #S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=Ftotal[j,k,m], S_in=selGroup[[m]]$removal[[l]])

            # changes added - selectivity
            if(is_multifleet) {
              # use Fleet 1 (assuming this is the most representative fleet) selectivity for all areas in movement calculations
              # I used this approach before too
              S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=Ftotal[j,k,m], S_in=selGroup[[m]][[1]]$removal[[l]])
            } else {
              # single fleet as original code
              S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=Ftotal[j,k,m], S_in=selGroup[[m]]$removal[[l]])
            }


            rows<-c((m-1)*dim(N[[l]])[1]+1,m*dim(N[[l]])[1])
            cols<-c(1,(dim(N[[l]])[1]*areas))
            P[rows[1]:rows[2],cols[1]:cols[2]]<- MoveMat(Surv_in=S, Move_in=TimeAreaObj@move, area_in=m)
          }
          Ntmp<-matrix(as.vector(N[[l]][,j,]), nrow=(dim(N[[l]])[1]*areas), ncol=1)
          N[[l]][,(j+1),]<-matrix(P%*%Ntmp, nrow=dim(N[[l]])[1], ncol=areas, byrow=FALSE)
        }
      }

      #Sampling - phase 1
      if(!is.null(StrategyObj)){
        dataObject<-c(list(j=j,
                           k=k,
                           is=is,
                           lh = lh,
                           areas = areas,
                           ageClasses = ageClasses,
                           N=N,
                           Z=Z,
                           catchNage=catchNage,
                           selGroup = selGroup,
                           selHist = selHist,
                           selPro = selPro,
                           SB=SB,
                           VB=VB,
                           RB=RB,
                           catchN=catchN,
                           catchB=catchB,
                           discN=discN,
                           discB=discB,
                           Ftotal=Ftotal,
                           SPR=SPR,
                           relSB=relSB,
                           recN=recN,
                           decisionData=decisionData,
                           decisionAnnual=decisionAnnual,
                           decisionLocal=decisionLocal
        ),
        inputObject
        )
        decisionData<-rbind(decisionData, do.call(get(StrategyObj@projectionName), list(phase=1, dataObject)))
      }
      }

      #stop ()


    #Optional exports for diagnostic mode.
    if(doDiagnostic & k==1) {
      Nexport = N
      catchNageExport = catchNage
      Zexport = Z
    }


    if(!is.null(hostName) & !is.null(waitName)){
      hostName$set(k/floor(TimeAreaObj@iterations)*100)
    }
  }
  if(!is.null(hostName) & !is.null(waitName)){
    waitName$hide()
  }

  # new addition: add multifleet data to the final output
  if(is_multifleet) {
    dynamics_multifleet <- list(
      Ftotal_by_fleet = Ftotal_by_fleet,
      catchB_by_fleet = catchB_by_fleet,
      catchN_by_fleet = catchN_by_fleet,
      discB_by_fleet = discB_by_fleet,
      discN_by_fleet = discN_by_fleet,
      fleet_proportions = fleet_proportions,
      nfleets = nfleets,

      # Adding new
      final_effort_proportions = final_effort_proportions,
      target_catch_proportions = target_catch_proportions,
      actual_catch_proportions = actual_catch_proportions,
      allocation_type = allocation_type
    )


  #save
  dynamics<-list(SB=SB, VB=VB, RB=RB, catchB=catchB, catchN=catchN, Ftotal=Ftotal,
                 discB=discB, discN=discN, SPR=SPR, relSB=relSB, recN=recN, ref = ref,
                 multifleet = dynamics_multifleet)
  } else {
    dynamics<-list(SB=SB, VB=VB, RB=RB, catchB=catchB, catchN=catchN, Ftotal=Ftotal,
                   discB=discB, discN=discN, SPR=SPR, relSB=relSB, recN=recN, ref = ref)
  }


  HCR<-list(decisionLocal=decisionLocal, decisionAnnual=decisionAnnual, decisionData=decisionData)
  return(list(dynamics=dynamics, HCR=HCR, iter=iter, N=Nexport, Z=Zexport, catchNage=catchNageExport))
}



#---------------------------------------
#Run the projection or MSE model
#---------------------------------------

#Roxygen header
#'Run the projection or MSE model
#'
#'Function for running projections or MSE
#'
#' @param LifeHistoryObj  A LifeHistory object. Required
#' @param TimeAreaObj A TimeArea object. Required
#' @param HistFisheryObj A Fishery object that characterizes the historical dynamics. Required as it is used in initial equilibrium and historical time dynamics (if applicable)
#' @param ProFisheryObj_list A Fishery object used in forward projection. Optional, only used when StrategyObj is supplied
#' @param StrategyObj A Strategy object. Optional
#' @param StochasticObj A Stochastic object. Optional
#' @param wd A working directly to save output. Required
#' @param fileName A file name for output. Required
#' @param seed A value used in base::set.seed function for producing consistent set of stochastic elements. Optional
#' @param doPlot Logical whether to produce diagnostic plots upon completing simulations. Default is FALSE (no plots)
#' @param customToCluster A character vector containing name or names of custom management strategies to export to the cluster (otherwise parallel processing will fail).
#' @param titleStrategy A title for management strategy being evaluated.
#' @param waitName When used within a shiny app, this function can update a host from the waiter package. See example.
#' @param hostName When used within a shiny app, this function can update a host from the waiter package. See example.
#' @importFrom grDevices dev.off png rainbow
#' @importFrom graphics mtext points
#' @importFrom snowfall sfInit sfLibrary sfLapply sfRemoveAll sfStop sfExport
#' @importFrom parallel detectCores
#' @importFrom methods is
#' @importFrom shinyWidgets updateProgressBar
#' @importFrom here here
#' @export


runProjection<-function(LifeHistoryObj, TimeAreaObj, HistFisheryObj, ProFisheryObj_list = NULL, StrategyObj = NULL, StochasticObj = NULL,MultifleetObj = NULL, IndexObj=NULL, CatchObsObj=NULL, LengthCompObj=NULL,
                        wd, fileName, seed = 1, doPlot = FALSE, doDiagnostic=F, customToCluster = NULL, titleStrategy = "No name", waitName=NULL, hostName=NULL){

  #-----------------------
  #Build inputObject
  #-----------------------
  TimeAreaObj@recArea <- TimeAreaObj@recArea / sum(TimeAreaObj@recArea) #Make sure this sums to 1

  #new addition: adding basic multifleet detection and validation
  is_multifleet <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1

  if(is_multifleet) {
    nfleets <- MultifleetObj@nfleets

    # basic validation - just check fleet proportions for now
    if(abs(sum(MultifleetObj@fleet_proportions) - 1.0) > 1e-6) {
      stop(paste("Fleet proportions must sum to 1.0. current sum:", sum(MultifleetObj@fleet_proportions)))
    }

    # if(length(MultifleetObj@fleet_selectivity_list) != nfleets) {
    #   stop(paste("fleet_selectivity_list must contain", nfleets, "Fishery objects"))
    # }

    #changed:
    if(length(MultifleetObj@fleet_selectivity_hist_list) != nfleets ||
       length(MultifleetObj@fleet_selectivity_proj_list) != nfleets) {
      stop(paste("fleet_selectivity_hist_list and fleet_selectivity_proj_list must each contain", nfleets, "Fishery objects"))
    }


    cat("Multifleet mode enabled with", nfleets, "fleets\n")
    cat("Fleet proportions:", paste(round(MultifleetObj@fleet_proportions, 3), collapse = ", "), "\n")

  } else {
    nfleets <- 1
    cat("Single fleet mode\n")
  }


  #------------------------------------------------
  #Build stochastic & uncertainty range parameters
  #------------------------------------------------
  set.seed(seed = seed)

  #Rec devs
  RdevMatrix<-recDev(LifeHistoryObj, TimeAreaObj, StochasticObj, StrategyObj)$Rmult

  #Initial depletion (SSB rel)
  Ddev<-bioDev(TimeAreaObj, StochasticObj)$Ddev

  #Historical cpue used only in projectionStrategy
  Cdev<-NULL
  if(is(StrategyObj, "Strategy") &&

     StrategyObj@projectionName == "projectionStrategy"
  ) Cdev<-cpueDev(TimeAreaObj, StrategyObj)$Cdev

  #Effort implementation error used only in projectionStrategy
  Edev<-NULL
  if(is(StrategyObj, "Strategy") &&
     StrategyObj@projectionName == "projectionStrategy"
  ) Edev<-effortImpErrorDev(TimeAreaObj, StrategyObj)$Edev

  #Life history parmeters
  LHdev<-lifehistoryDev(TimeAreaObj, StochasticObj)

  #Selectivity parameters
  Sdev<-selDev(TimeAreaObj, HistFisheryObj, ProFisheryObj_list, StochasticObj)

  #new addition: determine nfleets and call histEffortDev
  effective_nfleets <- if(is_multifleet) nfleets else 1
  cat("Creating histEffortDev with", effective_nfleets, "fleets\n")

  #Historical effort devs (adding multifleet)
  histEffortDev_result<-histEffortDev(TimeAreaObj, StochasticObj, effective_nfleets)

  # single fleet expects: histEffortDev[year, iteration, area] (3D)
  # multifleet expects:   histEffortDev[year, iteration, area, fleet] (4D)

  # Handle backward compatibility

  if(is_multifleet) {
    histEffortDev <- histEffortDev_result$Emult  # Use 4D array


    expected_dims <- c(1 + TimeAreaObj@historicalYears,
                       as.integer(floor(TimeAreaObj@iterations)),
                       TimeAreaObj@areas,
                       nfleets)

    cat("multifleet mode: using 4D histEffortDev array\n")
    cat("  expected dims:", paste(expected_dims, collapse = " x "), "\n")
    cat("  actual dims:  ", paste(dim(histEffortDev), collapse = " x "), "\n")


    # validation
    if(!all(dim(histEffortDev) == expected_dims)) {
      stop("multifleet histEffortDev dimension mismatch.\n",
           "  expected: ", paste(expected_dims, collapse = " x "), "\n",
           "  got:      ", paste(dim(histEffortDev), collapse = " x "))
    }


  } else {
    # For single fleet, check if 3D is available
    if(!is.null(histEffortDev_result$Emult_3D)) {
      histEffortDev <- histEffortDev_result$Emult_3D  # Use 3D array for backward compatibility
      cat("single fleet mode: using 3D histEffortDev array: dim =", paste(dim(histEffortDev), collapse = " x "), "\n")
      } else {
      histEffortDev <- histEffortDev_result$Emult[,,,1]  # Extract first fleet from 4D
      cat("single fleet mode: extracted 3D from 4D histEffortDev array: dim =", paste(dim(histEffortDev), collapse = " x "), "\n")
      }


  expected_dims <- c(1 + TimeAreaObj@historicalYears,
                     as.integer(floor(TimeAreaObj@iterations)),
                     TimeAreaObj@areas)

  cat("  expected dims:", paste(expected_dims, collapse = " x "), "\n")
  cat("  actual dims:  ", paste(dim(histEffortDev), collapse = " x "), "\n")


  # validate dimensions
  if(!all(dim(histEffortDev) == expected_dims)) {
    stop("single fleet histEffortDev dimension mismatch.\n",
         "  expected: ", paste(expected_dims, collapse = " x "), "\n",
         "  got:      ", paste(dim(histEffortDev), collapse = " x "))
  }
}

cat("histEffortDev validation passed!\n")






  #---------------------------------------
  #Initial checks that do not stop program
  #---------------------------------------

  print("
  #---------------
  #Initial checks
  #---------------
  ")

  #Check to see if uncertain initial bio created
  if(!is.null(StochasticObj)){
    if(length(StochasticObj@historicalBio) > 1) {
      print(paste("Uncertainty in initial biomass:", StochasticObj@historicalBio[1], "to", StochasticObj@historicalBio[2], TimeAreaObj@historicalBioType, "created."))
    } else {
      print(paste("Uncertainty in initial biomass:", "none"))
    }
  }

  #Check to see if uncertain life history specified and created
  if(!is.null(StochasticObj)){
    #Find LH params that are not null
    LHList<-names(LHdev[!unlist(lapply(LHdev, is.null))])
    if(NROW(LHList) > 0) {
      print(paste("Uncertainty in life history parameters:", LHList))
    } else {
      print(paste("Uncertainty in life history parameters:", "none"))
    }
  }

  #Check to see if uncertain fishery selectivity specified and created
  #Historical
  if(!is.null(StochasticObj)){
    #Find LH params that are not null
    selListHist<-names(Sdev$hist[!unlist(lapply(Sdev$hist, is.null))])
    if(NROW(selListHist) > 0) {
      print(paste("Uncertainty in historical fishery selectivity parameters:", selListHist))
    } else {
      print(paste("Uncertainty in historical fishery selectivity parameters:", "none"))
    }
    if(NROW(selListHist) > 0 & is.null(HistFisheryObj)) print("Uncertainty in historical fishery selectivity cannot be specified without also specifying HistFisheryObj")
  }

  #Projection
  if(!is.null(StochasticObj)){
    #Find params that are not null

    for(i in 1:TimeAreaObj@areas){
      selListPro<-names(Sdev$pro[[i]][!unlist(lapply(Sdev$pro[[i]], is.null))])
      if(NROW(selListPro) > 0) {
        print(paste("Area", i, "uncertainty in projection fishery selectivity parameters:", selListPro))
      } else {
        print(paste("Area", i, "uncertainty in projection fishery selectivity parameters:", "none"))
      }
      if(NROW(selListPro) > 0 & is.null(ProFisheryObj_list)) print("Uncertainty in projection fishery selectivity cannot be specified without also specifying ProFisheryObj")
    }
  }

  #Check to see if uncertain historical effort created
  if(!is.null(StochasticObj)){
    if(length(StochasticObj@histEffortSD) > 1) {
      print(paste("Uncertainty in historical effort:", StochasticObj@histEffortSD[1], "to", StochasticObj@histEffortSD[2], "created."))
    } else {
      print(paste("Uncertainty in historical effort:", "none"))
    }
  }

  #----------------------------------------------
  #Input checks that stop the program
  #----------------------------------------------
  proceedMSE<-TRUE

  #new additions for multfleet validations
  if(proceedMSE && is_multifleet) {
    # validate fleet selectivity objects
    # for(f in 1:nfleets) {
    #   if(is.null(MultifleetObj@fleet_selectivity_list[[f]])) {
    #     proceedMSE<-FALSE
    #     print(paste("Fleet", f, "selectivity object is missing"))
    #   }
    # }

    #changed
    for(f in 1:nfleets) {
      if(is.null(MultifleetObj@fleet_selectivity_hist_list[[f]]) ||
         is.null(MultifleetObj@fleet_selectivity_proj_list[[f]])) {
        proceedMSE<-FALSE
        print(paste("Fleet", f, "historical or projection selectivity object is missing"))
      }
    }


    # validate allocation type
    if(!MultifleetObj@allocation_type %in% c("effort", "catch")) {
      proceedMSE<-FALSE
      print("allocation_type must be 'effort' or 'catch'")
    }
  }


  #Is iterations specified correctly? (existing validation)
  if(proceedMSE && length(TimeAreaObj@iterations) == 0 ||
     TimeAreaObj@iterations < 1) {
      proceedMSE<-FALSE
      print("Iterations not specified correctly.")
  }

  #Can life history wrapper and sel wrapper be created for each iteration?
  if(proceedMSE){
    LHList<-names(LHdev[!unlist(lapply(LHdev, is.null))])
    selListHist<-names(Sdev$hist[!unlist(lapply(Sdev$hist, is.null))])
    selListPro<-lapply(1:TimeAreaObj@areas, function(x){
      names(Sdev$pro[[x]][!unlist(lapply(Sdev$pro[[x]], is.null))])
    })

    if(NROW(LHList) > 0 | NROW(selListHist) > 0 | NROW(unlist(selListPro)) > 0) {
      for(k in 1:floor(TimeAreaObj@iterations)){
        #LH
        LifeHistoryObj_TMP<-LifeHistoryObj
        if(NROW(LHList) > 0) {
          for(x in 1:NROW(LHList)) slot(LifeHistoryObj_TMP, LHList[x]) <- LHdev[[LHList[x]]][k]
        }

        #Hist sel
        HistFisheryObj_TMP<-HistFisheryObj
        if(NROW(selListHist) > 0){
          for(x in 1:NROW(selListHist)) slot(HistFisheryObj_TMP, selListHist[x]) <- Sdev$hist[[selListHist[x]]][k,]
        }

        #Pro sel
        ProFisheryObj_TMP<-lapply(1:TimeAreaObj@areas, function(x){
          TMP<-ProFisheryObj_list[[x]]
          if(NROW(selListPro[[x]]) > 0){
            for(y in 1:NROW(selListPro[[x]])) slot(TMP, selListPro[[x]][y]) <- Sdev$pro[[x]][[selListPro[[x]][y]]][k,]
          }
          TMP
        })

        #Setup
        lh<-LHwrapper(LifeHistoryObj_TMP, TimeAreaObj)
        selHist<-lapply(1:TimeAreaObj@areas, function(x){
          selWrapper(lh, TimeAreaObj, FisheryObj = HistFisheryObj_TMP, doPlot = FALSE)
        })
        selPro<-lapply(1:TimeAreaObj@areas, function(x){
          selWrapper(lh, TimeAreaObj, FisheryObj = ProFisheryObj_TMP[[x]], doPlot = FALSE)
        })
        if(is.null(lh)) {
          proceedMSE<-FALSE
          print(paste("Life history cannot be created. Check inputs. Stopped at interation", k))
        }
        for(x in 1:TimeAreaObj@areas){
          if(is.null(selHist[[x]])){
            proceedMSE<-FALSE
            print(paste("Historical selectivity cannot be created. Check inputs. Stopped at interation", k))
          }
        }
        for(x in 1:TimeAreaObj@areas){
          if(isTRUE(!is.null(StrategyObj) &  is.null(selPro[[x]]))){
            proceedMSE<-FALSE
            print(paste("Projection selectivity cannot be created. Check inputs. Stopped at interation", k))
          }
        }
      }
    } else {
      lh<-LHwrapper(LifeHistoryObj, TimeAreaObj)
      selHist<-lapply(1:TimeAreaObj@areas, function(x){
        selWrapper(lh, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = FALSE)
      })
      selPro<-lapply(1:TimeAreaObj@areas, function(x){
        selWrapper(lh, TimeAreaObj, FisheryObj = ProFisheryObj_list[[x]], doPlot = FALSE)
      })
      if(is.null(lh)) {
        proceedMSE<-FALSE
        print("Life history cannot be created. Check inputs.")
      }
      for(x in 1:TimeAreaObj@areas){
        if(is.null(selHist[[x]])){
          proceedMSE<-FALSE
          print("Historical selectivity cannot be created. Check inputs.")
        }
      }
      for(x in 1:TimeAreaObj@areas){
        if(isTRUE(!is.null(StrategyObj) &  is.null(selPro[[x]]))){
          proceedMSE<-FALSE
          print("Projection selectivity cannot be created. Check inputs.")
        }
      }
    }
  }

  #Rec devs failed
  if(proceedMSE && is.null(RdevMatrix)) {
    proceedMSE<-FALSE
    print("Inter-annual recruitment variation cannot be created. Check inputs.")
  }

  #Ddev
  if(proceedMSE && is.null(Ddev)) {
    proceedMSE<-FALSE
    print("Initial biomass variation cannot be created. Check inputs.")
  }

  #Cdev
  if(proceedMSE &&
     is(StrategyObj,"Strategy")  &&
     StrategyObj@projectionName == "projectionStrategy" &&
     is.null(Cdev)) {
    proceedMSE<-FALSE
    print("Initial CPUE variation cannot be created. Check inputs.")
  }

  #Edev
  if(proceedMSE &&
     is(StrategyObj,"Strategy")  &&
     StrategyObj@projectionName == "projectionStrategy" &&
     is.null(Edev)) {
    proceedMSE<-FALSE
    print("Effort implementation error cannot be created. Check inputs.")
  }

  #Historical effort devs failed
  if(proceedMSE && is.null(histEffortDev)) {
    proceedMSE<-FALSE
    print("Inter-annual historical effort variation cannot be created. Check inputs.")
  }

  #Number of areas not in agreement with dimensions of the move matrix.
  if(proceedMSE && isTRUE(TimeAreaObj@areas != dim(TimeAreaObj@move)[1] | TimeAreaObj@areas != dim(TimeAreaObj@move)[2])) {
    proceedMSE<-FALSE
    print("Number of areas not in agreement with dimensions of the move matrix. Check inputs.")
  }

  #Number of historical years does not match historical effort time series
  if(proceedMSE && isTRUE(TimeAreaObj@historicalYears > 0  && isTRUE(TimeAreaObj@areas != dim(TimeAreaObj@historicalEffort)[2] | TimeAreaObj@historicalYears != dim(TimeAreaObj@historicalEffort)[1]))) {
    proceedMSE<-FALSE
    print("Number of historical years does not match historical effort time series. Check inputs.")
  }

  #Historical and/or projection years must be at least 1.
  if(proceedMSE && isTRUE(TimeAreaObj@historicalYears + ifelse(is(StrategyObj, "Strategy")  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0) < 1)) {
    proceedMSE<-FALSE
    print("Historical and/or projection years must be at least 1. Check inputs.")
  }

  #Project strategy function missing
  if(proceedMSE && isTRUE(is(StrategyObj, "Strategy") &&
                             tryCatch({
                               get(StrategyObj@projectionName)
                               FALSE
                             }, error = function(c) TRUE)
  )) {
    proceedMSE<-FALSE
    print("Project strategy function missing. StrategyObj@projectionName must correspond to a named function.")
  }

  #When applying a management strategy, StrategyObj@projectionYears must be greater than 0
  if(proceedMSE &&
     isTRUE(is(StrategyObj, "Strategy") && length(StrategyObj@projectionYears) == 0) ||
     isTRUE(is(StrategyObj, "Strategy") && StrategyObj@projectionYears < 1)
  ){
    proceedMSE<-FALSE
    print("When applying a management strategy, StrategyObj@projectionYears must be greater than 0")
  }


  #---------------------------
  #Setup parallel processing
  #---------------------------

  #new addition: Only adding MultifleetObj = MultifleetObj

  #Test whether we can proceed to simulations
  if(
    isFALSE(proceedMSE)
  ) {
    warning("One or more components contain incomplete or erroneous information. Cannot proceed to simulation.")
    return(NULL)
  } else {

    ptm<-proc.time()
    #require(snowfall)
    #require(parallel)
    iterations <- floor(TimeAreaObj@iterations)

    if(detectCores() > 3 && iterations >= (detectCores() - 2) && is.null(waitName) && is.null(hostName)) {
      print("Running on multiple cores")
      cores<-min(iterations, (detectCores()-2))
      sfInit(parallel=T, cpus=cores)
      sfLibrary(fishSimGTG)
      if(!is.null(customToCluster)) sfExport(list = returnValue(customToCluster))
      input<-list()
      inputObject<-list()
      size<-floor(iterations/cores)
      for (i in 1:cores){
        input[[i]]<-c(size*(i-1)+1, ifelse(i==cores, iterations, size*i))
        inputObject[[i]]<-list(iter=c(size*(i-1)+1, ifelse(i==cores, iterations, size*i)),
                               RdevMatrix = RdevMatrix,
                               Ddev = Ddev,
                               Cdev = Cdev,
                               Edev = Edev,
                               LHdev = LHdev,
                               Sdev = Sdev,
                               histEffortDev = histEffortDev,
                               LifeHistoryObj = LifeHistoryObj,
                               TimeAreaObj = TimeAreaObj,
                               HistFisheryObj = HistFisheryObj,
                               ProFisheryObj_list = ProFisheryObj_list,
                               StrategyObj = StrategyObj,
                               StochasticObj = StochasticObj,
                               MultifleetObj = MultifleetObj,
                               IndexObj= IndexObj,
                               CatchObsObj= CatchObsObj,
                               LengthCompObj= LengthCompObj,
                               iterations=iterations,
                               waitName=waitName,
                               hostName=hostName,
                               doDiagnostic=doDiagnostic)
      }
      mseParallel<-sfLapply(inputObject, evalMSE)
      sfRemoveAll()
      sfStop()

      #-------------------------------
      #Ressemble from multiple cores
      #-------------------------------
      SB<-mseParallel[[1]]$dynamics$SB
      VB<-mseParallel[[1]]$dynamics$VB
      RB<-mseParallel[[1]]$dynamics$RB
      catchB<-mseParallel[[1]]$dynamics$catchB
      catchN<-mseParallel[[1]]$dynamics$catchN
      Ftotal<-mseParallel[[1]]$dynamics$Ftotal
      discB<-mseParallel[[1]]$dynamics$discB
      discN<-mseParallel[[1]]$dynamics$discN
      SPR<-mseParallel[[1]]$dynamics$SPR
      relSB<-mseParallel[[1]]$dynamics$relSB
      recN<-mseParallel[[1]]$dynamics$recN
      ref<-mseParallel[[1]]$dynamics$ref
      decisionAnnual<-mseParallel[[1]]$HCR$decisionAnnual
      decisionLocal<-mseParallel[[1]]$HCR$decisionLocal
      decisionData<-mseParallel[[1]]$HCR$decisionData

      # new addition: extract multifleet results if they exist
      multifleet_results <- NULL
      if(!is.null(mseParallel[[1]]$dynamics$multifleet)) {
        multifleet_results <- mseParallel[[1]]$dynamics$multifleet
        Ftotal_by_fleet <- multifleet_results$Ftotal_by_fleet
        catchB_by_fleet <- multifleet_results$catchB_by_fleet
        catchN_by_fleet <- multifleet_results$catchN_by_fleet
        discB_by_fleet <- multifleet_results$discB_by_fleet
        discN_by_fleet <- multifleet_results$discN_by_fleet
      }


      #Optional diagnostic outputs
      N<-mseParallel[[1]]$N
      catchNage<-mseParallel[[1]]$catchNage

      for (i in 2:cores){
        for(m in 1:TimeAreaObj@areas){
          SB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$SB[,input[[i]][1]:input[[i]][2],m]
          VB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$VB[,input[[i]][1]:input[[i]][2],m]
          RB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$RB[,input[[i]][1]:input[[i]][2],m]
          catchB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$catchB[,input[[i]][1]:input[[i]][2],m]
          catchN[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$catchN[,input[[i]][1]:input[[i]][2],m]
          Ftotal[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$Ftotal[,input[[i]][1]:input[[i]][2],m]
          discB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$discB[,input[[i]][1]:input[[i]][2],m]
          discN[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$discN[,input[[i]][1]:input[[i]][2],m]

          #new addition: reassemble multifleet arrays
          if(!is.null(multifleet_results)) {
            for(f in 1:multifleet_results$nfleets) {
              Ftotal_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$Ftotal_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
              catchB_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$catchB_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
              catchN_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$catchN_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
              discB_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$discB_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
              discN_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$discN_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
            }
          }
        }

          # continue with the rest of arraya (no modification needed)
        SPR[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$SPR[,input[[i]][1]:input[[i]][2]]
        relSB[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$relSB[,input[[i]][1]:input[[i]][2]]
        recN[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$recN[,input[[i]][1]:input[[i]][2]]
        ref[input[[i]][1]:input[[i]][2],]<-mseParallel[[i]]$dynamics$ref[input[[i]][1]:input[[i]][2],]
        decisionAnnual<-rbind(decisionAnnual, mseParallel[[i]]$HCR$decisionAnnual)
        decisionLocal<-rbind(decisionLocal, mseParallel[[i]]$HCR$decisionLocal)
        decisionData<-rbind(decisionData, mseParallel[[i]]$HCR$decisionData)
        }

        #new addition: reconstruct multifleet results structure

        if(!is.null(multifleet_results)) {
          multifleet_results$Ftotal_by_fleet <- Ftotal_by_fleet
          multifleet_results$catchB_by_fleet <- catchB_by_fleet
          multifleet_results$catchN_by_fleet <- catchN_by_fleet
          multifleet_results$discB_by_fleet <- discB_by_fleet
          multifleet_results$discN_by_fleet <- discN_by_fleet
        }




    } else {
      #single core processing
      mse<-evalMSE(inputObject=list(iter=c(1, iterations),
                                    RdevMatrix=RdevMatrix,
                                    Ddev=Ddev,
                                    Cdev = Cdev,
                                    Edev = Edev,
                                    LHdev = LHdev,
                                    Sdev = Sdev,
                                    histEffortDev = histEffortDev,
                                    LifeHistoryObj = LifeHistoryObj,
                                    TimeAreaObj = TimeAreaObj,
                                    HistFisheryObj = HistFisheryObj,
                                    ProFisheryObj_list = ProFisheryObj_list,
                                    StrategyObj = StrategyObj,
                                    StochasticObj = StochasticObj,
                                    MultifleetObj = MultifleetObj, #new addition
                                    IndexObj = IndexObj,
                                    CatchObsObj= CatchObsObj,
                                    LengthCompObj= LengthCompObj,
                                    iterations=iterations,
                                    waitName=waitName,
                                    hostName=hostName,
                                    doDiagnostic=doDiagnostic
                                    )
                   )

      SB<-mse$dynamics$SB
      VB<-mse$dynamics$VB
      RB<-mse$dynamics$RB
      catchB<-mse$dynamics$catchB
      catchN<-mse$dynamics$catchN
      Ftotal<-mse$dynamics$Ftotal
      discB<-mse$dynamics$discB
      discN<-mse$dynamics$discN
      SPR<-mse$dynamics$SPR
      relSB<-mse$dynamics$relSB
      recN<-mse$dynamics$recN
      ref<-mse$dynamics$ref
      decisionAnnual<-mse$HCR$decisionAnnual
      decisionLocal<-mse$HCR$decisionLocal
      decisionData<-mse$HCR$decisionData

      # new addition: extract multifleet results for single core mode
      multifleet_results <- NULL
      if(!is.null(mse$dynamics$multifleet)) {
        multifleet_results <- mse$dynamics$multifleet
      }

      #Optional diagnostic outputs
      N<-mse$N
      catchNage<-mse$catchNage
    }

    #---------------
    #Save results
    #---------------
    dynamics<-list(SB=SB, VB=VB, RB=RB, catchB=catchB, catchN=catchN, Ftotal=Ftotal, discB=discB, discN=discN, SPR=SPR, relSB=relSB, recN=recN, ref=ref, N=N, catchNage=catchNage)

    #new addditon: adding multifleet results if they exist
    if(!is.null(multifleet_results)) {
      dynamics$multifleet <- multifleet_results
    }

    HCR<-list(decisionLocal=decisionLocal, decisionAnnual=decisionAnnual, decisionData=decisionData)

    # adding MultifleetObj
    dt<-list(titleStrategy = titleStrategy, dynamics=dynamics, HCR=HCR, iterations=iterations, LifeHistoryObj=LifeHistoryObj, LHdev=LHdev, Sdev = Sdev, histEffortDev = histEffortDev, Ddev=Ddev, TimeAreaObj=TimeAreaObj, HistFisheryObj=HistFisheryObj, ProFisheryObj_list=ProFisheryObj_list,  StrategyObj= StrategyObj, StochasticObj=StochasticObj, MultifleetObj=MultifleetObj, IndexObj= IndexObj, CatchObsObj= CatchObsObj, LengthCompObj = LengthCompObj)
    saveRDS(dt, file=paste(wd, "/", fileName, ".rds", sep=""))

    #--------------------------------------------------------------------------------
    #Plot results (mostly for diagnostics, these are ugly - not publication quality)
    #--------------------------------------------------------------------------------
    if(doPlot) {

      rb<-rainbow(iterations)

      #Population level plots
      png(filename=paste(wd, "/", fileName, "_SPR.png",sep=""), width=4, height=4, units="in", res=300, bg="white", pointsize=12)
      par(mfrow=c(1,1), mar=c(4,4,3,1))
      plot(dt$dynamics$SPR[,1], type="l", las=1, ylab="", xlab = "Year", ylim=c(0,1), col=rb[1], main = "SPR")
      if(iterations > 1){
        for(k in 2:iterations){
          lines(dt$dynamics$SPR[,k], col=rb[k])
        }
      }
      dev.off()

      png(filename=paste(wd, "/", fileName, "_SBrel.png",sep=""), width=4, height=4, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(1,1), mar=c(4,4,3,1))
      plot(dt$dynamics$relSB[,1], type="l", las=1, ylab="", xlab = "Year", ylim=c(0,1), col=rb[1], main = "Relative spawning biomass")
      if(iterations > 1){
        for(k in 2:iterations){
          lines(dt$dynamics$relSB[,k], col=rb[k])
        }
      }
      dev.off()

      #Recruit over time
      png(filename=paste(wd, "/", fileName, "_recN.png",sep=""), width=4, height=4, units="in", res=300, bg="white", pointsize=12)
      par(mfrow=c(1,1), mar=c(4,4,3,1))
      plot(dt$dynamics$recN[,1], type="l", las=1, ylab="", xlab = "Year", ylim=c(min(dt$dynamics$recN),max(dt$dynamics$recN)), col=rb[1], main = "recruits N")
      if(iterations > 1){
        for(k in 2:iterations){
          lines(dt$dynamics$recN[,k], col=rb[k])
        }
      }
      dev.off()


      #S-R (need to reviw this plot)
      # png(filename=paste(wd, "/", fileName, "_SR.png",sep=""), width=4, height=4, units="in", res=300, bg="white", pointsize=12)
      # par(mfrow=c(1,1), mar=c(4,4,3,1))
      #
      # is<-solveD(lh, sel = selHist, doFit = FALSE, F_in = 0.01)
      # SRcurve<-t(sapply(seq(0, is$B0, length.out = 100), FUN=function(x){
      #   c(x/is$B0, recruit(LifeHistoryObj=dt$LifeHistoryObj, B0=is$B0, stock=x, forceR=FALSE, Rforced=0))
      # }))
      # plot(dt$dynamics$relSB[,1], dt$dynamics$recN[,1], type="b", las=1, ylim=c(min(dt$dynamics$recN),max(dt$dynamics$recN)), col=rb[1], ylab="Recruits", xlab = "Stock (rel SSB)", main = "Stock-recruit")
      # #text(dt$dynamics$relSB[,1], dt$dynamics$recN[,1], labels=1:NROW(dt$dynamics$recN[,1]))
      # if(iterations > 1){
      #   for(k in 2:iterations){
      #     lines(dt$dynamics$relSB[,k], dt$dynamics$recN[,k], type="b", col=rb[k])
      #     #text(dt$dynamics$relSB[,k], dt$dynamics$recN[,k], labels=1:NROW(dt$dynamics$recN[,k]))
      #   }
      # }
      # lines(SRcurve[,1], SRcurve[,2], type="l", col="black", las=1, ylab="Recruits", xlab = "Stock (rel SSB)", main = "Stock-recruit")
      #
      # dev.off()

      #----------------------
      # Area specific plots
      #---------------------

      #SSB
      png(filename=paste0(wd, "/", fileName, "_SB_Area.png"), width=9, height=3.5*ceiling(dt$TimeAreaObj@areas/2), units="in", res=96, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$SB[,1,m], type="l", las=1, ylab="", xlab = "Year", col=rb[1], ylim=c(min(dt$dynamics$SB[,,m]), max(dt$dynamics$SB[,,m])), main = "Spawning biomass")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$SB[,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()


      #Catch
      png(filename=paste0(wd, "/", fileName, "_catchB_Area.png"), width=9, height=3.5*ceiling(dt$TimeAreaObj@areas/2), units="in", res=96, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$catchB[,1,m], type="l", las=1, ylab="", xlab = "Year", col=rb[1], ylim=c(min(dt$dynamics$catchB[,,m]), max(dt$dynamics$catchB[,,m])), main = "catch in weight")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$catchB[,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()

      #F
      png(filename=paste0(wd, "/", fileName, "_F_Area.png"), width=9, height=3.5*ceiling(dt$TimeAreaObj@areas/2), units="in", res=96, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$Ftotal[,1,m], type="l", las=1, ylab="", xlab = "Year", col=rb[1], ylim=c(min(dt$dynamics$Ftotal[,,m]), max(dt$dynamics$Ftotal[,,m])), main = "Fishing mortality")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$Ftotal[,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()


      #rel change in SSB
      png(filename=paste0(wd, "/", fileName, "_SBchange_Area.png"), width=9, height=3.5*ceiling(dt$TimeAreaObj@areas/2), units="in", res=96, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$SB[,1,m]/dt$dynamics$SB[1,1,m], type="l", las=1, ylab="", col=rb[1], ylim=c(0,2), main = "Relative spawning biomass")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$SB[,k,m]/dt$dynamics$SB[1,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()

      # new addition: add multifleet-specific plots if multifleet results exist
      if(!is.null(multifleet_results)) {
        # Fleet-specific F plots
        png(filename=paste0(wd, "/", fileName, "_F_by_Fleet.png"), width=12, height=8, units="in", res=96, bg="white", pointsize=12)
        par(mfrow=c(ceiling(multifleet_results$nfleets/2), 2), mar=c(4,4,3,1))
        for(f in 1:multifleet_results$nfleets) {
          for(m in 1:dt$TimeAreaObj@areas) {
            if(f == 1 && m == 1) {
              plot(multifleet_results$Ftotal_by_fleet[,1,m,f], type="l", las=1, ylab="F", xlab = "Year",
                   col=rb[1], main = paste("Fleet", f, "Area", m))
            } else {
              plot(multifleet_results$Ftotal_by_fleet[,1,m,f], type="l", las=1, ylab="F", xlab = "Year",
                   col=rb[1], main = paste("Fleet", f, "Area", m))
            }
            if(iterations > 1){
              for(k in 2:iterations){
                lines(multifleet_results$Ftotal_by_fleet[,k,m,f], col=rb[k])
              }
            }
          }
        }
        dev.off()

        # Fleet-specific catch plots
        png(filename=paste0(wd, "/", fileName, "_Catch_by_Fleet.png"), width=12, height=8, units="in", res=96, bg="white", pointsize=12)
        par(mfrow=c(ceiling(multifleet_results$nfleets/2), 2), mar=c(4,4,3,1))
        for(f in 1:multifleet_results$nfleets) {
          for(m in 1:dt$TimeAreaObj@areas) {
            plot(multifleet_results$catchB_by_fleet[,1,m,f], type="l", las=1, ylab="Catch", xlab = "Year",
                 col=rb[1], main = paste("Fleet", f, "Area", m))
            if(iterations > 1){
              for(k in 2:iterations){
                lines(multifleet_results$catchB_by_fleet[,k,m,f], col=rb[k])
              }
            }
          }
        }
        dev.off()
      }
    }
    print("Simulation time in minutes: ")
    print((proc.time()-ptm)/60)
  }
}



#---------------------------------------
#Read in data from existing projection or MSE model
#---------------------------------------

#Roxygen header
#'Read in data from existing projection or MSE model
#'
#'Function for running projections or MSE
#'
#' @param wd A working directly where output is saved. Required
#' @param fileName A file name previously used to create output. Required
#' @export

readProjection<-function(wd, fileName){
  readRDS(file=paste(wd, "/", fileName, ".rds", sep=""))
}
