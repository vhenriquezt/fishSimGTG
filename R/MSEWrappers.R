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

  #--------------------
  #1. Unpack dataObject
  #--------------------
  TimeAreaObj <- StrategyObj <- LifeHistoryObj <- HistFisheryObj <- ProFisheryObj_list <- iterations <- iter <- Ddev <- Edev <- LHdev <- Sdev <- Cdev <- Edev <- histEffortDev <- RdevMatrix <- doDiagnostic <- MultifleetObj <- NULL
  for(r in 1:NROW(inputObject)) assign(names(inputObject)[r], inputObject[[r]])

  #The code has two paths: single and multifleet
  is_multifleet <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1

  #---------------
  #2. Arrays setup
  #---------------

  controlRuleYear<-c(FALSE, rep(FALSE,(TimeAreaObj@historicalYears)), rep(TRUE, ifelse(is(StrategyObj, "Strategy")  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)))
  years <- 1 + TimeAreaObj@historicalYears + ifelse(is(StrategyObj, "Strategy")  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)
  areas <- TimeAreaObj@areas

  SB<-array(dim=c(years, iterations, areas))
  catchN<-array(dim=c(years, iterations, areas))
  catchB<-array(dim=c(years, iterations, areas))
  discN<-array(dim=c(years, iterations, areas))
  discB<-array(dim=c(years, iterations, areas))
  SPR<-array(dim=c(years, iterations))
  relSB<-array(dim=c(years, iterations))
  recN<-array(dim=c(years, iterations))

  #Single fleet
  if(!is_multifleet) {
    VB<-array(dim=c(years, iterations, areas))
    RB<-array(dim=c(years, iterations, areas))
    Ftotal<-array(dim=c(years, iterations, areas))
  }

  #Multifleet 4D arrays [years, iterations, areas, fleets]
  if(is_multifleet) {

    nfleets <- MultifleetObj@nfleets
    fleet_proportions <- MultifleetObj@fleet_proportions
    cat("Multifleet mode detected with", nfleets, "fleets\n")
    Ftotal_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    catchB_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    catchN_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    discB_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    discN_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    RB_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    VB_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
    Ftotal_by_fleet[] <- NA
    catchB_by_fleet[] <- NA
    catchN_by_fleet[] <- NA
    discB_by_fleet[] <- NA
    discN_by_fleet[] <- NA
    final_effort_proportions <- array(dim=c(iterations, nfleets))
    target_catch_proportions <- array(dim=c(iterations, nfleets))
    actual_catch_proportions <- array(dim=c(iterations, nfleets))
    allocation_type <- array(dim=iterations)
  }

  #Diagnostic mode
  Nexport<-NULL
  catchNageExport<-NULL
  Zexport<-NULL

  #Recording of management strategy details
  decisionData<-data.frame()
  decisionAnnual<-data.frame()
  decisionLocal<-data.frame()

  #Setup capture of benchmarks
  ref<-array(dim = c(iterations, 10))

  #-------------------------------------------
  #3. Deteministic LH and Sel, if present (saves computation time, if present)
  #-------------------------------------------

  #LH uncertainty list
  LHList<-names(LHdev[!unlist(lapply(LHdev, is.null))])

  #Selectivity uncertainty list
  #Single fleet
  if(!is_multifleet){
    selListHist<-names(Sdev$hist[!unlist(lapply(Sdev$hist, is.null))])
    selListPro<-lapply(1:TimeAreaObj@areas, function(x){
      names(Sdev$pro[[x]][!unlist(lapply(Sdev$pro[[x]], is.null))])
    })
  }
  #Multi fleet
  if(is_multifleet){
    selListHist <- character(0)  # empty character vector
    selListPro <- replicate(TimeAreaObj@areas, character(0), simplify = FALSE)
  }

  #Is deterministic, so save time by make calculations only once.
  if(NROW(LHList) == 0 & NROW(selListHist) == 0 & NROW(unlist(selListPro)) == 0){
    lh<-LHwrapper(LifeHistoryObj, TimeAreaObj)
    ageClasses <- lh$ageClasses
    if(!is.null(lh) & lh$LifeHistory@Steep < 0.21) lh$LifeHistory@Steep <- 0.21
    if(!is.null(lh) & lh$LifeHistory@Steep > 1) lh$LifeHistory@Steep <- 1

    #Single fleet
    if(!is_multifleet){
      selHist<-lapply(1:TimeAreaObj@areas, function(x){
        selWrapper(lh, TimeAreaObj, FisheryObj = HistFisheryObj)
      })
      selPro<-lapply(1:TimeAreaObj@areas, function(x){
        selWrapper(lh, TimeAreaObj, FisheryObj = ProFisheryObj_list[[x]])
      })
      refCalc<-gtgYPRWrapper_Fonly(lh=lh, sel=selHist[[1]])
      for(k in iter[1]:iter[2]) ref[k, ]<-as.matrix(refCalc$sim)[1,]
      colnames(ref)<-names(refCalc$sim)
    }

    #Multi fleet
    if(is_multifleet){
      #Bill Edit
      selHistBlock <- lapply(1:TimeAreaObj@areas, function(x){
        lapply(1:nfleets, function(f) {
          bl<-NROW(MultifleetObj@fleet_block_hist_list[[f]])
          lapply(1:bl, function(b) {
            selWrapper(lh, TimeAreaObj,
                      FisheryObj = MultifleetObj@fleet_selectivity_hist_list[[f]][[b]])
          })
        })
      })
      selPro <- lapply(1:TimeAreaObj@areas, function(x){
        lapply(1:nfleets, function(f) {
          if(is(StrategyObj, "Strategy")){
            selWrapper(lh, TimeAreaObj,
                       FisheryObj = MultifleetObj@fleet_selectivity_proj_list[[x]][[f]])
          } else {
            NULL
          }
        })
      })

      #Note: I am using area 1 and fleet 1 for this calculation (it could be changed)
      #Bill Edit
      indX<-sapply(1:NROW(MultifleetObj@fleet_block_hist_list[[1]]), function(x){1 %in% MultifleetObj@fleet_block_hist_list[[1]][[x]]})
      refCalc<-gtgYPRWrapper_Fonly(lh=lh, sel=selHistBlock[[1]][[1]][[which(indX)[1]]])
    }
  }


  #------------------------------
  #4. Run simulator of k iterations
  #------------------------------

  #Used when fishSimGTG runs in a shiny app
  if(!is.null(hostName) & !is.null(waitName)){
    waitName$show()
  }

  #step through iterations k
  for(k in iter[1]:iter[2]){

    #-------------------------------------------------------------------
    #5. Setup iteration-specific life history & selectivity (if present)
    #-------------------------------------------------------------------
    if(NROW(LHList) > 0 | NROW(selListHist) > 0 | NROW(unlist(selListPro)) > 0){

      LifeHistoryObj_TMP<-LifeHistoryObj
      if(NROW(LHList) > 0){
        for(x in 1:NROW(LHList)) slot(LifeHistoryObj_TMP, LHList[x]) <- LHdev[[LHList[x]]][k]
      }
      lh<-LHwrapper(LifeHistoryObj_TMP, TimeAreaObj)
      ageClasses <- lh$ageClasses
      if(!is.null(lh) & lh$LifeHistory@Steep < 0.21) lh$LifeHistory@Steep <- 0.21
      if(!is.null(lh) & lh$LifeHistory@Steep > 1) lh$LifeHistory@Steep <- 1

      #Single fleet
      if(!is_multifleet){
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
        #Hist sel
        selHist<-lapply(1:TimeAreaObj@areas, function(x){
          selWrapper(lh, TimeAreaObj, FisheryObj = HistFisheryObj_TMP, doPlot = FALSE)
        })
        #Pro sel
        selPro<-lapply(1:TimeAreaObj@areas, function(x){
          selWrapper(lh, TimeAreaObj, FisheryObj = ProFisheryObj_TMP[[x]], doPlot = FALSE)
        })
        refCalc<-gtgYPRWrapper_Fonly(lh=lh, sel=selHist[[1]])
      }

      #Multi fleet
      if(is_multifleet) {
        #Fleet-specific fishery objects directly from MultifleetObj (no stochasticity)
        #Bill Edit
        selHistBlock <- lapply(1:TimeAreaObj@areas, function(x){
          lapply(1:nfleets, function(f) {
            bl<-NROW(MultifleetObj@fleet_block_hist_list[[f]])
            lapply(1:bl, function(b) {
              selWrapper(lh, TimeAreaObj,
                         FisheryObj = MultifleetObj@fleet_selectivity_hist_list[[f]][[b]])
            })
          })
        })

        selPro<-lapply(1:TimeAreaObj@areas, function(area){
          lapply(1:nfleets, function(f) {
            if(is(StrategyObj, "Strategy")){
              selWrapper(lh, TimeAreaObj,
                         FisheryObj = MultifleetObj@fleet_selectivity_proj_list[[x]][[f]])
            } else {
              NULL
            }
          })
        })
        #Bill Edit
        indX<-sapply(1:NROW(MultifleetObj@fleet_block_hist_list[[1]]), function(x){1 %in% MultifleetObj@fleet_block_hist_list[[1]][[x]]})
        refCalc<-gtgYPRWrapper_Fonly(lh=lh, sel=selHistBlock[[1]][[1]][[which(indX)[1]]])
      }

      ref[k, ]<-as.matrix(refCalc$sim)[1,]
      colnames(ref)<-names(refCalc$sim)
    }


    #-----------------------------------
    #6. Initial equilibrium calculations
    #-----------------------------------

    #Single feet
    if(!is_multifleet) {
      is<-solveD(lh, sel = selHist[[1]], doFit = TRUE, D_type = TimeAreaObj@historicalBioType, D_in = Ddev[k])
    }

    #Multi fleet
    if(is_multifleet) {

      #Initial equilibrium (before movement between areas)
      #Bill Edit
      selHist<-lapply(1:TimeAreaObj@areas, function(area){
        lapply(1:nfleets, function(f) {
          indX<-sapply(1:NROW(MultifleetObj@fleet_block_hist_list[[f]]), function(x){1 %in% MultifleetObj@fleet_block_hist_list[[f]][[x]]})
          selHistBlock[[area]][[f]][[which(indX)[1]]]
        })
      })

      is <- solveD_multifleet2(lh = lh, sel_list = selHist[[1]], doFit = TRUE, D_type = TimeAreaObj@historicalBioType,
                               D_in = Ddev[k], fleet_proportions = fleet_proportions,
                               allocation_type = MultifleetObj@allocation_type)

      #Store the final proportions to report after runProjection()
      if(is$allocation_type == "effort") final_effort_proportions[k,] <- is$final_effort_proportions
      if(is$allocation_type == "catch"){
        target_catch_proportions[k,] <- is$target_catch_proportions
        actual_catch_proportions[k,] <- is$actual_catch_proportions
      }
      allocation_type[k] <- is$allocation_type

      #debugging
      #CRITICAL: Verify values immediately after equilibrium calculation
      cat("=== CRITICAL DEBUG in evalMSE (iteration", k, ") ===\n")
      cat("IMMEDIATELY after solveD_multifleet2() call:\n")
      cat("  is$final_effort_proportions:", if(!is.null(is$final_effort_proportions)) round(is$final_effort_proportions, 4) else "NULL", "\n")
      cat("  is$actual_catch_proportions:", if(!is.null(is$actual_catch_proportions)) round(is$actual_catch_proportions, 4) else "NULL", "\n")
      cat("  is$target_catch_proportions:", if(!is.null(is$target_catch_proportions)) round(is$target_catch_proportions, 4) else "NULL", "\n")
      cat("  is$allocation_type:", is$allocation_type, "\n")

      if(!is.null(is$final_effort_proportions) && !is.null(is$actual_catch_proportions)) {
        cat("  Are effort/catch identical in evalMSE?", identical(is$final_effort_proportions, is$actual_catch_proportions), "\n")
        cat("  Max difference:", max(abs(is$final_effort_proportions - is$actual_catch_proportions)), "\n")
      }
      cat("================================================\n")
    }

    #Burn-in to calibrate N by area, noting effect of movement (this is for area distribution)
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
        if(j< yrsTmp){
          P<-matrix(nrow=ageClasses*areas, ncol=ageClasses*areas)
          for(m in 1:areas){
            #Single fleet
            if(!is_multifleet) {
              S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=is$Feq, S_in=selHist[[m]]$removal[[l]] )
            }
            #Multi fleet
            if(is_multifleet) {
              #calculate total Z from all fleet contributions
              #Z = M + sum_across_fleets(F_fleet * selectivity_fleet)
              # NEVER combine selectivities - each fleet contributes independently
              total_fishing_mortality <- sapply(1:ageClasses, function(age) {
                sum(sapply(1:nfleets, function(ff) {
                  is$F_by_fleet[ff] * selHist[[m]][[ff]]$removal[[l]][age]
                }))
              })

              #Note that fishing mortality accounted for in total_fishing_mortality vector, thus F_in = 1 is just a dummy input
              S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=1, S_in=total_fishing_mortality )
            }
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

    #-----------------------------------------------------------------
    #7. Initiate lists and arrays for iteration k and populate j = 1
    #-----------------------------------------------------------------

    ##########
    #Abundance initial conditions in iteration k
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

    #############
    #Catches and discards
    #Catch-at-age list by gtg and Z list (holds arrays by ageClasses, years, areas) in iteration k
    catchNage<-list()
    Z<-list()
    for(l in 1:lh$gtg){
      catchNage[[l]]<-array(dim=c(ageClasses, years, areas))
      Z[[l]]<-array(dim=c(ageClasses, years, areas))
    }

    #Single fleet
    if(!is_multifleet){
      for(m in 1:areas){

        #By-area fishing mortality
        Ftotal[1,k,m] <- is$Feq

        #Calculate catch matrices by gtg and age
        for(l in 1:lh$gtg){
          Z[[l]][,1,m] <- Ftotal[1,k,m]*selHist[[m]]$removal[[l]] + lh$LifeHistory@M
          catchNage[[l]][,1,m] <- Ftotal[1,k,m]*selHist[[m]]$keep[[l]]/(Z[[l]][,1,m])*(1-exp(-Z[[l]][,1,m]))*N[[l]][,1,m]
        }

        #By-area catches and discards
        catchN[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(catchNage[[x]][,1,m])))
        catchB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*catchNage[[x]][,1,m])))
        discN[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal[1,k,m]*selHist[[m]]$discard[[x]]/(Ftotal[1,k,m]*selHist[[m]]$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[1,k,m]*selHist[[m]]$removal[[x]]-lh$LifeHistory@M))*N[[x]][,1,m])))
        discB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal[1,k,m]*selHist[[m]]$discard[[x]]/(Ftotal[1,k,m]*selHist[[m]]$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[1,k,m]*selHist[[m]]$removal[[x]]-lh$LifeHistory@M))*N[[x]][,1,m])))
      }
    }

    #Multi fleet
    if(is_multifleet){
      #Initialize fleet-specific catch arrays for iteration 1
      #Structure: catchNage_by_fleet[[fleet]][[gtg]][age, year, area]
      catchNage_by_fleet <- list()
      for(f in 1:nfleets) {
        catchNage_by_fleet[[f]] <- list()
        for(l in 1:lh$gtg) {
          catchNage_by_fleet[[f]][[l]] <- array(dim=c(ageClasses, years, areas))
        }
      }
      for(m in 1:areas){
        #By-area and each fleet fishing mortality
        Ftotal_by_fleet[1,k,m,] <- is$F_by_fleet
        #Calculate catch matrices by gtg and age
        for(l in 1:lh$gtg){
          #calculate total Z from all fleet contributions
          #Z = M + sum_across_fleets(F_fleet * selectivity_fleet)
          # NEVER combine selectivities - each fleet contributes independently
          total_fishing_mortality <- sapply(1:ageClasses, function(age) {
            sum(sapply(1:nfleets, function(ff) {
              Ftotal_by_fleet[1,k,m,ff] * selHist[[m]][[ff]]$removal[[l]][age]
            }))
          })
          #total mortality shared by all fleets: Z = M + total_fishing_mortality
          Z[[l]][,1,m] <- total_fishing_mortality + lh$LifeHistory@M
          for(f in 1:nfleets){
            #fleet-specific catch (using using Baranov equation with shared Z)
            catchNage_by_fleet[[f]][[l]][,1,m] <- Ftotal_by_fleet[1,k,m,f] * selHist[[m]][[f]]$keep[[l]] /
              Z[[l]][,1,m] * (1-exp(-Z[[l]][,1,m])) * N[[l]][,1,m]
          }
        }

        #By-area and fleet catches and discards
        for(f in 1:nfleets){
          # fleet totals: sum across GTGs and ages for each fleet
          catchN_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(catchNage_by_fleet[[f]][[x]][,1,m])))
          catchB_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*catchNage_by_fleet[[f]][[x]][,1,m])))
          discN_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal_by_fleet[1,k,m,f]*selHist[[m]][[f]]$discard[[x]]/(Z[[x]][,1,m])*(1-exp(-Z[[x]][,1,m]))*N[[x]][,1,m])))
          discB_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal_by_fleet[1,k,m,f]*selHist[[m]][[f]]$discard[[x]]/(Z[[x]][,1,m])*(1-exp(-Z[[x]][,1,m]))*N[[x]][,1,m])))
        }

        #By-area catches and discards
        catchN[1,k,m] <- sum(catchN_by_fleet[1,k,m,1:nfleets])
        catchB[1,k,m] <- sum(catchB_by_fleet[1,k,m,1:nfleets])
        discN[1,k,m] <- sum(discN_by_fleet[1,k,m,1:nfleets])
        discB[1,k,m] <- sum(discB_by_fleet[1,k,m,1:nfleets])

        #Total catchNage across fleets: Sum fleet-specific catches to create total catch-at-age
        #create total catchNage across fleets for backward compatibility
        for(l in 1:lh$gtg){
          catchNage[[l]][,1,m] <- rowSums(sapply(1:nfleets, function(f) catchNage_by_fleet[[f]][[l]][,1,m]))
        }
      }
    }

    ###########
    #Biomass & recruitment
    for(m in 1:areas) {
      SB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]][,1,m]*lh$mat[[x]]*lh$W[[x]])[2:ageClasses])))
    }
    SPR[1,k]<-(sum(SB[1,k,])/is$Req)/(is$B0/lh$LifeHistory@R0)
    relSB[1,k]<-sum(SB[1,k,])/is$B0
    recN[1,k]<-is$Req

    #Single fleet
    if(!is_multifleet){
      for(m in 1:areas){
        #By-area biomass
        VB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist[[m]]$vul[[x]]*lh$W[[x]])))
        RB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist[[m]]$keep[[x]]*lh$W[[x]])))
      }
    }

    #Multi fleet
    if(is_multifleet){
      #Bill edit: VB and RB should be by-fleet (not for selected fleet)
      #added here:
      for(m in 1:areas){
        for(f in 1:nfleets) {
          VB_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist[[m]][[f]]$vul[[x]]*lh$W[[x]])))
          RB_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist[[m]][[f]]$keep[[x]]*lh$W[[x]])))
        }
      }
    }

    #--------------------
    #8. Time dynamics
    #--------------------
    for (j in 2:years){

      #---------------------------------------------------------------------------
      #9. Selgroup correct sel group based on historical or projection time period
      #---------------------------------------------------------------------------
      #Bill Edit
      if(is_multifleet){
        selHist<-lapply(1:TimeAreaObj@areas, function(area){
          lapply(1:nfleets, function(f) {
            indX<-sapply(1:NROW(MultifleetObj@fleet_block_hist_list[[f]]), function(x){(j-1) %in% MultifleetObj@fleet_block_hist_list[[f]][[x]]})
            selHistBlock[[area]][[f]][[which(indX)[1]]]
          })
        })
      }

      if(controlRuleYear[j]) selGroup <- selPro
      if(!controlRuleYear[j]) selGroup <- selHist

      #-------------------------------------------------------------
      #10. Annual regulation decisions - known as phase 2 in any MP
      #-------------------------------------------------------------

      #Single fleet
      if(!is_multifleet){
        dataObject<-c(list(j=j,
                           k=k,
                           is=is,
                           lh = lh,
                           areas = areas,
                           ageClasses = ageClasses,
                           is_multifleet = is_multifleet,
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
      }

      #Multi fleet
      if(is_multifleet){
        dataObject<-c(list(j=j,
                           k=k,
                           is=is,
                           lh = lh,
                           areas = areas,
                           ageClasses = ageClasses,
                           is_multifleet = is_multifleet,
                           nfleets = nfleets,
                           N=N,
                           Z=Z,
                           catchNage=catchNage,
                           selGroup = selGroup,
                           selHist = selHist,
                           selPro = selPro,
                           SB=SB,
                           VB_by_fleet=VB_by_fleet,
                           RB_by_fleet=RB_by_fleet,
                           catchN=catchN,
                           catchB=catchB,
                           discN=discN,
                           discB=discB,
                           discN_by_fleet = discN_by_fleet,
                           catchNage_by_fleet = catchNage_by_fleet,
                           Ftotal_by_fleet=Ftotal_by_fleet,
                           catchB_by_fleet = catchB_by_fleet,
                           catchN_by_fleet = catchN_by_fleet,
                           discB_by_fleet = discB_by_fleet,
                           final_effort_proportions = final_effort_proportions,
                           actual_catch_proportions = actual_catch_proportions, #Bill edit: would we need to include target_catch_proportions?
                           allocation_type = allocation_type,
                           SPR=SPR,
                           relSB=relSB,
                           recN=recN,
                           decisionData=decisionData,
                           decisionAnnual=decisionAnnual,
                           decisionLocal=decisionLocal
        ),
        inputObject
        )
      }

      if(controlRuleYear[j]) decisionAnnual<-rbind(decisionAnnual, do.call(get(StrategyObj@projectionName), list(phase=2, dataObject)))

      #-------------------------------------------------------------
      #11. HCR and F calculation - known as phase 2 in any MP
      #-------------------------------------------------------------

      #Single fleet
      if(!is_multifleet){
        dataObject<-c(list(j=j,
                           k=k,
                           is=is,
                           lh = lh,
                           areas = areas,
                           ageClasses = ageClasses,
                           is_multifleet = is_multifleet,
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
      }

      #Multi fleet
      if(is_multifleet){
        dataObject<-c(list(j=j,
                           k=k,
                           is=is,
                           lh = lh,
                           areas = areas,
                           ageClasses = ageClasses,
                           is_multifleet = is_multifleet,
                           nfleets = nfleets,
                           N=N,
                           Z=Z,
                           catchNage=catchNage,
                           selGroup = selGroup,
                           selHist = selHist,
                           selPro = selPro,
                           SB=SB,
                           VB_by_fleet=VB_by_fleet,
                           RB_by_fleet=RB_by_fleet,
                           catchN=catchN,
                           catchB=catchB,
                           discN=discN,
                           discB=discB,
                           Ftotal_by_fleet=Ftotal_by_fleet,
                           catchB_by_fleet = catchB_by_fleet,
                           catchN_by_fleet = catchN_by_fleet,
                           discB_by_fleet = discB_by_fleet,
                           discN_by_fleet = discN_by_fleet,
                           catchNage_by_fleet = catchNage_by_fleet,
                           final_effort_proportions = final_effort_proportions,
                           actual_catch_proportions = actual_catch_proportions, #Bill edit: would we need to include target_catch_proportions?
                           allocation_type = allocation_type,
                           SPR=SPR,
                           relSB=relSB,
                           recN=recN,
                           decisionData=decisionData,
                           decisionAnnual=decisionAnnual,
                           decisionLocal=decisionLocal
        ),
        inputObject
        )
      }

      if(controlRuleYear[j]) { decisionLocal<-rbind(decisionLocal, do.call(get(StrategyObj@projectionName), list(phase=3, dataObject)))
      } else { decisionLocal<-rbind(decisionLocal, do.call(fixedStrategy, list(phase=3, dataObject)))}

      #-------------------------------------------------------------
      #12. Update arrays for current time step
      #-------------------------------------------------------------

      #Bill edit. I've completed refactoring of fixedStrategy() and thus modified the code below.

      #SB and recruits
      for(m in 1:areas) SB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]][,j,m]*lh$mat[[x]]*lh$W[[x]])[2:ageClasses])))
      Rtmp<-recruit(LifeHistoryObj = lh$LifeHistory, B0=is$B0, stock=sum(SB[j,k,]))
      SPR[j,k]<-(sum(SB[j,k,])/Rtmp)/(is$B0/lh$LifeHistory@R0)
      relSB[j,k]<-sum(SB[j,k,])/is$B0
      recN[j,k]<-Rtmp*RdevMatrix[j,k]
      for (l in 1:lh$gtg) N[[l]][1,j,]<- Rtmp*lh$recProb[l]*TimeAreaObj@recArea*RdevMatrix[j,k]

      #Single fleet
      if(!is_multifleet){
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


      #Multi fleet
      if(is_multifleet){

        for(m in 1:areas){
          for(f in 1:nfleets) {
            VB_by_fleet[j,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup[[m]][[f]]$vul[[x]]*lh$W[[x]])))
            RB_by_fleet[j,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup[[m]][[f]]$keep[[x]]*lh$W[[x]])))
          }
        }

        for(m in 1:areas){

          #By-area and each fleet fishing mortality
          for(f in 1:nfleets) {
            #Fleet-specific F from 4D array [year, iteration, area, fleet]
            xRow<-which(decisionLocal$year==j & decisionLocal$iteration==k & decisionLocal$area==m & decisionLocal$fleet==f)
            Ftotal_by_fleet[j,k,m,f] <- decisionLocal$Flocal[xRow]
          }

          #Calculate catch matrices by gtg and age
          for(l in 1:lh$gtg){
            #calculate total Z from all fleet contributions
            #Z = M + sum_across_fleets(F_fleet * selectivity_fleet)
            # NEVER combine selectivities - each fleet contributes independently
            total_fishing_mortality <- sapply(1:ageClasses, function(age) {
              sum(sapply(1:nfleets, function(ff) {
                Ftotal_by_fleet[j,k,m,ff] * selGroup[[m]][[ff]]$removal[[l]][age]
              }))
            })

            #total mortality shared by all fleets: Z = M + total_fishing_mortality
            Z[[l]][,j,m] <- total_fishing_mortality + lh$LifeHistory@M

            #fleet-specific catch (using using Baranov equation with shared Z)
            for(f in 1:nfleets) {
              catchNage_by_fleet[[f]][[l]][,j,m] <- Ftotal_by_fleet[j,k,m,f] * selGroup[[m]][[f]]$keep[[l]] /
                Z[[l]][,j,m] * (1-exp(-Z[[l]][,j,m])) * N[[l]][,j,m]
            }
          }

          # fleet totals: sum across GTGs and ages for each fleet
          for(f in 1:nfleets) {
            catchN_by_fleet[j,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(catchNage_by_fleet[[f]][[x]][,j,m])))
            catchB_by_fleet[j,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*catchNage_by_fleet[[f]][[x]][,j,m])))
            discN_by_fleet[j,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal_by_fleet[j,k,m,f]*selGroup[[m]][[f]]$discard[[x]]/(Z[[x]][,j,m])*(1-exp(-Z[[x]][,j,m]))*N[[x]][,j,m])))
            discB_by_fleet[1,k,m,f] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal_by_fleet[j,k,m,f]*selGroup[[m]][[f]]$discard[[x]]/(Z[[x]][,j,m])*(1-exp(-Z[[x]][,j,m]))*N[[x]][,j,m])))
          }

          catchN[j,k,m] <- sum(catchN_by_fleet[j,k,m,1:nfleets])
          catchB[j,k,m] <- sum(catchB_by_fleet[j,k,m,1:nfleets])
          discN[j,k,m] <- sum(discN_by_fleet[j,k,m,1:nfleets])
          discB[j,k,m] <- sum(discB_by_fleet[j,k,m,1:nfleets])

          #create total catchNage across fleets for backward compatibility
          for(l in 1:lh$gtg){
            catchNage[[l]][,j,m] <- rowSums(sapply(1:nfleets, function(f) catchNage_by_fleet[[f]][[l]][,j,m]))
          }
        }
      }

      #Next year abundance, move through each gtg (no changes needed)
      for (l in 1:lh$gtg){
        if(j<years){
          P<-matrix(nrow=ageClasses*areas, ncol=ageClasses*areas)
          for(m in 1:areas){
            #Single fleet
            if(!is_multifleet) {
              S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=Ftotal[j,k,m], S_in=selGroup[[m]]$removal[[l]])
            }
            #Multi fleet
            if(is_multifleet) {
              #calculate total Z from all fleet contributions
              #Z = M + sum_across_fleets(F_fleet * selectivity_fleet)
              # NEVER combine selectivities - each fleet contributes independently
              total_fishing_mortality <- sapply(1:ageClasses, function(age) {
                sum(sapply(1:nfleets, function(ff) {
                  Ftotal_by_fleet[j,k,m,ff] * selGroup[[m]][[ff]]$removal[[l]][age]
                }))
              })
              #Note that fishing mortality accounted for in total_fishing_mortality vector, thus F_in = 1 is just a dummy input
              S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=1, S_in=total_fishing_mortality )
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
        #Single fleet
        if(!is_multifleet){
          dataObject<-c(list(j=j,
                             k=k,
                             is=is,
                             lh = lh,
                             areas = areas,
                             ageClasses = ageClasses,
                             is_multifleet = is_multifleet,
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
        }
        #Multi fleet
        if(is_multifleet){
          dataObject<-c(list(j=j,
                             k=k,
                             is=is,
                             lh = lh,
                             areas = areas,
                             ageClasses = ageClasses,
                             is_multifleet = is_multifleet,
                             nfleets = nfleets,
                             N=N,
                             Z=Z,
                             catchNage=catchNage,
                             selGroup = selGroup,
                             selHist = selHist,
                             selPro = selPro,
                             SB=SB,
                             VB_by_fleet=VB_by_fleet,
                             RB_by_fleet=RB_by_fleet,
                             catchN=catchN,
                             catchB=catchB,
                             discN=discN,
                             discB=discB,
                             Ftotal_by_fleet=Ftotal_by_fleet,
                             catchB_by_fleet = catchB_by_fleet,
                             catchN_by_fleet = catchN_by_fleet,
                             discB_by_fleet = discB_by_fleet,
                             discN_by_fleet = discN_by_fleet,
                             catchNage_by_fleet = catchNage_by_fleet,
                             final_effort_proportions = final_effort_proportions,
                             actual_catch_proportions = actual_catch_proportions, #Bill edit: would we need to include target_catch_proportions?
                             allocation_type = allocation_type,
                             SPR=SPR,
                             relSB=relSB,
                             recN=recN,
                             decisionData=decisionData,
                             decisionAnnual=decisionAnnual,
                             decisionLocal=decisionLocal
          ),
          inputObject
          )
        }
        decisionData<-rbind(decisionData, do.call(get(StrategyObj@projectionName), list(phase=1, dataObject)))
      }
    } #end j loop

    #Optional exports for diagnostic mode.
    if(doDiagnostic & k==1) {
      Nexport = N
      catchNageExport = catchNage
      Zexport = Z
    }

    if(!is.null(hostName) & !is.null(waitName)){
      hostName$set(k/floor(TimeAreaObj@iterations)*100)
    }
  } #End k loop

  if(!is.null(hostName) & !is.null(waitName)){
    waitName$hide()
  }

  #Save multi fleet data
  if(is_multifleet) {
    dynamics_multifleet <- list(
      Ftotal_by_fleet = Ftotal_by_fleet,
      catchB_by_fleet = catchB_by_fleet,
      catchN_by_fleet = catchN_by_fleet,
      discB_by_fleet = discB_by_fleet,
      discN_by_fleet = discN_by_fleet,
      RB_by_fleet = RB_by_fleet,
      VB_by_fleet = RB_by_fleet,
      final_effort_proportions = final_effort_proportions,
      target_catch_proportions = target_catch_proportions,
      actual_catch_proportions = actual_catch_proportions,
      allocation_type = allocation_type
    )

    #save
    dynamics<-list(SB=SB, catchB=catchB, catchN=catchN, discB=discB, discN=discN, SPR=SPR, relSB=relSB, recN=recN, ref = ref,
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


runProjection<-function(LifeHistoryObj, TimeAreaObj, HistFisheryObj, ProFisheryObj_list = NULL, StrategyObj = NULL, StochasticObj = NULL, MultifleetObj = NULL, IndexObj=NULL, CatchObsObj=NULL, LengthCompObj=NULL,
                        wd, fileName, seed = 1, doPlot = FALSE, doDiagnostic=F, customToCluster = NULL, titleStrategy = "No name", waitName=NULL, hostName=NULL){

  #-----------------------
  #Build inputObject
  #-----------------------
  TimeAreaObj@recArea <- TimeAreaObj@recArea / sum(TimeAreaObj@recArea) #Make sure this sums to 1

  #Multifleet detection and validation
  is_multifleet <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1

  if(is_multifleet) {
    #nfleets parameter
    nfleets <- MultifleetObj@nfleets
    #override these to avoid using them
    HistFisheryObj <- NULL
    ProFisheryObj_list <- NULL
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
  if(is_multifleet) { #Empty (no) selectivity stochasticity for multifleet
    Sdev <- list(hist = list(), pro = replicate(TimeAreaObj@areas, list(), simplify = FALSE))
  } else {
    #for single fleet: keep all selectivity stochasticity
    Sdev<-selDev(TimeAreaObj, HistFisheryObj, ProFisheryObj_list, StochasticObj)
  }

  #Historical effort devs
  if(is_multifleet) {
    #Using 4D histEffortDev array
    histEffortDev<-histEffortDev(TimeAreaObj, StochasticObj, is_multifleet = is_multifleet, nfleets = nfleets)$Emult
  } else {
    #Using 3D histEffortDev array
    histEffortDev<-histEffortDev(TimeAreaObj, StochasticObj)$Emult
  }

  #---------------------------------------
  #Initial checks that do not stop program
  #---------------------------------------

  print("
  #---------------
  #Initial checks
  #---------------
  ")

  #Single or multifleet mode
  if(is_multifleet) {
    print(paste("Multifleet mode enabled with", nfleets, "fleets."))
    print("Note: Selectivity stochasticity disabled for multifleet mode.")
    print("HistFisheryObj and ProFisheryObj_list are not needed and will be ignored.")
    print("Fleet selectivities are defined in MultifleetObj@fleet_selectivity_hist_list and MultifleetObj@fleet_selectivity_proj_list.")
  } else {
    print("Single fleet mode.")
  }

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
  if(!is.null(StochasticObj)){
    if(is_multifleet) {
      print("Uncertainty in historical fishery selectivity parameters: disabled for multifleet mode")
      print("Uncertainty in projection fishery selectivity parameters: disabled for multifleet mode")
    } else {
      #Historical
      selListHist<-names(Sdev$hist[!unlist(lapply(Sdev$hist, is.null))])
      if(NROW(selListHist) > 0) {
        print(paste("Uncertainty in historical fishery selectivity parameters:", selListHist))
      } else {
        print(paste("Uncertainty in historical fishery selectivity parameters:", "none"))
      }
      #Projection
      for(i in 1:TimeAreaObj@areas){
        selListPro<-names(Sdev$pro[[i]][!unlist(lapply(Sdev$pro[[i]], is.null))])
        if(NROW(selListPro) > 0) {
          print(paste("Area", i, "uncertainty in projection fishery selectivity parameters:", selListPro))
        } else {
          print(paste("Area", i, "uncertainty in projection fishery selectivity parameters:", "none"))
        }
      }
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

  #Is iterations specified correctly?
  if(length(TimeAreaObj@iterations) == 0 || TimeAreaObj@iterations < 1) {
      stop("Iterations not specified correctly.")
  }

  #Number of areas not in agreement with dimensions of the move matrix.
  if(isTRUE(TimeAreaObj@areas != dim(TimeAreaObj@move)[1] | TimeAreaObj@areas != dim(TimeAreaObj@move)[2])) {
    stop("Number of areas not in agreement with dimensions of the move matrix. Check inputs.")
  }

  #Number of historical years does not match historical effort time series
  if(is_multifleet) {
    if(isTRUE(TimeAreaObj@historicalYears > 0  && isTRUE(TimeAreaObj@areas != dim(MultifleetObj@fleet_historicalEffort)[2] | TimeAreaObj@historicalYears != dim(MultifleetObj@fleet_historicalEffort)[1]))) {
      stop("Multi fleet: number of historical years does not match historical effort time series. Check inputs.")
    }
  } else {
    if(isTRUE(TimeAreaObj@historicalYears > 0  && isTRUE(TimeAreaObj@areas != dim(TimeAreaObj@historicalEffort)[2] | TimeAreaObj@historicalYears != dim(TimeAreaObj@historicalEffort)[1]))) {
      stop("Number of historical years does not match historical effort time series. Check inputs.")
    }
  }

  #Historical and/or projection years must be at least 1.
  if(isTRUE(TimeAreaObj@historicalYears + ifelse(is(StrategyObj, "Strategy")  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0) < 1)) {
    stop("Historical and/or projection years must be at least 1. Check inputs.")
  }

  #Multifleet issues
  if(is_multifleet) {

    #Check fleet proportions sum to 1
    if(abs(sum(MultifleetObj@fleet_proportions) - 1.0) > 1e-6) {
      print("Fleet proportions:", paste(round(MultifleetObj@fleet_proportions, 3), collapse = ", "), "\n")
      stop(paste("Multi fleet: fleet proportions must sum to 1.0. current sum:", sum(MultifleetObj@fleet_proportions)))
    }

    #Check fleet allocation type
    if(!MultifleetObj@allocation_type %in% c("effort", "catch")) {
      stop("Multi fleet: fleet allocation_type must be 'effort' or 'catch'")
    }

    #Check for correct number of fleets in historical effort
    if(dim(MultifleetObj@fleet_historicalEffort)[3] !=  nfleets) {
      stop("Multi fleet: incorrect number of fleets in historical effort. Check inputs.")
    }

    #Historical sel - Check for correct number of fleets
    if(length(MultifleetObj@fleet_selectivity_hist_list) != nfleets) {
      stop(paste("Multi fleet: fleet_selectivity_hist_list must contain", nfleets, "Fishery objects"))
    }

    #Historical sel - check that list dimensions match blocks and fishery objects
    for(f in 1:nfleets) {
      if(NROW(MultifleetObj@fleet_block_hist_list[[f]]) != NROW(MultifleetObj@fleet_selectivity_hist_list[[f]])){
        stop(paste("Fleet", f, "historical: mismatch between number of blocks and Fishery objects"))
      }
    }

    #Historical sel - Check that each fleet slot contains a Fishery object
    #Bill Edit
    for(f in 1:nfleets) {
      if(!any(sapply(1:NROW(MultifleetObj@fleet_block_hist_list[[f]]), function(x){is(MultifleetObj@fleet_selectivity_hist_list[[f]][[x]], "Fishery")}))) {
        stop(paste("Fleet", f, "historical Fishery object is missing"))
      }
    }

    #Historical sel - check to see if blocks contain each year in historical and only once
    for(f in 1:nfleets) {
      tmp<-MultifleetObj@fleet_block_hist_list[[f]]
      #Contains entries of length equal to TimeAreaObj@historicalYears
      if(isFALSE(NROW(unique(unlist(tmp))) == TimeAreaObj@historicalYears & min(unlist(tmp)) == 1 & max(unlist(tmp)) == TimeAreaObj@historicalYears)){
        stop(paste("Fleet", f, "Problem with historical Fishery blocks"))
      }
    }

    #Projection sel - check for correct number of areas
    if(is(StrategyObj, "Strategy") && length(MultifleetObj@fleet_selectivity_proj_list) != TimeAreaObj@areas) {
      stop(paste("Multi fleet: fleet_selectivity_proj_list[[areas]] must contain", TimeAreaObj@areas, "areas"))
    }

    #Projection sel - check for correct number of Fishery by area
    if(is(StrategyObj, "Strategy") &&
       any(sapply(1:TimeAreaObj@areas, function(x){length(MultifleetObj@fleet_selectivity_proj_list[[x]]) != nfleets}))) {
       stop(paste("Multi fleet: fleet_selectivity_proj_list must contain[[areas]][[nfleets]]", nfleets, "Fishery objects"))
    }

    #Projection sel - Check that each fleet slot contains a Fishery object
    for(f in 1:nfleets) {
      if(is(StrategyObj, "Strategy") &&
         any(sapply(1:TimeAreaObj@areas, function(x){!is(MultifleetObj@fleet_selectivity_proj_list[[x]][[f]], 'Fishery')}))) {
         stop(paste("Multi fleet: fleet_selectivity_proj_list, fleet", f, "Fishery object is missing"))
      }
    }
  }

  #Can life history and selectivity wrappers be created for each iteration?
  #LH uncertainty list
  LHList<-names(LHdev[!unlist(lapply(LHdev, is.null))])

  #Single fleet
  if(!is_multifleet) {
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
          stop(paste("Life history cannot be created. Check inputs. Stopped at interation", k))
        }
        for(x in 1:TimeAreaObj@areas){
          if(is.null(selHist[[x]])){
            stop(paste("Historical selectivity cannot be created. Check inputs. Stopped at interation", k))
          }
        }
        for(x in 1:TimeAreaObj@areas){
          if(isTRUE(!is.null(StrategyObj) &  is.null(selPro[[x]]))){
            stop(paste("Projection selectivity cannot be created. Check inputs. Stopped at interation", k))
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
        stop("Life history cannot be created. Check inputs.")
      }
      for(x in 1:TimeAreaObj@areas){
        if(is.null(selHist[[x]])){
          stop("Historical selectivity cannot be created. Check inputs.")
        }
      }
      for(x in 1:TimeAreaObj@areas){
        if(isTRUE(!is.null(StrategyObj) &  is.null(selPro[[x]]))){
          stop("Projection selectivity cannot be created. Check inputs.")
        }
      }
    }
  }

  #Multi fleet
  if(is_multifleet) {
    if(NROW(LHList) > 0){
      for(k in 1:floor(TimeAreaObj@iterations)){
        #LH
        LifeHistoryObj_TMP<-LifeHistoryObj
        if(NROW(LHList) > 0) {
          for(x in 1:NROW(LHList)) slot(LifeHistoryObj_TMP, LHList[x]) <- LHdev[[LHList[x]]][k]
        }
        #Setup
        lh<-LHwrapper(LifeHistoryObj_TMP, TimeAreaObj)
        selHistBlock <- lapply(1:TimeAreaObj@areas, function(x){
          lapply(1:nfleets, function(f) {
            bl<-NROW(MultifleetObj@fleet_block_hist_list[[f]])
            lapply(1:bl, function(b) {
              selWrapper(lh, TimeAreaObj,
                         FisheryObj = MultifleetObj@fleet_selectivity_hist_list[[f]][[b]])
            })
          })
        })

        selPro <- lapply(1:TimeAreaObj@areas, function(area){
          lapply(1:nfleets, function(f) {
            if(is(StrategyObj, "Strategy")) {
              selWrapper(lh, TimeAreaObj,
                         FisheryObj = MultifleetObj@fleet_selectivity_proj_list[[area]][[f]])
            } else {
              NULL
            }
          })
        })
        if(is.null(lh)) {
          stop("Multi fleet: Life history cannot be created. Check inputs. Stopped at interation", k)
        }
        #Bill Edit
        for(x in 1:TimeAreaObj@areas){
          for(f in 1:nfleets){
            for(b in 1:NROW(MultifleetObj@fleet_block_hist_list[[f]])){
              if(is.null(selHistBlock[[x]][[f]][[b]])){
                stop("Multi fleet: Historical selectivity cannot be created. Check inputs.")
              }
            }
          }
        }
        for(x in 1:TimeAreaObj@areas){
          for(f in 1:nfleets){
            if(is(StrategyObj, "Strategy") &&  is.null(selPro[[x]][[f]])){
              stop(paste("Multi fleet: Projection selectivity cannot be created. Check inputs. Stopped at interation", k))
            }
          }
        }
      }
    } else {
      lh<-LHwrapper(LifeHistoryObj, TimeAreaObj)
      #Bill Edit
      selHistBlock <- lapply(1:TimeAreaObj@areas, function(x){
        lapply(1:nfleets, function(f) {
          bl<-NROW(MultifleetObj@fleet_block_hist_list[[f]])
          lapply(1:bl, function(b) {
            selWrapper(lh, TimeAreaObj,
                       FisheryObj = MultifleetObj@fleet_selectivity_hist_list[[f]][[b]])
          })
        })
      })
      selPro <- lapply(1:TimeAreaObj@areas, function(area){
        lapply(1:nfleets, function(f) {
          #check if projection sel is available (do we have startegy and proj sel exist for area/fleet?)
          if(is(StrategyObj, "Strategy")) {
            selWrapper(lh, TimeAreaObj,
                       FisheryObj = MultifleetObj@fleet_selectivity_proj_list[[area]][[f]])
          } else {
            NULL
          }
        })
      })
      if(is.null(lh)) {
        stop("Multi fleet: Life history cannot be created. Check inputs.")
      }
      #Bill Edit
      for(x in 1:TimeAreaObj@areas){
        for(f in 1:nfleets){
          for(b in 1:NROW(MultifleetObj@fleet_block_hist_list[[f]])){
            if(is.null(selHistBlock[[x]][[f]][[b]])){
              stop("Multi fleet: Historical selectivity cannot be created. Check inputs.")
            }
          }
        }
      }
      for(x in 1:TimeAreaObj@areas){
        for(f in 1:nfleets){
          if(is(StrategyObj, "Strategy") &&  is.null(selPro[[x]][[f]])){
            stop("Projection selectivity cannot be created. Check inputs.")
          }
        }
      }
    }
  }

  #Rec devs failed
  if(is.null(RdevMatrix)) {
    stop("Inter-annual recruitment variation cannot be created. Check inputs.")
  }

  #Ddev
  if(is.null(Ddev)) {
    stop("Initial biomass variation cannot be created. Check inputs.")
  }

  #Cdev
  if(is(StrategyObj,"Strategy")  &&
     StrategyObj@projectionName == "projectionStrategy" &&
     is.null(Cdev)) {
    stop("Initial CPUE variation cannot be created. Check inputs.")
  }

  #Edev
  if(is(StrategyObj,"Strategy")  &&
     StrategyObj@projectionName == "projectionStrategy" &&
     is.null(Edev)) {
    stop("Effort implementation error cannot be created. Check inputs.")
  }

  #Historical effort devs failed to be created
  if(is.null(histEffortDev)) {
    stop("Inter-annual historical effort variation cannot be created. Check inputs.")
  }

  #Historical effort devs have incorrect dimensions
  if(is_multifleet) {
    expected_dims <- c(1 + TimeAreaObj@historicalYears,
                       as.integer(floor(TimeAreaObj@iterations)),
                       TimeAreaObj@areas,
                       nfleets)

    expected_dims2 <- c(1 + dim(MultifleetObj@fleet_historicalEffort)[1],
                        as.integer(floor(TimeAreaObj@iterations)),
                        dim(MultifleetObj@fleet_historicalEffort)[2],
                        dim(MultifleetObj@fleet_historicalEffort)[3])
  } else {
    expected_dims <- c(1 + TimeAreaObj@historicalYears,
                       as.integer(floor(TimeAreaObj@iterations)),
                       TimeAreaObj@areas)

    expected_dims2 <- c(1 + dim(TimeAreaObj@historicalEffort)[1],
                        as.integer(floor(TimeAreaObj@iterations)),
                        dim(TimeAreaObj@historicalEffort)[2])
  }

  if(!all(dim(histEffortDev) == expected_dims) || !all(dim(histEffortDev) == expected_dims2)) {
    stop("Fleet histEffortDev dimension mismatch.\n",
         "  expected: ", paste(expected_dims, collapse = " x "), "\n",
         "  got:      ", paste(dim(histEffortDev), collapse = " x "))
  }


  #Management strategy function missing
  if(isTRUE(is(StrategyObj, "Strategy") &&
                             tryCatch({
                               get(StrategyObj@projectionName)
                               FALSE
                             }, error = function(c) TRUE)
  )) {
    stop("Project strategy function missing. StrategyObj@projectionName must correspond to a named function.")
  }

  #When applying a management strategy, StrategyObj@projectionYears must be greater than 0
  if(isTRUE(is(StrategyObj, "Strategy") && length(StrategyObj@projectionYears) == 0) ||
     isTRUE(is(StrategyObj, "Strategy") && StrategyObj@projectionYears < 1)
  ){
    stop("When applying a management strategy, StrategyObj@projectionYears must be greater than 0")
  }


  #---------------------------
  #Setup parallel processing
  #---------------------------

  ptm<-proc.time()
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

    #Single fleet
    if(!is_multifleet){
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

      #Optional diagnostic outputs
      N<-mseParallel[[1]]$N
      Z<-mseParallel[[1]]$Z
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
        }
        SPR[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$SPR[,input[[i]][1]:input[[i]][2]]
        relSB[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$relSB[,input[[i]][1]:input[[i]][2]]
        recN[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$recN[,input[[i]][1]:input[[i]][2]]
        ref[input[[i]][1]:input[[i]][2],]<-mseParallel[[i]]$dynamics$ref[input[[i]][1]:input[[i]][2],]
        decisionAnnual<-rbind(decisionAnnual, mseParallel[[i]]$HCR$decisionAnnual)
        decisionLocal<-rbind(decisionLocal, mseParallel[[i]]$HCR$decisionLocal)
        decisionData<-rbind(decisionData, mseParallel[[i]]$HCR$decisionData)
      }

    }

    #Multi fleet
    if(is_multifleet){
      SB<-mseParallel[[1]]$dynamics$SB
      catchB<-mseParallel[[1]]$dynamics$catchB
      catchN<-mseParallel[[1]]$dynamics$catchN
      discB<-mseParallel[[1]]$dynamics$discB
      discN<-mseParallel[[1]]$dynamics$discN
      SPR<-mseParallel[[1]]$dynamics$SPR
      relSB<-mseParallel[[1]]$dynamics$relSB
      recN<-mseParallel[[1]]$dynamics$recN
      ref<-mseParallel[[1]]$dynamics$ref
      decisionAnnual<-mseParallel[[1]]$HCR$decisionAnnual
      decisionLocal<-mseParallel[[1]]$HCR$decisionLocal
      decisionData<-mseParallel[[1]]$HCR$decisionData

      Ftotal_by_fleet <- mseParallel[[1]]$dynamics$multifleet$Ftotal_by_fleet
      catchB_by_fleet <- mseParallel[[1]]$dynamics$multifleet$catchB_by_fleet
      catchN_by_fleet <- mseParallel[[1]]$dynamics$multifleet$catchN_by_fleet
      discB_by_fleet <- mseParallel[[1]]$dynamics$multifleet$discB_by_fleet
      discN_by_fleet <- mseParallel[[1]]$dynamics$multifleet$discN_by_fleet
      RB_by_fleet <- mseParallel[[1]]$dynamics$multifleet$RB_by_fleet
      VB_by_fleet <- mseParallel[[1]]$dynamics$multifleet$RB_by_fleet

      final_effort_proportions <- mseParallel[[1]]$dynamics$multifleet$final_effort_proportions
      target_catch_proportions <- mseParallel[[1]]$dynamics$multifleet$target_catch_proportions
      actual_catch_proportions <- mseParallel[[1]]$dynamics$multifleet$actual_catch_proportions
      allocation_type <- mseParallel[[1]]$dynamics$multifleet$allocation_type

      for (i in 2:cores){
        for(m in 1:TimeAreaObj@areas){
          SB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$SB[,input[[i]][1]:input[[i]][2],m]
          catchB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$catchB[,input[[i]][1]:input[[i]][2],m]
          catchN[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$catchN[,input[[i]][1]:input[[i]][2],m]
          discB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$discB[,input[[i]][1]:input[[i]][2],m]
          discN[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$discN[,input[[i]][1]:input[[i]][2],m]
          for(f in 1:nfleets) {
            Ftotal_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$Ftotal_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
            catchB_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$catchB_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
            catchN_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$catchN_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
            discB_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$discB_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
            discN_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$discN_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
            RB_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$RB_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
            VB_by_fleet[,input[[i]][1]:input[[i]][2],m,f] <- mseParallel[[i]]$dynamics$multifleet$VB_by_fleet[,input[[i]][1]:input[[i]][2],m,f]
          }
        }
        SPR[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$SPR[,input[[i]][1]:input[[i]][2]]
        relSB[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$relSB[,input[[i]][1]:input[[i]][2]]
        recN[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$recN[,input[[i]][1]:input[[i]][2]]
        ref[input[[i]][1]:input[[i]][2],]<-mseParallel[[i]]$dynamics$ref[input[[i]][1]:input[[i]][2],]
        decisionAnnual<-rbind(decisionAnnual, mseParallel[[i]]$HCR$decisionAnnual)
        decisionLocal<-rbind(decisionLocal, mseParallel[[i]]$HCR$decisionLocal)
        decisionData<-rbind(decisionData, mseParallel[[i]]$HCR$decisionData)

        final_effort_proportions[input[[i]][1]:input[[i]][2],] <- mseParallel[[i]]$dynamics$multifleet$final_effort_proportions[input[[i]][1]:input[[i]][2],]
        target_catch_proportions[input[[i]][1]:input[[i]][2],] <- mseParallel[[i]]$dynamics$multifleet$target_catch_proportions[input[[i]][1]:input[[i]][2],]
        actual_catch_proportions[input[[i]][1]:input[[i]][2],] <- mseParallel[[i]]$dynamics$multifleet$actual_catch_proportions[input[[i]][1]:input[[i]][2],]

        allocation_type[input[[i]][1]:input[[i]][2]] <- mseParallel[[i]]$dynamics$multifleet$allocation_type[input[[i]][1]:input[[i]][2]]
      }
    }

    #Optional diagnostic outputs
    N<-mseParallel[[1]]$N
    Z<-mseParallel[[1]]$Z
    catchNage<-mseParallel[[1]]$catchNage


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

    #Single fleet
    if(!is_multifleet){
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
    }

    #Multi fleet
    if(is_multifleet){
      SB<-mse$dynamics$SB
      catchB<-mse$dynamics$catchB
      catchN<-mse$dynamics$catchN
      discB<-mse$dynamics$discB
      discN<-mse$dynamics$discN
      SPR<-mse$dynamics$SPR
      relSB<-mse$dynamics$relSB
      recN<-mse$dynamics$recN
      ref<-mse$dynamics$ref
      decisionAnnual<-mse$HCR$decisionAnnual
      decisionLocal<-mse$HCR$decisionLocal
      decisionData<-mse$HCR$decisionData

      Ftotal_by_fleet <- mse$dynamics$multifleet$Ftotal_by_fleet
      catchB_by_fleet <- mse$dynamics$multifleet$catchB_by_fleet
      catchN_by_fleet <- mse$dynamics$multifleet$catchN_by_fleet
      discB_by_fleet <- mse$dynamics$multifleet$discB_by_fleet
      discN_by_fleet <- mse$dynamics$multifleet$discN_by_fleet
      RB_by_fleet <- mse$dynamics$multifleet$RB_by_fleet
      VB_by_fleet <- mse$dynamics$multifleet$RB_by_fleet

      final_effort_proportions <- mse$dynamics$multifleet$final_effort_proportions
      target_catch_proportions <- mse$dynamics$multifleet$target_catch_proportions
      actual_catch_proportions <- mse$dynamics$multifleet$actual_catch_proportions
      allocation_type <- mse$dynamics$multifleet$allocation_type
    }

    #Optional diagnostic outputs
    N<-mse$N
    Z<-mse$Z
    catchNage<-mse$catchNage
  }

  #---------------
  #Save results
  #---------------

  #Single fleet
  if(!is_multifleet){
    dynamics<-list(SB=SB, VB=VB, RB=RB, catchB=catchB, catchN=catchN, Ftotal=Ftotal, discB=discB, discN=discN, SPR=SPR, relSB=relSB, recN=recN, ref=ref, N=N, Z=Z, catchNage=catchNage)
  }

  #Single fleet
  if(is_multifleet){
    dynamics<-list(SB=SB, catchB=catchB, catchN=catchN, discB=discB, discN=discN, SPR=SPR, relSB=relSB, recN=recN, ref=ref, N=N, Z=Z, catchNage=catchNage)
    dynamics$multifleet <- list(
      Ftotal_by_fleet = Ftotal_by_fleet,
      catchB_by_fleet = catchB_by_fleet,
      catchN_by_fleet = catchN_by_fleet,
      discB_by_fleet = discB_by_fleet,
      discN_by_fleet = discN_by_fleet,
      RB_by_fleet = RB_by_fleet,
      VB_by_fleet = VB_by_fleet,
      final_effort_proportions = final_effort_proportions,
      target_catch_proportions = target_catch_proportions,
      actual_catch_proportions = actual_catch_proportions,
      allocation_type = allocation_type
    )
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


    #F
    #Single fleet
    if(!is_multifleet) {
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
    }

    #Multi fleet
    if(is_multifleet) {
      # Fleet-specific F plots
      png(filename=paste0(wd, "/", fileName, "_F_by_Fleet.png"), width=12, height=8, units="in", res=96, bg="white", pointsize=12)
      par(mfrow=c(ceiling(nfleets/2)*dt$TimeAreaObj@areas, 2), mar=c(4,4,3,1))
      for(f in 1:nfleets) {
        for(m in 1:dt$TimeAreaObj@areas) {
          if(f == 1 && m == 1) {
            plot(dt$dynamics$multifleet$Ftotal_by_fleet[,1,m,f], type="l", las=1, ylab="F", xlab = "Year",
                 col=rb[1], main = paste("Fleet", f, "Area", m))
          } else {
            plot(dt$dynamics$multifleet$Ftotal_by_fleet[,1,m,f], type="l", las=1, ylab="F", xlab = "Year",
                 col=rb[1], main = paste("Fleet", f, "Area", m))
          }
          if(iterations > 1){
            for(k in 2:iterations){
              lines(dt$dynamics$multifleet$Ftotal_by_fleet[,k,m,f], col=rb[k])
            }
          }
        }
      }
      dev.off()

      # Fleet-specific catch plots
      png(filename=paste0(wd, "/", fileName, "_Catch_by_Fleet.png"), width=12, height=8, units="in", res=96, bg="white", pointsize=12)
      par(mfrow=c(ceiling(nfleets/2)*dt$TimeAreaObj@areas, 2), mar=c(4,4,3,1))
      for(f in 1:nfleets) {
        for(m in 1:dt$TimeAreaObj@areas) {
          plot(dt$dynamics$multifleet$catchB_by_fleet[,1,m,f], type="l", las=1, ylab="Catch", xlab = "Year",
               col=rb[1], main = paste("Fleet", f, "Area", m))
          if(iterations > 1){
            for(k in 2:iterations){
              lines(dt$dynamics$multifleet$catchB_by_fleet[,k,m,f], col=rb[k])
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
