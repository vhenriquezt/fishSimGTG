


#---------------
#Example
#---------------
devtools::load_all()
#library(fishSimGTG)
library(here)

#----------------------------
#Create a LifeHistory object
#----------------------------

#---Populate LifeHistory object
#---Contains the life history parameters
LifeHistoryObj <- new("LifeHistory")
LifeHistoryObj@title<-"Hawaiian Uhu - Parrotfish"
LifeHistoryObj@speciesName<-"Chlorurus perspicillatus"
LifeHistoryObj@Linf<-53.2
LifeHistoryObj@K<-0.225
LifeHistoryObj@t0<- -1.48
LifeHistoryObj@L50<-35
LifeHistoryObj@L95delta<-5.25
LifeHistoryObj@M<-0.16
LifeHistoryObj@L_type<-"FL"
LifeHistoryObj@L_units<-"cm"
LifeHistoryObj@LW_A<-0.0136
LifeHistoryObj@LW_B<-3.109
LifeHistoryObj@Steep<-0.6
LifeHistoryObj@isHermaph<-TRUE
LifeHistoryObj@H50<-46.2
LifeHistoryObj@H95delta<-11.8
LifeHistoryObj@recSD<-0 #Run with no rec var'n to see deterministic trends
LifeHistoryObj@R0 <- 1000
LifeHistoryObj@recRho <- 0

#---Populate a TimeArea object
#---Contains basic inputs about time and space needed to establish simulation bounds
#---The effort matrix is set as multipliers of initial equilibrium effort
TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Example"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 3
TimeAreaObj@historicalYears = 10
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(1:1, nrow = 10, ncol = 2, byrow = FALSE)


#-----------------------------
#Setup fishery characteristics
#-----------------------------

#---Pupulate a Fishery object
#---Contains selectivity, retention and discard characteristics
#---Not sure how to set this up? Type ?selWrapper
HistFisheryObj<-new("Fishery")
HistFisheryObj@title<-"Example"
HistFisheryObj@vulType<-"explogFlex"
HistFisheryObj@vulParams<-c(0.05, 24.20,  0.80, 24.00)
HistFisheryObj@retType<-"full"
HistFisheryObj@retMax <- 1
HistFisheryObj@Dmort <- 0.2

#To simply display to the console
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj)
selWrapper(lh = lhOut, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = TRUE)



#---Pupulate a Fishery object
#---Contains selectivity, retention and discard characteristics
#---Not sure how to set this up? Type ?selWrapper
HistFisheryObj<-new("Fishery")
HistFisheryObj@title<-"Example"
HistFisheryObj@vulType<-"gillnetMasterNormal"
HistFisheryObj@vulParams<-c(4.19,0.71,5.08)
HistFisheryObj@retType<-"logistic"
HistFisheryObj@retParams <- c(25,0.1)
HistFisheryObj@retMax <- 1
HistFisheryObj@Dmort <- 0.2

#To simply display to the console
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj)
selWrapper(lh = lhOut, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = TRUE)


HistFisheryObj<-new("Fishery")
HistFisheryObj@title<-"Example"
HistFisheryObj@vulType<-"logistic"
HistFisheryObj@vulParams<-c(20,1)
HistFisheryObj@retType<-"full"
#HistFisheryObj@retParams <- c(25,0.1)
HistFisheryObj@retMax <- 0.8
HistFisheryObj@Dmort <- 0

#To simply display to the console
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj)
selWrapper(lh = lhOut, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = TRUE)


ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Example"
ProFisheryObj@vulType<-"explog"
ProFisheryObj@vulParams<-c(0.15,40,0.5)
ProFisheryObj@retType<-"logistic"
ProFisheryObj@retParams <- c(30.1, 0.1)
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

ProFisheryObj2<-new("Fishery")
ProFisheryObj2@title<-"Example"
ProFisheryObj2@vulType<-"explog"
ProFisheryObj2@vulParams<-c(0.15,40,0.5)
ProFisheryObj2@retType<-"logistic"
ProFisheryObj2@retParams <- c(30.1, 0.1)
ProFisheryObj2@retMax <- 1
ProFisheryObj2@Dmort <- 0

ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj)


#-----------------------------
#Setup stochastic object
#-----------------------------
StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.3, 0.6)

proFisheryVul_list<-list()
proFisheryVul_list[[1]]<-matrix(c(0.1,0.2,35,45,0.4,0.5), nrow=2, byrow=FALSE)
proFisheryVul_list[[2]]<-matrix(c(0.1,0.2,35,45,0.4,0.5), nrow=2, byrow=FALSE)
StochasticObj@proFisheryVul_list<-proFisheryVul_list

proFisheryRet_list<-list()
proFisheryRet_list[[1]]<-matrix(c(30.1,30.1,0.1,0.1), nrow=2, byrow=FALSE)
proFisheryRet_list[[2]]<-matrix(c(30.1,30.1,0.1,0.1), nrow=2, byrow=FALSE)
StochasticObj@proFisheryRet_list<-proFisheryRet_list

# proFisheryDmort_list<-list()
# proFisheryDmort_list[[1]]<-matrix(c(0.2,0.2), nrow=2, byrow=FALSE)
# proFisheryDmort_list[[2]]<-matrix(c(0.8,0.9), nrow=2, byrow=FALSE)
# StochasticObj@proFisheryDmort_list<-proFisheryDmort_list


#------------------------
#Setup a Strategy object
#------------------------

StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE), CPUE = c(1,2), CPUEtype = "retN")



#----------------
#Run projection
#----------------

#---At this stage, it should be clear that we've set a size limit (even if a useless size limit) via ProFisheryObj and left all other aspects of the fishery constant
#---We will now run the projection

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              ProFisheryObj_list = ProFisheryObj_list,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here(),
              fileName = "Test1",
              doPlot = TRUE,
              titleStrategy = "Test1"
)


#---No report-ready plots are produced. However, see Batch_Projection_example for built report-ready plotting functionality
#---To create custom plots, start by reading in all the data generated by the simulation.

X<-readProjection(wd = here(),
                  fileName = "Test1"
)

