
devtools::load_all()

#----------------------------
#Create a LifeHistory object
#----------------------------

#---Populate LifeHistory object
#---Contains the life history parameters
# Populate the life history object
LifeHistoryObj <- new("LifeHistory")
LifeHistoryObj@title<-"Hawaiian Uhu - Parrotfish"
LifeHistoryObj@speciesName<-"Chlorurus perspicillatus"
LifeHistoryObj@Linf<-53.2
LifeHistoryObj@K<-0.225
LifeHistoryObj@t0<- -1.48
LifeHistoryObj@L50<-35
LifeHistoryObj@L95delta<-5
LifeHistoryObj@M<-0.16
LifeHistoryObj@L_type<-"FL"
LifeHistoryObj@L_units<-"cm"
LifeHistoryObj@LW_A<-0.0136
LifeHistoryObj@LW_B<-3.109
LifeHistoryObj@Steep<-0.6
LifeHistoryObj@isHermaph<-TRUE
LifeHistoryObj@H50<-46.2
LifeHistoryObj@H95delta<-11.8
LifeHistoryObj@recSD<-0
LifeHistoryObj@recRho<-0
LifeHistoryObj@R0<-10000

#---Populate a TimeArea object
#---Contains basic inputs about time and space needed to establish simulation bounds
#---The effort matrix is set as multipliers of initial equilibrium fishing mortality
TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Example"
TimeAreaObj@gtg = 13
TimeAreaObj@gtgCV = 0.3
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 3
TimeAreaObj@historicalYears = 50
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(1:1, nrow = 50, ncol = 2, byrow = FALSE)

# populate the historical fishery object
HistFisheryObj<-new("Fishery")
HistFisheryObj@title<-"Test"
HistFisheryObj@vulType<-"logistic"
HistFisheryObj@vulParams<-c(40, 1)
HistFisheryObj@retType<-"full"
HistFisheryObj@retMax <- 1
HistFisheryObj@Dmort <- 0

StrategyObj <- new("Strategy")
StrategyObj@title <- "Testing"           # add a descriptive name
StrategyObj@projectionYears <- 10                  #future projection
StrategyObj@projectionName <- "constantCatchStrategy" # Which function to use - next: the parameters of projectionStrategy
StrategyObj@projectionParams <- list(aveYrs = 5)

# Create a new Catchobs class
catch_obs1 <- new("CatchObs",
                  catchID = "Fishery Catch",
                  title = "Fishery catch observations",
                  areas = c(1, 2),
                  catchYears = 1:60,
                  reporting_rates = rep(1, 60),
                  obs_CVs = matrix(0:0, nrow = 60, ncol=2)
)

# Now run the simulation
runProjection(
  LifeHistoryObj = LifeHistoryObj,
  TimeAreaObj = TimeAreaObj,
  HistFisheryObj = HistFisheryObj,
  ProFisheryObj_list = list(HistFisheryObj, HistFisheryObj),
  StrategyObj = StrategyObj,
  CatchObsObj = catch_obs1,
  wd = here(),
  fileName = "Test_constantCatch",
  doPlot = TRUE
)

out1<-readProjection(
  wd = here(),
  fileName = "Test_constantCatch"
)
out1$HCR$decisionData
out1$HCR$decisionAnnual

out1$HCR$decisionData$CPUE_1
out1$HCR$decisionData$observed_catch
out1$dynamics$Ftotal
