

#---------------
#Example
#---------------
devtools::load_all()
library(fishSimGTG)
#library(here)

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
TimeAreaObj@iterations = 20
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

#---Visualize life history. Does everything make sense?
#---Optional, create a plot of life history that is useful for reports.

#To save to file (for reports?)
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj, wd = here(), imageName = "LifeHistory", dpi = 300, doPlot = TRUE)

#Note that LHwrapper returns all the details of the life history
lhOut


#---Visualize fishery vulnerability. Does everything make sense?
#---Optional, create a plot of life history that is useful for reports.
#To simply display to the console
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj)
selWrapper(lh = lhOut, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = TRUE)

#To save to file (for reports?)
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj)
selWrapper(lh = lhOut, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = TRUE, wd = here(), imageName = "Vulnerability", dpi = 300)


#--------------------------------------------------
#At this stage, we can simulate historical dynamics
#--------------------------------------------------

#---Notice that the item TimeAreaObj@historicalEffort was set for constant effort for 10 years
#---The output is constant SSB and catch because we specified the same effort conditions to occur into the future
#---The runProjection function is the master function for conducting simulation modeling
#---The runProjection funtion produces many low quality plots. These are meant for the analyst to verify the results.
runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              wd = here::here(),
              fileName = "HistoricalDy",
              doPlot = TRUE,
              titleStrategy = "HistoricalDy"
)

#---What if we change the historical effort trend?
#---An increasing trend that doubles effort after 10 years
TimeAreaObj@historicalEffort<-matrix(seq(1,3,length.out = 50), nrow = 50, ncol = 2, byrow = FALSE)

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              wd = here(),
              fileName = "Test2",
              doPlot = TRUE,
              titleStrategy = "Test2"
)

#---Now, lets add uncertainty in the initial depletion
#---Current, we have the following parameters specificed TimeAreaObj@historicalBio & TimeAreaObj@historicalBioType, which specify a constant SSB/SSB0 for all iterations
#---But when we use a Stochastic object, it overrides this constant value
StochasticObj<-new("Stochastic")
#StochasticObj@historicalBio = c(0.3, 0.6)
StochasticObj@histEffortSD = c(0.1, 0.3)

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              StochasticObj = StochasticObj,
              wd = here(),
              fileName = "Test3",
              doPlot = TRUE,
              titleStrategy = "Test3"
)

#---Finally, lets introduce inter-annual recruitment variation
#---Recruitment variation is specified in the LifeHistory object.
LifeHistoryObj@recSD<-0.6
runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              StochasticObj = StochasticObj,
              wd = here(),
              fileName = "Test4",
              doPlot = TRUE,
              titleStrategy = "Test4"
)
