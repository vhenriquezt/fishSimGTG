# initial testings
 rm(list=ls())
 #devtools::install()
 devtools::load_all()
# #devtools::document()
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

lh <- LHwrapper(LifeHistoryObj = lh, TimeAreaObj = ta)

# Create multifleet object
fishery1 <- new("Fishery")
fishery1@vulType <- "logistic"
fishery1@vulParams <- c(9.0, 1.5)
fishery1@retType <- "full"
fishery1@retMax <- 1
fishery1@Dmort <- 0

fishery2 <- new("Fishery")
fishery2@vulType <- "logistic"
fishery2@vulParams <- c(11.0, 2.0)
fishery2@retType <- "full"
fishery2@retMax <- 1
fishery2@Dmort <- 0

MultifleetObj <- new("Multifleet")
MultifleetObj@nfleets <- 2
MultifleetObj@fleet_proportions <- c(0.6, 0.4)
MultifleetObj@allocation_type <- "catch"
MultifleetObj@fleet_selectivity_list <- list(fishery1, fishery2)

fleet_proportions <- MultifleetObj@fleet_proportions

#selectivity list
hist_sel_list <- lapply(1:MultifleetObj@nfleets, function(f) {
  selWrapper(lh, ta, FisheryObj = MultifleetObj@fleet_selectivity_list[[f]], doPlot = FALSE)
})

#fleet_proportions <- MultifleetObj@fleet_proportions
Ddev <- c(0.4)  # Target depletion
k <- 1

#test
is <- solveD_multifleet2(lh = lh,
sel_list = hist_sel_list,
doFit = TRUE,
D_type = ta@historicalBioType,
D_in = Ddev[k],
#D_type = "relB",
#D_in = 0.4,
fleet_proportions = fleet_proportions,
allocation_type = MultifleetObj@allocation_type)

# Check results
print(is$Feq)  # Total F
print(is$F_by_fleet)  # Fleet-specific F
print(is$D)  # Achieved depletion
print(is$catchB_by_fleet)  # Fleet-specific catches



# COMPARING 1 FLEET APPRAOCH (solveD_multifleet2)
MultifleetObj <- new("Multifleet")
MultifleetObj@nfleets <- 1
MultifleetObj@fleet_proportions <- c(1)
MultifleetObj@allocation_type <- "catch"
MultifleetObj@fleet_selectivity_list <- list(fishery1)

fleet_proportions <- MultifleetObj@fleet_proportions

#selectivity list
hist_sel_list <- lapply(1:MultifleetObj@nfleets, function(f) {
  selWrapper(lh, ta, FisheryObj = MultifleetObj@fleet_selectivity_list[[f]], doPlot = FALSE)
})

#fleet_proportions <- MultifleetObj@fleet_proportions
Ddev <- c(0.4)  # Target depletion
k <- 1

#test
is <- solveD_multifleet2(lh = lh,
                         sel_list = hist_sel_list,
                         doFit = TRUE,
                         D_type = ta@historicalBioType,
                         D_in = Ddev[k],
                         #D_type = "relB",
                         #D_in = 0.4,
                         fleet_proportions = fleet_proportions,
                         allocation_type = MultifleetObj@allocation_type)

# Check results
print(is$Feq)  # Total F
print(is$F_by_fleet)  # Fleet-specific F
print(is$D)  # Achieved depletion
print(is$catchB_by_fleet)  # Fleet-specific catches


# COMPARING 1 FLEET APPRAOCH (solveD)
sel1 <- selWrapper(lh, ta, fishery1, doPlot = FALSE)
is2 <- solveD(lh =lh,
              sel1,
              doFit = TRUE,
              D_type = "relB",
              D_in = 0.4)
print(is2$Feq)
