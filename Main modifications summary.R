General summary of changes to covert single fleet to multifleet approach

1. S4 class structure

1) Adding S4 multifleet class 

setClass("Multifleet",
         representation(
           title = "character",
           nfleets = "numeric",
           fleet_proportions = "numeric",        # Proportional allocation
           allocation_type = "character",        # "effort" or "catch"
           allocation_type = "character",        
           fleet_selectivity_hist_list = "list", # Historical selectivity per fleet
           fleet_selectivity_proj_list = "list"  # Projection selectivity per fleet
         )
)

# uses lists of Fishery objects
# historical vs projection selectivity
# flexible allocation - both effort-based and catch-based
# Backward compatibility - single fleet is just nfleets = 1


2. General math framework

1) single fleet:
# Simple mortality calculation
Z = M + F × selectivity
2) multifleet:
# Shared mortality
# each fleet contributes independently to total mortality
# fleets interact through shared mortality rather than averaging selectivities.
Z = M + sum(F_fleet × selectivity_fleet)

3. Data structure

1) single fleet:
3D arrays: [years, iterations, areas]

# 3D arrays only
SB <- array(dim=c(years, iterations, areas))
VB <- array(dim=c(years, iterations, areas))
catchB <- array(dim=c(years, iterations, areas))
histEffortDev <- array(dim=c(years, iterations, areas))  # 3D


2) Multifleet:
4D arrays: [years, iterations, areas, fleets] #for fleet-specific tracking
arrays[years, iterations, areas]              #maintains 3D arrays for backward compatibility

# 3D arrays maintained for backward compatibility
SB <- array(dim=c(years, iterations, areas))
VB <- array(dim=c(years, iterations, areas))

# NEW: 4D arrays for fleet-specific tracking
Ftotal_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
catchB_by_fleet <- array(dim=c(years, iterations, areas, nfleets))
histEffortDev <- array(dim=c(years, iterations, areas, nfleets))  # 4D



4. Selectivity structure

1) Single fleet:
Simple selectivity: selHist[[area]]

# Simple nested structure
selHist[[area]]$vul[[gtg]]
selPro[[area]]$keep[[gtg]]


2) Multifleet: 

Selectivity
Nested selectivity: selHist[[area]][[fleet]]

# Double-nested structure for fleet-specific selectivity
selHist[[area]][[fleet]]$vul[[gtg]]
selPro[[area]][[fleet]]$keep[[gtg]]

# Separate historical and projection selectivity lists
MultifleetObj@fleet_selectivity_hist_list[[f]]
MultifleetObj@fleet_selectivity_proj_list[[f]]


5. Equilibrium calcs.

1) Single Fleet: solveD()

# Simple equilibrium
is <- solveD(lh, sel = selHist[[1]], doFit = TRUE, ...)

2) Multifleet: solveD_multifleet2()
# Multi-fleet equilibrium with iterative allocation
is <- solveD_multifleet2(lh = lh, sel_list = hist_sel_list, 
                        fleet_proportions = fleet_proportions,
                        allocation_type = "effort" or "catch")

# Returns: 
# - Feq (total F)
# - F_by_fleet (vector of fleet-specific F)
# - final_effort_proportions, target_catch_proportions, actual_catch_proportions


When allocation_type= "catch" 
#Iterative process to find effort proportions that achieve target catch proportions.
#We specify the target catch proportions, and the algorithm finds the effort 
#proportions needed to achieve those catch splits.

fleet_proportions <- c(0.7, 0.3)  # want Fleet 1 to catch 70%, Fleet 2 to catch 30%
allocation_type <- "catch"
# Final: effort split [0.67, 0.33] achieves target catch split [0.7, 0.3]


allocation_type = "effort"
#Specify how fishing effort (F values) should be split among fleets.
fleet_proportions <- c(0.6, 0.4)  # Fleet 1 gets 60% of effort, Fleet 2 gets 40%

# If total equilibrium F = 0.5
# Simple, direct allocation
F_fleet1 = 0.5 × 0.6 = 0.3
F_fleet2 = 0.5 × 0.4 = 0.2


6. Flet detection, mode switching, and backward compatibility (handle both single and multi)

1) mode swithcing

is_multifleet <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1

if(is_multifleet) {
    # Multifleet calculations
    nfleets <- MultifleetObj@nfleets
} else {
    # Single fleet path (unchanged)
    nfleets <- 1
}


7. Historical Effort Deviations

1) Single Fleet:

histEffortDev <- function(TimeAreaObj, StochasticObj) {
    # Creates 3D array [years, iterations, areas]
    Emult <- array(dim=c(years, iterations, areas))
}

2) Multifleet:

histEffortDev <- function(TimeAreaObj, StochasticObj, nfleets = 1) {
    # Creates 4D array [years, iterations, areas, nfleets]
    Emult <- array(dim=c(years, iterations, areas, nfleets))
    
    # Backward compatibility handling
    if(nfleets == 1) {
    # return 3D version for single fleet
    }
}

8. Management Strategy 
# Now this section includes an Area loop + Fleet loop and fleet specific data  (histEffortDev[j,k,m,f]) # 4D array (added fleet dimension))
# the returning structure was modified: return(year, iteration, area, fleet, F)

# The function now returns multiple F values per area (one per fleet + one total), 
# So evalMSE() uses:

# total F for population mortality/survival
# individual F for fleet-specific catch calculations

1) single fleet:

fixedStrategy <- function(phase, dataObject) {
    # Simple F calculation
    Flocal <- rbind(Flocal, c(j, k, m, 
        TimeAreaObj@historicalEffort[yr,m] * is$Feq * histEffortDev[j,k,m]))
}

# OLD (Single Fleet):
return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3], Flocal=Flocal[,4]))
#           4 columns: year, iteration, area, F



2) multifleet:
fixedStrategy <- function(phase, dataObject) {
    # Fleet detection
    is_multifleet <- !is.null(MultifleetObj) && MultifleetObj@nfleets >= 1
    
    if(is_multifleet) {
        # Fleet-specific F calculations
        for(f in 1:nfleets) {
            F_fleet <- F_eq_by_fleet[f] * TimeAreaObj@historicalEffort[yr,m] * 
                      histEffortDev[yr + 1, k, m, f]  # 4D indexing
        }
        # Returns fleet-specific and total F
    } else {
        # Original single fleet logic unchanged
    }
}

# NEW (Multifleet):
colnames(Flocal) <- c("year", "iteration", "area", "fleet", "Flocal")
return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3],
            fleet=Flocal[,4], Flocal=Flocal[,5]))
#           5 columns: year, iteration, area, fleet, F


#very similar integration in evalMSE() as single fleet
#Same function calls
#Same data structures (little enhanced)
#Same extraction logic (just one extra filter)
#The total F value (fleet=0) is identical to what single fleet would return


9. Population Dynamics - Time Loop
# Main differences: 
#how mortality and catches are calculated
#data structure changes (ues 3D + 4D arrays)
#nested selectivity: selGroup[[m]][[f]]
#fleet-specific catch tracking: catchNage_by_fleet[[f]][[l]]

1) single fleet:

# Simple mortality and catch calculations
Z[[l]][,j,m] <- Ftotal[j,k,m]*selGroup[[m]]$removal[[l]] + lh$LifeHistory@M
catchNage[[l]][,j,m] <- Ftotal[j,k,m]*selGroup[[m]]$keep[[l]]/
    (Z[[l]][,j,m])*(1-exp(-Z[[l]][,j,m]))*N[[l]][,j,m]


2) multifleet:
# Combined mortality from all fleets
total_fishing_mortality <- sapply(1:ageClasses, function(age) {
    sum(sapply(1:nfleets, function(f) {
        F_by_fleet_current[f] * selGroup[[m]][[f]]$removal[[l]][age]
    }))
})
Z[[l]][,j,m] <- total_fishing_mortality + lh$LifeHistory@M

# Fleet-specific catches using shared Z 
for(f in 1:nfleets) {
    catchNage_by_fleet[[f]][[l]][,j,m] <- 
        F_by_fleet_current[f] * selGroup[[m]][[f]]$keep[[l]] /
        Z[[l]][,j,m] * (1-exp(-Z[[l]][,j,m])) * N[[l]][,j,m]
}


10) Results structure

1) single fleet:

dynamics <- list(SB=SB, VB=VB, RB=RB, catchB=catchB, catchN=catchN, 
                Ftotal=Ftotal, ...)


2) multifleet:
# All original results plus fleet-specific tracking
dynamics <- list(
    # Original single fleet results (backward compatible)
    SB=SB, VB=VB, RB=RB, catchB=catchB, catchN=catchN, Ftotal=Ftotal,
    
    # NEW: Multifleet-specific results
    multifleet = list(
        Ftotal_by_fleet = Ftotal_by_fleet,
        catchB_by_fleet = catchB_by_fleet,
        catchN_by_fleet = catchN_by_fleet,
        discB_by_fleet = discB_by_fleet,
        discN_by_fleet = discN_by_fleet,
        fleet_proportions = fleet_proportions,
        final_effort_proportions = final_effort_proportions,
        target_catch_proportions = target_catch_proportions,
        actual_catch_proportions = actual_catch_proportions,
        allocation_type = allocation_type
    )
)


11. Some validations (there are more along the code)

# validation checks
sum(F_by_fleet) == Ftotal  # Fleet F values sum to total F
sum(catchB_by_fleet) == catchB  # Fleet catches sum to total catch
abs(actual_catch_proportions - target_catch_proportions) < 0.01  # 1% tolerance


12. Basic MP (the cistomized MP)
# added the fleet column to simpleMP_multi to maintain data structure 
# consistency with the multifleet fixedStrategy() format.