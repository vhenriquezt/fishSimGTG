Details about Selectivity in multifleet approach

1) The selectivity wrapper (selWrapper()) function remains unchanged in the multifleet implementation. 
What changes is how it is called and how the results are organized.

Multifleet Selectivity Structure:

# HISTORICAL SELECTIVITY
selHist <- lapply(1:TimeAreaObj@areas, function(area) {
    lapply(1:nfleets, function(f) {
        selWrapper(lh, TimeAreaObj,
                   FisheryObj = MultifleetObj@fleet_selectivity_hist_list[[f]],
                   doPlot = FALSE)
    })
})

# PROJECTION SELECTIVITY
selPro <- lapply(1:TimeAreaObj@areas, function(area) {
    lapply(1:nfleets, function(f) {
        selWrapper(lh, TimeAreaObj,
                   FisheryObj = MultifleetObj@fleet_selectivity_proj_list[[f]],
                   doPlot = FALSE)
    })
})

Resulting structure:
selHist[[area]][[fleet]]$vul[[gtg]][age]     # 4-level nesting
selHist[[area]][[fleet]]$keep[[gtg]][age]    
selHist[[area]][[fleet]]$removal[[gtg]][age]

set up:
# define fleet selectivities
fleet1_hist <- new("Fishery")
fleet1_hist@vulType <- "logistic"
fleet1_hist@vulParams <- c(8.0, 1.0)

fleet2_hist <- new("Fishery")
fleet2_hist@vulType <- "logistic"  
fleet2_hist@vulParams <- c(12.0, 2.0) 

MultifleetObj@fleet_selectivity_hist_list <- list(fleet1_hist, fleet2_hist)



The model loop is across areas and fleets:
for(area in 1:2) {
    selHist[[area]] <- list()  # initialize fleet list for this area
    
    # inner loop: Fleets
    for(fleet in 1:2) {
        # each fleet gets its own selectivity calculation
        selHist[[area]][[fleet]] <- selWrapper(
            lh, 
            TimeAreaObj,
            FisheryObj = MultifleetObj@fleet_selectivity_hist_list[[fleet]],
            doPlot = FALSE
        )
    }
}

Result:
# Area 1, Fleet 1 selectivity
selHist[[1]][[1]]  
# Area 1, Fleet 2 selectivity  
selHist[[1]][[2]] 
# Area 2, Fleet 1 selectivity (same as Area 1 since using same FisheryObj)
selHist[[2]][[1]] 
# Area 2, Fleet 2 selectivity
selHist[[2]][[2]] 

 In mortality calculations:
 #for each age, sum contributions from all fleets
total_fishing_mortality <- sapply(1:ageClasses, function(age) {
    sum(sapply(1:nfleets, function(f) {
        F_by_fleet[f] * selHist[[area]][[f]]$removal[[gtg]][age]
    }))
})

# Example with actual numbers:
# Age 5, GTG 1, Area 1
# Fleet 1: F=0.2, selectivity=0.3 → contribution = 0.06
# Fleet 2: F=0.1, selectivity=0.8 → contribution = 0.08
# Total fishing mortality at age 5 = 0.06 + 0.08 = 0.14

Fleet-Specific Catch:
# Each fleet catch uses its own selectivity but shared Z
for(f in 1:nfleets) {
    catchNage_by_fleet[[f]][[gtg]][age,year,area] <- 
        F_by_fleet[f] * selHist[[area]][[f]]$keep[[gtg]][age] /  
        Z_total[age] * (1-exp(-Z_total[age])) * N[[gtg]][age,year,area]
}


Some assumptions/decision for simplicity:

1) In the current implementation, each fleet uses the same selectivity pattern across all areas:
# Fleet 1 uses fleet1_hist in ALL areas
selHist[[1]][[1]] = selWrapper(..., fleet1_hist)  # Area 1, Fleet 1
selHist[[2]][[1]] = selWrapper(..., fleet1_hist)  # Area 2, Fleet 1 (same)

This assumes fleets operate consistently across areas. 
To have area-specific selectivity per fleet, I need to create something like this:

MultifleetObj@fleet_selectivity_hist_list <- list(
    list(fleet1_area1, fleet1_area2),  # Fleet 1 different in each area
    list(fleet2_area1, fleet2_area2)   # Fleet 2 different in each area
)

#For blanck grouper this approach should be enough 

2) Movement and benchmark clculations (for simplicity use fleet 1 selectivity)
The code uses Fleet 1 as "representative" for certain calculations:

# movement matrix uses Fleet 1 selectivity
S <- SurvMat(ageClasses, M_in=M, F_in=Ftotal[j,k,m], 
             S_in=selGroup[[m]][[1]]$removal[[l]])  # Fleet 1

# benchmark calculations use Fleet 1
refCalc <- gtgYPRWrapper_Fonly(lh=lh, sel=selHist[[1]][[1]])

In summary:
Each fleet maintains independent selectivity: No averaging or combining
Selectivity affects both F and catch: through the removal and keep curves
All fleets share the same Z: but contribute differently based on their selectivity
Selectivity determines fleet interactions: fleets targeting different sizes compete less