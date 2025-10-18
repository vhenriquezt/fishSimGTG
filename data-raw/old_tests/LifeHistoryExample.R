

devtools::load_all()
library(LBSPR)

#Life history
LifeHistoryExample<-new("LifeHistory")
LifeHistoryExample@Linf<-100
LifeHistoryExample@L50<-66
LifeHistoryExample@L95delta<-1
LifeHistoryExample@MK<-1.5
LifeHistoryExample@LW_A<-0.01
LifeHistoryExample@LW_B<-3
LifeHistoryExample@title<-"Example fish"
LifeHistoryExample@shortDescription<-"Simulated life history of a fish based on B-H invariants"
LifeHistoryExample@speciesName<-"Example fish"
LifeHistoryExample@L_type<-"TL"
LifeHistoryExample@L_units<-"cm"
LifeHistoryExample@Walpha_units<-"g"
LifeHistoryExample@K<-0.2
LifeHistoryExample@M<-0.3
LifeHistoryExample@t0<-0
LifeHistoryExample@Tmax<- floor(-log(0.01)/0.3)
LifeHistoryExample@Steep<-0.99
LifeHistoryExample@R0<-1000
LifeHistoryExample@recSD<-0.6
LifeHistoryExample@recRho<-0
LifeHistoryExample@isHermaph<-FALSE

usethis::use_data(LifeHistoryExample, overwrite = TRUE)

#Use life history in YPR sim examples
lbsprSimExample<-lbsprSimWrapper(LifeHistoryExample)
usethis::use_data(lbsprSimExample, overwrite = TRUE)

#lbsprSimExample<-lbsprSimWrapper(LifeHistoryExample)
#usethis::use_data(lbsprSimExample, overwrite = TRUE)

#New update
gtgSimExample<-gtgYPRWrapper(LifeHistoryExample)
usethis::use_data(gtgSimExample, overwrite = TRUE)



