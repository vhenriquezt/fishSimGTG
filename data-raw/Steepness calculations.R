
#Steepnes re-scaling of beta distribution

#Regular
a <- 3.89
b <- 1.52

a/(a+b)
(a - 1)/(a + b - 2)

X<-rbeta(100000,a,b)



#Scaled
ascaled <- 2.6
bscaled <- 1.42
upper <- 1
lower <- 0.21

ascaled/(ascaled+bscaled)*(upper-lower)+lower
(ascaled - 1)/(ascaled + bscaled - 2)*(upper-lower)+lower

Xscaled <- rbeta(100000,ascaled,bscaled)*(upper-lower)+lower

par(mfrow = c(2,1))
hist(X)
hist(Xscaled)

mean(X)
mean(Xscaled)

quantile(X, probs = c(0.025, 0.5, 0.75, 0.975))
quantile(Xscaled, probs = c(0.025, 0.5, 0.75, 0.975))
