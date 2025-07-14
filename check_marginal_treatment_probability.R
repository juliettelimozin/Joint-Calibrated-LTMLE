# check marginal treatment probabilities

source("/home/juliette/Calibrated-weights-sequential-trial-emulation/dgm_2nd_simulation_biometrics.R")

simdata<-DATA_GEN(ns = 1000000,conf = 0.2)
                  
p0_1 <- dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]/dim(simdata[simdata$t == 0,])[1]
 
p1_1 <- dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]
p1_0 <- dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 0,])[1]

p2_1 <- dim(simdata[simdata$t == 2 & simdata$CA == 3,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]
p2_0 <- dim(simdata[simdata$t == 2 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]

p1_1
p1_0
p2_1
p2_0

simdata<-DATA_GEN(ns = 1000000,conf = 4, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -4.2, treat_prev_d1_2 = -3.8, treat_prev_d0_2 = -4.1)
p1_1 <- dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]
p1_0 <- dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 0,])[1]
p2_1 <- dim(simdata[simdata$t == 2 & simdata$CA == 3,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]
p2_0 <- dim(simdata[simdata$t == 2 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]

p1_1
p1_0
p2_1
p2_0

simdata<-DATA_GEN(ns = 1000000,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)
p1_1 <- dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]
p1_0 <- dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 0,])[1]
p2_1 <- dim(simdata[simdata$t == 2 & simdata$CA == 3,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]
p2_0 <- dim(simdata[simdata$t == 2 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]

p1_1
p1_0
p2_1
p2_0

