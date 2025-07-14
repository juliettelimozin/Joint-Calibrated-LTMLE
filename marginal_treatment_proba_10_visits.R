# check marginal treatment probabilities

source("/home/juliette/Calibrated-weights-sequential-trial-emulation/dgm_10_visits.R")

simdata<-DATA_GEN_TEN(ns = 500000,conf = 0.2)

p0_1 <- dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]/dim(simdata[simdata$t == 0,])[1]
p0_0 <- dim(simdata[simdata$t == 0 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 0,])[1]

p1_1 <- dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]
p1_0 <- dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 0,])[1]

p2_1 <- dim(simdata[simdata$t == 2 & simdata$CA == 3,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]
p2_0 <- dim(simdata[simdata$t == 2 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]

p3_1 <- dim(simdata[simdata$t == 3 & simdata$CA == 4,])[1]/dim(simdata[simdata$t == 2 & simdata$CA == 3,])[1]
p3_0 <- dim(simdata[simdata$t == 3 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 2 & simdata$CA == 0,])[1]

p4_1 <- dim(simdata[simdata$t == 4 & simdata$CA == 5,])[1]/dim(simdata[simdata$t == 3 & simdata$CA == 4,])[1]
p4_0 <- dim(simdata[simdata$t == 4 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 3 & simdata$CA == 0,])[1]

p5_1 <- dim(simdata[simdata$t == 5 & simdata$CA == 6,])[1]/dim(simdata[simdata$t == 4 & simdata$CA == 5,])[1]
p5_0 <- dim(simdata[simdata$t == 5 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 4 & simdata$CA == 0,])[1]

p6_1 <- dim(simdata[simdata$t == 6 & simdata$CA == 7,])[1]/dim(simdata[simdata$t == 5 & simdata$CA == 6,])[1]
p6_0 <- dim(simdata[simdata$t == 6 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 5 & simdata$CA == 0,])[1]

p7_1 <- dim(simdata[simdata$t == 7 & simdata$CA == 8,])[1]/dim(simdata[simdata$t == 6 & simdata$CA == 7,])[1]
p7_0 <- dim(simdata[simdata$t == 7 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 6 & simdata$CA == 0,])[1]

p8_1 <- dim(simdata[simdata$t == 8 & simdata$CA == 9,])[1]/dim(simdata[simdata$t == 7 & simdata$CA == 8,])[1]
p8_0 <- dim(simdata[simdata$t == 8 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 7 & simdata$CA == 0,])[1]

p9_1 <- dim(simdata[simdata$t == 9 & simdata$CA == 10,])[1]/dim(simdata[simdata$t == 8 & simdata$CA == 9,])[1]
p9_0 <- dim(simdata[simdata$t == 9 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 8 & simdata$CA == 0,])[1]

p1_1
p1_0
p2_1
p2_0
p3_1
p3_0
p4_1
p4_0
p5_1
p5_0
p6_1
p6_0
p7_1
p7_0
p8_1
p8_0
p9_1
p9_0

p0_1*p1_1*p2_1*p3_1*p4_1*p5_1*p6_1*p7_1*p8_1*p9_1
p0_0*p1_0*p2_0*p3_0*p4_0*p5_0*p6_0*p7_0*p8_0*p9_0

ggplot(simdata[simdata$t == 9 & (simdata$CA == 10 | simdata$CA == 0),], aes(x = X3, color = factor(A), fill = factor(A))) +
  geom_density(position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), name = "Treatment") +
  scale_color_manual(values = c("red", "blue"), name = "Treatment") +
  labs(title = paste("Density of X by Treatment, conf =",x[i]), x = "X", y = "Density") +
  theme_minimal()

simdata<-DATA_GEN_TEN(ns = 500000,conf = 2.5, 
                  treat_prev_d1 = c(2,-0.5,-3,-5.5,-8,-10.5,-13,-15.5,-18), 
                  treat_prev_d0 = c(-4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55))


p0_1 <- dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]/dim(simdata[simdata$t == 0,])[1]
p0_0 <- dim(simdata[simdata$t == 0 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 0,])[1]

p1_1 <- dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]
p1_0 <- dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 0,])[1]

p2_1 <- dim(simdata[simdata$t == 2 & simdata$CA == 3,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 2,])[1]
p2_0 <- dim(simdata[simdata$t == 2 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 1 & simdata$CA == 0,])[1]

p3_1 <- dim(simdata[simdata$t == 3 & simdata$CA == 4,])[1]/dim(simdata[simdata$t == 2 & simdata$CA == 3,])[1]
p3_0 <- dim(simdata[simdata$t == 3 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 2 & simdata$CA == 0,])[1]

p4_1 <- dim(simdata[simdata$t == 4 & simdata$CA == 5,])[1]/dim(simdata[simdata$t == 3 & simdata$CA == 4,])[1]
p4_0 <- dim(simdata[simdata$t == 4 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 3 & simdata$CA == 0,])[1]

p5_1 <- dim(simdata[simdata$t == 5 & simdata$CA == 6,])[1]/dim(simdata[simdata$t == 4 & simdata$CA == 5,])[1]
p5_0 <- dim(simdata[simdata$t == 5 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 4 & simdata$CA == 0,])[1]

p6_1 <- dim(simdata[simdata$t == 6 & simdata$CA == 7,])[1]/dim(simdata[simdata$t == 5 & simdata$CA == 6,])[1]
p6_0 <- dim(simdata[simdata$t == 6 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 5 & simdata$CA == 0,])[1]

p7_1 <- dim(simdata[simdata$t == 7 & simdata$CA == 8,])[1]/dim(simdata[simdata$t == 6 & simdata$CA == 7,])[1]
p7_0 <- dim(simdata[simdata$t == 7 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 6 & simdata$CA == 0,])[1]

p8_1 <- dim(simdata[simdata$t == 8 & simdata$CA == 9,])[1]/dim(simdata[simdata$t == 7 & simdata$CA == 8,])[1]
p8_0 <- dim(simdata[simdata$t == 8 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 7 & simdata$CA == 0,])[1]

p9_1 <- dim(simdata[simdata$t == 9 & simdata$CA == 10,])[1]/dim(simdata[simdata$t == 8 & simdata$CA == 9,])[1]
p9_0 <- dim(simdata[simdata$t == 9 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 8 & simdata$CA == 0,])[1]

p1_1
p1_0
p2_1
p2_0
p3_1
p3_0
p4_1
p4_0
p5_1
p5_0
p6_1
p6_0
p7_1
p7_0
p8_1
p8_0
p9_1
p9_0

p0_1*p1_1*p2_1*p3_1*p4_1*p5_1*p6_1*p7_1*p8_1*p9_1
p0_0*p1_0*p2_0*p3_0*p4_0*p5_0*p6_0*p7_0*p8_0*p9_0

ggplot(simdata[simdata$t == 9 & (simdata$CA == 10 | simdata$CA == 0),], aes(x = X3, color = factor(A), fill = factor(A))) +
  geom_density(position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), name = "Treatment") +
  scale_color_manual(values = c("red", "blue"), name = "Treatment") +
  labs(title = paste("Density of X by Treatment, conf =",x[i]), x = "X", y = "Density") +
  theme_minimal()
