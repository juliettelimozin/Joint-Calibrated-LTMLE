library(cobalt)
library(ggplot2)
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/dgm_2nd_simulation_biometrics.R")

#x <- seq(0.5,4.0, 0.5)
x <- c(-5,-0.2)
X3_dist <- y_mean <- range_y <-always_treat <- never_treat <- initiators <- rep(1, length(x))

for (i in 1:length(x)){
  simdata<-DATA_GEN(ns = 200, nv = 3, conf =x[i])
  X3_dist[i] <- abs(mean(simdata[simdata$t == 2 & simdata$CA == 0,]$X3) - mean(simdata[simdata$t == 2 & simdata$CA == 3,]$X3))
  y_mean[i]<- mean(simdata[simdata$t == 2,]$Y)
  range_y[i] <- max(simdata$Y) - min(simdata$Y)
  initiators[i] <- dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]
  always_treat[i] <- dim(simdata[simdata$t == 2 & simdata$CA == 3,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 1,])[1]
  never_treat[i] <- dim(simdata[simdata$t == 2 & simdata$CA == 0,])[1]/dim(simdata[simdata$t == 0 & simdata$CA == 0,])[1]
  
  print(ggplot(simdata[simdata$t == 2 & (simdata$CA == 3 | simdata$CA == 0),], aes(x = X3, color = factor(A), fill = factor(A))) +
          geom_density(position = "identity", alpha = 0.5) +
          scale_fill_manual(values = c("red", "blue"), name = "Treatment") +
          scale_color_manual(values = c("red", "blue"), name = "Treatment") +
          labs(title = paste("Density of X by Treatment, conf =",x[i]), x = "X", y = "Density") +
          theme_minimal())
  print(ggplot(simdata[simdata$t == 2 & (simdata$CA == 3 | simdata$CA == 0),], aes(x = Y, color = factor(A), fill = factor(A))) +
          geom_density(position = "identity", alpha = 0.5) +
          scale_fill_manual(values = c("red", "blue"), name = "Treatment") +
          scale_color_manual(values = c("red", "blue"), name = "Treatment") +
          labs(title = paste("Density of Y by Treatment, conf =",x[i]), x = "Y", y = "Density") +
          theme_minimal())
  }


   


# Specify formula: treatment ~ covariates

