library(cobalt)
library(ggplot2)
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/dgm_2nd_simulation_biometrics.R")

x <- seq(0,0.7, 0.05)

X3_dist <- y_mean <- range_y <- rep(1, length(x))

for (i in 1:length(x)){
  simdata<-DATA_GEN(ns = 200, treat_prev = 0, conf =x[i])
  X3_dist[i] <- abs(mean(simdata[simdata$t == 4 & simdata$CA == 0,]$X3) - mean(simdata[simdata$t == 4 & simdata$CA == 5,]$X3))
  y_mean[i]<- mean(simdata[simdata$t == 4,]$Y)
  range_y[i] <- max(simdata$Y) - min(simdata$Y)
  print(ggplot(simdata[simdata$t == 4 & (simdata$CA == 5 | simdata$CA == 0),], aes(x = X3, color = factor(A), fill = factor(A))) +
          geom_histogram(position = "identity", alpha = 0.5, bins = 100) +
          scale_fill_manual(values = c("red", "blue"), name = "Treatment") +
          scale_color_manual(values = c("red", "blue"), name = "Treatment") +
          labs(title = paste("Density of X by Treatment, conf =",x[i]), x = "X", y = "Density") +
          theme_minimal())
  }

   


# Specify formula: treatment ~ covariates

