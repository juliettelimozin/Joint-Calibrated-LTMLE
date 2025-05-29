rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))

GetDataForAnalysis <- function(end.time) {
  load("LongTMLE/data/dataCd4")
  uniquePatients <- unique(z$patient)
  n <- length(uniquePatients)
  L <- C1 <- C2 <- A1 <- matrix(NA, nrow=n, ncol=end.time)
  Y <- matrix(NA, nrow=n, ncol=end.time + 1)
  ID <- W1.factor <- W3 <- W4 <- end.time.patient <- rep(NA, n)
  
  for (p in 1:n) {
    index <- which(z$patient == uniquePatients[p])
    ID[p] <- GetBaseline(z$patient[index])
    W1.factor[p] <- GetBaseline(z$bas_age_gp[index])
    W3[p] <- GetBaseline(z$gender[index])
    W4[p] <- GetBaseline(z$adv_stage[index])
    for (i in index) {
      t <- z$qrt[i]
      if (t <= (end.time + 1)) {
        Y[p, t] <- z$death[i]
      }
      if (t <= end.time && !z$death[i]) {  #don't record other info in quarter of death
        L[p, t] <- sqrt(z$cd4[i])
        C1[p, t] <- 1 - z$cens[i]  #code uses 1 as uncensored (data has 1 as censored)
        C2[p, t] <- 1 - z$lost_yn2[i] #code uses 1 as uncensored (data has 1 as censored)
        A1[p, t] <- z$sl_yn[i]
      }
    }
    end.time.patient[p] <- max(z$qrt[index])
  }

  stopifnot(setequal(W1.factor, 1:3))
  W1.factor <- factor(W1.factor)
  W1 <- as.numeric(W1.factor == 1) #could leave these as factors, but will convert to binaries to match simulation 
  W2 <- as.numeric(W1.factor == 2)
  
  #condition on being alive after qrt 1
  index <- Y[, 1] == 0
  ID <- ID[index]
  W1 <- W1[index]
  W2 <- W2[index]
  W3 <- W3[index]
  W4 <- W4[index]
  L <- L[index, ]
  C1 <- C1[index, ]
  C2 <- C2[index, ]
  A1 <- A1[index, ]
  Y <- Y[index, ]
  
  for (t in 2:end.time) {
    A1[A1[, t-1] == 1, t] <- 1 #if A1 was 1 last period, it stays 1
  }
  for (t in 2:(end.time + 1)) {
    Y[Y[, t-1] == 1, t] <- 1 #if Y was 1 last period, it stays 1
  }
  
  return(data.frame(ID, CreateDataFrame.matchdata(W1, W2, W3, W4, L, C1, C2, A1, Y, end.time)))
}

GetBaseline <- function(x) {
  u <- unique(x)
  if (length(u) != 1) stop ("length(u) != 1")
  return(u)
}

GenerateData.rct <- function(sim.opts, abar=NULL) {
  # W A1 A2 A3 ... Y
  est.psi0 <- !is.null(abar) 
  end.time <- sim.opts$end.time
  n <- sim.opts$n
  
  W <- rnorm(n)
  if (est.psi0) { 
    mean.A <- mean(abar)
  } else {
    if (sim.opts$confounder) {
      p <- ltmle:::Bound(sim.opts$prob.a + W, bounds=c(0.38, 0.62)) 
    } else {
      p <- sim.opts$prob.a
    }
    g0 <- matrix(p, nrow=n, ncol=end.time) #note: this is recycling p in confounder case
    A <- matrix(rbinom(n * end.time, size=1, prob=g0), nrow=n, ncol=end.time)  
    mean.A <- rowMeans(A)
  }
  if (sim.opts$linear.Y) {
    Y <- rexpit(W - mean.A)
  } else {
    Y <- rexpit(W - 4*mean.A^2)
  }
  if (est.psi0) return(mean(Y))
  data <- data.frame(W, A=A, Y)
  if (end.time == 1) names(data)[2] <- "A.1" #if end.time>1, this is automatic
  
  return(data)
}

GenerateData.matchdata <- function(sim.opts, abar=NULL) {
  #match data coefficients
  A1_CONSTANT_TERM <- -5
  C1_CONSTANT_TERM <- 2
  L_COEFF <- 1.5
  PATIENT_SEEN_CONSTANT <- 0.4
  
  Y_A_COEFF <- 0.9
  Y_CONSTANT_TERM <- -5.8
  Y_L_COEFF <- 0.7
  MIN_L <- -4
  MAX_L <- 4
  
  CalcL <- function(W1=NULL, W2=NULL, W3=NULL, W4, PrevL=NULL, A1=NULL) {
    n <- length(W4)
    if (is.null(W1)) {
      L <- rnorm(n) - W4  #this is for t==1
    } else {
      L <- 0.1*W1 - 0.1*W2 - 0.1*W3 - 0.5*W4 + 0.9*PrevL + A1 + rnorm(n)
    }
    L <- ltmle:::Bound(L, c(MIN_L, MAX_L))
    return(L)
  }
  
  CalcY <- function(W1, W2, W3, W4, L, A1) {
    Y <- rexpit(Y_CONSTANT_TERM - 0.1*W1 - 0.1*W2 + 0.1*W3 - 0.2*W4 - Y_L_COEFF*L - Y_A_COEFF*A1)
    return(Y)
  }
  est.psi0 <- !is.null(abar)
  end.time <- sim.opts$end.time
  n <- sim.opts$n
  
  #Data generating process:
  # W1 W2 W3 W4 Y0=0 L0 C1_0 C2_0 A0 Y1 L1 C1_1 C2_1 A1 ... Y9 L9 C1_9 C2_9 A9 Y10
  
  W1 <- rbinom(n, 1, 0.3)  #age<30
  W2 <- (1-W1) * rbinom(n, 1, 0.5) #age 30-39
  W3 <- rbinom(n, 1, 0.5) #gender
  W4 <- rbinom(n, 1, 0.3) #advanced stage  
  
  L <- observedL <- C1 <- A1 <- C2 <- patient.seen <- matrix(NA, nrow=n, ncol=end.time)
  Y <- matrix(NA, nrow=n, ncol=end.time+1)
  
  uncensored.alive <- rep(TRUE, nrow=n)
  alive <- rep(TRUE, nrow=n)
  #store in matrix as 1:end.time but then record names as 0:end.time-1 (except Y1:Yend.time)
  #the coding for Y is slightly inefficient but I think it's more readable
  for (t in 1:(end.time + 1)) {
    if (t==1) {
      Y[, 1] <- 0
      patient.seen[, 1] <- rep(1, n) 
      L[, 1] <- CalcL(W4=W4)  
      observedL[, 1] <- L[, 1]
    } else {
      Y[uncensored.alive, t] <- CalcY(W1[uncensored.alive], W2[uncensored.alive], W3[uncensored.alive], W4[uncensored.alive], L[uncensored.alive, t-1], A1[uncensored.alive, t-1])
      Y[Y[, t-1]==1, t] <- 1 #if Y was 1 last period, it stays 1 (this doesn't explicitly condition on uncensored.alive, but if you died, you must have been uncensored; we set Y=1 even after death )
      uncensored.alive <- uncensored.alive & (!Y[,t])  
      if (t == (end.time + 1)) break #don't calculate after end.time+1 for other variables
      
      patient.seen[uncensored.alive, t] <- rexpit(PATIENT_SEEN_CONSTANT + (0.1*W1 - 0.2*W2 + 0.3*W3 + 0.1*W4)[uncensored.alive] - 0.1*observedL[uncensored.alive, t-1] + 0.2*A1[uncensored.alive, t-1])
      L[uncensored.alive, t] <- CalcL(W1=W1[uncensored.alive], W2=W2[uncensored.alive], W3=W3[uncensored.alive], W4=W4[uncensored.alive], PrevL=L[uncensored.alive, t-1], A1=A1[uncensored.alive, t-1])
      
      observedL[uncensored.alive & patient.seen[, t], t] <- L[uncensored.alive & patient.seen[,t], t] #if patient.seen, L is observed
      observedL[uncensored.alive & !patient.seen[, t], t] <- observedL[uncensored.alive & !patient.seen[, t], t-1] #if patient not seen, use L from last period
    }
    
    if (est.psi0) {
      C1[uncensored.alive,t] <- 1  #For estimating PSI0 (all uncensored)
    } else {
      C1[uncensored.alive, t] <- rexpit(C1_CONSTANT_TERM + (0.1*W1 + 0.2*W2 + 0.1*W3 + 0.1*W4 + 0.1*L[, 1])[uncensored.alive]) #C1 depends only on baseline L, not Lt
      uncensored.alive <- uncensored.alive & C1[,t]    
    }  
    
    if (est.psi0) {
      C2[uncensored.alive, t] <- 1  #For estimating PSI0 (all uncensored)
    } else {
      if (t <= 3) { #there's currently no risk of C2 at t=3 since patient is always seen at t=1
        C2[uncensored.alive, t] <- 1
      } else {
        C2[uncensored.alive, t] <- as.double(apply(patient.seen[uncensored.alive, (t-2):t]==1, 1, any)) #C2 is 0 if patient has not been seen for last 3 periods
      }
      uncensored.alive <- uncensored.alive & C2[, t]
    }
    
    if (est.psi0) {
      A1[uncensored.alive, t] <- abar[t]
    } else {
      A1[uncensored.alive, t] <- rexpit(A1_CONSTANT_TERM + (0.1*W1 + 0.1*W2 + 0.2*W3 + 0.2*W4)[uncensored.alive] - L_COEFF*L[uncensored.alive, t] + rnorm(sum(uncensored.alive)))
      
      A1[uncensored.alive & !patient.seen[, t], t] <- 0 #if patient not seen this period and was not on treatment, treatment can't start (if was on treatment, we will set to stay on treatment below)
      
      if (t > 1) A1[uncensored.alive & A1[, t-1]==1, t] <- 1 #if A1 was 1 last period, it stays 1  
    }
  } 
  
  if (est.psi0) {
    psi0 <- mean(Y[, end.time + 1])   
    return(psi0)
  }
  
  return(CreateDataFrame.matchdata(W1, W2, W3, W4, L=observedL, C1, C2, A1, Y, end.time))
}

CreateDataFrame.matchdata <- function(W1, W2, W3, W4, L, C1, C2, A1, Y, end.time) {
  stopifnot(all(Y[, 1] == 0)) #Y0 should be all 0
  d <- data.frame(W1, W2, W3, W4)
  for (t in 1:(end.time + 1)) {
    if (t > 1) {
      d <- data.frame(d, Y[, t])   #leave Y0 out since it's all 0
      names(d)[ncol(d)] <- paste0("Y", t-1)
    }
    if (t == (end.time + 1)) {
      break  #record Y from 1 to end.time, all others from 0 to end.time-1
    }
    d <- data.frame(d, L[, t], BinaryToCensoring(is.uncensored=C1[, t]), BinaryToCensoring(is.uncensored=C2[, t]), A1[, t])
    names(d)[ncol(d) - 3] <- paste0("L", t-1)
    names(d)[ncol(d) - 2] <- paste0("C1_", t-1)
    names(d)[ncol(d) - 1] <- paste0("C2_", t-1)
    names(d)[ncol(d)] <- paste0("A1_", t-1)
  }
  return(d)
}

GetForm.rct <- function(sim.opts, nodes) {
  if (sim.opts$gcorrect && sim.opts$confounder) {
    gform <- paste(nodes$names[nodes$AC], "~ W")
  } else {
    gform <- paste(nodes$names[nodes$AC], "~ 1")
  }
  
  if (sim.opts$Qcorrect) {
    sumA.string <- paste0("A.", 1:sim.opts$end.time, collapse=" + ")
    if (sim.opts$linear.Y) {
      Qform <- paste0("Q.kplus1 ~ W + I(", sumA.string, ")")
    } else {
      Qform <- paste0("Q.kplus1 ~ W + I((", sumA.string, ") ^ 2)")
    }
  } else {
    Qform <- "Q.kplus1 ~ 1"
  }
  names(Qform) <- "Y"
  return(list(g=gform, Q=Qform))
}

GetForm.matchdata <- function(sim.opts, nodes) {
  end.time <- sim.opts$end.time
  Wstr <- " ~ W1 + W2 + W3 + W4"
  if (sim.opts$gcorrect) {
    gform <- NULL
    for (t in 0:(end.time - 1)) {
      C1form <- paste0("C1_", t, Wstr, " + L0")
      Aform <- paste0("A1_", t, Wstr, " + L", t)
      if (t <= 2) {
        C2form <- paste0("C2_", t, " ~ 1")
      } else {
        C2form <- paste0("C2_", t, Wstr, " + L", t - 1, " + A1_", t - 1)
      }
      gform <- c(gform, C1form, C2form, Aform)
    }
  } else {
    gform <- paste0(nodes$names[nodes$AC], " ~ 1")
  }
  
  if (sim.opts$Qcorrect) {
    Qform <- paste0("Q.kplus1", Wstr, " + L", 0:(end.time - 1), " + A1_", 0:(end.time - 1))
  } else {
    Qform <- rep("Q.kplus1 ~ 1", length(nodes$Y)) 
  }
  names(Qform) <- nodes$names[nodes$Y] 
  return(list(g=gform, Q=Qform))
}

GetRegimens.matchdata <- function(sim.opts) {
  SwitchTime <- function(x) {
    padded.x <- c(x, 1)
    switch.time <- min(which(padded.x==1)) - 1 #if never switch, switch time = end.time
    return(switch.time)
  }
  
  n <- sim.opts$n
  end.time.set <- 1:sim.opts$end.time
  end.time <- sim.opts$end.time
  regimens <- array(NaN, dim=c(n, end.time, end.time + 1)) #patients x time x regimens
  for (t in 1:(end.time + 1)) {
    abar <- rep(0, end.time)
    if (t <= end.time) abar[t:end.time] <- 1  #abar[end.time+1] is all zeros (never switch)
    regimens[, , t] <- matrix(abar, nrow=n, ncol=end.time, byrow=TRUE)
  }
  
  num.regimens <- dim(regimens)[3]
  summary.measures <- array(dim=c(num.regimens, 2, length(end.time.set))) #regimens x summary.measures x final.ynodes
  for (i in 1:length(end.time.set)) {
    for (j in 1:num.regimens) {
      summary.measures[j, , i] <- cbind(SwitchTime(regimens[1, 1:end.time.set[i], j]), end.time.set[i])
    }
  }
  colnames(summary.measures) <- c("switch.time", "TIME")
  
  return(list(regimens=regimens, summary.measures=summary.measures))
}

GetRegimens.rct <- function(sim.opts) {
  n <- sim.opts$n
  end.time <- sim.opts$end.time
  
  library("combinat")
  reg <- hcube(c(rep(2, end.time))) - 1
  num.regimens <- nrow(reg)
  summary.measures <- array(rowSums(reg), dim=c(num.regimens, 1, 1)) #regimens x summary.measures x final.ynodes
  colnames(summary.measures) <- "sumA"
  
  regimens <- array(NaN, dim=c(n, end.time, nrow(reg))) #patients x time x regimens
  for (i in 1:nrow(reg)) {
    abar <- reg[i, ]
    regimens[, , i] <- matrix(abar, nrow=n, ncol=end.time, byrow=TRUE)
  }
  
  return(list(regimens=regimens, summary.measures=summary.measures))
}

GetSmallRegimens <- function(sim.opts) {
  #For static regimens, we just need the regimens for n=1 (saves time)
  sim.opts$n <- 1
  return(GetRegimens(sim.opts))
}

GetTrueBeta <- function(sim.opts) { 
  if (sim.opts$calc.true.beta) {
    EstTrueBeta <- function(i) {
      if (class(sim.opts)=="matchdata") cat("In EstTrueBeta i = ", i,"\n")
      set.seed(1000 + i + sim.opts$base.seed)
      plist <- EstPsi0(sim.opts)
      #psi0 and weights are matricies num.regimens x num.final.Ynodes with NA for duplicates
      beta <- ltmle:::FitMSM(plist$psi0, GetSmallRegimens(sim.opts)$summary.measures, sim.opts$working.msm, IC=NULL, plist$weights)
      return(beta)
    }
    sim.opts$n <- sim.opts$n.for.true.beta
    beta.set <- t(simplify2array(mclapply(1:sim.opts$niter.for.true.beta, EstTrueBeta, mc.cores=detectCores())))
    print(summary(beta.set))
    beta <- colMeans(beta.set)
  } else {
    s <- GetSmallRegimens(sim.opts)$summary.measures
    num.betas <- ncol(model.matrix(as.formula(sim.opts$working.msm), data=data.frame(Y=1, ltmle:::drop3(s[, , 1, drop=FALSE]))))
    beta <- rep(NA, num.betas)
  }
  return(beta)
}

EstPsi0.matchdata <- function(sim.opts) {
  end.time.set <- 1:sim.opts$end.time
  num.end.times <- length(end.time.set)
  num.regimens <- dim(GetSmallRegimens(sim.opts)$regimens)[3] #number of regimens at max time
  psi0 <- weights <- matrix(nrow=num.regimens, ncol=num.end.times)
  
  regimens <- GetRegimens(sim.opts)$regimens 
  data.full <- GenerateData(sim.opts) #data where A is not set to abar
  Anodes <- grep("^A", names(data.full))
  Ynodes <- grep("^Y", names(data.full))
  Cnodes <- grep("^C", names(data.full))

  sim.opts2 <- sim.opts
  for (j in 1:num.end.times) {
    end.time <- end.time.set[j] 
    sim.opts2$end.time <- end.time 
    final.Ynode <- Ynodes[end.time]
    num.regimens <- dim(GetSmallRegimens(sim.opts2)$regimens)[3] #number of regimens at time = end.time.set[j]  (to avoid duplicates)
    for (i in 1:num.regimens) {
      abar <- regimens[1, , i]  #regimens are the same for all patients
      psi0[i, j] <- GenerateData(sim.opts2, abar=abar)
      weights[i, j] <- ltmle:::ComputeGA(data.full, Anodes, Cnodes, ltmle:::drop3(regimens[, , i, drop=F]), final.Ynode, weight.msm=sim.opts2$use.weights)
    }
  }
  return(list(psi0=psi0, weights=weights))
}

EstPsi0.rct <- function(sim.opts) { 
  regimens <- GetSmallRegimens(sim.opts)$regimens
  num.regimens <- dim(regimens)[3]
  
  psi0 <- numeric(num.regimens)
  weights <- rep(1, num.regimens)
  for (i in 1:num.regimens) {
    abar <- regimens[1, , i]
    psi0[i] <- GenerateData(sim.opts, abar)
    if (sim.opts$use.weights) {
      weights[i] <- sim.opts$prob.a^(sum(abar == 1)) * (1 - sim.opts$prob.a)^(sum(abar == 0)) 
    }
  }
  return(list(psi0=as.matrix(psi0), weights=as.matrix(weights)))
}

GetForm <- function(sim.opts, nodes) UseMethod("GetForm")
GetRegimens <- function(sim.opts) UseMethod("GetRegimens")
GenerateData <- function(sim.opts, abar=NULL) UseMethod("GenerateData")
EstPsi0 <- function(sim.opts) UseMethod("EstPsi0")

CreateSimOpts <- function(sim.class, n, end.time, niter, gcorrect, Qcorrect, use.weights, name, calc.true.beta=T, base.seed=0, linear.Y=NA, confounder=NA, prob.a=NA, n.for.true.beta=10000, niter.for.true.beta=ifelse(sim.class=="rct",10000, 500)) {
  obj <- list(n=n, end.time=end.time, niter=niter, gcorrect=gcorrect, Qcorrect=Qcorrect, use.weights=use.weights, name=name, calc.true.beta=calc.true.beta, base.seed=base.seed, n.for.true.beta=n.for.true.beta, niter.for.true.beta=niter.for.true.beta)
  
  rct.args <- c(linear.Y, confounder, prob.a)
  if (sim.class == "rct") {
    obj$linear.Y <- linear.Y
    obj$confounder <- confounder
    obj$prob.a <- prob.a
    obj$standard.dynamic.iptw <- TRUE #also use standard dynamic (this was false in older verions, now true in both cases)
    obj$working.msm <- "Y ~ sumA"
    obj$lower.gbound <- 0.001
    obj$det.g.fun <- NULL
    if (any(is.na(rct.args))) stop("missing argument for rct")
  } else if (sim.class == "matchdata") {
    obj$standard.dynamic.iptw <- TRUE #also use standard dynamic 
    obj$working.msm <- "Y ~ TIME + I(pmax(TIME - switch.time, 0))"
    obj$lower.gbound <- 0.01
    obj$det.g.fun <- MaintainTreatment
    if (!all(is.na(rct.args))) stop("extra argument for matchdata")
  } else {
    stop("bad sim.class")
  }
  
  class(obj) <- sim.class
  return(obj)
}

GetNodes <- function(data) {
  node.names <- names(data)
  Anodes <- grep("^A", node.names)
  Ynodes <- grep("^Y", node.names)
  Cnodes <- grep("^C", node.names)
  Lnodes <- grep("^L", node.names)[-1] #ok for rct and match.data, but not in general
  return(list(A=Anodes, Y=Ynodes, L=Lnodes, C=Cnodes, names=node.names, AC=sort(c(Anodes, Cnodes))))
}

OutputLatex <- function(r, caption, col.index) {
  library("Hmisc")
  colnames(r$bias) <- colnames(r$beta) <- colnames(r$inCI) <- paste0("$\\hat{\\beta}_", 0:(ncol(r$bias)-1),"$")
  
  bias <- colMeans(r$bias)
  v <- apply(r$beta, c(2,3), var)
  bse <- bias / sqrt(v)
  mse <- colMeans(r$bias ^ 2)
  coverage <- colMeans(r$inCI)
  
  tabl <-  rbind(bias, bse, v, mse, coverage)
  tabl <- tabl[, col.index]
  rgroup <- c("Bias", "Bias / SE", "Variance", "MSE", "Coverage")
  num.betas <- ncol(r$bias)
  latex(tabl, file="", title="", dec=4, caption=caption, n.rgroup=rep(num.betas, 5), rgroup=rgroup)  #print latex code to screen
  #print(latex(tabl, file="test.tex", append=F, title="", dec=4, caption=caption, n.rgroup=rep(num.betas, 5), rgroup=rgroup)) #show nicely formatted output
  #rownames(tabl) <- rep(rgroup, each=num.betas)
  #print(tabl)
  invisible(tabl)
}

OutputLatexBootstrap <- function(result, bootstrap.est) {
  library("Hmisc")
  estimators <- colnames(bootstrap.est)
  
  tabl <- NULL
  for (j in 1:4) {
    tabl <- rbind(tabl, (t(rbind(est=result$rBeta[[j]], "95\\% CI (IC)"=result$rCI[[j]][, 1], "95\\% CI (IC)"=result$rCI[[j]][, 2], "95\\% CI (boot)"=apply(bootstrap.est[,j,], 1, quantile, 0.025),  "95\\% CI (boot)"=apply(bootstrap.est[,j,], 1, quantile, 0.975)))))
  }
  rownames(tabl)[c(1,4,7,10)] <- "$\\hat{\\beta}_0$"
  rownames(tabl)[c(2,5,8,11)] <- "$\\hat{\\beta}_1$"
  rownames(tabl)[c(3,6,9,12)] <- "$\\hat{\\beta}_2$"
  latex(tabl, file="", append=F, title="", dec=4, caption="Data Analysis", n.rgroup=rep(3, length(estimators)), rgroup=estimators) #print latex code to screen
  #print( latex(tabl, file="test.tex", append=F, title="", dec=4, caption="Data Analysis", n.rgroup=rep(3, length(estimators)), rgroup=estimators)) #show nicely formatted output
  invisible(tabl)
}

GetConfInt <- function(r, estimator="tmle") {
  return(cbind(as.numeric(summary(r, estimator=estimator)$cmat[, "CI 2.5%"]), as.numeric(summary(r, estimator=estimator)$cmat[, "CI 97.5%"])))
}

StandardDynamicIPTW <- function(data, nodes, working.msm, regimens, summary.measures, final.Ynodes, cum.g) {
  library(geepack)
  num.final.Ynodes <- length(final.Ynodes)
  Y.vec <- X.mat <- weight.vec <- id <- NULL
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- final.Ynodes[j]
    is.duplicate <- duplicated(regimens[, nodes$A < final.Ynode, , drop=FALSE], MARGIN=3)
    for (i in which(! is.duplicate)) {
      abar <- ltmle:::drop3(regimens[, , i, drop=F])
      abar <- abar[, nodes$A < final.Ynode, drop=FALSE]
      
      uncensored <- ltmle:::IsUncensored(data, nodes$C, final.Ynode)
      intervention.match <- ltmle:::InterventionMatch(data, abar, nodes$A, final.Ynode)  
      index <- uncensored & intervention.match
      
      col.index <- which.max(nodes$AC[nodes$AC < final.Ynode]) #this is right (fixes old bug)
      
      Y <- data[index, final.Ynode]
      g <- cum.g[index, col.index, i] 
      X <- ltmle:::repmat(summary.measures[i, , j, drop=FALSE], sum(index), 1)
      
      gA <- sum(index) / nrow(data)
      weight <- gA / g
      
      Y.vec <- c(Y.vec, Y)
      X.mat <- rbind(X.mat, X)
      weight.vec <- c(weight.vec, weight) 
      id <- c(id, which(index))
    }
  }
  sort.index <- sort.int(id, index.return=TRUE)$ix
  Y.vec <- Y.vec[sort.index]
  X.mat <- X.mat[sort.index, , drop=FALSE]
  weight.vec <- weight.vec[sort.index]
  id <- id[sort.index]
  colnames(X.mat) <- colnames(summary.measures)
  m <- geeglm(formula(working.msm), family="binomial", data=data.frame(Y=Y.vec, X.mat, weight.vec, id), weights=weight.vec, id=id)
  beta <- coef(m)
  std.dev <- summary(m)$coef[,"Std.err"]
  CI <- ltmle:::GetCI(beta, std.dev)
  #PlotCurve(m, end.time=max(X.mat[,"TIME"]))
  return(list(beta=beta, CI=CI, msm=m, IC=matrix(NA, nrow(data), length(beta))))
}

PlotCurve <- function(coef.msm, end.time) {
  TIME <- 1:end.time
  y <- matrix(nrow=end.time+1, ncol=end.time)
  for (switch.time in 0:end.time) { 
    X <- cbind(1, TIME, pmax(TIME - switch.time, 0))
    y[switch.time+1, ] <- 1 - plogis(X %*% coef.msm)
    cat("switch time = ", switch.time, "\n")
    print(data.frame(TIME, y=y[switch.time+1, ]))
  }
  plot(x=c(1, end.time), y=c(min(y), 1), type="n", xlab="time", ylab="Survival Probability")
  for (switch.time in 1:(end.time+1)) { 
    lines(TIME, y[switch.time, ])
  }
  invisible(y)
}


TestMSM <- function(sim.opts) {
  niter <- sim.opts$niter
  rList <- mclapply(1:niter, RunIter, mc.cores=detectCores(), sim.opts=sim.opts)
  set.seed(1000 + sim.opts$base.seed)
  true.beta <- GetTrueBeta(sim.opts)
  estimators <- names(rList[[1]][[1]])
  num.estimators <- length(estimators)
  bias <- beta <- inCI <- array(dim=c(niter, length(true.beta), num.estimators), dimnames=list(NULL, names(true.beta), estimators))
  IC <- array(dim=c(niter, sim.opts$n, length(true.beta), num.estimators))
  CI <- array(dim = c(niter, length(true.beta), 2, num.estimators), dimnames = list(NULL, names(true.beta), c("lower", "upper"), estimators))
  maxIC <- matrix(nrow=niter, ncol=num.estimators)
  for (i in 1:niter) {
    if (! is.null(rList[[i]])) {
      for (j in 1:num.estimators) {
        beta[i, , j] <- rList[[i]]$rBeta[[j]]
        bias[i, , j] <- rList[[i]]$rBeta[[j]] - true.beta
        inCI[i, , j] <- true.beta > rList[[i]]$rCI[[j]][, 1] & true.beta < rList[[i]]$rCI[[j]][, 2]
        CI[i, , , j] <- rList[[i]]$rCI[[j]]
        maxIC[i, j] <- max(abs(colSums(rList[[i]]$rIC[[j]])))
        IC[i, , , j] <- rList[[i]]$rIC[[j]]
      }
    }
  }
  return(list(bias=bias, beta=beta, inCI=inCI, maxIC=maxIC, true.beta=true.beta, CI=CI, IC=IC)) 
}

RunIter <- function(i, sim.opts) {
  cat("i = ", i, "\n")
  try.object <- try({
    set.seed(i + 10000000 + sim.opts$base.seed)
    data <- GenerateData(sim.opts)
    est <- GetEst(data, sim.opts)
  })
  if (class(try.object) == "try-error") {
    log.file <- paste0("LongTMLE/logs/error-file-",format(Sys.time(), "%b%d%Y%H.%M.%S"))
    msg <- paste("error in RunIter, i = ", i, "\n\n try.object=\n", try.object, "log.file=", log.file, "\n")
    cat(msg)
    save(list=ls(), file=log.file)
    return(NULL)
  }
  return(est)
}

GetEst <- function(data, sim.opts) {
  l <- GetRegimens(sim.opts)  
  regimens <- l$regimens
  summary.measures <- l$summary.measures
  det.g.fun <- sim.opts$det.g.fun
  gbounds <- c(sim.opts$lower.gbound, 1)  
  nodes <- GetNodes(data)
  form <- GetForm(sim.opts=sim.opts, nodes=nodes)
  final.Ynodes <- nodes$Y

  r1 <- ltmleMSM(data, Anodes=nodes$A, Cnodes=nodes$C, Lnodes=nodes$L, Ynodes=nodes$Y, Qform=form$Q, gform=form$g, gbounds=gbounds, deterministic.g.function=det.g.fun, SL.library=NULL, regimes=regimens, working.msm=sim.opts$working.msm, summary.measures=summary.measures, final.Ynodes=final.Ynodes, pooledMSM=TRUE, weight.msm=sim.opts$use.weights, estimate.time=FALSE, mhte.iptw=TRUE, iptw.only=F, survivalOutcome=T) #pooledMSM=T
  r2 <- ltmleMSM(data, Anodes=nodes$A, Cnodes=nodes$C, Lnodes=nodes$L, Ynodes=nodes$Y, Qform=form$Q, gform=form$g, gbounds=gbounds, deterministic.g.function=det.g.fun, SL.library=NULL, regimes=regimens, working.msm=sim.opts$working.msm, summary.measures=summary.measures, final.Ynodes=final.Ynodes, pooledMSM=FALSE, weight.msm=sim.opts$use.weights, estimate.time=FALSE, mhte.iptw=TRUE, iptw.only=F,survivalOutcome=T) #pooledMSM=F
  
  if (sim.opts$standard.dynamic.iptw) {
    r3 <- StandardDynamicIPTW(data, nodes, sim.opts$working.msm, regimens, summary.measures, final.Ynodes, cum.g=r2$cum.g)
  } else {
    num.beta <- length(r1$beta)
    n <- nrow(data)
    r3 <- list(beta=rep(NA, num.beta), CI=matrix(nrow=num.beta, ncol=2), IC=matrix(nrow=n, ncol=num.beta))
  }
  
  rBeta <- list(r1$beta, r2$beta, r3$beta, r2$beta.iptw)
  rCI <- list(GetConfInt(r1), GetConfInt(r2), r3$CI, GetConfInt(r2, estimator="iptw"))
  rIC <- list(r1$IC, r2$IC, r3$IC, r2$IC.iptw)
  names(rBeta) <- names(rCI) <- names(rIC) <- c("Pooled TMLE", "Stratified TMLE", "Standard IPTW", "Stratified IPTW") 
  return(list(rBeta=rBeta, rCI=rCI, rIC=rIC))
}

Bootstrap <- function(d, id.column, num.bootstraps, sim.opts) {
  RunOne <- function(i) {
    set.seed(i + 10000000 + sim.opts$base.seed)
    sample.id <- sample(unique(d[, id.column]), replace=TRUE)
    index <- d[, id.column] %in% sample.id
    cat(" i = ", i)
    d.indexed <- d[index, ]
    sim.opts.indexed <- sim.opts
    sim.opts.indexed$n <- nrow(d.indexed)
    
    try.object <- try({
      est <- GetEst(data=d.indexed, sim.opts=sim.opts.indexed)
      outputs <- simplify2array(est$rBeta)
      save(outputs, file="LongTMLE/logs/bootstrap-int.RData")
    })
    if (class(try.object) == "try-error") {
      cat("error in RunOne:\n")
      cat(try.object)
    }
    return(outputs)
  }
  l <- mclapply(1:num.bootstraps, RunOne, mc.cores=detectCores())
  return(simplify2array(l))
}