rm(list = ls())
library("deSolve")
library("DEoptim")
library("stats")
library("stats4")
library("Hmisc")
# library("bbmle")
# library("rgl")
# library("akima")
# library("RColorBrewer")

## ----------------------------------------------------------------------------
## The model
## ----------------------------------------------------------------------------

Threestates = function(t,state,parameters) {
  with(as.list(c(state,parameters)),{
    dU = - pBeta*U*V
    dI = pBeta*U*V - pDelta*I
    dV = pP*I - pC*V
    list(c(dU,dI,dV))
  })
}

## ----------------------------------------------------------------------------
## Objective function least square
## ----------------------------------------------------------------------------

# Obj <- function(x) {
# 	pBeta 	<- x[1]
#   	pDelta 	<- x[2]
#   	pP 		<- x[3]
#   	pC 		<- x[4]
#   	paras 	<- c(pBeta = pBeta, pDelta = pDelta, pP = pP, pC = pC)
#     out     <- ode(y = state, times = times, func = Threestates, parms = paras)
#     select  <- out[round(out[,1],2) %in% round(cDa[,1],2), ]
#     RMS     <- sum( (cDa[,2] - log10(select[,4]) )^2 ) / dim(select)[2]
#   	return(RMS)
# }

## ----------------------------------------------------------------------------
## Select reference parameters
## ----------------------------------------------------------------------------
## Draw the desire data
## Extract data points
## Fit the model to the data
## Arbitrarily adjust afterward
## ----------------------------------------------------------------------------

# data.draw <- read.table("data.draw.csv", header = T, sep = ",")
# cDa <- data.draw
# fit <- DEoptim(Obj, lower, upper, DEoptim.control(trace = 1) )
# summary(fit)

times   <- seq(0, 12, by = 0.01)
state   <- c(U = 10^6, I = 0, V = 10)

lower   = log10(c(1e-7, 1e-2, 1e+0, 1e-1))
upper   = log10(c(1e-3, 1e+2, 1e+2, 1e+2))
# lower   = log10(c(1e-7, 0.001,  1, 1e-1))
# upper   = log10(c(1e-3, 10,   100, 1e+2))

paras   <- c(pBeta = 1e-05, pDelta = 1.6, pP = 5, pC = 3.7)

# The "correct" data
Data0   <- ode(state, times, Threestates, paras)
Data0   <- as.data.frame(Data0)


## ----------------------------------------------------------------------------
## Data generate
## ----------------------------------------------------------------------------
## Generate using adjusted parameter set
## Adding noise
## Randomly select time points and number of data points
## Adding noise as random, assume log-normal distribution of
## viral load measurements
## ----------------------------------------------------------------------------

ChooseData <- function(timepoints, n.datapoints, std, plot=FALSE) {
    ## ------------------------------------------------------
    ## std is the desire standard deviation
    ## set.seed(2015)
    ## extract only data in chosen time points
    ## Using as character to solve the floating point issue when using %in%
    ## ------------------------------------------------------
    sel <- Data0[which(as.character(Data0[, "time"]) %in% as.character(timepoints)), "V"]
    tmp <- NULL
    for (i in 1:length(sel)) {
        tmp[[i]] <- log10(sel[i]) + rnorm(100, 0, std) # Normal on log scale
        tmp[[i]] <- sample(tmp[[i]], n.datapoints, replace = TRUE)
    }
    # export from list to data frame and give column names
    tmp <- cbind(rep(timepoints, each = n.datapoints), unlist(tmp))
    tmp <- data.frame(tmp)
    colnames(tmp) <- c("time","V")
    if (plot) print(PlotData(tmp))
    return(logV = tmp)
}

## ----------------------------------------------------------------------------
## Function that take the data resulted from ChooseData and plot it
## ----------------------------------------------------------------------------

PlotData <- function(data, origin = TRUE) {
    plot(times, log10(Data0[, "V"]), type = "n", ylim = c(0, 8),
    ylab = expression(paste(log[10], "(Viral load)") ),
    xlab = "Day post infection (dpi)", las = 1)
    if (origin) lines(times, log10(Data0[, "V"]), lwd = 1, lty = 2, col = 1)
    points(data[,"time"], data[,"V"], col = "darkgray", pch = 16)
    abline(h = 0.5, lty = 2, col = "gray", lwd = 2)
    text(4, 0.2, "Undetectable level")
}

## ----------------------------------------------------------------------------
## RMSE Objective function for Least Square Optimization
## ----------------------------------------------------------------------------

RMSE <- function(x) {
    pBeta   <- x[1]
    pDelta  <- x[2]
    pP      <- x[3]
    pC      <- x[4]
    paras   <- c(pBeta = pBeta, pDelta = pDelta, pP = pP, pC = pC)
    out     <- ode(y = state, times = times, func = Threestates, parms = paras)
    # Select data match the time points
    sel  <- out[out[,1] %in% round(cDa[,1],2), ]
    # calculate the RMS
    tmp1    <- NULL
    ii = 1 # help to prevent time is not interger -> cannot use as index
    for ( i in unique(cDa[,1]) ) {
        tmp1[[ii]] <- (cDa[cDa$time == i, "V"] - log10(sel[sel[,1] == i, 4]))^2
        ii = ii + 1
    }
    RMSE    <- sqrt(mean(unlist(tmp1), na.rm = TRUE)) # index != time points
    return(RMSE)
}

## ----------------------------------------------------------------------------
## Weighted objective function for bootstrap
## ----------------------------------------------------------------------------

WMSE <- function(x) {
    pBeta   <- x[1]
    pDelta  <- x[2]
    pP      <- x[3]
    pC      <- x[4]
    paras   <- c(pBeta = pBeta, pDelta = pDelta, pP = pP, pC = pC)
    out     <- ode(y = state, times = times, func = Threestates, parms = paras)
    # Select data match the time points
    sel  <- out[out[,1] %in% round(cDa[,1],2), ]
    # calculate the RMS
    tmp1    <- NULL
    ii = 1 # help to prevent time is not interger -> cannot use as index
    for ( i in unique(cDa[,1]) ) {
        tmp1[[ii]] <- W * (cDa[cDa$time == i, "V"] - log10(sel[sel[,1] == i, 4]))^2 # note that i indicates a vector
        ii = ii + 1
    }
    RMSE    <- sqrt(mean(unlist(tmp1), na.rm = TRUE)) # index != time points
    return(RMSE)
}

## ----------------------------------------------------------------------------
## Objective function for Least ABSOLUTE
## ----------------------------------------------------------------------------

MAE <- function(x) {
    pBeta   <- x[1]
    pDelta  <- x[2]
    pP      <- x[3]
    pC      <- x[4]
    paras   <- c(pBeta = pBeta, pDelta = pDelta, pP = pP, pC = pC)
    out     <- ode(y = state, times = times, func = Threestates, parms = paras)
    # Select data match the time points
    sel  <- out[out[,1] %in% round(cDa[,1],2), ]
    # calculate the RMS
    tmp1    <- NULL
    ii = 1 # help to prevent time is not interger -> cannot use as index
    for ( i in unique(cDa[,1]) ) {
        tmp1[[ii]] <- abs(cDa[cDa$time == i, "V"] - log10(sel[sel[,1] == i, 4]))
        ii = ii + 1
    }
    MAE    <- mean(unlist(tmp1), na.rm = TRUE)
    return(MAE)
}

## ----------------------------------------------------------------------------
## MLE for DE
## ----------------------------------------------------------------------------

MLE <- function(x) {
    pBeta   <- x[1]
    pDelta  <- x[2]
    pP      <- x[3]
    pC      <- x[4]
    sig     <- x[5]
    paras   <- 10^c(pBeta = pBeta, pDelta = pDelta, pP = pP, pC = pC)
    # Avoid errors in computation
    tryCatch( {
        out   <- ode(state, times, Threestates, paras)
        sel   <- out[as.character(out[,1]) %in% as.character(cDa$time), ]
        # make sure viral load does not go negative
        if (sum(sel[, "V"] < 0) > 0) { 
            mll <-  1e+8
            return(mll)
        } else {
            # Number of replicates
            nm  <- dim(cDa)[1]/length(unique(cDa$time))
            # Residuals
            x   <- cDa$V - rep(log10(sel[, "V"]), each = nm)
            # Normal density, return in log prob.
            ll  <- dnorm(x, mean = 0, sd = sig, log = TRUE)
            mll <- -sum(ll)
            cat(".")
            return(mll)
        }}, error = function(e) {
            mll <- 1e+8
            return(mll)
        }
    )
}

W <- NA
DEoptions <- DEoptim.control(parallelType = 1, packages = c("deSolve"),
    parVar = c("Threestates","state", "lower", "upper", "times", "cDa","W"),
    trace = 1, itermax = 10000, steptol = 50, reltol = 1e-8, F = 0.8, CR = 0.9)


## ----------------------------------------------------------------------------
## Predictive error function loops through number of data points/time points
## ----------------------------------------------------------------------------
## Took time scheme, standard deviation as input
## Generate data
## Obtain parameter and save
## Calculate predictive error and save
## ----------------------------------------------------------------------------

W <- NA
DEoptions <- DEoptim.control(parallelType = 1, packages = c("deSolve"),
    parVar = c("Threestates","state", "lower", "upper", "times", "cDa","W"),
    trace = 1, itermax = 10000, steptol = 50, reltol = 1e-8, F = 0.8, CR = 0.9)

PER <- function(t.scheme, max.n, std, robust = FALSE) {
    # ------------------------------------------------------
    # t.scheme  : data collection scheme
    # max.n     : maximum number of data points/time points
    # std       : desired standard deviation on log10 scale
    # robust    : true will use MAE instead of RMSE
    # ------------------------------------------------------
    n.data <- Pred.Err <- Est.Par <- Pred <- NULL
    for (i in 1:max.n) {
        cDa <<- ChooseData(t.scheme, i, std)
        if(robust) {
            ppp <- DEoptim(MAE, lower, upper, DEoptions)$optim$bestmem
        } else {
            ppp <- DEoptim(RMSE, lower, upper, DEoptions)$optim$bestmem
        }
        ppp <- c(pBeta=ppp[[1]], pDelta=ppp[[2]], pP=ppp[[3]], pC=ppp[[4]])
        pre <- ode(y = state, times = times, func = Threestates, parms = ppp)
        sel <- pre[pre[,1] %in% t.scheme, 4]
        ref <- Data0[Data0[,1] %in% t.scheme,4]
        PEr <- sqrt(mean((sel - ref)^2))
        # Storing data
        n.data          <- c(n.data, i)
        Pred.Err        <- c(Pred.Err, PEr)
        Est.Par[[i]]    <- ppp
        Pred[[i]]       <- pre
    }
    return(list(Dat = data.frame(n.data, Pred.Err), Est.Par = Est.Par, Pred = Pred))
}

MatPro <- function(input.data, methods = 'RMSE', PL = FALSE, boot = FALSE, n.boot = 1000) {
    require(deSolve)
    require(DEoptim)
    ## ------------------------------------------------------------------------
    startt <- proc.time()
    cDa <<- input.data
    ## ------------------------------------------------------------------------
    message("Optimization starts")
    ## ------------------------------------------------------------------------
    method <- get(methods, mode = "function")
    DEargs <- list(method, lower, upper, DEoptions)
    if (methods == 'MLE') DEargs <- list(method, c(lower, 1e-8), c(upper, 2), DEoptions)
    fit <- do.call("DEoptim", DEargs)
    ppp <- fit$optim$bestmem
    ppp <- c(pBeta=ppp[[1]], pDelta=ppp[[2]], pP=ppp[[3]], pC=ppp[[4]])
    obj <- fit$optim$bestval
    ## ------------------------------------------------------------------------
    message("Optimization finished!", "\n", ppp, "\n")
    ## ------------------------------------------------------------------------
    ## Profile likelihood
    ## ------------------------------------------------------------------------
    if(PL) {
        message("Profiling starts", "\n")
        pro.ll <-  NULL
        for (v in 1:length(ppp)) {
            # Creating parameter sequence
            tmpl    <- seq(lower[v], ppp[[v]], length.out = 100)
            tmpl    <- tmpl[order(tmpl, decreasing = TRUE)[cumsum(1:13)]]
            tmpr    <- seq(ppp[[v]], upper[v], length.out = 100)
            tmpr    <- tmpr[cumsum(1:13)]
            par.seq <- sort(unique(c(lower[v], tmpl,ppp[[v]],tmpr, upper[v])))
            ppl     <- NULL
            # Optim and save obj and coresponding var
            for (p in par.seq) {
                DEargs  <- list(method, replace(lower, v, p), replace(upper, v, p), DEoptions)
                if (methods == 'MLE') DEargs <- list(method, c(replace(lower, v, p), 1e-8), c(replace(upper, v, p), 2), DEoptions)
                fit     <- do.call("DEoptim", DEargs)
                ppl     <- c(ppl, fit$optim$bestval)
            }
            pro.ll[[v]] <- cbind(par.seq, ppl)
        }
    } else {
        pro.ll <-  NULL
    }
    # if (boot == TRUE & method %in% "MAE")
        # stop("Cannot use weighted bootstrap with MAE")
    if (boot == TRUE) {
        message("Bootstrap starts", "\n")
        # weighted bootstrap
        Weights <- matrix(NA, dim(cDa)[1], n.boot)
        set.seed(2015)
        for (i in 1:n.boot) {
          Weights[,i] <- rexp(dim(cDa)[1])
        }
        Vmem <- matrix(NA, n.boot, 4)
        for (i in 1:n.boot) {
            W <<- Weights[,i]
            DEargs <- list(WMSE, lower, upper, DEoptions)
            fit <- do.call("DEoptim", DEargs)
            Vmem[i,1] <- fit$optim$bestmem[[1]]
            Vmem[i,2] <- fit$optim$bestmem[[2]]
            Vmem[i,3] <- fit$optim$bestmem[[3]]
            Vmem[i,4] <- fit$optim$bestmem[[4]]
        }
        Vmem <- as.data.frame(Vmem)
        colnames(Vmem) <- c("beta", "delta", "p", "c")
    } else {
        Vmem = NULL
    }
    # Save running time
    endt <- proc.time() - startt
    return(list(Est.Pars = ppp, Profile = pro.ll, Boot = Vmem, runtime = endt))
}

## ----------------------------------------------------------------------------
## MLE for optim
## ----------------------------------------------------------------------------
minusLLL <- function(pBeta, pDelta, pP, pC, mu, sig) {
    tryCatch( {
        paras <- c(pBeta = 10^(pBeta), pDelta = 10^(pDelta), pP = 10^(pP), pC = 10^(pC))
        out   <- ode(state, times, Threestates, paras)
        sel   <- out[as.character(out[,1]) %in% as.character(cDa[,1]), ]
        if (sum(sel[, "V"] < 0) > 0) { # check if viral load goes negative
            mll <-  1e+8
            return(mll)
        } else {
            # Residuals
            nm  <- dim(cDa)[1]/length(unique(cDa$time))
            x   <- cDa$V - rep(log10(sel[, "V"]), each = nm)
            ll  <- dnorm(x, mean = mu, sd = sig, log = TRUE)
            mll <- -sum(ll)
            # cat(".")
            return(mll)
        }}, error = function(e) {
            mll <- 1e+8
            return(mll)
        })
}

updateFit <- function(fitlnm) {
    # Update the fit with better starting value
    fit.old <- fitlnm
    # bestpar <- coef(fit.old)
    startv  <- coef(fit.old)
    # startv  <- list(pBeta = 10^(bestpar[[1]]), pDelta = 10^(bestpar[[2]]), pP = 10^(bestpar[[3]]), pC = 10^(bestpar[[4]]), sig = bestpar[[6]])
    BFGScontrols[[1]] <- startv
    fit.new <-  do.call("mle2", c(minusLLL, BFGScontrols))
    if ( logLik(fit.new) > logLik(fit.old) ) {
        message("Result is improved! Updated!")
        fit <- fit.new
    } else {
        message("Result is NOT improved! Keep the first estimates.")
        fit <- fit.old
    }
    return(fit = fit)
}

plotODEandData <- function(data, ODeoutput) {
    par(mfrow = c(1,2))
    plot(data[,1], 10^(data[, 2]), type= "p", main = "Raw data")
    lines(ODeoutput[,1], (ODeoutput[,4]), col = 2)
    plot(data[,1], (data[, 2]), type= "p", main = "Log10")
    lines(ODeoutput[,1], log10(ODeoutput[,4]), col = 2)
}

# mmle <- function(cDa, update=FALSE, plot=FALSE) {
mmle <- function(cDa, update=FALSE) {
    fitlnm  <-  do.call("mle2", c(minusLLL, BFGScontrols))
    if (update) {
        message("Updating parameters...")
        fitlnm <- updateFit(fitlnm)
    }
    # bestpar <- coef(fitlnm)
    # out     <- ode(y = state, times = times, func = Threestates, 
        # parms = c(pBeta = 10^(bestpar[[1]]), pDelta = 10^(bestpar[[2]]), 
        # pP = 10^(bestpar[[3]]), pC = 10^(bestpar[[4]])))
    # if (plot) 
    #     plotODEandData(cDa, out)
    # return(list(ll = logLik(fitlnm), Est = 10^bestpar[1:4], sig = bestpar[6], hessian = fitlnm@details$hessian,  Pred = out[,c(1,4)]))
    # return(list(ll = logLik(fitlnm), Est = 10^bestpar[1:4], sig = bestpar[6], hessian = fitlnm@details$hessian))
    return(fit = fitlnm)
}

BFGScontrols <- list(
    start = list(pBeta = lower[1], pDelta = lower[2], pP = lower[3], pC = lower[4], sig = 1e-8),
    # start = list(pBeta = -6, pDelta = 0, pP = 1, pC = 1, sig = 0.5),
    fixed = c(mu = 0), method = "L-BFGS-B",
    lower = c(pBeta = lower[1], pDelta = lower[2], pP = lower[3], pC = lower[4], sig = 1e-8),
    upper = c(pBeta = upper[1], pDelta = upper[2], pP = upper[3], pC = upper[4], sig = 2),
    control = list(maxit = 10000, factr = 1e+8))

Krofile <- function(input.data, fit, indexvar) {
    startt <- proc.time()
    cDa    <- input.data
    message("Profiling starts...", "\n")
    pro.ll <-  NULL
    ppp <- coef(fit)[1:4]
    pll <- logLik(fit)
    for (v in indexvar) {
        cat(v, fill = TRUE)
        # Creating parameter sequence
        tmpl    <- seq(lower[v], ppp[[v]], length.out = 100)
        tmpl    <- tmpl[order(tmpl, decreasing = TRUE)[cumsum(1:13)]]
        tmpr    <- seq(ppp[[v]], upper[v], length.out = 100)
        tmpr    <- tmpr[cumsum(1:13)]
        par.seq <- sort(c(lower[v], tmpl, ppp[[v]], tmpr, upper[v]))
        message("Parameters sequence", "\n")
        cat(par.seq, "\n")
        ppl     <- NULL
        BFGScontrols.new <- BFGScontrols
        for (p in par.seq) {
            fixv <- c(0, p)
            names(fixv) <- c("mu", names(paras[v]))
            BFGScontrols.new$fixed <- fixv
            BFGScontrols.new$start <- BFGScontrols$start[-v]
            BFGScontrols.new$lower <- BFGScontrols$lower[-v]
            BFGScontrols.new$upper <- BFGScontrols$upper[-v]
            fit.i     <- do.call("mle2", c(minusLLL, BFGScontrols.new))
            BFGScontrols.new$start <- as.list(coef(fit.i))
            ppl     <- c(ppl, logLik(fit.i))
            fit.ii     <- do.call("mle2", c(minusLLL, BFGScontrols.new))
            if (logLik(fit.ii) > logLik(fit.i)) {
                message("Result is improved! Updated!")
                ppl     <- c(ppl, logLik(fit.ii))
            } else {
                message("Result is NOT improved! Keep the first estimates.")
                ppl     <- c(ppl, logLik(fit.i))
            }
            cat(p, "\t", ppl, "\n")
        }
        pro.ll[[v]] <- cbind(par.seq, ppl)
    }
    # Save running time
    endt <- proc.time() - startt
    return(list(Est.Pars = ppp, Profile = pro.ll, runtime = endt))
}

## ----------------------------------------------------------------------------
## Time points schemes
## ----------------------------------------------------------------------------

t3  <- c( round(seq(0, 12*24, by = 3)/24, 2) ) [-1]
t6  <- c( round(seq(0, 12*24, by = 6)/24, 2) ) [-1]
t8  <- c( round(seq(0, 12*24, by = 8)/24, 2) ) [-1]
t12 <- c( round(seq(0, 12*24, by = 12)/24, 2) ) [-1]
t24 <- c( round(seq(0, 12*24, by = 24)/24, 2) ) [-1]
tn1 <- c(1, 2, 3, 5, 7, 9)
tn2 <- c(round(seq(8, 72, by = 8)/24, 2), seq(4, 12, by = 1))

# tnames <- c("t3", "t6", "t8", "t12", "t24", "tn1", "tn2")
tnames <- c("t3", "t6", "t8", "t12", "t24")
ndata  <- unique(c(seq(1, 15, by = 2), seq(15, 30, by = 5)))
# varseq <- unique(c(0.01, seq(0, 0.5, by = 0.05)[-1], seq(0.5, 1, by = 0.1)))
varseq <- unique(c(0.01, seq(0, 0.01, by = 0.05)[-1], seq(0.1, 1, by = 0.1)))
N <- length(tnames)*length(varseq)*length(ndata)

# -----------------------------------------------------------------------------
# Generate 1232 settings
# -----------------------------------------------------------------------------
# set.seed(123)
# Datas <- NULL
# ii    <- 1
# for (t in 1:length(tnames)) {
#     for (n in 1:length(ndata)) {
#         for (v in 1:length(varseq)) {
#             Datas[[ii]] <- ChooseData(get(tnames[t]), ndata[n], varseq[v])
#             Datas[[ii]]$v = v
#             Datas[[ii]]$n = n
#             Datas[[ii]]$t = t
#             ii <- ii + 1
#         }
#     }
# }

# Split the job for computers and run
# -----------------------------------

# ivec <- split(1:N, ceiling(seq_along(1:N)/142))
# computer <- c("simm-linux2", "simm-linux5", "simm-linux6", "simm-linux9", "simm-linux11", "Dell")
# ivec <- ivec[[which(computer == Sys.info()["nodename"])]]

# starT <- proc.time()
# for (i in ivec) {
#     cat(i, "\n", "---", "\n")
#     cDa         <- as.data.frame(Datas[[i]])
#     Fits[[i]]   <- mmle(cDa, update = TRUE)
# }
# (Runtime <- (proc.time() - starT) / 60) # in minutes

# dput(Fits, paste(Sys.info()["nodename"],"mle.txt", sep = "-"))

# -----------------------------------------------------------------------------
# Best case
# -----------------------------------------------------------------------------

pnames <- c("beta", "delta", "p", "c")

# # cDa <- ChooseData(t3, 30, 0.01)

# # tmp <- mmle(cDa, TRUE)
# bestpar <- coef(tmp)[1:4]
# # PLfit <- Krofile(cDa, tmp, 1:4)

# PLfit.old <- PLfit

# PLfit$Profile[[1]] <- rbind(PLfit$Profile[[1]], cbind(bestpar[1], logLik(tmp2)))
# PLfit$Profile[[2]] <- rbind(PLfit$Profile[[2]], cbind(bestpar[2], logLik(tmp2)))
# PLfit$Profile[[3]] <- rbind(PLfit$Profile[[3]], cbind(bestpar[3], logLik(tmp2)))
# PLfit$Profile[[4]] <- rbind(PLfit$Profile[[4]], cbind(bestpar[4], logLik(tmp2)))

# par(mfrow = c(1,4))
# for (i in 1:4) {
#     dff = 10
#     if (i == 2) dff = 21
#     plot(smooth.spline(PLfit$Profile[[i]], df = dff), type = "l", col = "skyblue", lwd = 3, ylab = "Log likelihood", xlab=parse(text = pnames[i]), cex.lab=1.5)
#     abline(v=log10(paras[i]), lty = 2, col = 2)
#     points(bestpar[[i]], min(PLfit$Profile[[i]]), col = "darkgray", pch = 17, cex = 1.7)
# }
# size.onefour <- dev.size()
# dev.copy(png, "bestcase.png", width = size.onefour[1], height = size.onefour[2], res = 300, units = "in")
# dev.off()

# -----------------------------------------------------------------------------
# Undetectable
# -----------------------------------------------------------------------------

# cDa <- ChooseData(t3, 30, 0.01)
# cDa2 <- cDa[10^(cDa$V) >= 50,]

# tmp2 <- mmle(cDa2, TRUE)
# bestpar2 <- coef(tmp2)[1:4]
# # PLfit2 <- Krofile(cDa2, tmp2, 1:4)
# PLfit2.old <- PLfit2

# PLfit2$Profile[[1]] <- rbind(PLfit2$Profile[[1]], cbind(bestpar2[1], logLik(tmp2)))
# PLfit2$Profile[[2]] <- rbind(PLfit2$Profile[[2]], cbind(bestpar2[2], logLik(tmp2)))
# PLfit2$Profile[[3]] <- rbind(PLfit2$Profile[[3]], cbind(bestpar2[3], logLik(tmp2)))
# PLfit2$Profile[[4]] <- rbind(PLfit2$Profile[[4]], cbind(bestpar2[4], logLik(tmp2)))

# par(mfrow = c(1,4))
# for (i in 1:4) {
#     dff = 12
#     if (i == 2) dff = 21
#     # plot(PLfit2$Profile[[i]], type = "b", col = "skyblue", lwd = 3, ylab = "Log likelihood", xlab=parse(text = pnames[i]), cex.lab=1.5)
#     plot(smooth.spline(PLfit2$Profile[[i]], df = dff), type = "l", col = "skyblue", lwd = 3, ylab = "Log likelihood", xlab=parse(text = pnames[i]), cex.lab=1.5)
#     abline(v=log10(paras[i]), lty = 2, col = 2)
#     points(bestpar2[[i]], min(PLfit2$Profile[[i]]), col = "darkgray", pch = 17, cex = 1.7)
# }
# size.onefour <- dev.size()
# dev.copy(png, "undetectable.png", width = size.onefour[1], height = size.onefour[2], res = 300, units = "in")
# dev.off()

# -----------------------------------------------------------------------------
# Viral load effect
# -----------------------------------------------------------------------------

computer <- c("simm-linux2", "simm-linux5", "simm-linux6", "simm-linux9", "simm-linux11", "Dell")

cDa <- ChooseData(t3, 30, 0.01)
cDa <- cDa[10^(cDa$V) >= 50,]

Vv <- c(5, 15, 20, 25, 40, 50)
Vv <- Vv[which(computer == Sys.info()["nodename"])]

state   <- c(U = 10^6, I = 0, V = Vv)

PLDE <- MatPro(cDa, methods = 'MLE', PL = TRUE)

dput(PLfit, paste(Sys.info()["nodename"],"DEllV.txt", sep = "-"))
