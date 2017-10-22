## Information from the submitted job
cmdArgs <- commandArgs(trailingOnly = TRUE)

nw <- as.integer(cmdArgs[1])
fracid <- as.integer(cmdArgs[2])
id <- as.integer(cmdArgs[3])

cvs <- rep(1:10, each = 4)
fixs <- rep(rep(1:2, 2), times = 10)
rans <- rep(rep(1:2, each = 2), times = 10)
tmp <- cbind(cv = cvs, fix = fixs, ran = rans)
wid <- cbind(rbind(tmp, tmp), nn = rep(1:2, each = 40))

cid <- wid[id, 1]
fid <- wid[id, 2]
rid <- wid[id, 3]
nid <- wid[id, 4]

trainPart <- vector("list", nw)
if (nid == 2) {
    fnms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e4.rds")
    for (ll in 1:nw) {
        nms <- paste0("dem_mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e4_k_", nw, "_part_", ll, ".rds")
        train <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/sub/", nms))
        group <- train$group
        grpLbl <- sort(unique(group))
        ngroup <- length(grpLbl)
        ranefList0 <- list()
        fixefList0 <- list()
        rmatList0 <- list()
        ylist0 <- list()
        grpIdx0 <- list()
        for (gg in 1:ngroup) {
            grpIdx0[[gg]] <- which(group == grpLbl[gg])
            ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
            fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
            ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
            rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
        }
        trainPart[[ll]] <- list(z =  ranefList0,
                                x = fixefList0,
                                y = ylist0,
                                r = rmatList0)
    }
} else {
    fnms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e5.rds")
    for (ll in 1:nw) {
        nms <- paste0("dem_mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e5_k_", nw,"_part_", ll, ".rds")
        train <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/sub/", nms))

        group <- train$group
        grpLbl <- sort(unique(group))
        ngroup <- length(grpLbl)
        ranefList0 <- list()
        fixefList0 <- list()
        rmatList0 <- list()
        ylist0 <- list()
        grpIdx0 <- list()
        for (gg in 1:ngroup) {
            grpIdx0[[gg]] <- which(group == grpLbl[gg])
            ranefList0[[gg]] <- train$z[grpIdx0[[gg]], , drop = FALSE]
            fixefList0[[gg]] <- train$x[grpIdx0[[gg]], , drop = FALSE]
            ylist0[[gg]] <- train$y[grpIdx0[[gg]]]
            rmatList0[[gg]] <- diag(1, length(grpIdx0[[gg]]))
        }
        trainPart[[ll]] <- list(z =  ranefList0,
                                x = fixefList0,
                                y = ylist0,
                                r = rmatList0)
    }
}

## no. of individuals
nsample0 <- sum(sapply(sapply(lapply(trainPart, function(X) lapply(X$x, function(Y) nrow(Y))), unlist), length))
## no. of observations
nobs0 <- sum(sapply(sapply(lapply(trainPart, function(X) lapply(X$x, function(Y) nrow(Y))), unlist), sum))
cat("( m,  n ): (", nsample0, ", ", nobs0, ")\n")

library(Rmpi)
# mpi.spawn.Rslaves(nslaves = nw) has been already executed in .Rprofile
cat("Workers reporting:\n")
mpi.remote.exec(mpi.comm.rank())

source("~/dem/mixef/code/dem_estep_sync.R")
mpi.remote.exec(rm(list = ls()))
mpi.remote.exec(ls())
mpi.bcast.Robj2slave(recvData)
mpi.bcast.Robj2slave(distEstep)
mpi.remote.exec(ls())

mpi.bcast.cmd(recvData())

for (ww in 1:nw) {
    ranefListPart <- trainPart[[ww]]$z
    fixefListPart <- trainPart[[ww]]$x
    rmatListPart <- trainPart[[ww]]$r
    ylistPart <- trainPart[[ww]]$y

    mpi.send.Robj(list(ranefListPart = ranefListPart), ww, 1)
    mpi.send.Robj(list(fixefListPart = fixefListPart), ww, 2)
    mpi.send.Robj(list(rmatListPart = rmatListPart), ww, 3)
    mpi.send.Robj(list(ylistPart = ylistPart), ww, 4)

    cat("send data to worker: ", ww, "\n")
}

cat("snapshot of workers GlobalEnv:\n")
mpi.remote.exec(head(ranefListPart[[1]])[ , 1:5])
mpi.remote.exec(head(fixefListPart[[1]])[ , 1:5])
mpi.remote.exec(head(rmatListPart[[1]])[ , 1:5])
mpi.remote.exec(sum(sapply(ylistPart, sum)))

source("~/dem/mixef/code/dem_mstep_sync.R")

## MPI setup
nworkers <- nw # same as mpi.comm.size() - 1
workerTasks <- rep(1, nworkers)
returnWorkers <- 0
closedWorkers <- 0

## EM specs
niter <- 1000
emats <- vector("list", nworkers)
bmats <- vector("list", nworkers)
fvecs <- vector("list", nworkers)
quads <- numeric(nworkers)
logLiks <- numeric(nworkers)
logLikVec <- numeric(niter)
logLikVecDist <- numeric(niter)

frac <- fracid / 100
nactv <- floor(frac * nworkers)

## EM parameter initialization
library(MCMCpack)

fixefs <- c(10, 20)
ranefs <- c(3, 6)

nranef0 <- ranefs[rid]; nfixef0 <- fixefs[fid]
muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- -(2 + nfixef0); sig0 <- 0
eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
dmat0 <- riwish(2 * nranef0, diag(1, nranef0))
errVar0 <- runif(1, 1, 10)
emat0 <- diag(1, nfixef0)

parEst <- list(dmat = dmat0, errVar = errVar0, fixCovDivErr = emat0, fixMean = muBeta0)
workerTrack <- matrix(0, nworkers, niter) ## track

## Call the distEstep on all the workers to get them ready to undertake tasks
mpi.bcast.cmd(distEstep())

## record keeping
## 1. workerTask[ww] = 0 if work has been assigned to worker 'ww'
##                   = 1 if worker 'ww' is waiting/idle
## 2. workerTrack[ww, iter] = 0 if worker 'ww' didn't return
##                              sufficient statistics at iteration 'iter'
##                          = 1 if worker 'ww' returned sufficient
##                              statistics at iteration 'iter'

## run some initializations
for (ww in 1:nworkers) {
    msg <- mpi.recv.Robj(ww, 0)
    cat("initial msg recvd from worker: ", ww, "\n")
    mpi.send.Robj(parEst, ww, 1)
    workerTasks[ww] <- 0 # task has been assigned
    cat("DEM setp started at worker: ", ww, "\n")
}

for (ww in 1:nworkers) {
    msg <- mpi.recv.Robj(ww, 1)
    workerTasks[ww] <- 1 # worker is idle
    ## sufficient stats for updating beta estimate
    emats[[ww]] <- msg$emat
    fvecs[[ww]] <- msg$fvec
    cat("recvd beta suff. stat. from worker: ", ww, "\n")
}

if (all(sum(sapply(emats, function(x) !(is.null(x)))))) {
    fixPost <- masterUpdateFixPost(emats, fvecs, parEst$errVar, muBeta0, sigBetaInv0)
    parEst$fixCovDivErr <- fixPost$fixCovDivErr
    parEst$fixMean <- fixPost$fixMean
    for (ww in 1:nworkers) {
        mpi.send.Robj(parEst, ww, 2)
        workerTasks[ww] <- 0 # task has been assigned
    }
    cat("sent parameter estimates to workers \n")
} else {
    stop("problem in receiving tag=1 \n")
}

for (ww in 1:nworkers) {
    msg <- mpi.recv.Robj(ww, 2)
    workerTasks[ww] <- 1 # worker is idle
    bmats[[ww]] <- msg$bmat
    cat("recvd covariance suff. stat. from worker: ", ww, "\n")
}

if (all(sum(sapply(bmats, function(x) !(is.null(x)))))) {
    parEst$dmat <- masterUpdateDmat(bmats, parEst$errVar, nsample0, eta0, tmat0)
    for (ww in 1:nworkers) {
        mpi.send.Robj(parEst, ww, 1)
        workerTasks[ww] <- 0 # task has been assigned
    }
    cat("sent updated Dmat to workers \n")
} else {
    stop("problem in receiving tag=2 \n")
}

for (ww in 1:nworkers) {
    msg <- mpi.recv.Robj(ww, 1)
    workerTasks[ww] <- 1 # worker is idle
    ## sufficient stats for updating beta estimate
    emats[[ww]] <- msg$emat
    fvecs[[ww]] <- msg$fvec
    cat("recvd beta suff. stat. from worker to update err. var.: ", ww, "\n")
}

if (all(sum(sapply(emats, function(x) !(is.null(x)))))) {
    fixPost <- masterUpdateFixPost(emats, fvecs, parEst$errVar, muBeta0, sigBetaInv0)
    parEst$fixCovDivErr <- fixPost$fixCovDivErr
    parEst$fixMean <- fixPost$fixMean
    for (ww in 1:nworkers) {
        mpi.send.Robj(parEst, ww, 3)
        workerTasks[ww] <- 0 # task has been assigned
    }
    cat("sent parameter estimates to workers to finish an EM iteration \n")
} else {
    stop("problem in receiving tag=1, before updating error variance \n")
}

for (ww in 1:nworkers) {
    msg <- mpi.recv.Robj(ww, 3)
    workerTasks[ww] <- 1 # worker is idle
    ## sufficient stats for updating err estimate
    quads[ww] <- msg$quad
}

## update error variance and finish first round of EM
parEst$errVar <- masterUpdateErrVar(quads, parEst$fixMean, nobs0, muBeta0, sigBetaInv0, sig0, nu0)
## update fixed effects posterior covariance matrix to account for the
## correct error variance
fixPost$fixCov <- parEst$errVar * parEst$fixCovDivErr

## send the current parameter estimates to all the workers for
## estimating the log likelihood at the end of every iteration
for (ww in 1:nworkers) {
    mpi.send.Robj(parEst, ww, 4)
    workerTasks[ww] <- 0 # task has been assigned
}

for (ww in 1:nworkers) {
    msg <- mpi.recv.Robj(ww, 4)
    workerTasks[ww] <- 1 # worker is idle
    ## log likelihood contribution of worker "ww"
    logLiks[ww] <- msg$logLik
}

## begin DEM iterations now!
## rcrdDmat <- vector("list", niter)
## rcrdErrVar <- numeric(niter)
## rcrdFixPost <- vector("list", niter)
llk0 <- 1e7
startTime <- proc.time()
for (its in 1:niter) {
    if (its %% 5 == 0) cat("DEM iteration: ", its, "\n")

    ## RANDOM EFFECTS COVARIANCE
    ## 1. send to all workers
    for (ww in 1:nworkers) {
        mpi.send.Robj(parEst, ww, 2)
        workerTasks[ww] <- 0 # task has been assigned
    }
    ## active & inactive workers
    actv <- sort(sample(1:nworkers, nactv, replace = FALSE))
    inactv <- setdiff(1:nworkers, actv)
    workerTrack[actv, its] <- 1
    ## 2. recv from active set and ignore the inactive set
    for (ww in actv) {
        msg <- mpi.recv.Robj(ww, 2)
        workerTasks[ww] <- 1 # worker is idle
        bmats[[ww]] <- msg$bmat
    }
    for (ww in inactv) {
        ## ignore message
        msg <- mpi.recv.Robj(ww, 2)
    }
    ## update
    parEst$dmat <- masterUpdateDmat(bmats, parEst$errVar, nsample0, eta0, tmat0)
    ## rcrdDmat[[its]] <- parEst$dmat

    ## FIXED EFFECTS MEAN
    ## 1. send to all workers
    for (ww in 1:nworkers) {
        mpi.send.Robj(parEst, ww, 1)
        workerTasks[ww] <- 0 # task has been assigned
    }
    ## 2. recv from active set and ignore the inactive set
    for (ww in actv) {
        msg <- mpi.recv.Robj(ww, 1)
        workerTasks[ww] <- 1 # worker is idle
        ## sufficient stats for updating beta estimate
        emats[[ww]] <- msg$emat
        fvecs[[ww]] <- msg$fvec
    }
    for (ww in inactv) {
        ## ignore message
        msg <- mpi.recv.Robj(ww, 1)
    }
    ## update
    fixPost <- masterUpdateFixPost(emats, fvecs, parEst$errVar, muBeta0, sigBetaInv0)
    parEst$fixCovDivErr <- fixPost$fixCovDivErr
    parEst$fixMean <- fixPost$fixMean

    ## ERROR VARIANCE
    ## 1. send to all workers
    for (ww in 1:nworkers) {
        mpi.send.Robj(parEst, ww, 3)
        workerTasks[ww] <- 0 # task has been assigned
    }
    ## 2. recv from active set and ignore the inactive set
    for (ww in actv) {
        msg <- mpi.recv.Robj(ww, 3)
        workerTasks[ww] <- 1 # worker is idle
        ## sufficient stats for updating beta estimate
        quads[ww] <- msg$quad
    }
    for (ww in inactv) {
        ## ignore message
        msg <- mpi.recv.Robj(ww, 3)
    }
    ## update
    parEst$errVar <- masterUpdateErrVar(quads, parEst$fixMean, nobs0, muBeta0, sigBetaInv0, sig0, nu0)
    fixPost$fixCov <- parEst$errVar * parEst$fixCovDivErr
    ## rcrdErrVar[its] <- parEst$errVar
    ## rcrdFixPost[[its]] <- fixPost[1:2]

    ## LOG LIKELIHOOD
    ## 1. send parameter estimates to all workers
    for (ww in 1:nworkers) {
        mpi.send.Robj(parEst, ww, 4)
        workerTasks[ww] <- 0 # task has been assigned
    }
    ## 2. recv from ALL the workers
    for (ww in 1:nworkers) {
        msg <- mpi.recv.Robj(ww, 4)
        workerTasks[ww] <- 1 # worker is idle
        ## log likelihood contribution of worker "ww"
        logLiks[ww] <- msg$logLik
    }
    ## calc log-likelihood
    logLikVec[its] <- Reduce("+", logLiks)

    if (abs(llk0 - logLikVec[its]) > 1e-7) {
        llk0 <- logLikVec[its]
    } else {
        break()
    }
}
endTime <- proc.time()
demTime <- endTime - startTime

finalCnt <- its
## rcrdDmat <- rcrdDmat[1:finalCnt]
## rcrdErrVar <- rcrdErrVar[1:finalCnt]
## rcrdFixPost <- rcrdFixPost[1:finalCnt]

res <- list(
    pars = parEst,
    track = workerTrack[ , 1:finalCnt],
    logLik = logLikVec[1:finalCnt],
    niters = finalCnt,
    time = demTime
)

if (nid == 2) {
    rnms <- paste0("cv_", cid, "_p_", fid, "_q_", rid, "_n_1e4_k_", nw, "_frac_", fracid, ".rds")
} else {
    rnms <- paste0("cv_", cid, "_p_", fid, "_q_", rid, "_n_1e5_k_", nw, "_frac_", fracid, ".rds")
}

fname <- paste0("/Shared/ssrivastva/dem/mixef/result/dem/dem_res_", rnms)
saveRDS(res, fname)

for (ww in 1:nworkers) {
    cat("terminating job on worker: ", ww, "\n")
    mpi.send.Robj(0, ww, 666)
    closedWorkers <- closedWorkers + 1
}

workerSumm <- list()
for (ww in 1:nworkers) {
    workerSumm[[ww]] <- mpi.recv.Robj(ww, 666)
}

cat("Closing workers \n")
mpi.close.Rslaves()
mpi.quit()
