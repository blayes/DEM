updateFixRanMeanVar <- function (ylist, fixefList, ranefList, rmatList, dmat, errVar, muBeta0, sigBetaInv0) {
    nsample <- length(fixefList)

    umatList <- vector("list", nsample)
    xtuxList <- vector("list", nsample)
    xtuyList <- vector("list", nsample)
    for (ii in 1:nsample) {
        umatInv <- rmatList[[ii]] + tcrossprod(ranefList[[ii]] %*% dmat, ranefList[[ii]])
        umatList[[ii]] <- chol2inv(chol(umatInv))
        xtuxList[[ii]] <- crossprod(fixefList[[ii]], umatList[[ii]] %*% fixefList[[ii]])
        xtuyList[[ii]] <- crossprod(fixefList[[ii]], umatList[[ii]] %*% ylist[[ii]])
    }
    sumXtux <- Reduce("+", xtuxList)
    sumXtuy <- Reduce("+", xtuyList)
    tmp1 <- chol2inv(chol(sigBetaInv0 + sumXtux))
    postFixVar <- errVar * tmp1
    postFixMean <- tmp1 %*% (crossprod(sigBetaInv0, muBeta0) + sumXtuy)

    postRanVarList <- vector("list", nsample)
    postRanMeanList <- vector("list", nsample)
    for (ii in 1:nsample) {
        tmp2 <- umatList[[ii]] %*% fixefList[[ii]]
        tmp3 <- crossprod(fixefList[[ii]], umatList[[ii]])
        pmat <- umatList[[ii]] - tmp2 %*%  tmp1 %*% tmp3
        postRanVarList[[ii]] <- errVar * (dmat - dmat %*% crossprod(ranefList[[ii]], pmat %*% ranefList[[ii]]) %*% dmat)
        postRanMeanList[[ii]] <- dmat %*% crossprod(ranefList[[ii]], umatList[[ii]] %*% (ylist[[ii]] - fixefList[[ii]] %*% postFixMean))
    }

    list(
        fixMean = postFixMean,
        fixCov = postFixVar,
        ranMeanList = postRanMeanList,
        ranCovList = postRanVarList,
        umatList = umatList
    )
}

updateDmat <- function (ranMeanList, ranCovList, errVar, eta0, tmat0) {
    nsample <- length(ranMeanList)
    nranef <- ncol(ranCovList[[1]])

    bbtList <- vector("list", nsample)
    for (ii in 1:nsample) {
        bbtList[[ii]] <- tcrossprod(ranMeanList[[ii]], ranMeanList[[ii]]) + ranCovList[[ii]]
    }
    sumBbt <- Reduce("+", bbtList)

    (tmat0 + sumBbt) / (errVar * (nsample + eta0 + nranef + 1))
}

updateErrVar <- function (ylist, fixefList, ranefList, rmatList, dmat, muBeta0, sigBetaInv0, sig0, nu0) {
    nsample <- length(fixefList)
    nobs <- sum(sapply(ranefList, nrow))

    umatList <- vector("list", nsample)
    xtuxList <- vector("list", nsample)
    xtuyList <- vector("list", nsample)
    for (ii in 1:nsample) {
        umatInv <- rmatList[[ii]] + tcrossprod(ranefList[[ii]] %*% dmat, ranefList[[ii]])
        umatList[[ii]] <- chol2inv(chol(umatInv))
        xtuxList[[ii]] <- crossprod(fixefList[[ii]], umatList[[ii]] %*% fixefList[[ii]])
        xtuyList[[ii]] <- crossprod(fixefList[[ii]], umatList[[ii]] %*% ylist[[ii]])
    }
    sumXtux <- Reduce("+", xtuxList)
    sumXtuy <- Reduce("+", xtuyList)
    tmp1 <- chol2inv(chol(sigBetaInv0 + sumXtux))
    postFixMean <- tmp1 %*% (crossprod(sigBetaInv0, muBeta0) + sumXtuy)

    trm1 <- 1 / (nobs + nu0 + 2)
    trm2 <- sig0 * nu0
    trm3 <- 0
    for (ii in 1:nsample) {
        rr <- ylist[[ii]] - fixefList[[ii]] %*% postFixMean
        trm3 <- trm3 + drop(crossprod(rr, umatList[[ii]] %*% rr))
    }
    trm4 <- drop(crossprod((muBeta0 - postFixMean), sigBetaInv0 %*% (muBeta0 - postFixMean)))

    list(errVar = trm1 * (trm2 + trm3 + trm4),
         fixMean = postFixMean
         )
}

logLike <- function (ylist, fixefList, ranefList, rmatList, fixMean, dmat, errVar, muBeta0, sigBetaInv0, sig0, nu0, eta0, tmat0) {
    nsample <- length(ylist)
    nfixef <- ncol(fixefList[[1]])
    nranef <- ncol(ranefList[[1]])

    umatList <- vector("list", nsample)
    umatInvList <- vector("list", nsample)
    xtuxList <- vector("list", nsample)
    quadVec <- numeric(nsample)
    for (ii in 1:nsample) {
        umatInvList[[ii]] <- rmatList[[ii]] + tcrossprod(ranefList[[ii]] %*% dmat, ranefList[[ii]])
        umatList[[ii]] <- chol2inv(chol(umatInvList[[ii]]))
        xtuxList[[ii]] <- crossprod(fixefList[[ii]], umatList[[ii]] %*% fixefList[[ii]])
        rr <- ylist[[ii]] - fixefList[[ii]] %*% fixMean
        quadVec[ii] <- drop(crossprod(rr, umatList[[ii]] %*% rr))
    }

    sumXtux <- Reduce("+", xtuxList)

    trm1 <- -0.5 * sum(log(eigen((sumXtux + sigBetaInv0) / errVar, only.values = TRUE)$values))
    trm2 <- -0.5 * Reduce("+", lapply(umatInvList, function(x) sum(log(eigen(x * errVar, only.values = TRUE)$values))))
    trm3 <- -0.5 * nfixef * log(errVar)
    trm4 <- -0.5 * (sum(quadVec) + drop(crossprod((muBeta0 - fixMean), sigBetaInv0 %*% (muBeta0 - fixMean)))) / errVar

    ptrm1 <- -0.5 * (nu0 + 2) * log(errVar) -  (nu0 * sig0) / (2 * errVar)
    ptrm2 <- -0.5 * (eta0 + nranef + 1) * sum(log(eigen(errVar * dmat, only.values = TRUE)$values))
    ptrm3 <- -0.5 * sum(chol2inv(chol(dmat)) * tmat0) / errVar

    trm1 + trm2 + trm3 + trm4 + ptrm1 + ptrm2 + ptrm3
}

mixefEm <- function (ylist, fixefList, ranefList, rmatList, niter = 1000, dmat0, errVar0, muBeta0, sigBetaInv0, sig0, nu0, eta0, tmat0) {
    library(MCMCpack)

    dmat <- dmat0
    errVar <- errVar0
    llk0 <- 1e6
    llk1 <- 1e10

    logLikVec <- numeric(niter)
    cnt <- 0

    startTime <- proc.time()
    while(abs(llk1 - llk0) > 1e-7 && (cnt < niter)) {
        cnt <- cnt + 1
        if (cnt %% 10 == 0) cat("mixef iter: ", cnt, "\n")

        llk0 <- llk1
        resFixRan <- updateFixRanMeanVar(ylist, fixefList, ranefList, rmatList, dmat, errVar, muBeta0, sigBetaInv0)
        dmat <- updateDmat(resFixRan$ranMeanList, resFixRan$ranCovList, errVar, eta0, tmat0)
        resErrVar <- updateErrVar(ylist, fixefList, ranefList, rmatList, dmat, muBeta0, sigBetaInv0, sig0, nu0)
        errVar <- resErrVar$errVar
        fixMean <- as.numeric(resErrVar$fixMean)
        logLikVec[cnt] <- llk1 <- logLike(ylist, fixefList, ranefList, rmatList, fixMean, dmat, errVar, muBeta0, sigBetaInv0, sig0, nu0, eta0, tmat0)
    }
    endTime <- proc.time()

    list(
        errVar = errVar,
        dmat = dmat,
        fixMean = fixMean,
        logLik = logLikVec[1:cnt],
        iter = cnt,
        time = endTime - startTime
    )
}
