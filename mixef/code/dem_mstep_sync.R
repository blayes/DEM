masterUpdateFixPost <- function(ematList, fvecList, errVar, muBeta0, sigBetaInv0) {
    emat <- Reduce("+", ematList)
    fvec <- Reduce("+", fvecList)

    tmp1 <- chol2inv(chol(emat + sigBetaInv0))

    fixMean <- drop(tmp1 %*% (fvec + sigBetaInv0 %*% muBeta0))
    fixCov <- errVar * tmp1

    list(
        fixMean = fixMean,
        fixCov = fixCov,
        fixCovDivErr = tmp1
    )
}

masterUpdateDmat <- function(bmatList, errVar, nsample, eta0, tmat0) {
    nranef <- nrow(bmatList[[1]])

    sumBmat <- Reduce("+", bmatList)

    (sumBmat + tmat0) / (errVar * (nsample + eta0 + nranef + 1))
}

masterUpdateErrVar <- function(quadVec, fixMean, nobs, muBeta0, sigBetaInv0, sig0, nu0) {

    trm1 <- 1 / (nobs + nu0 + 2)
    trm2 <- sig0 * nu0
    trm3 <- sum(quadVec)
    trm4 <- drop(crossprod((muBeta0 - fixMean), sigBetaInv0 %*% (muBeta0 - fixMean)))

    trm1 * (trm2 + trm3 + trm4)
}
