recvData <- function () {
    msg <- mpi.recv.Robj(0, tag = 1)
    ranefListPart <<- msg$ranefListPart; rm(msg)

    msg <- mpi.recv.Robj(0, tag = 2)
    fixefListPart <<- msg$fixefListPart; rm(msg)

    msg <- mpi.recv.Robj(0, tag = 3)
    rmatListPart <<- msg$rmatListPart; rm(msg)

    msg <- mpi.recv.Robj(0, tag = 4)
    ylistPart <<- msg$ylistPart; rm(msg)
}

distEstep <- function () {
    ## update contribution of worker towards the calculation of
    ## sufficient statistics related to the fixed effects
    distUpdateFixMeanVar <<- function (ylist, fixefList, ranefList, rmatList, dmat) {
        nsample <- length(fixefList)

        umatInvList <- vector("list", nsample)
        umatList <- vector("list", nsample)
        xtuxList <- vector("list", nsample)
        xtuyList <- vector("list", nsample)
        for (ii in 1:nsample) {
            umatInvList[[ii]] <- rmatList[[ii]] + tcrossprod(ranefList[[ii]] %*% dmat, ranefList[[ii]])
            umatList[[ii]] <- chol2inv(chol(umatInvList[[ii]]))
            xtuxList[[ii]] <- crossprod(fixefList[[ii]], umatList[[ii]] %*% fixefList[[ii]])
            xtuyList[[ii]] <- crossprod(fixefList[[ii]], umatList[[ii]] %*% ylist[[ii]])
        }
        sumXtux <- Reduce("+", xtuxList)
        sumXtuy <- Reduce("+", xtuyList)

        ## send the first two to the master and keep the last two in
        ## worker memory; they are used in the next three functions
        list (
            emat = sumXtux,
            fvec = sumXtuy,
            umatList = umatList,
            umatInvList = umatInvList
        )
    }

    ## update contribution of worker towards the calculation of sufficient
    ## statistics related to the random effects
    distUpdateRanVar <<- function (ylist, fixefList, ranefList, rmatList, dmat, errVar,
                                   fixMean, fixCovDivErr, umatList) {
        nsample <- length(fixefList)

        postRanVarList <- vector("list", nsample)
        postRanMeanList <- vector("list", nsample)
        for (ii in 1:nsample) {
            tmp2 <- umatList[[ii]] %*% fixefList[[ii]]
            tmp3 <- crossprod(fixefList[[ii]], umatList[[ii]])
            pmat <- umatList[[ii]] - tmp2 %*%  fixCovDivErr %*% tmp3
            postRanVarList[[ii]] <- errVar * (dmat - dmat %*% crossprod(ranefList[[ii]], pmat %*% ranefList[[ii]]) %*% dmat)
            postRanMeanList[[ii]] <- drop(dmat %*% crossprod(ranefList[[ii]], umatList[[ii]] %*% (ylist[[ii]] - fixefList[[ii]] %*% fixMean)))
        }

        bbtList <- vector("list", nsample)
        for (ii in 1:nsample) {
            bbtList[[ii]] <- tcrossprod(postRanMeanList[[ii]], postRanMeanList[[ii]]) + postRanVarList[[ii]]
        }
        sumBbt <- Reduce("+", bbtList)

        ## contribution towards the sufficient statistics in the
        ## estimation of covariance of random matrix of the
        ## random effects
        sumBbt
    }

    ## update sufficient statistics contribution of worker in updating
    ## the error variance
    distUpdateErrVar <<- function(ylist, fixefList, fixMean, umatList) {
        nsample <- length(fixefList)

        trm3 <- 0
        for (ii in 1:nsample) {
            rr <- ylist[[ii]] - fixefList[[ii]] %*% fixMean
            trm3 <- trm3 + drop(crossprod(rr, umatList[[ii]] %*% rr))
        }

        trm3
    }

    ## log likelihood of a single observation in a mixed-effects model
    logLike <<- function (yvec, xmat, zmat, rmat, dmat, sigmasq, fixMean) {
        mu <- drop(xmat %*% fixMean)
        Si <- (zmat %*% tcrossprod(dmat, zmat) + rmat) * sigmasq
        yvec <- yvec - mu
        ll <-  -log(det(Si)) / 2 - sum(drop(crossprod(yvec, solve(Si))) * yvec) / 2

        ll
    }

    ## calculate the log-likelihood contribution of worker to overall log-likelihood
    distLogLik <<- function (ylist, fixefList, ranefList, rmatList, dmat, errVar, fixMean) {
        nsample <- length(fixefList)
        ll <- 0
        for (ii in 1:nsample) {
            ll <- ll + logLike(ylist[[ii]], fixefList[[ii]], ranefList[[ii]],
                              rmatList[[ii]], dmat, errVar, fixMean)
        }

        ll
    }

    ## Tags of messages ***sent to*** the master:
    ##     0 = ready to start EM;
    ##     1 = sufficient statistics related to the fixed effects;
    ##     2 = sufficient statistics related to covariance matrix of
    ##         the random effects;
    ##     3 = sufficient statistics related to error variance;
    ##     4 = contribution to the log likelihood.
    ##   666 = finished!
    ## Tags of messages ***received from*** the master:
    ##     1 = covariance matrix of random effects to update mean of
    ##         fixed effects;
    ##     2 = covariance matrix of random effects, error variance,
    ##         beta, and (error-variance)^{-1} * cov-beta required to
    ##         update sufficient statistics for calculating the
    ##         covariance matrix of random effects;
    ##     3 = updated beta (depending on the covariance matrix)
    ##         required to update the error variance.
    ##     4 = parameter estimates required to the calculate the
    ##         contribution to the log likelihood.
    ##   666 = finish!

    ## initalize global variables
    cnt <<- 0 # keep track of iterations
    done <<- 0 # done with EM?
    ## sufficient statistics to be updated using the
    ## distUpdateFixMeanVar function
    umats <<- umatInvs <<- list()
    task <<- taskInfo <<- tag <<- NULL # task assigned by the master
    infoList <<- list() # keep track of tasks -- not used currently
                        # due to memory issues

    ## Signal that master that worker is ready for the job
    mpi.send.Robj(0, dest = 0, tag = 0)

    while (done == 0) {
        ## Receive a task from the master
        task <<- mpi.recv.Robj(0, tag = mpi.any.tag())
        taskInfo <<- mpi.get.sourcetag()
        tag <<- taskInfo[2]
        ## update for record keeping
        ## infoList <<- c(infoList, task)

        if (tag == 1) {
            resBeta <- distUpdateFixMeanVar(ylistPart, fixefListPart, ranefListPart, rmatListPart, task$dmat)
            ## update the global values of sufficient statistics
            umats <<- resBeta$umatList
            umatInvs <<- resBeta$umatInvList

            task <<- taskInfo <<- tag <<- NULL # clean up

            mpi.send.Robj(obj = list(emat = resBeta$emat,
                                     fvec = resBeta$fvec),
                          dest = 0, tag = 1)
        } else if (tag == 2) {
            resRan <- distUpdateRanVar(ylistPart, fixefListPart, ranefListPart, rmatListPart, task$dmat,
                                       task$errVar, task$fixMean, task$fixCovDivErr, umats)

            task <<- taskInfo <<- tag <<- NULL # clean up

            mpi.send.Robj(obj = list(bmat = resRan), dest = 0, tag = 2)
        } else if (tag == 3) {
            ## complete EM iteration by updating the sufficient
            ## statistics related to the error variance
            resErrVar <- distUpdateErrVar(ylistPart, fixefListPart, task$fixMean, umats)

            task <<- taskInfo <<- tag <<- NULL # clean up

            mpi.send.Robj(obj = list(quad = resErrVar), dest = 0, tag = 3)
        } else if (tag == 4) {
            cnt <<- cnt + 1 # done with one round of EM
            ## calculate the contribution of the worker to the log-likelihood
            resLogLik <- distLogLik(ylistPart, fixefListPart, ranefListPart, rmatListPart,
                                    task$dmat, task$errVar, task$fixMean)

            mpi.send.Robj(obj = list(logLik = resLogLik), dest = 0, tag = 4)
        } else {
            done <<- 1
            ## major clean up
            umats <<- umatInvs <<- list()
            task <<- taskInfo <<- tag <<- NULL
            mpi.send.Robj(obj = list(iters = cnt, infoList = infoList),
                          dest = 0, tag = 666)
            cnt <<- 0
            infoList <<- list()
        }
    }
}
