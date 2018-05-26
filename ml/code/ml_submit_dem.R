cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("lmer.R")

    train <- readRDS(paste0("../data/ml_train_", id, ".rds"))
    xmat <- as.matrix(train$x)
    zmat <- as.matrix(train$z)
    yvec <- as.numeric(train$y)
    group <- as.integer(train$group)

    res <- mixefLmer(train$y, train$x, train$z, train$group, family = "gaussian")
    fname <- paste0("/Shared/ssrivastva/dem/ml/result/lmer/ml_lmer_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("vandyk00.R")
    library(MCMCpack)

    train <- readRDS(paste0("../data/ml_train_", id, ".rds"))

    group <- train$group
    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)

    ranefList0 <- list()
    fixefList0 <- list()
    rmatList0 <- list()
    ylist0 <- list()
    grpIdx0 <- list()
    cat("formating data ...", "\n")
    for (ii in 1:ngroup) {
        grpIdx0[[ii]] <- which(group == grpLbl[ii])
        ranefList0[[ii]] <- train$z[grpIdx0[[ii]], , drop = FALSE]
        fixefList0[[ii]] <- train$x[grpIdx0[[ii]], , drop = FALSE]
        ylist0[[ii]] <- train$y[grpIdx0[[ii]]]
        rmatList0[[ii]] <- diag(1, length(grpIdx0[[ii]]))
    }
    cat("done", "\n")

    nranef0 <- ncol(train$z); nfixef0 <- ncol(train$x)
    muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- -(2 + nfixef0); sig0 <- 0
    eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
    dmat0 <- diag(1, nranef0)
    errVar0 <- 10

    rm(train)
    gc()
    res <- mixefEm(ylist0, fixefList0, ranefList0, rmatList0, niter = 1000, dmat0, errVar0, muBeta0, sigBetaInv0, sig0, nu0, eta0, tmat0)
    fname <- paste0("/Shared/ssrivastva/dem/ml/result/em/ml_em_", id, ".rds")
    saveRDS(res, fname)
} else {
    source("lmer.R")

    train <- readRDS(paste0("../data/dem_ml_train_cv_", id, ".rds"))
    trainPart <- vector("list", 20)
    for (ll in 1:20) {
        xmat <- as.matrix(train[[ll]]$x)
        zmat <- as.matrix(train[[ll]]$z)
        yvec <- as.numeric(train[[ll]]$y)
        group <- as.integer(train[[ll]]$group)

        trainPart[[ll]] <- mixefLmer(yvec, xmat, zmat, group, family = "gaussian")
    }

    fname <- paste0("/Shared/ssrivastva/dem/ml/result/lmer/meta_lmer_", id, ".rds")
    saveRDS(trainPart, fname)
}
