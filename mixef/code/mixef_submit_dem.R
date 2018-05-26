cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

nfixef <- c(10, 20)
nranef <- c(3, 6)

if (mtd == 0) {
    source("simulate_data.R")
    cvs <- rep(1:10, each = 4)
    fixs <- rep(rep(1:2, 2), times = 10)
    rans <- rep(rep(1:2, each = 2), times = 10)

    cid <- cvs[id]
    fid <- fixs[id]
    rid <- rans[id]

    ngroup <- 1e5
    nobs <- 1e7
    nms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e5.rds")

    set.seed(12345)
    repData <- genData(ngroup, nobs, nfixef[fid], nranef[rid])
    fname <- paste0("/Shared/ssrivastva/dem/mixef/data/full/", nms)
    saveRDS(repData, fname)

    ngroup <- 1e4
    nobs <- 1e6
    nms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e4.rds")

    set.seed(12345)
    repData <- genData(ngroup, nobs, nfixef[fid], nranef[rid])
    fname <- paste0("/Shared/ssrivastva/dem/mixef/data/full/", nms)
    saveRDS(repData, fname)

} else if (mtd == 1) {
    source("lmer.R")

    cvs <- rep(1:10, each = 4)
    fixs <- rep(rep(1:2, 2), times = 10)
    rans <- rep(rep(1:2, each = 2), times = 10)
    tmp <- cbind(cv = cvs, fix = fixs, ran = rans)
    wid <- cbind(rbind(tmp, tmp), nn = rep(1:2, each = 40))

    cid <- wid[id, 1]
    fid <- wid[id, 2]
    rid <- wid[id, 3]
    nid <- wid[id, 4]

    if (nid == 1) {
        nms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e5.rds")
        train <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/full/", nms))

        res <- mixefLmer(train$y, train$x, train$z, train$group, family = "gaussian")
        fname <- paste0("/Shared/ssrivastva/dem/mixef/result/lmer/lmer_res_", nms)
        saveRDS(res, fname)
    } else {
        nms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e4.rds")
        train <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/full/", nms))

        res <- mixefLmer(train$y, train$x, train$z, train$group, family = "gaussian")
        fname <- paste0("/Shared/ssrivastva/dem/mixef/result/lmer/lmer_res_", nms)
        saveRDS(res, fname)
    }
} else if (mtd == 2) {
    source("vandyk00.R")
    library(MCMCpack)

    cvs <- rep(1:10, each = 4)
    fixs <- rep(rep(1:2, 2), times = 10)
    rans <- rep(rep(1:2, each = 2), times = 10)
    tmp <- cbind(cv = cvs, fix = fixs, ran = rans)
    wid <- cbind(rbind(tmp, tmp), nn = rep(1:2, each = 40))

    cid <- wid[id, 1]
    fid <- wid[id, 2]
    rid <- wid[id, 3]
    nid <- wid[id, 4]

    if (nid == 1) {
        nms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e5.rds")
        train <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/full/", nms))

        group <- train$group
        grpLbl <- sort(unique(group))
        ngroup <- length(grpLbl)

        ranefList0 <- list()
        fixefList0 <- list()
        rmatList0 <- list()
        ylist0 <- list()
        grpIdx0 <- list()
        for (ii in 1:ngroup) {
            grpIdx0[[ii]] <- which(group == grpLbl[ii])
            ranefList0[[ii]] <- train$z[grpIdx0[[ii]], , drop = FALSE]
            fixefList0[[ii]] <- train$x[grpIdx0[[ii]], , drop = FALSE]
            ylist0[[ii]] <- train$y[grpIdx0[[ii]]]
            rmatList0[[ii]] <- diag(1, length(grpIdx0[[ii]]))
        }

        nranef0 <- nranef[rid]; nfixef0 <- nfixef[fid]
        muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- -(2 + nfixef0); sig0 <- 0
        eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
        dmat0 <- diag(1, nranef0)
        errVar0 <- 10

        res <- mixefEm(ylist0, fixefList0, ranefList0, rmatList0, niter = 1000, dmat0, errVar0, muBeta0, sigBetaInv0, sig0, nu0, eta0, tmat0)
        fname <- paste0("/Shared/ssrivastva/dem/mixef/result/em/em_res_", nms)
        saveRDS(res, fname)
    } else {
        nms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e4.rds")
        train <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/full/", nms))

        group <- train$group
        grpLbl <- sort(unique(group))
        ngroup <- length(grpLbl)

        ranefList0 <- list()
        fixefList0 <- list()
        rmatList0 <- list()
        ylist0 <- list()
        grpIdx0 <- list()
        for (ii in 1:ngroup) {
            grpIdx0[[ii]] <- which(group == grpLbl[ii])
            ranefList0[[ii]] <- train$z[grpIdx0[[ii]], , drop = FALSE]
            fixefList0[[ii]] <- train$x[grpIdx0[[ii]], , drop = FALSE]
            ylist0[[ii]] <- train$y[grpIdx0[[ii]]]
            rmatList0[[ii]] <- diag(1, length(grpIdx0[[ii]]))
        }

        nranef0 <- nranef[rid]; nfixef0 <- nfixef[fid]
        muBeta0 <- rep(0, nfixef0); sigBetaInv0 <- diag(0, nfixef0); nu0 <- -(2 + nfixef0); sig0 <- 0
        eta0 <- -(nranef0 + 1); tmat0 <- diag(0, nranef0);
        dmat0 <- diag(1, nranef0)
        errVar0 <- 10

        res <- mixefEm(ylist0, fixefList0, ranefList0, rmatList0, niter = 1000, dmat0, errVar0, muBeta0, sigBetaInv0, sig0, nu0, eta0, tmat0)
        fname <- paste0("/Shared/ssrivastva/dem/mixef/result/em/em_res_", nms)
        saveRDS(res, fname)
    }
} else {
    source("lmer.R")

    cvs <- rep(1:10, each = 4)
    fixs <- rep(rep(1:2, 2), times = 10)
    rans <- rep(rep(1:2, each = 2), times = 10)
    tmp <- cbind(cv = cvs, fix = fixs, ran = rans)
    wid <- cbind(rbind(tmp, tmp), nn = rep(1:2, each = 40))

    cid <- wid[id, 1]
    fid <- wid[id, 2]
    rid <- wid[id, 3]
    nid <- wid[id, 4]

    for (nw in c(10, 20)) {
        trainPart <- vector("list", nw)
        if (nid == 1) {
            fnms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e5.rds")
            for (ll in 1:nw) {
                cat("sub: ", ll, "\n")
                nms <- paste0("dem_mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e5_k_", nw, "_part_", ll, ".rds")
                train <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/sub/", nms))
                xmat <- as.matrix(train$x)
                zmat <- as.matrix(train$z)
                yvec <- as.numeric(train$y)
                group <- as.integer(train$group)
                trainPart[[ll]] <- mixefLmer(yvec, xmat, zmat, group, family = "gaussian")
            }
            nms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e5_nsub_", nw, ".rds")
            fname <- paste0("/Shared/ssrivastva/dem/mixef/result/lmer/meta_lmer_res_", nms)
            saveRDS(trainPart, fname)
        } else {
            fnms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e4.rds")
            for (ll in 1:nw) {
                cat("sub: ", ll, "\n")
                nms <- paste0("dem_mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e4_k_", nw,"_part_", ll, ".rds")
                train <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/sub/", nms))
                xmat <- as.matrix(train$x)
                zmat <- as.matrix(train$z)
                yvec <- as.numeric(train$y)
                group <- as.integer(train$group)
                trainPart[[ll]] <- mixefLmer(yvec, xmat, zmat, group, family = "gaussian")
            }
            nms <- paste0("mix_cv_", cid, "_p_", fid, "_q_", rid, "_n_1e4_nsub_", nw, ".rds")
            fname <- paste0("/Shared/ssrivastva/dem/mixef/result/lmer/meta_lmer_res_", nms)
            saveRDS(trainPart, fname)
        }
    }
}
