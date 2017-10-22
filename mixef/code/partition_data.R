##  Partition the data when m = 1e4, n = 1e6
rm(list=ls())

partitionData1 <- function(dat, npart, loc) {
    grpSplit <- split(1:nrow(dat$x), dat$group)
    partsIdx <- sample(1:npart, length(grpSplit), replace = TRUE)
    partData <- list()
    for (ll in 1:npart) {
        grpIdx <- which(partsIdx == ll)
        idx <- unlist(grpSplit[grpIdx])
        partData$nobs <- length(idx)
        partData$x <- dat$x[idx, ]
        partData$y <- dat$y[idx]
        partData$z <- dat$z[idx, ]
        partData$group <- dat$group[idx]
        partData$idx <- idx
        nms <- paste0("dem_mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e4", "_k_", npart, "_part_", ll, ".rds")
        saveRDS(partData, paste0(loc, "/", nms))
    }
}

nparts <- c(10, 20)

ngroup <- 1e4
nobs <- 1e6

set.seed(12345)

for (cc in 1:10) {
    for (pp in 1:2) {
        for (qq in 1:2) {
            cat("cc ", cc, " pp ", pp, " qq ", qq, "\n")
            nms <- paste0("mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e4.rds")
            repData <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/full/", nms))
            for (kk in 1:2) {
                cat("(cc, pp, qq, kk): (", cc, ", ", pp, ", ", qq, ", ", kk, ")\n")
                partitionData1(repData, nparts[kk], "/Shared/ssrivastva/dem/mixef/data/sub")
            }
        }
    }
}

## Partition the data when m = 1e5, n = 1e7

rm(list=ls())

partitionData2 <- function(dat, npart, loc) {
    grpSplit <- split(1:nrow(dat$x), dat$group)
    partsIdx <- sample(1:npart, length(grpSplit), replace = TRUE)
    partData <- list()
    for (ll in 1:npart) {
        grpIdx <- which(partsIdx == ll)
        idx <- unlist(grpSplit[grpIdx])
        partData$nobs <- length(idx)
        partData$x <- dat$x[idx, ]
        partData$y <- dat$y[idx]
        partData$z <- dat$z[idx, ]
        partData$group <- dat$group[idx]
        partData$idx <- idx
        nms <- paste0("dem_mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_5e4", "_k_", npart, "_part_", ll, ".rds")
        saveRDS(partData, paste0(loc, "/", nms))
    }
}

nparts <- c(10, 20)

ngroup <- 1e4
nobs <- 1e6

set.seed(12345)

for (cc in 1:10) {
    for (pp in 1:2) {
        for (qq in 1:2) {
            cat("cc ", cc, " pp ", pp, " qq ", qq, "\n")
            nms <- paste0("mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_5e4.rds")
            repData <- readRDS(paste0("/Shared/ssrivastva/dem/mixef/data/full/", nms))
            for (kk in 1:2) {
                cat("(cc, pp, qq, kk): (", cc, ", ", pp, ", ", qq, ", ", kk, ")\n")
                partitionData2(repData, nparts[kk], "/Shared/ssrivastva/dem/mixef/data/sub")
            }
        }
    }
}
