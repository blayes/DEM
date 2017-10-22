rm(list = ls())
setwd("~/dem/ml/code/")
fdata <- readRDS("../data/ml_dem.rds")

set.seed(12345)

ncv <- 10
cvgrps <- sample(1:ncv, length(unique(fdata$group)), replace = TRUE)
grpSplit <- split(1:nrow(fdata$x), fdata$group)

testIdx <- list()
trainIdx <- list()
train <- list()
test <- list()
for (ii in 1:ncv) {
    idx <- (cvgrps == ii)
    testIdx <- unlist(grpSplit[idx])
    test <- list(x = fdata$x[unlist(grpSplit[idx]), ],
                       z = fdata$z[unlist(grpSplit[idx]), ],
                       y = fdata$y[unlist(grpSplit[idx])],
                       group = fdata$group[unlist(grpSplit[idx])]
                       )
    trainIdx <- setdiff(1:nrow(fdata$x), unlist(grpSplit[idx]))
    train <- list(x = fdata$x[-unlist(grpSplit[idx]), ],
                        z = fdata$z[-unlist(grpSplit[idx]), ],
                        y = fdata$y[-unlist(grpSplit[idx])],
                        group = fdata$group[-unlist(grpSplit[idx])]
                        )
    saveRDS(train, paste0("../data/ml_train_", ii, ".rds"))
    saveRDS(test, paste0("../data/ml_test_", ii, ".rds"))

    saveRDS(testIdx, paste0("../data/test_idx_", ii, ".rds"))
    saveRDS(trainIdx, paste0("../data/train_idx_", ii, ".rds"))
    cat("done with ", ii, "\n")
}

rm(list = ls())

ncv <- 10
npart <- 20
parts <- vector("list", npart)
names(parts) <- paste0("k", 1:npart)

set.seed(12345)
for (cc in 1:ncv) {
    train <- readRDS(paste0("../data/ml_train_", cc, ".rds"))
    lst <- train
    grpSplit <- split(1:nrow(lst$x), lst$group)
    partsIdx <- sample(1:npart, length(grpSplit), replace = TRUE)
    for (ll in 1:npart) {
        grpIdx <- which(partsIdx == ll)
        idx <- unlist(grpSplit[grpIdx])
        parts[[ll]]$nobs <- length(idx)
        parts[[ll]]$x <- lst$x[idx, ]
        parts[[ll]]$y <- lst$y[idx]
        parts[[ll]]$z <- lst$z[idx, ]
        parts[[ll]]$group <- lst$group[idx]
        parts[[ll]]$idx <- idx
    }
    saveRDS(parts, paste0("../data/dem_ml_train_cv_", cc, ".rds"))
}
