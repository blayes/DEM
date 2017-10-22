## See mbest package at http://ptrckprry.com/code/ Perry (2017) in JRSS-B.
genData <- function (ngroup, nobs, nfixef, nranef) {
    library(matrixStats)
    library(Matrix)

    ## fixed effects coefficients
    fixef <- rep(c(-2, 2), length = nfixef)
    if (nranef == 3) {
        ranefCorr <- matrix(c(1, -0.4, 0.3,
                              -0.4, 1, 0.001,
                              0.3, 0.001, 1),
                            nranef, nranef)
    } else {
        ranefCorr <- as.matrix(bdiag(rep(list(matrix(c(1, -0.4, 0.3,
                                                       -0.4, 1, 0.001,
                                                       0.3, 0.001, 1),
                                                     3, 3)), 2)))
    }
    ranefCov <- outer(sqrt(1:nranef), sqrt(1:nranef)) * ranefCorr
    ranefCovSqrt <- chol(ranefCov)

                                        # generate coefficients
    u <- matrix(rnorm(ngroup * nranef), ngroup, nranef)
    ranef <- u %*% ranefCovSqrt

    ## generate group
    suppressWarnings({ # ignore warning about using Walker's alias method
        group <- sample.int(ngroup, nobs, replace=TRUE)
    })

    ## generate feature  matrices with Pr(x[i,j] = +1) = P(x[i,j] = -1) = 1/2,
    x <- matrix(sample(c(-1, +1), nobs * nfixef, replace=TRUE), nobs, nfixef)
    z <- matrix(sample(c(-1, +1), nobs * nranef, replace=TRUE), nobs, nranef)

    ## compute linear predictors and generate observations
    mu <- drop(x %*% fixef) + rowSums(z * ranef[group, ])
    y <- rnorm(nobs, mean=mu, sd=1)

    list(ngroup = ngroup, nobs = nobs,
         fixef = fixef,
         ranefCov = ranefCov,
         ranefCovSqrt = ranefCovSqrt,
         group = group, x = x, z = z, y = y)
}
