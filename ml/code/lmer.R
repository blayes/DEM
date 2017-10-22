mixefLmer <- function (yvec, xmat, zmat, group, family = "gaussian") {
    library(lme4)

    if (family == "gaussian") {
        control <- lmerControl(calc.derivs = FALSE,
                               check.nobs.vs.nlev = "ignore",
                               check.nlev.gtr.1 = "ignore",
                               check.nobs.vs.nRE = "ignore",
                               check.rankX = "stop.deficient",
                               check.scaleX = "warning")
    } else {
        ## implement the glmer part later
        control <- glmerControl(calc.derivs = FALSE,
                                check.nobs.vs.nlev = "ignore",
                                check.nlev.gtr.1 = "ignore",
                                check.nobs.vs.nRE = "ignore",
                                check.rankX = "stop.deficient",
                                check.scaleX = "warning")

    }

    startTime <- proc.time()
    if (family == "gaussian") {
        res <- lmer(yvec ~ xmat - 1  + (zmat - 1 | group), control = control)
    } else {
        res <- glmer(yvec ~ xmat - 1 + (zmat - 1 | group), control = control, family = family)
    }
    endTime <- proc.time()

    ranefCov <- VarCorr(res)[["group"]]
    ranefCorr <- attr(ranefCov, "correlation")
    ranefSd <- attr(ranefCov, "stddev")
    ranef <- as.matrix(ranef(res)[["group"]])
    attr(ranefCov, "stddev") <- NULL
    attr(ranefCov, "correlation") <- NULL

    fixef <- fixef(res)
    fixefCov <- as.matrix(vcov(res))

    dispersion <- (sigma(res))^2
    ypred <- fitted.values(res)

    converged <- is.null(res@optinfo$conv$lme4)
    niter <- res@optinfo$feval

    list(fixef = fixef,
         fixefCov = fixefCov,
         ranef = ranef,
         ranefCov = ranefCov,
         ranefCorr = ranefCorr,
         ranefSd = ranefSd,
         dispersion = dispersion,
         fittedValues = ypred,
         converged = converged,
         niter = niter,
         time = endTime - startTime)
}
