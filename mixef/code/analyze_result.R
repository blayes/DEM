rm(list = ls())
library(matrixStats)
library(xtable)
library(RColorBrewer)
colors <- brewer.pal(6, "Set1")

setwd("/Shared/ssrivastva/dem/mixef/result")

fracid <- c(30, 50, 70)
lmerList <- vector("list", 10)
metaList <- vector("list", 10)
demList <- vector("list", 10)
emList <- vector("list", 10)
iemList <- vector("list", 10)
for (cc in 1:10) {
    lmerList[[cc]] <- vector("list", 2)
    metaList[[cc]] <- vector("list", 2)
    demList[[cc]] <- vector("list", 2)
    iemList[[cc]] <- vector("list", 2)
    emList[[cc]] <- vector("list", 2)
    for (pp in 1:2) {
        lmerList[[cc]][[pp]] <- vector("list", 2)
        demList[[cc]][[pp]] <- vector("list", 2)
        metaList[[cc]][[pp]] <- vector("list", 2)
        iemList[[cc]][[pp]] <- vector("list", 2)
        emList[[cc]][[pp]] <- vector("list", 2)
        for (qq in 1:2) {
            lmerList[[cc]][[pp]][[qq]] <- vector("list", 2)
            metaList[[cc]][[pp]][[qq]] <- vector("list", 2)
            demList[[cc]][[pp]][[qq]] <- vector("list", 2)
            iemList[[cc]][[pp]][[qq]] <- vector("list", 2)
            emList[[cc]][[pp]][[qq]] <- vector("list", 2)
            for (nn in 1:2) {
                if (nn == 1) {
                    lnms <- paste0("lmer/lmer_res_mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e4.rds")
                    enms <- paste0("em/em_res_mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e4.rds")
                } else {
                    lnms <- paste0("lmer/lmer_res_mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e5.rds")
                    enms <- paste0("em/em_res_mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e5.rds")
                }
                lmerList[[cc]][[pp]][[qq]][[nn]] <- readRDS(lnms)
                emList[[cc]][[pp]][[qq]][[nn]] <- readRDS(enms)
                demList[[cc]][[pp]][[qq]][[nn]] <- vector("list", 2)
                iemList[[cc]][[pp]][[qq]][[nn]] <- vector("list", 2)
                metaList[[cc]][[pp]][[qq]][[nn]] <- vector("list", 2)
                for (kk in 1:2) {
                    if (nn == 1) {
                        mnms <- paste0("lmer/meta_lmer_res_mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e4_nsub_", c(10, 20)[kk], ".rds")
                        inms <- paste0("iem/iem_res_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e4_k_", c(10, 20)[kk], "_iem.rds")
                    } else {
                        mnms <- paste0("lmer/meta_lmer_res_mix_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e5_nsub_", c(10, 20)[kk], ".rds")
                        inms <- paste0("iem/iem_res_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e5_k_", c(10, 20)[kk], "_iem.rds")
                    }
                    iemList[[cc]][[pp]][[qq]][[nn]][[kk]] <- readRDS(inms)
                    metaList[[cc]][[pp]][[qq]][[nn]][[kk]] <- readRDS(mnms)
                    demList[[cc]][[pp]][[qq]][[nn]][[kk]] <- vector("list", 3)
                    for (ff in 1:3) {
                        if (nn == 1) {
                            dnms <- paste0("dem/dem_res_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e4_k_", c(10, 20)[kk], "_frac_", fracid[ff], ".rds")
                        } else {
                            dnms <- paste0("dem/dem_res_cv_", cc, "_p_", pp, "_q_", qq, "_n_1e5_k_", c(10, 20)[kk], "_frac_", fracid[ff], ".rds")
                        }
                        demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]] <- readRDS(dnms)
                    }
                }
            }
        }
    }
}

## test the number of times the monotone ascent of the likelihood was
## violated in DEM and IEM
decList <- demList
idecList <- iemList

for (cc in 1:10) {
    for (pp in 1:2) {
        for (qq in 1:2) {
            for (nn in 1:2) {
                for (kk in 1:2) {
                    for (ff in 1:3) {
                        decList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]] <- sum(diff(demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]$logLik) < 0)
                    }
                    idecList[[cc]][[pp]][[qq]][[nn]][[kk]] <- sum(diff(iemList[[cc]][[pp]][[qq]][[nn]][[kk]]$logLik) < 0)
                }
            }
        }
    }
}

decArr <- array(NA, dim = c(2, 2, 2, 2, 3, 10), dimnames = list(c("p1", "p2"),
                                                               c("q1", "q2"),
                                                               c("n1", "n2"),
                                                               c("k1", "k2"),
                                                               c("f30", "f50", "f70"),
                                                               paste0("cv", 1:10)))

idecArr <- array(NA, dim = c(2, 2, 2, 2, 10), dimnames = list(c("p1", "p2"),
                                                               c("q1", "q2"),
                                                               c("n1", "n2"),
                                                               c("k1", "k2"),
                                                               paste0("cv", 1:10)))
for (pp in 1:2) {
    for (qq in 1:2) {
        for (nn in 1:2) {
            for (kk in 1:2) {
                for (ff in 1:3) {
                    for (cc in 1:10) {
                        decArr[pp, qq, nn, kk, ff, cc] <- decList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]
                        idecArr[pp, qq, nn, kk, cc] <- idecList[[cc]][[pp]][[qq]][[nn]][[kk]]
                    }
                }
            }
        }
    }
}

## comparison of parameter estimates, iterations, and run-time
varList <- list()
covList <- list()
fixList <- list()
errList <- list()
timeList <- list()
timeList <- list()
iterList <- list()
trackList <- list()
logLikList <- list()
for (pp in 1:2) {
    varList[[pp]] <- vector("list", 2)
    covList[[pp]] <- vector("list", 2)
    fixList[[pp]] <- vector("list", 2)
    errList[[pp]] <- vector("list", 2)
    timeList[[pp]] <- vector("list", 2)
    iterList[[pp]] <- vector("list", 2)
    trackList[[pp]] <- vector("list", 2)
    logLikList[[pp]] <- vector("list", 2)
    for (qq in 1:2) {
        varList[[pp]][[qq]] <- vector("list", 2)
        covList[[pp]][[qq]] <- vector("list", 2)
        fixList[[pp]][[qq]] <- vector("list", 2)
        errList[[pp]][[qq]] <- vector("list", 2)
        timeList[[pp]][[qq]] <- vector("list", 2)
        iterList[[pp]][[qq]] <- vector("list", 2)
        trackList[[pp]][[qq]] <- vector("list", 2)
        logLikList[[pp]][[qq]] <- vector("list", 2)
        for (nn in 1:2) {
            varList[[pp]][[qq]][[nn]] <- vector("list", 2)
            covList[[pp]][[qq]][[nn]] <- vector("list", 2)
            fixList[[pp]][[qq]][[nn]] <- vector("list", 2)
            errList[[pp]][[qq]][[nn]] <- vector("list", 2)
            timeList[[pp]][[qq]][[nn]] <- vector("list", 2)
            iterList[[pp]][[qq]][[nn]] <- vector("list", 2)
            trackList[[pp]][[qq]][[nn]] <- vector("list", 2)
            logLikList[[pp]][[qq]][[nn]] <- vector("list", 2)
            for (kk in 1:2) {
                varList[[pp]][[qq]][[nn]][[kk]] <- vector("list", 3)
                covList[[pp]][[qq]][[nn]][[kk]] <- vector("list", 3)
                fixList[[pp]][[qq]][[nn]][[kk]] <- vector("list", 3)
                errList[[pp]][[qq]][[nn]][[kk]] <- vector("list", 3)
                timeList[[pp]][[qq]][[nn]][[kk]] <- vector("list", 2)
                iterList[[pp]][[qq]][[nn]][[kk]] <- vector("list", 2)
                trackList[[pp]][[qq]][[nn]][[kk]] <- vector("list", 2)
                logLikList[[pp]][[qq]][[nn]][[kk]] <- vector("list", 2)
                for (ff in 1:3) {
                    varList[[pp]][[qq]][[nn]][[kk]][[ff]] <- vector("list", 10)
                    covList[[pp]][[qq]][[nn]][[kk]][[ff]] <- vector("list", 10)
                    fixList[[pp]][[qq]][[nn]][[kk]][[ff]] <- vector("list", 10)
                    errList[[pp]][[qq]][[nn]][[kk]][[ff]] <- vector("list", 10)
                    timeList[[pp]][[qq]][[nn]][[kk]][[ff]] <- vector("list", 10)
                    iterList[[pp]][[qq]][[nn]][[kk]][[ff]] <- vector("list", 10)
                    trackList[[pp]][[qq]][[nn]][[kk]][[ff]] <- vector("list", 10)
                    logLikList[[pp]][[qq]][[nn]][[kk]][[ff]] <- vector("list", 10)
                    for (cc in 1:10) {
                        emCov <- emList[[cc]][[pp]][[qq]][[nn]]$dmat * emList[[cc]][[pp]][[qq]][[nn]]$errVar
                        emFix <- emList[[cc]][[pp]][[qq]][[nn]]$fixMean
                        emErr <- emList[[cc]][[pp]][[qq]][[nn]]$errVar
                        demCov <- demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]$pars$dmat * demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]$pars$errVar
                        demFix <- demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]$pars$fixMean
                        demErr <- demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]$pars$errVar
                        iemCov <- iemList[[cc]][[pp]][[qq]][[nn]][[kk]]$pars$dmat * iemList[[cc]][[pp]][[qq]][[nn]][[kk]]$pars$errVar
                        iemFix <- iemList[[cc]][[pp]][[qq]][[nn]][[kk]]$pars$fixMean
                        iemErr <- iemList[[cc]][[pp]][[qq]][[nn]][[kk]]$pars$errVar

                        lmerCov <- lmerList[[cc]][[pp]][[qq]][[nn]]$ranefCov
                        lmerFix <- lmerList[[cc]][[pp]][[qq]][[nn]]$fixef
                        lmerErr <- lmerList[[cc]][[pp]][[qq]][[nn]]$dispersion

                        metaFix <- rowMeans(do.call(cbind, lapply(metaList[[cc]][[pp]][[qq]][[nn]][[kk]], function(x) x$fixef)))
                        metaCov <- matrix(rowMeans(do.call(cbind, lapply(metaList[[cc]][[pp]][[qq]][[nn]][[kk]],
                                                                         function(x) as.numeric(x$ranefCov)))),
                                          nrow(emCov), nrow(emCov))
                        metaErr <- rowMeans(do.call(cbind, lapply(metaList[[cc]][[pp]][[qq]][[nn]][[kk]], function(x) x$dispersion)))

                        fixList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]] <- rbind(lmerFix, metaFix, emFix, iemFix, demFix)
                        rownames(fixList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]]) <- c("lme4", "meta", "em", "iem", paste0("dem", fracid[ff]))
                        varList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]] <- rbind(diag(lmerCov), diag(metaCov), diag(emCov), diag(iemCov), diag(demCov))
                        rownames(varList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]]) <- c("lme4", "meta", "em", "iem", paste0("dem", fracid[ff]))
                        covList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]] <- rbind(lmerCov[lower.tri(lmerCov)],
                                                                             metaCov[lower.tri(metaCov)],
                                                                             emCov[lower.tri(emCov)],
                                                                             iemCov[lower.tri(iemCov)],
                                                                             demCov[lower.tri(demCov)])
                        rownames(covList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]]) <- c("lme4", "meta", "em", "iem", paste0("dem", fracid[ff]))
                        errList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]] <- c(lmerErr, metaErr, emErr, iemErr, demErr)
                        names(errList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]]) <- c("lme4", "meta", "em", "iem", paste0("dem", fracid[ff]))

                        timeList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]] <- c(iemList[[cc]][[pp]][[qq]][[nn]][[kk]]$time[3],
                                                                    demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]$time[3]) / emList[[cc]][[pp]][[qq]][[nn]]$time[3]
                        names(timeList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]]) <- c("iem", paste0("dem", fracid[ff]))

                        iterList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]] <- c(iemList[[cc]][[pp]][[qq]][[nn]][[kk]]$niter,
                                                                    demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]$niter) / emList[[cc]][[pp]][[qq]][[nn]]$iter
                        names(iterList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]]) <- c("iem", paste0("dem", fracid[ff]))

                        trackList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]] <- rbind(rowMeans(iemList[[cc]][[pp]][[qq]][[nn]][[kk]]$track),
                                                                         rowMeans(demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]$track))
                        rownames(trackList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]]) <- c("iem", paste0("dem", fracid[ff]))
                        logLikList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]] <- c(max(iemList[[cc]][[pp]][[qq]][[nn]][[kk]]$logLik),
                                                                            max(as.numeric(demList[[cc]][[pp]][[qq]][[nn]][[kk]][[ff]]$logLik))) / max(emList[[cc]][[pp]][[qq]][[nn]]$logLik)
                        names(logLikList[[pp]][[qq]][[nn]][[kk]][[ff]][[cc]]) <- c("iem", paste0("dem", fracid[ff]))
                    }
                }
            }
        }
    }
}

varArr <- array(NA, dim = c(4, 2, 2, 2, 2, 3, 2),
                dimnames = list(c("lmer", "meta", "iem", "dem"),
                                paste0("p", 1:2),
                                paste0("q", 1:2),
                                paste0("n", 1:2),
                                paste0("k", 1:2),
                                paste0("f", 1:3),
                                c("mean", "sd")))
fixArr <- covArr <- errArr <- varArr

for (pp in 1:2) {
    for (qq in 1:2) {
        for (nn in 1:2) {
            for (kk in 1:2) {
                for (ff in 1:3) {
                    varArr[1, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[1, ]) - sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[1, ]) - sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                    varArr[2, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[2, ]) - sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[2, ]) - sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                    varArr[3, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[4, ]) - sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[4, ]) - sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                    varArr[4, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[5, ]) - sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[5, ]) - sapply(varList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                }
            }
        }
    }
}

for (pp in 1:2) {
    for (qq in 1:2) {
        for (nn in 1:2) {
            for (kk in 1:2) {
                for (ff in 1:3) {
                    covArr[1, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[1, ]) - sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[1, ]) - sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                    covArr[2, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[2, ]) - sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[2, ]) - sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                    covArr[3, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[4, ]) - sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[4, ]) - sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                    covArr[4, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[5, ]) - sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[5, ]) - sapply(covList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                }
            }
        }
    }
}

for (pp in 1:2) {
    for (qq in 1:2) {
        for (nn in 1:2) {
            for (kk in 1:2) {
                for (ff in 1:3) {
                    fixArr[1, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[1, ]) - sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[1, ]) - sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                    fixArr[2, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[2, ]) - sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[2, ]) - sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                    fixArr[3, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[4, ]) - sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[4, ]) - sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                    fixArr[4, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                         function(x) x[5, ]) - sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]],
                                                                                                      function(x) x[3, ]))^2)),
                                                       sd(sqrt((sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[5, ]) - sapply(fixList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3, ]))^2)))
                }
            }
        }
    }
}

for (pp in 1:2) {
    for (qq in 1:2) {
        for (nn in 1:2) {
            for (kk in 1:2) {
                for (ff in 1:3) {
                    errArr[1, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[1]) - sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3]))^2)),
                                                         sd(sqrt((sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[1]) - sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3]))^2)))

                    errArr[2, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[2]) - sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3]))^2)),
                                                         sd(sqrt((sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[2]) - sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3]))^2)))

                    errArr[3, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[4]) - sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3]))^2)),
                                                       sd(sqrt((sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[4]) - sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3]))^2)))
                    errArr[4, pp, qq, nn, kk, ff, ] <- c(sqrt(mean((sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[5]) - sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3]))^2)),
                                                       sd(sqrt((sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[5]) - sapply(errList[[pp]][[qq]][[nn]][[kk]][[ff]], function(x) x[3]))^2)))
                }
            }
        }
    }
}

fixTbl <- rbind(c(fixArr[1, 1, 1, 1, "k1", 1, "mean"], fixArr[2, 1, 1, 1, "k1", 1, "mean"], fixArr[3, 1, 1, 1, "k1", 1, "mean"], fixArr[4, 1, 1, 1, "k1", , "mean"]),
                c(fixArr[1, 2, 1, 1, "k1", 1, "mean"], fixArr[2, 2, 1, 1, "k1", 1, "mean"], fixArr[3, 2, 1, 1, "k1", 1, "mean"], fixArr[4, 2, 1, 1, "k1", , "mean"]),
                c(fixArr[1, 1, 2, 1, "k1", 1, "mean"], fixArr[2, 1, 2, 1, "k1", 1, "mean"], fixArr[3, 1, 2, 1, "k1", 1, "mean"], fixArr[4, 1, 2, 1, "k1", , "mean"]),
                c(fixArr[1, 2, 2, 1, "k1", 1, "mean"], fixArr[2, 2, 2, 1, "k1", 1, "mean"], fixArr[3, 2, 2, 1, "k1", 1, "mean"], fixArr[4, 2, 2, 1, "k1", , "mean"]),
                c(fixArr[1, 1, 1, 2, "k1", 1, "mean"], fixArr[2, 1, 1, 2, "k1", 1, "mean"], fixArr[3, 1, 1, 2, "k1", 1, "mean"], fixArr[4, 1, 1, 2, "k1", , "mean"]),
                c(fixArr[1, 2, 1, 2, "k1", 1, "mean"], fixArr[2, 2, 1, 2, "k1", 1, "mean"], fixArr[3, 2, 1, 2, "k1", 1, "mean"], fixArr[4, 2, 1, 2, "k1", , "mean"]),
                c(fixArr[1, 1, 2, 2, "k1", 1, "mean"], fixArr[2, 1, 2, 2, "k1", 1, "mean"], fixArr[3, 1, 2, 2, "k1", 1, "mean"], fixArr[4, 1, 2, 2, "k1", , "mean"]),
                c(fixArr[1, 2, 2, 2, "k1", 1, "mean"], fixArr[2, 2, 2, 2, "k1", 1, "mean"], fixArr[3, 2, 2, 2, "k1", 1, "mean"], fixArr[4, 2, 2, 2, "k1", , "mean"]),
                c(fixArr[1, 1, 1, 1, "k2", 1, "mean"], fixArr[2, 1, 1, 1, "k2", 1, "mean"], fixArr[3, 1, 1, 1, "k2", 1, "mean"], fixArr[4, 1, 1, 1, "k2", , "mean"]),
                c(fixArr[1, 2, 1, 1, "k2", 1, "mean"], fixArr[2, 2, 1, 1, "k2", 1, "mean"], fixArr[3, 2, 1, 1, "k2", 1, "mean"], fixArr[4, 2, 1, 1, "k2", , "mean"]),
                c(fixArr[1, 1, 2, 1, "k2", 1, "mean"], fixArr[2, 1, 2, 1, "k2", 1, "mean"], fixArr[3, 1, 2, 1, "k2", 1, "mean"], fixArr[4, 1, 2, 1, "k2", , "mean"]),
                c(fixArr[1, 2, 2, 1, "k2", 1, "mean"], fixArr[2, 2, 2, 1, "k2", 1, "mean"], fixArr[3, 2, 2, 1, "k2", 1, "mean"], fixArr[4, 2, 2, 1, "k2", , "mean"]),
                c(fixArr[1, 1, 1, 2, "k2", 1, "mean"], fixArr[2, 1, 1, 2, "k2", 1, "mean"], fixArr[3, 1, 1, 2, "k2", 1, "mean"], fixArr[4, 1, 1, 2, "k2", , "mean"]),
                c(fixArr[1, 2, 1, 2, "k2", 1, "mean"], fixArr[2, 2, 1, 2, "k2", 1, "mean"], fixArr[3, 2, 1, 2, "k2", 1, "mean"], fixArr[4, 2, 1, 2, "k2", , "mean"]),
                c(fixArr[1, 1, 2, 2, "k2", 1, "mean"], fixArr[2, 1, 2, 2, "k2", 1, "mean"], fixArr[3, 1, 2, 2, "k2", 1, "mean"], fixArr[4, 1, 2, 2, "k2", , "mean"]),
                c(fixArr[1, 2, 2, 2, "k2", 1, "mean"], fixArr[2, 2, 2, 2, "k2", 1, "mean"], fixArr[3, 2, 2, 2, "k2", 1, "mean"], fixArr[4, 2, 2, 2, "k2", , "mean"])
                )
colnames(fixTbl) <- c("lmer", "meta", "iem", "dem30", "dem50", "dem70")
rownames(fixTbl) <- c("p1q1n1k1", "p2q1n1k1", "p1q2n1k1", "p2q2n1k1", "p1q1n2k1", "p2q1n2k1", "p1q2n2k1", "p2q2n2k1",
                      "p1q1n1k2", "p2q1n1k2", "p1q2n1k2", "p2q2n1k2", "p1q1n2k2", "p2q1n2k2", "p1q2n2k2", "p2q2n2k2"
                      )

xtable(t(format(round(fixTbl[1:8, ], 4), nsmall = 4)))
xtable(t(format(round(fixTbl[9:16, ], 4), nsmall = 4)))

varTbl <- rbind(c(varArr[1, 1, 1, 1, "k1", 1, "mean"], varArr[2, 1, 1, 1, "k1", 1, "mean"], varArr[3, 1, 1, 1, "k1", 1, "mean"], varArr[4, 1, 1, 1, "k1", , "mean"]),
                c(varArr[1, 2, 1, 1, "k1", 1, "mean"], varArr[2, 2, 1, 1, "k1", 1, "mean"], varArr[3, 2, 1, 1, "k1", 1, "mean"], varArr[4, 2, 1, 1, "k1", , "mean"]),
                c(varArr[1, 1, 2, 1, "k1", 1, "mean"], varArr[2, 1, 2, 1, "k1", 1, "mean"], varArr[3, 1, 2, 1, "k1", 1, "mean"], varArr[4, 1, 2, 1, "k1", , "mean"]),
                c(varArr[1, 2, 2, 1, "k1", 1, "mean"], varArr[2, 2, 2, 1, "k1", 1, "mean"], varArr[3, 2, 2, 1, "k1", 1, "mean"], varArr[4, 2, 2, 1, "k1", , "mean"]),
                c(varArr[1, 1, 1, 2, "k1", 1, "mean"], varArr[2, 1, 1, 2, "k1", 1, "mean"], varArr[3, 1, 1, 2, "k1", 1, "mean"], varArr[4, 1, 1, 2, "k1", , "mean"]),
                c(varArr[1, 2, 1, 2, "k1", 1, "mean"], varArr[2, 2, 1, 2, "k1", 1, "mean"], varArr[3, 2, 1, 2, "k1", 1, "mean"], varArr[4, 2, 1, 2, "k1", , "mean"]),
                c(varArr[1, 1, 2, 2, "k1", 1, "mean"], varArr[2, 1, 2, 2, "k1", 1, "mean"], varArr[3, 1, 2, 2, "k1", 1, "mean"], varArr[4, 1, 2, 2, "k1", , "mean"]),
                c(varArr[1, 2, 2, 2, "k1", 1, "mean"], varArr[2, 2, 2, 2, "k1", 1, "mean"], varArr[3, 2, 2, 2, "k1", 1, "mean"], varArr[4, 2, 2, 2, "k1", , "mean"]),
                c(varArr[1, 1, 1, 1, "k2", 1, "mean"], varArr[2, 1, 1, 1, "k2", 1, "mean"], varArr[3, 1, 1, 1, "k2", 1, "mean"], varArr[4, 1, 1, 1, "k2", , "mean"]),
                c(varArr[1, 2, 1, 1, "k2", 1, "mean"], varArr[2, 2, 1, 1, "k2", 1, "mean"], varArr[3, 2, 1, 1, "k2", 1, "mean"], varArr[4, 2, 1, 1, "k2", , "mean"]),
                c(varArr[1, 1, 2, 1, "k2", 1, "mean"], varArr[2, 1, 2, 1, "k2", 1, "mean"], varArr[3, 1, 2, 1, "k2", 1, "mean"], varArr[4, 1, 2, 1, "k2", , "mean"]),
                c(varArr[1, 2, 2, 1, "k2", 1, "mean"], varArr[2, 2, 2, 1, "k2", 1, "mean"], varArr[3, 2, 2, 1, "k2", 1, "mean"], varArr[4, 2, 2, 1, "k2", , "mean"]),
                c(varArr[1, 1, 1, 2, "k2", 1, "mean"], varArr[2, 1, 1, 2, "k2", 1, "mean"], varArr[3, 1, 1, 2, "k2", 1, "mean"], varArr[4, 1, 1, 2, "k2", , "mean"]),
                c(varArr[1, 2, 1, 2, "k2", 1, "mean"], varArr[2, 2, 1, 2, "k2", 1, "mean"], varArr[3, 2, 1, 2, "k2", 1, "mean"], varArr[4, 2, 1, 2, "k2", , "mean"]),
                c(varArr[1, 1, 2, 2, "k2", 1, "mean"], varArr[2, 1, 2, 2, "k2", 1, "mean"], varArr[3, 1, 2, 2, "k2", 1, "mean"], varArr[4, 1, 2, 2, "k2", , "mean"]),
                c(varArr[1, 2, 2, 2, "k2", 1, "mean"], varArr[2, 2, 2, 2, "k2", 1, "mean"], varArr[3, 2, 2, 2, "k2", 1, "mean"], varArr[4, 2, 2, 2, "k2", , "mean"])
                )
colnames(varTbl) <- c("lmer", "meta", "iem", "dem30", "dem50", "dem70")
rownames(varTbl) <- c("p1q1n1k1", "p2q1n1k1", "p1q2n1k1", "p2q2n1k1", "p1q1n2k1", "p2q1n2k1", "p1q2n2k1", "p2q2n2k1",
                      "p1q1n1k2", "p2q1n1k2", "p1q2n1k2", "p2q2n1k2", "p1q1n2k2", "p2q1n2k2", "p1q2n2k2", "p2q2n2k2"
                      )

xtable(t(format(round(varTbl[1:8, ], 4), nsmall = 4)))
xtable(t(format(round(varTbl[9:16, ], 4), nsmall = 4)))

covTbl <- rbind(c(covArr[1, 1, 1, 1, "k1", 1, "mean"], covArr[2, 1, 1, 1, "k1", 1, "mean"], covArr[3, 1, 1, 1, "k1", 1, "mean"], covArr[4, 1, 1, 1, "k1", , "mean"]),
                c(covArr[1, 2, 1, 1, "k1", 1, "mean"], covArr[2, 2, 1, 1, "k1", 1, "mean"], covArr[3, 2, 1, 1, "k1", 1, "mean"], covArr[4, 2, 1, 1, "k1", , "mean"]),
                c(covArr[1, 1, 2, 1, "k1", 1, "mean"], covArr[2, 1, 2, 1, "k1", 1, "mean"], covArr[3, 1, 2, 1, "k1", 1, "mean"], covArr[4, 1, 2, 1, "k1", , "mean"]),
                c(covArr[1, 2, 2, 1, "k1", 1, "mean"], covArr[2, 2, 2, 1, "k1", 1, "mean"], covArr[3, 2, 2, 1, "k1", 1, "mean"], covArr[4, 2, 2, 1, "k1", , "mean"]),
                c(covArr[1, 1, 1, 2, "k1", 1, "mean"], covArr[2, 1, 1, 2, "k1", 1, "mean"], covArr[3, 1, 1, 2, "k1", 1, "mean"], covArr[4, 1, 1, 2, "k1", , "mean"]),
                c(covArr[1, 2, 1, 2, "k1", 1, "mean"], covArr[2, 2, 1, 2, "k1", 1, "mean"], covArr[3, 2, 1, 2, "k1", 1, "mean"], covArr[4, 2, 1, 2, "k1", , "mean"]),
                c(covArr[1, 1, 2, 2, "k1", 1, "mean"], covArr[2, 1, 2, 2, "k1", 1, "mean"], covArr[3, 1, 2, 2, "k1", 1, "mean"], covArr[4, 1, 2, 2, "k1", , "mean"]),
                c(covArr[1, 2, 2, 2, "k1", 1, "mean"], covArr[2, 2, 2, 2, "k1", 1, "mean"], covArr[3, 2, 2, 2, "k1", 1, "mean"], covArr[4, 2, 2, 2, "k1", , "mean"]),
                c(covArr[1, 1, 1, 1, "k2", 1, "mean"], covArr[2, 1, 1, 1, "k2", 1, "mean"], covArr[3, 1, 1, 1, "k2", 1, "mean"], covArr[4, 1, 1, 1, "k2", , "mean"]),
                c(covArr[1, 2, 1, 1, "k2", 1, "mean"], covArr[2, 2, 1, 1, "k2", 1, "mean"], covArr[3, 2, 1, 1, "k2", 1, "mean"], covArr[4, 2, 1, 1, "k2", , "mean"]),
                c(covArr[1, 1, 2, 1, "k2", 1, "mean"], covArr[2, 1, 2, 1, "k2", 1, "mean"], covArr[3, 1, 2, 1, "k2", 1, "mean"], covArr[4, 1, 2, 1, "k2", , "mean"]),
                c(covArr[1, 2, 2, 1, "k2", 1, "mean"], covArr[2, 2, 2, 1, "k2", 1, "mean"], covArr[3, 2, 2, 1, "k2", 1, "mean"], covArr[4, 2, 2, 1, "k2", , "mean"]),
                c(covArr[1, 1, 1, 2, "k2", 1, "mean"], covArr[2, 1, 1, 2, "k2", 1, "mean"], covArr[3, 1, 1, 2, "k2", 1, "mean"], covArr[4, 1, 1, 2, "k2", , "mean"]),
                c(covArr[1, 2, 1, 2, "k2", 1, "mean"], covArr[2, 2, 1, 2, "k2", 1, "mean"], covArr[3, 2, 1, 2, "k2", 1, "mean"], covArr[4, 2, 1, 2, "k2", , "mean"]),
                c(covArr[1, 1, 2, 2, "k2", 1, "mean"], covArr[2, 1, 2, 2, "k2", 1, "mean"], covArr[3, 1, 2, 2, "k2", 1, "mean"], covArr[4, 1, 2, 2, "k2", , "mean"]),
                c(covArr[1, 2, 2, 2, "k2", 1, "mean"], covArr[2, 2, 2, 2, "k2", 1, "mean"], covArr[3, 2, 2, 2, "k2", 1, "mean"], covArr[4, 2, 2, 2, "k2", , "mean"])
                )
colnames(covTbl) <- c("lmer", "meta", "iem", "dem30", "dem50", "dem70")
rownames(covTbl) <- c("p1q1n1k1", "p2q1n1k1", "p1q2n1k1", "p2q2n1k1", "p1q1n2k1", "p2q1n2k1", "p1q2n2k1", "p2q2n2k1",
                      "p1q1n1k2", "p2q1n1k2", "p1q2n1k2", "p2q2n1k2", "p1q1n2k2", "p2q1n2k2", "p1q2n2k2", "p2q2n2k2"
                      )

xtable(t(format(round(covTbl[1:8, ], 4), nsmall = 4)))
xtable(t(format(round(covTbl[9:16, ], 4), nsmall = 4)))

errTbl <- rbind(c(errArr[1, 1, 1, 1, "k1", 1, "mean"], errArr[2, 1, 1, 1, "k1", 1, "mean"], errArr[3, 1, 1, 1, "k1", 1, "mean"], errArr[4, 1, 1, 1, "k1", , "mean"]),
                c(errArr[1, 2, 1, 1, "k1", 1, "mean"], errArr[2, 2, 1, 1, "k1", 1, "mean"], errArr[3, 2, 1, 1, "k1", 1, "mean"], errArr[4, 2, 1, 1, "k1", , "mean"]),
                c(errArr[1, 1, 2, 1, "k1", 1, "mean"], errArr[2, 1, 2, 1, "k1", 1, "mean"], errArr[3, 1, 2, 1, "k1", 1, "mean"], errArr[4, 1, 2, 1, "k1", , "mean"]),
                c(errArr[1, 2, 2, 1, "k1", 1, "mean"], errArr[2, 2, 2, 1, "k1", 1, "mean"], errArr[3, 2, 2, 1, "k1", 1, "mean"], errArr[4, 2, 2, 1, "k1", , "mean"]),
                c(errArr[1, 1, 1, 2, "k1", 1, "mean"], errArr[2, 1, 1, 2, "k1", 1, "mean"], errArr[3, 1, 1, 2, "k1", 1, "mean"], errArr[4, 1, 1, 2, "k1", , "mean"]),
                c(errArr[1, 2, 1, 2, "k1", 1, "mean"], errArr[2, 2, 1, 2, "k1", 1, "mean"], errArr[3, 2, 1, 2, "k1", 1, "mean"], errArr[4, 2, 1, 2, "k1", , "mean"]),
                c(errArr[1, 1, 2, 2, "k1", 1, "mean"], errArr[2, 1, 2, 2, "k1", 1, "mean"], errArr[3, 1, 2, 2, "k1", 1, "mean"], errArr[4, 1, 2, 2, "k1", , "mean"]),
                c(errArr[1, 2, 2, 2, "k1", 1, "mean"], errArr[2, 2, 2, 2, "k1", 1, "mean"], errArr[3, 2, 2, 2, "k1", 1, "mean"], errArr[4, 2, 2, 2, "k1", , "mean"]),
                c(errArr[1, 1, 1, 1, "k2", 1, "mean"], errArr[2, 1, 1, 1, "k2", 1, "mean"], errArr[3, 1, 1, 1, "k2", 1, "mean"], errArr[4, 1, 1, 1, "k2", , "mean"]),
                c(errArr[1, 2, 1, 1, "k2", 1, "mean"], errArr[2, 2, 1, 1, "k2", 1, "mean"], errArr[3, 2, 1, 1, "k2", 1, "mean"], errArr[4, 2, 1, 1, "k2", , "mean"]),
                c(errArr[1, 1, 2, 1, "k2", 1, "mean"], errArr[2, 1, 2, 1, "k2", 1, "mean"], errArr[3, 1, 2, 1, "k2", 1, "mean"], errArr[4, 1, 2, 1, "k2", , "mean"]),
                c(errArr[1, 2, 2, 1, "k2", 1, "mean"], errArr[2, 2, 2, 1, "k2", 1, "mean"], errArr[3, 2, 2, 1, "k2", 1, "mean"], errArr[4, 2, 2, 1, "k2", , "mean"]),
                c(errArr[1, 1, 1, 2, "k2", 1, "mean"], errArr[2, 1, 1, 2, "k2", 1, "mean"], errArr[3, 1, 1, 2, "k2", 1, "mean"], errArr[4, 1, 1, 2, "k2", , "mean"]),
                c(errArr[1, 2, 1, 2, "k2", 1, "mean"], errArr[2, 2, 1, 2, "k2", 1, "mean"], errArr[3, 2, 1, 2, "k2", 1, "mean"], errArr[4, 2, 1, 2, "k2", , "mean"]),
                c(errArr[1, 1, 2, 2, "k2", 1, "mean"], errArr[2, 1, 2, 2, "k2", 1, "mean"], errArr[3, 1, 2, 2, "k2", 1, "mean"], errArr[4, 1, 2, 2, "k2", , "mean"]),
                c(errArr[1, 2, 2, 2, "k2", 1, "mean"], errArr[2, 2, 2, 2, "k2", 1, "mean"], errArr[3, 2, 2, 2, "k2", 1, "mean"], errArr[4, 2, 2, 2, "k2", , "mean"])
                )
colnames(errTbl) <- c("lmer", "meta", "iem", "dem30", "dem50", "dem70")
rownames(errTbl) <- c("p1q1n1k1", "p2q1n1k1", "p1q2n1k1", "p2q2n1k1", "p1q1n2k1", "p2q1n2k1", "p1q2n2k1", "p2q2n2k1",
                      "p1q1n1k2", "p2q1n1k2", "p1q2n1k2", "p2q2n1k2", "p1q1n2k2", "p2q1n2k2", "p1q2n2k2", "p2q2n2k2"
                      )

xtable(t(format(round(errTbl[1:8, ], 4), nsmall = 4)))
xtable(t(format(round(errTbl[9:16, ], 4), nsmall = 4)))

resLogLik <- list()
resNiter <- list()
resTrack <- list()
resTime <- list()
for (pp in 1:2) {
    resLogLik[[pp]] <- list()
    resNiter[[pp]] <- list()
    resTrack[[pp]] <- list()
    resTime[[pp]] <- list()
    for (qq in 1:2) {
        resLogLik[[pp]][[qq]] <- list()
        resNiter[[pp]][[qq]] <- list()
        resTrack[[pp]][[qq]] <- list()
        resTime[[pp]][[qq]] <- list()
        for (nn in 1:2) {
            resLogLik[[pp]][[qq]][[nn]] <- list()
            resNiter[[pp]][[qq]][[nn]] <- list()
            resTrack[[pp]][[qq]][[nn]] <- list()
            resTime[[pp]][[qq]][[nn]] <- list()
            for (kk in 1:2) {
                resLogLik[[pp]][[qq]][[nn]][[kk]] <- cbind(sapply(logLikList[[pp]][[qq]][[nn]][[kk]][[1]], function(x) x[1]),
                                                           do.call(cbind, lapply(logLikList[[pp]][[qq]][[nn]][[kk]], function(x) sapply(x, function(y) y[2]))))
                resNiter[[pp]][[qq]][[nn]][[kk]] <- cbind(sapply(iterList[[pp]][[qq]][[nn]][[kk]][[1]], function(x) x[1]),
                                                          do.call(cbind, lapply(iterList[[pp]][[qq]][[nn]][[kk]], function(x) sapply(x, function(y) y[2]))))
                resTrack[[pp]][[qq]][[nn]][[kk]] <- cbind(as.numeric(sapply(trackList[[pp]][[qq]][[nn]][[kk]][[1]], function(x) x[1, ])),
                                                          do.call(cbind, lapply(trackList[[pp]][[qq]][[nn]][[kk]], function(x) as.numeric(sapply(x, function(y) y[2, ])))))
                resTime[[pp]][[qq]][[nn]][[kk]] <- cbind(sapply(timeList[[pp]][[qq]][[nn]][[kk]][[1]], function(x) x[1]),
                                                         do.call(cbind, lapply(timeList[[pp]][[qq]][[nn]][[kk]], function(x) sapply(x, function(y) y[2]))))
                colnames(resLogLik[[pp]][[qq]][[nn]][[kk]]) <- c("IEM", "DEM-30", "DEM-50", "DEM70")
                colnames(resNiter[[pp]][[qq]][[nn]][[kk]]) <- c("IEM", "DEM-30", "DEM-50", "DEM70")
                colnames(resTrack[[pp]][[qq]][[nn]][[kk]]) <- c("IEM", "DEM-30", "DEM-50", "DEM70")
                colnames(resTime[[pp]][[qq]][[nn]][[kk]]) <- c("IEM", "DEM-30", "DEM-50", "DEM70")
            }
        }
    }
}

arrLogLik <- array(NA, dim = c(4, 2, 2, 2, 2),
                   dimnames = list(c("iem", "dem30", "dem50", "dem70"),
                                   c("p1", "p2"),
                                   c("q1", "q2"),
                                   c("n1", "n2"),
                                   c("k1", "k2")))

arrNiter <- arrTime <- arrTrack <- arrLogLik

for (pp in 1:2) {
    for (qq in 1:2) {
        for (nn in 1:2) {
            for (kk in 1:2) {
                arrLogLik[, pp, qq, nn, kk] <- paste0(format(round(colMeans(resLogLik[[pp]][[qq]][[nn]][[kk]]), 2), nsmall = 2), " (",
                                                      format(round(colSds(resLogLik[[pp]][[qq]][[nn]][[kk]]), 2), nsmall = 2), ")")
                arrNiter[, pp, qq, nn, kk] <- paste0(format(round(colMeans(resNiter[[pp]][[qq]][[nn]][[kk]]), 2), nsmall = 2), " (",
                                                     format(round(colSds(resNiter[[pp]][[qq]][[nn]][[kk]]), 2), nsmall = 2), ")")
                arrTime[, pp, qq, nn, kk] <- paste0(format(round(colMeans(resTime[[pp]][[qq]][[nn]][[kk]]), 2), nsmall = 2), " (",
                                                    format(round(colSds(resTime[[pp]][[qq]][[nn]][[kk]]), 2), nsmall = 2), ")")
                arrTrack[, pp, qq, nn, kk] <- paste0(format(round(colMeans(resTrack[[pp]][[qq]][[nn]][[kk]]), 2), nsmall = 2), " (",
                                                       format(round(colSds(resTrack[[pp]][[qq]][[nn]][[kk]]), 2), nsmall = 2), ")")
            }
        }
    }
}

tblLogLik <- cbind(arrLogLik[ , 1, 1, 1, 1], arrLogLik[ , 2, 1, 1, 1],
                   arrLogLik[ , 1, 2, 1, 1], arrLogLik[ , 2, 2, 1, 1],
                   arrLogLik[ , 1, 1, 2, 1], arrLogLik[ , 2, 1, 2, 1],
                   arrLogLik[ , 1, 2, 2, 1], arrLogLik[ , 2, 2, 2, 1],
                   arrLogLik[ , 1, 1, 1, 2], arrLogLik[ , 2, 1, 1, 2],
                   arrLogLik[ , 1, 2, 1, 2], arrLogLik[ , 2, 2, 1, 2],
                   arrLogLik[ , 1, 1, 2, 2], arrLogLik[ , 2, 1, 2, 2],
                   arrLogLik[ , 1, 2, 2, 2], arrLogLik[ , 2, 2, 2, 2]
                   )
rownames(tblLogLik) <- c("iem", "dem30", "dem50", "dem70")
colnames(tblLogLik) <- c("p1q1n1k1", "p2q1n1k1", "p1q2n1k1", "p2q2n1k1", "p1q1n2k1", "p2q1n2k1", "p1q2n2k1", "p2q2n2k1",
                         "p1q1n1k2", "p2q1n1k2", "p1q2n1k2", "p2q2n1k2", "p1q1n2k2", "p2q1n2k2", "p1q2n2k2", "p2q2n2k2")

tblNiter <- cbind(arrNiter[ , 1, 1, 1, 1], arrNiter[ , 2, 1, 1, 1],
                   arrNiter[ , 1, 2, 1, 1], arrNiter[ , 2, 2, 1, 1],
                   arrNiter[ , 1, 1, 2, 1], arrNiter[ , 2, 1, 2, 1],
                   arrNiter[ , 1, 2, 2, 1], arrNiter[ , 2, 2, 2, 1],
                   arrNiter[ , 1, 1, 1, 2], arrNiter[ , 2, 1, 1, 2],
                   arrNiter[ , 1, 2, 1, 2], arrNiter[ , 2, 2, 1, 2],
                   arrNiter[ , 1, 1, 2, 2], arrNiter[ , 2, 1, 2, 2],
                   arrNiter[ , 1, 2, 2, 2], arrNiter[ , 2, 2, 2, 2]
                   )
rownames(tblNiter) <- c("iem", "dem30", "dem50", "dem70")
colnames(tblNiter) <- c("p1q1n1k1", "p2q1n1k1", "p1q2n1k1", "p2q2n1k1", "p1q1n2k1", "p2q1n2k1", "p1q2n2k1", "p2q2n2k1",
                         "p1q1n1k2", "p2q1n1k2", "p1q2n1k2", "p2q2n1k2", "p1q1n2k2", "p2q1n2k2", "p1q2n2k2", "p2q2n2k2")

tblTime <- cbind(arrTime[ , 1, 1, 1, 1], arrTime[ , 2, 1, 1, 1],
                   arrTime[ , 1, 2, 1, 1], arrTime[ , 2, 2, 1, 1],
                   arrTime[ , 1, 1, 2, 1], arrTime[ , 2, 1, 2, 1],
                   arrTime[ , 1, 2, 2, 1], arrTime[ , 2, 2, 2, 1],
                   arrTime[ , 1, 1, 1, 2], arrTime[ , 2, 1, 1, 2],
                   arrTime[ , 1, 2, 1, 2], arrTime[ , 2, 2, 1, 2],
                   arrTime[ , 1, 1, 2, 2], arrTime[ , 2, 1, 2, 2],
                   arrTime[ , 1, 2, 2, 2], arrTime[ , 2, 2, 2, 2]
                   )
rownames(tblTime) <- c("iem", "dem30", "dem50", "dem70")
colnames(tblTime) <- c("p1q1n1k1", "p2q1n1k1", "p1q2n1k1", "p2q2n1k1", "p1q1n2k1", "p2q1n2k1", "p1q2n2k1", "p2q2n2k1",
                         "p1q1n1k2", "p2q1n1k2", "p1q2n1k2", "p2q2n1k2", "p1q1n2k2", "p2q1n2k2", "p1q2n2k2", "p2q2n2k2")

tblTrack <- cbind(arrTrack[ , 1, 1, 1, 1], arrTrack[ , 2, 1, 1, 1],
                   arrTrack[ , 1, 2, 1, 1], arrTrack[ , 2, 2, 1, 1],
                   arrTrack[ , 1, 1, 2, 1], arrTrack[ , 2, 1, 2, 1],
                   arrTrack[ , 1, 2, 2, 1], arrTrack[ , 2, 2, 2, 1],
                   arrTrack[ , 1, 1, 1, 2], arrTrack[ , 2, 1, 1, 2],
                   arrTrack[ , 1, 2, 1, 2], arrTrack[ , 2, 2, 1, 2],
                   arrTrack[ , 1, 1, 2, 2], arrTrack[ , 2, 1, 2, 2],
                   arrTrack[ , 1, 2, 2, 2], arrTrack[ , 2, 2, 2, 2]
                   )
rownames(tblTrack) <- c("iem", "dem30", "dem50", "dem70")
colnames(tblTrack) <- c("p1q1n1k1", "p2q1n1k1", "p1q2n1k1", "p2q2n1k1", "p1q1n2k1", "p2q1n2k1", "p1q2n2k1", "p2q2n2k1",
                         "p1q1n1k2", "p2q1n1k2", "p1q2n1k2", "p2q2n1k2", "p1q1n2k2", "p2q1n2k2", "p1q2n2k2", "p2q2n2k2")

xtable(tblLogLik[ , 1:8])
xtable(tblLogLik[ , 9:16])

xtable(tblNiter[ , 1:8])
xtable(tblNiter[ , 9:16])

xtable(tblTrack[ , 1:8])
xtable(tblTrack[ , 9:16])

nfix <- c("p = 10", "p = 20")
nran <- c("q = 3", "q = 6")
npart <- c("K = 10", "K = 20")

pdf("~/dem/mixef/result/img/mixef_iters.pdf", 30, 10)
par(mfrow = c(2, 8))
par(mar = c(0, 0, 0, 0), oma = c(4, 6.5, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (nn in 1:2) {
    for (kk in 1:2) {
        for (qq in 1:2) {
            for (pp in 1:2) {
                boxplot((resNiter[[pp]][[qq]][[nn]][[kk]]), ylim = c(0, range(unlist(resNiter))[2]), ylab = NA, axes = FALSE, lwd = 2,
                        boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1, staplelwd = 3, medlwd = 4, outcex = 1.5)
                abline (h = 1, lty = "dotted", lwd = 5)
                if (nn == 1) {
                    mtext(npart[kk], side = 3, line = -3.2, cex = 2)
                    mtext(expression(m == 10^4 ~ ", " ~ n == 10^6), side = 3, line = -6, cex = 2)
                    mtext(paste0(nfix[pp], ", ", nran[qq]), side = 3, line = -8.5, cex = 2)
                } else {
                    mtext(npart[kk], side = 3, line = -3.2, cex = 2)
                    mtext(expression(m == 10^5 ~ ", " ~ n == 10^7), side = 3, line = -6, cex = 2)
                    mtext(paste0(nfix[pp], ", ", nran[qq]), side = 3, line = -8.5, cex = 2)
                }
                box(col = "grey40", lwd = 4)
                if (nn == 2) {
                    axis(side = 1, tck = -.01, labels = NA)
                    mtext(c("IEM", "30%", "50%", "70%"), at = 1:4, side = 1, line = 1.5, cex = 2)
                }
                mtext(expression("No. of distributed EM iterations"  %/%  "No. of EM iterations"), side = 2, outer = TRUE, cex = 2, line = 3.5)
                if ((nn == 1 & kk == 1 & pp == 1 & qq == 1) | (nn == 2 & kk == 1 & pp == 1 & qq == 1)) {
                    axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
                    mtext(format(axTicks(2)), at = axTicks(2), side = 2, line = 0.3, las = 1, cex = 2)
                }
            }
        }
    }
}
dev.off()

pdf("~/dem/mixef/result/img/mixef_time.pdf", 30, 10)
par(mfrow = c(2, 8))
par(mar = c(0, 0, 0, 0), oma = c(4, 6.5, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (nn in 1:2) {
    for (kk in 1:2) {
        for (qq in 1:2) {
            for (pp in 1:2) {
                if (nn == 1) {
                    boxplot((resTime[[pp]][[qq]][[nn]][[kk]]), ylim = c(0, 14), ylab = NA, axes = FALSE, lwd = 2,
                            boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1,
                            staplelwd = 3, medlwd = 4, outcex = 1.5)
                    mtext(npart[kk], side = 3, line = -3.2, cex = 2)
                    mtext(expression(m == 10^4 ~ ", " ~ n == 10^6), side = 3, line = -6, cex = 2)
                    mtext(paste0(nfix[pp], ", ", nran[qq]), side = 3, line = -8.5, cex = 2)
                } else {
                    boxplot((resTime[[pp]][[qq]][[nn]][[kk]]), ylim = c(0, 6), ylab = NA, axes = FALSE, lwd = 2,
                            boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1,
                            staplelwd = 3, medlwd = 4, outcex = 1.5)
                    mtext(npart[kk], side = 3, line = -3.2, cex = 2)
                    mtext(expression(m == 10^5 ~ ", " ~ n == 10^7), side = 3, line = -6, cex = 2)
                    mtext(paste0(nfix[pp], ", ", nran[qq]), side = 3, line = -8.5, cex = 2)
                }
                abline (h = 1, lty = "dotted", lwd = 5)
                box(col = "grey40", lwd = 4)
                if (nn == 2) {
                    axis(side = 1, tck = -.01, labels = NA)
                    mtext(c("IEM", "30%", "50%", "70%"), at = 1:4, side = 1, line = 1.5, cex = 2)
                }
                mtext(expression("Distributed EM time"  %/%  "EM time"), side = 2, outer = TRUE, cex = 2, line = 3.5)
                if ((nn == 1 & kk == 1 & pp == 1 & qq == 1) | (nn == 2 & kk == 1 & pp == 1 & qq == 1)) {
                    axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
                    mtext(format(axTicks(2)), at = axTicks(2), side = 2, line = 0.3, las = 1, cex = 2)
                }
            }
        }
    }
}
dev.off()

pdf("~/dem/mixef/result/img/mixef_track.pdf", 30, 10)
par(mfrow = c(2, 8))
par(mar = c(0, 0, 0, 0), oma = c(4, 6.5, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (nn in 1:2) {
    for (kk in 1:2) {
        for (qq in 1:2) {
            for (pp in 1:2) {
                boxplot((resTrack[[pp]][[qq]][[nn]][[kk]]), ylim = c(0, 1.15), ylab = NA, axes = FALSE, lwd = 2,
                        boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1, staplelwd = 3, medlwd = 4, outcex = 1.5)
                if (nn == 1) {
                    mtext(npart[kk], side = 3, line = -3.2, cex = 2)
                    mtext(expression(m == 10^4 ~ ", " ~ n == 10^6), side = 3, line = -6, cex = 2)
                    mtext(paste0(nfix[pp], ", ", nran[qq]), side = 3, line = -8.5, cex = 2)
                } else {
                    mtext(npart[kk], side = 3, line = -3.2, cex = 2)
                    mtext(expression(m == 10^5 ~ ", " ~ n == 10^7), side = 3, line = -6, cex = 2)
                    mtext(paste0(nfix[pp], ", ", nran[qq]), side = 3, line = -8.5, cex = 2)
                }
                box(col = "grey40", lwd = 4)
                if (nn == 2) {
                    axis(side = 1, tck = -.01, labels = NA)
                    mtext(c("IEM", "30%", "50%", "70%"), at = 1:4, side = 1, line = 1.5, cex = 2)
                }
                mtext(expression("Empirical estimate of " ~ gamma), side = 2, outer = TRUE, cex = 2, line = 3.8)
                if ((nn == 1 & kk == 1 & pp == 1 & qq == 1) | (nn == 2 & kk == 1 & pp == 1 & qq == 1)) {
                    axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
                    mtext(format(axTicks(2)), at = axTicks(2), side = 2, line = 0.3, las = 1, cex = 2)
                }
            }
        }
    }
}
dev.off()
