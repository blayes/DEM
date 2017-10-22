rm(list = ls())
library(matrixStats)
library(xtable)
library(RColorBrewer)
colors <- brewer.pal(6, "Set1")

setwd("/Shared/ssrivastva/dem/ml/result")

lmerList <- list()
demList <- list()
emList <- list()
iemList <- list()
metaList <- list()

fracid <- c(30, 50, 70)
for (cc in 1:10) {
    lmerList[[cc]] <- readRDS(paste0("lmer/ml_lmer_", cc, ".rds"))
    emList[[cc]] <- readRDS(paste0("em/ml_em_", cc, ".rds"))
    metaList[[cc]] <- readRDS(paste0("lmer/meta_lmer_", cc, ".rds"))
    iemList[[cc]] <- readRDS(paste0("iem/ml_iem_cv_", cc, ".rds"))
    demList[[cc]] <- vector("list", 3)
    for (ff in 1:3) {
        demList[[cc]][[ff]] <- readRDS(paste0("dem/ml_dem_cv_", cc, "_frac_", fracid[ff], ".rds"))
    }
}

decLogLikMat <- matrix(NA, 4, 10)
decLogLikMat[1, ] <- sapply(iemList, function(x) sum(diff(x$logLik) < 0))
decLogLikMat[2, ] <- sapply(demList, function(x) sum(diff(x[[1]]$logLik) < 0))
decLogLikMat[3, ] <- sapply(demList, function(x) sum(diff(x[[2]]$logLik) < 0))
decLogLikMat[4, ] <- sapply(demList, function(x) sum(diff(x[[3]]$logLik) < 0))

metaLmerList <- list()
for (cc in 1:10) {
    metaLmerList[[cc]] <- list()
    metaLmerList[[cc]]$fixef <- rowMeans(do.call(cbind, lapply(metaList[[cc]], function(x) x$fixef)))
    metaLmerList[[cc]]$ranefCov <- matrix(rowMeans(do.call(cbind, lapply(metaList[[cc]], function(x) as.numeric(x$ranefCov)))),
                                          length(metaLmerList[[cc]]$fixef), length(metaLmerList[[cc]]$fixef))
    metaLmerList[[cc]]$dispersion <- rowMeans(do.call(cbind, lapply(metaList[[cc]], function(x) x$dispersion)))
}

ranArr <- array(NA, dim = c(7, 21, 10),
                dimnames = list(c("lmer", "meta", "em", "iem", paste0("dem", fracid)),
                                c(paste0("sig", 1:6), paste0("cov", 1:15)),
                                paste0("cv", 1:10)))
fixArr <- array(NA, dim = c(7, 6, 10),
                dimnames = list(c("lmer", "meta", "em", "iem", paste0("dem", fracid)),
                                paste0("fix", 1:6),
                                paste0("cv", 1:10)))
errVarMat <- matrix(NA, 7, 10,
                dimnames = list(c("lmer", "meta", "em", "iem", paste0("dem", fracid)),
                                paste0("cv", 1:10)))

logLikMat <- matrix(NA, 4, 10,
                dimnames = list(c("iem", paste0("dem", fracid)),
                                paste0("cv", 1:10)))

niterMat <- matrix(NA, 4, 10,
                   dimnames = list(c("iem", paste0("dem", fracid)),
                                   paste0("cv", 1:10)))

timeMat <- matrix(NA, 4, 10,
                   dimnames = list(c("iem", paste0("dem", fracid)),
                                   paste0("cv", 1:10)))

trackArr <- array(NA, dim = c(4, 20, 10),
                  dimnames = list(c("iem", paste0("dem", fracid)),
                                  paste0("mach", 1:20),
                                  paste0("cv", 1:10)))
for (cc in 1:10) {
    ranArr[1, , cc] <- c(diag(lmerList[[cc]]$ranefCov),
                         lmerList[[cc]]$ranefCov[lower.tri(lmerList[[cc]]$ranefCov)]
                         )
    ranArr[2, , cc] <- c(diag(metaLmerList[[cc]]$ranefCov),
                         metaLmerList[[cc]]$ranefCov[lower.tri(metaLmerList[[cc]]$ranefCov)]
                         )
    ranArr[3, , cc] <- c(diag(emList[[cc]]$errVar * emList[[cc]]$dmat),
                         emList[[cc]]$errVar * emList[[cc]]$dmat[lower.tri(emList[[cc]]$dmat)]
                         )
    ranArr[4, , cc] <- c(diag(iemList[[cc]]$pars$errVar * iemList[[cc]]$pars$dmat),
                         iemList[[cc]]$pars$errVar * iemList[[cc]]$pars$dmat[lower.tri(iemList[[cc]]$pars$dmat)]
                         )
    ranArr[5, , cc] <- c(diag(demList[[cc]][[1]]$pars$errVar * demList[[cc]][[1]]$pars$dmat),
                         demList[[cc]][[1]]$pars$errVar * demList[[cc]][[1]]$pars$dmat[lower.tri(demList[[cc]][[1]]$pars$dmat)]
                         )
    ranArr[6, , cc] <- c(diag(demList[[cc]][[2]]$pars$errVar * demList[[cc]][[2]]$pars$dmat),
                         demList[[cc]][[2]]$pars$errVar * demList[[cc]][[2]]$pars$dmat[lower.tri(demList[[cc]][[2]]$pars$dmat)]
                         )
    ranArr[7, , cc] <- c(diag(demList[[cc]][[3]]$pars$errVar * demList[[cc]][[3]]$pars$dmat),
                         demList[[cc]][[3]]$pars$errVar * demList[[cc]][[3]]$pars$dmat[lower.tri(demList[[cc]][[3]]$pars$dmat)]
                         )
}

for (cc in 1:10) {
    fixArr[1, , cc] <- lmerList[[cc]]$fixef
    fixArr[2, , cc] <- metaLmerList[[cc]]$fixef
    fixArr[3, , cc] <- emList[[cc]]$fixMean
    fixArr[4, , cc] <- iemList[[cc]]$pars$fixMean
    fixArr[5, , cc] <- demList[[cc]][[1]]$pars$fixMean
    fixArr[6, , cc] <- demList[[cc]][[2]]$pars$fixMean
    fixArr[7, , cc] <- demList[[cc]][[3]]$pars$fixMean
}

for (cc in 1:10) {
    errVarMat[1, cc] <- lmerList[[cc]]$dispersion
    errVarMat[2, cc] <- metaLmerList[[cc]]$dispersion
    errVarMat[3, cc] <- emList[[cc]]$errVar
    errVarMat[4, cc] <- iemList[[cc]]$pars$errVar
    errVarMat[5, cc] <- demList[[cc]][[1]]$pars$errVar
    errVarMat[6, cc] <- demList[[cc]][[2]]$pars$errVar
    errVarMat[7, cc] <- demList[[cc]][[3]]$pars$errVar
}

for (cc in 1:10) {
    trackArr[1, , cc] <- rowMeans(iemList[[cc]]$track)
    trackArr[2, , cc] <- rowMeans(demList[[cc]][[1]]$track)
    trackArr[3, , cc] <- rowMeans(demList[[cc]][[2]]$track)
    trackArr[4, , cc] <- rowMeans(demList[[cc]][[3]]$track)
}

for (cc in 1:10) {
    logLikMat[1, cc] <- max(iemList[[cc]]$logLik) / (max(emList[[cc]]$logLik))
    logLikMat[2, cc] <- max(demList[[cc]][[1]]$logLik) / (max(emList[[cc]]$logLik))
    logLikMat[3, cc] <- max(demList[[cc]][[2]]$logLik) / (max(emList[[cc]]$logLik))
    logLikMat[4, cc] <- max(demList[[cc]][[3]]$logLik) / (max(emList[[cc]]$logLik))
}

for (cc in 1:10) {
    niterMat[1, cc] <- iemList[[cc]]$niter / emList[[cc]]$iter
    niterMat[2, cc] <- demList[[cc]][[1]]$niter / emList[[cc]]$iter
    niterMat[3, cc] <- demList[[cc]][[2]]$niter / emList[[cc]]$iter
    niterMat[4, cc] <- demList[[cc]][[3]]$niter / emList[[cc]]$iter
}

for (cc in 1:10) {
    timeMat[1, cc] <- iemList[[cc]]$time[3] / emList[[cc]]$time[3]
    timeMat[2, cc] <- demList[[cc]][[1]]$time[3] / (emList[[cc]]$time[3])
    timeMat[3, cc] <- demList[[cc]][[2]]$time[3] / (emList[[cc]]$time[3])
    timeMat[4, cc] <- demList[[cc]][[3]]$time[3] / (emList[[cc]]$time[3])
}

trk <- rbind(as.numeric(trackArr[1, , ]), as.numeric(trackArr[2, , ]), as.numeric(trackArr[3, , ]), as.numeric(trackArr[4, , ]))
rownames(trk) <- c("iem", "dem30", "dem50", "dem70")


emat <- errVarMat[-3, ] - matrix(t(errVarMat[3, ]), nrow = 6, ncol = 10)

farr <- fixArr[-3, , ]
for (cc in 1:10) {
    farr[ , , cc] <- fixArr[-3, , cc] - matrix(fixArr[3, , cc], nrow = 6, ncol = 6, byrow = TRUE)
}

rarr <- ranArr[-3, , ]
for (cc in 1:10) {
    rarr[ , , cc] <- ranArr[-3, , cc] - matrix(ranArr[3, , cc], nrow = 6, ncol = 21, byrow = TRUE)
}

merrs <- round(sqrt(rowMeans(emat^2)), 4)
serrs <- round(rowSds(abs(emat)), 4)

mfixs <- round(sqrt(apply(farr^2, 1:2, mean)), 4)
sfixs <- round(apply(abs(farr), 1:2, sd), 4)

mrans <- round(sqrt(apply(rarr^2, 1:2, mean)), 4)
srans <- round(apply(abs(rarr), 1:2, sd), 4)

res <- cbind(format(mfixs, nsmall = 4),
             format(mrans, nsmall = 4),
             format(merrs, nsmall = 4))
colnames(res) <- c(colnames(mfixs),
                   colnames(mrans),
                   "err")
xtable(res)

mloglik <- rowMeans(logLikMat)
sloglik <- rowSds(logLikMat)
xtable(t(paste0(format(round(mloglik, 5), nsmall = 5), " (", format(round(sloglik, 5), nsmall = 5), ")")))

pdf("~/dem/ml/result/img/iters.pdf", 6, 9)
par(cex=1)
par(mar = c(0, 0, 0, 0), oma = c(3, 4, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
boxplot(t(niterMat), ylim = c(0, 14), ylab = NA, axes = FALSE, lwd = 2,
        boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1, staplelwd = 3, medlwd = 4, outcex = 1.5)
abline (h = 1, lty = "dotted", lwd = 5)
box(col = "grey40", lwd = 4)
axis(side = 1, tck = -.01, labels = NA)
mtext(c("IEM", "30%", "50%", "70%"), at = 1:4, side = 1, line = 1.5, cex = 2)
axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
mtext(format(axTicks(2)), at = axTicks(2), side = 2, line = 0.3, las = 1, cex = 2)
mtext(expression("No. of distributed EM iterations"  %/%  "No. of EM iterations"), side = 2, outer = TRUE, cex = 2, line = 2.2)
dev.off()


pdf("~/dem/ml/result/img/time.pdf", 6, 9)
par(cex = 1)
par(mar = c(0, 0, 0, 0), oma = c(3, 4, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
boxplot(t(timeMat), ylim = c(0, 4), ylab = NA, axes = FALSE, lwd = 2,
        boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1,
        staplelwd = 3, medlwd = 4, outcex = 1.5)
abline (h = 1, lty = "dotted", lwd = 5)
box(col = "grey40", lwd = 4)
axis(side = 1, tck = -.01, labels = NA)
mtext(c("IEM", "30%", "50%", "70%"), at = 1:4, side = 1, line = 1.5, cex = 2)
axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
mtext(format(axTicks(2)), at = axTicks(2), side = 2, line = 0.3, las = 1, cex = 2)
mtext(expression("Distributed EM time"  %/%  "EM time"), side = 2, outer = TRUE, cex = 2, line = 2.2)
dev.off()

pdf("~/dem/ml/result/img/track.pdf", 6, 9)
par(cex = 1)
par(mar = c(0, 0, 0, 0), oma = c(3, 4.8, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
boxplot(t(trk), ylim = c(0, 1), ylab = NA, axes = FALSE, lwd = 2,
        boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1, staplelwd = 3, medlwd = 4, outcex = 1.5)
box(col = "grey40", lwd = 4)
axis(side = 1, tck = -.01, labels = NA)
mtext(c("IEM", "30%", "50%", "70%"), at = 1:4, side = 1, line = 1.5, cex = 2)
axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
mtext(format(axTicks(2)), at = axTicks(2), side = 2, line = 0.3, las = 1, cex = 2)
mtext(expression("Empirical estimate of " ~ gamma), side = 2, outer = TRUE, cex = 2, line = 2.8)
dev.off()

logLikVals <- list()
for (cc in 1:10) {
    logLikVals[[cc]] <- list()
    logLikVals[[cc]][[1]] <- iemList[[cc]]$logLik
    logLikVals[[cc]][[2]] <- emList[[cc]]$logLik
    logLikVals[[cc]][[3]] <- demList[[cc]][[1]]$logLik
    logLikVals[[cc]][[4]] <- demList[[cc]][[2]]$logLik
    logLikVals[[cc]][[5]] <- demList[[cc]][[3]]$logLik
}

cc <- 1
pdf("~/dem/ml/result/img/loglik.pdf", 20, 10)
par(cex = 1)
par(mar = c(0, 0, 0, 0), oma = c(5, 6.2, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
plot(logLikVals[[cc]][[3]] / 1e5, ylab = NA, axes = FALSE, lwd = 4, type = "b", lty = "solid", pch = 0, col = colors[1], cex = 0.5)
lines(logLikVals[[cc]][[2]] / 1e5, lwd = 5, type = "b", lty = "solid", pch = 1, col = colors[2], cex = 1)
lines(logLikVals[[cc]][[4]] / 1e5, lwd = 3, type = "b", lty = "solid", pch = 2, col = colors[3], cex = 1.5)
lines(logLikVals[[cc]][[5]] / 1e5, lwd = 3, type = "b", lty = "solid", pch = 5, col = colors[4], cex = 2)
box(col = "grey40", lwd = 4)
axis(side = 1, tck = -.01, labels = NA)
axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
mtext(format(axTicks(2)), at = axTicks(2), side = 2, line = 0.5, las = 1, cex = 2)
axis(side = 1, tck = -0.01, labels = NA, lwd = 3)
mtext(format(axTicks(1)), at = axTicks(1), side = 1, line = 1.5, cex = 2)
mtext(expression("Iteration"), side = 1, outer = TRUE, cex = 2, line = 2)
mtext(expression("Log likelihood"  %*%  0.00001), side = 2, outer = TRUE, cex = 2, line = 4.2)
legend("bottom",
       c(expression("DEM with " ~ gamma  == 1.0 ~ "   "), expression("DEM with " ~ gamma == 0.7 ~ "   "), expression("DEM with " ~ gamma == 0.5 ~ "   "), expression("DEM with " ~ gamma == 0.3 ~ "   ")),
       lty = rep("solid", 4),
       pch = c(1, 5, 2, 0),
       col = colors[c(2, 4, 3, 1)],
       lwd = 5, bty = "n", cex = 2, ncol = 4)
dev.off()
