library(limma)
library(hexbin)

##########################################################
## load data and design
##########################################################

## load design data.frame
prDes <- readRDS("../data/GSE4051_design.rds")
str(prDes) # 39 obs. of  4 variables

## load gene expression data.frame
prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)
## 'data.frame':  29949 obs. of  39 variables:

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

##########################################################
## helper functions for finding and featuring genes
##########################################################

source("98_lectSupp-prepareDataFunction.r")
source("99_lectSupp-plotItFunction.r")

##########################################################
## source in-house functions for post-processing mlm objects
##########################################################

## code written by Rick White
source("80_anova-mlm.r")

## was helpful when printing stuff to screen to insert into lecture slides
options(width = 110)

##########################################################
## limma
##########################################################
jDesMat <- model.matrix(~ gType * devStage, prDes)

## make some stuff to show in lecture slides
## ridiculous machination to print a version to screen with small
## variable names
str(jDesMat)
cbind(colnames(jDesMat))
foo <- jDesMat
dimnames(jDesMat)
colnames(foo) <- paste0("X", sprintf("%02d", seq_len(ncol(jDesMat))))
cbind(prDes, foo)

## fitting two-way ANOVA to all probesets at once
system.time(jFit <- lmFit(prDat, jDesMat))
str(jFit, max.level = 3)
summary(jFit)

## can we get hits yet?
topTable(jFit)
topTable(jFit, coef = 2)
topTable(jFit, coef = "gTypeNrlKO")
## NO

ebFit <- eBayes(jFit)
str(ebFit, max.level = 3)
summary(ebFit)

## can we get hits now?
topTable(ebFit)
## YES

##########################################################
## lm for multivariate regression
##########################################################
## prepare to use homegrown mlm code
## responses must be a matrix, one *row* per response
prMat <- t(as.matrix(prDat))
twAnova <- lm(prMat ~ gType * devStage, prDes)
twSumm <- summary(twAnova)
twAnovaNoInter <- lm(prMat ~ gType + devStage, prDes)
twNiSumm <- summary(twAnovaNoInter)
dsOnly <- lm(prMat ~ devStage, prDes)
dsSumm <- summary(dsOnly)

##########################################################
## testing gTypeNrlKO
##########################################################

## this is equivalent: topTable(ebFit, coef = 2)
## but much more self-documenting; USE NAMES!
## you will be glad you did when you revisit/read your code
topTable(ebFit, coef = "gTypeNrlKO")

hits <- topTable(ebFit, coef = "gTypeNrlKO")

stripplotIt(hitDat <- prepareData(head(rownames(hits), 6)),
            scales = list(rot = c(45, 0)))

dev.print(pdf,
          "../figs/limma/stripplot-NrlKO.E16-hits.pdf",
          width = 10, height = 6)

##########################################################
## testing for interaction
##########################################################
colnames(coef(ebFit))                 # parameters in the linear model
grep(":", colnames(coef(ebFit)))      # isolate the interaction terms
hits <- topTable(ebFit, coef = grep(":", colnames(coef(ebFit))))

stripplotIt(hitDat <- prepareData(head(rownames(hits), 6)),
            scales = list(rot = c(45, 0)))

dev.print(pdf,
          "../figs/limma/stripplot-interaction-hits.pdf",
          width = 10, height = 6)

(huh <- rownames(hits)[1])
stripplotIt(huhDat <- prepareData(huh))

huhFitBig <- lm(gExp ~ gType * devStage, huhDat)
huhFitSmall <- lm(gExp ~ gType + devStage, huhDat)
anova(huhFitSmall, huhFitBig)
summary(huhFitBig)

hits <- topTable(ebFit, coef = grep(":", colnames(coef(ebFit))),
                 n = Inf, sort = "none")

foo <- anova(twAnova, twAnovaNoInter)[ , "M1:Pval"]
str(foo)

## compare the p-values for testing whether all interaction terms are zero:
## classical stats vs. limma's emp Bayes stats
xyplot(hits$P.Value ~ foo, xlab = "interaction p-value, lm",
       ylab = "interaction p-value, limma", aspect = 1, xbins = 50,
       panel = function(x, y, ...) {
           panel.hexbinplot(x, y, ...)
           panel.abline(a = 0, b = 1, col = "orange")
       })

dev.print(pdf,
          "../figs/limma/smoothScatter-interactionPvalsLimmaVsLm.pdf",
          width = 7, height = 7)

densityplot(~ foo + hits$P.Value, plot.points = FALSE, n = 300, xlab = "p-value")

dev.print(pdf,
          "../figs/limma/densityplot-interactionPvalsLimmaVsLm.pdf",
          width = 7, height = 7)

table(hits$P.Value > foo)
## FALSE  TRUE
## 17702 12247

table(hits$P.Value > foo)/nrow(prDat)
##     FALSE      TRUE
## 0.5910715 0.4089285

sum(foo < 0.01) # 1874
sum(hits$P.Value < 0.01) # 1888

##########################################################
## does genotype matter?
##########################################################
colnames(coef(ebFit))                 # parameters in the linear model
grep("gType", colnames(coef(ebFit)))  # isolate the genotype terms
hits3 <- topTable(ebFit,
                  coef = grep("gType", colnames(coef(ebFit))),
                  n = Inf)

stripplotIt(hitDat <- prepareData(c(head(rownames(hits3), 3),
                                    tail(rownames(hits3), 3))),
            scales = list(rot = c(45, 0)))

dev.print(pdf,
          "../figs/limma/stripplot-gType-hits.pdf",
          width = 10, height = 6)

huh <- twSumm$Full[ , "Pval"]
foo <- anova(twAnova, dsOnly)[ , "M1:Pval"]
str(foo)

hits3 <- topTable(ebFit,
                  coef = grep("gType", colnames(coef(ebFit))),
                  n = Inf, sort = "none")

xyplot(foo ~ hits3$P.Value)
densityplot(~ foo + hits3$P.Value, plot.points = FALSE, n = 300, xlab = "p-value")

## need this for multiple testing lecture to talk about multiple
## testing across contrasts
topTable(ebFit, coef = grep("gType", colnames(coef(ebFit))))

hits3 <- topTable(ebFit,
                  coef = grep("gType", colnames(coef(ebFit))))

##########################################################
## looking at variance
##########################################################

## 'raw' residual variance
str(jFit$sigma)

## df associated with the prior
ebFit$df.prior
## [1] 3.102649

## location associated with the prior
ebFit$s2.prior
## [1] 0.06166446

## comparing the prior mean to the mean of the 'raw' residual variances
exp(mean(log(jFit$sigma ^ 2)))
## [1] 0.08498974

## "Specifically, d0 and s20 are estimated by equating empirical to expected
## values for the first two moments of log s2g." So why can't I reproduce
## ebFit$s2.prior above? Because the methods of moments estimation is not quite
## simple as it sounds. In fact, d0 is obtained first, then plugged into another
## expression to yield s20. The result is described like so: "This estimate for
## s20 is usually somewhat less than the mean of the s2g in recognition of the
## skewness of the F-distribution." Which is compatible with what I see here.
## All quotes are from Smyth 2004 http://www.statsci.org/smyth/pubs/ebayes.pdf

## post-eb residual variance
str(ebFit$s2.post)

summary(jFit$sigma)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## 0.0704  0.1993  0.2723  0.3305  0.4219  1.4260

cbind(raw = summary(jFit$sigma ^ 2), eb = summary(ebFit$s2.post))
#              raw      eb
# Min.    0.004956 0.01044
# 1st Qu. 0.039710 0.04183
# Median  0.074170 0.07296
# Mean    0.140700 0.13300
# 3rd Qu. 0.178000 0.16670
# Max.    2.033000 1.84200

table(ebFit$s2.post > jFit$sigma^2)
## FALSE  TRUE
## 17129 12820

table(ebFit$s2.post > jFit$sigma^2)/nrow(prDat)
##    FALSE     TRUE
## 0.571939 0.428061

xyplot(ebFit$s2.post ~ jFit$sigma ^ 2,
       aspect = 1,
       prepanel = function(x, y, ...) {
           rng <- range(x, y, finite = TRUE)
           ## comment next line out to get all vars!
           #rng = c(0, 0.1)
           list(xlim = rng, ylim = rng)
       },
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#          panel.smoothScatter(x, y, ...)
           panel.abline(a = 0, b = 1, col = "grey40")
           with(ebFit,
                panel.abline(a = (df.prior/median(df.total)) * s2.prior,
                             b = (df.residual/df.total),
                             col = "orange"))
      })

dev.print(pdf,
          #"../figs/limma/xyplot-rawAndModeratedResidVarZoom.pdf",
          "../figs/limma/xyplot-rawAndModeratedResidVar.pdf",
          width = 7, height = 7)
