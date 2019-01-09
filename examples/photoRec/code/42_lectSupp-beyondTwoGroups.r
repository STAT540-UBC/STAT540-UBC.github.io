## TO DO (?): refactor using ggplot2, plyr

## assume working directory is where this file lives

library(lattice)

## load design data.frame
prDes <- readRDS("../data/GSE4051_design.rds")
str(prDes)

## load gene expression data.frame
prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)
## 'data.frame':  29949 obs. of  39 variables:

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

## in 04_dea-no-empirical-bayes.r, I've already done differential expression
## analysis by devStage for all probesets; load those results; use them to pick
## genes to feature here; also good for sanity checks
deDevStage <- read.table("../results/dea-devStage-mlm-stats.txt")
str(deDevStage)
# 'data.frame':  29949 obs. of  4 variables:
# $ F    : num  1.393 0.542 0.706 1.63 5.515 ...
# $ pVal : num  0.25722 0.70624 0.59329 0.18935 0.00156 ...
# $ sigma: num  0.151 0.603 0.398 0.168 0.108 ...
# $ rank : num  24743 29149 28541 23024 5201 ...


## used for selecting genes to feature as hits
## (tryMe <- which(deDevStage$rank == 1))
## miniDat <- as.vector(t(prDat[tryMe, ]))
## miniDat <- data.frame(gene = tryMe, gExp = miniDat)
## miniDat <- data.frame(prDes, miniDat)
## stripplot(gExp ~ devStage, miniDat)
## rownames(prDat)[tryMe]

## 17483 gradual increase at first, then huge at 4_weeks
## 3966, 24029 steady increase
## 9345 in between those two scenarios
## 23043 big gains, almost all done by P10

## I've selected this as our hit
(theHit <- which(rownames(prDat) == "1440645_at")) # 17843
## and this as our boring gene
(theBore <- which(rownames(prDat) == "1443184_at")) # 18898

(keepers <- data.frame(row = c(theBore, theHit),
                       probesetID = I(rownames(prDat)[c(theBore, theHit)])))
## WOW! forget to protect probesetID with I() and you get totally the
## wrong data; great example of Jenny's Factor Law!
miniDat <- as.vector(t(prDat[keepers$probesetID, ]))
miniDat <- data.frame(gene = rep(c("theBore", "theHit"), each = nrow(prDes)),
                      gExp = miniDat)
str(miniDat)
miniDat <- data.frame(prDes, miniDat)
str(miniDat)


stripplot(gExp ~ devStage | gene, miniDat,
          grid = TRUE, type = c('p', 'a'),
          jitter.data = TRUE)

dev.print(pdf, "../figs/devStage/stripplot.pdf", width = 7, height = 4)

stripplot(gExp ~ devStage | gene, miniDat,
          subset = gene == "theBore", ylim = c(6.5, 8.5),
          grid = TRUE, type = c('p', 'a'),
          jitter.data = TRUE)

dev.print(pdf, "../figs/devStage/stripplot-theBore.pdf", width = 7, height = 4)

stripplot(gExp ~ devStage | gene, miniDat,
          subset = gene == "theHit",
          grid = TRUE, type = c('p', 'a'),
          jitter.data = TRUE)

dev.print(pdf, "../figs/devStage/stripplot-theHit.pdf", width = 7, height = 4)

with(miniDat,
     tapply(gExp, list(devStage, gene), mean))

(theBoreAvgs <- with(subset(miniDat, gene == "theBore"),
                     tapply(gExp, devStage, mean)))
data.frame(cellMeans = theBoreAvgs,
           txEffects = theBoreAvgs - theBoreAvgs[1])

boreFit <- lm(gExp ~ devStage, miniDat, gene == "theBore")
summary(boreFit)

(theHitAvgs <- with(subset(miniDat, gene == "theHit"),
                     tapply(gExp, devStage, mean)))
data.frame(cellMeans = theHitAvgs,
           txEffects = theHitAvgs - theHitAvgs[1])

hitFit <- lm(gExp ~ devStage, miniDat, gene == "theHit")
summary(hitFit)
summary(hitFit)$coef

with(subset(miniDat, gene == "theHit"),
     model.matrix(~ devStage))

hitFitCellMeans <- lm(gExp ~ 0 + devStage, miniDat, gene == "theHit")
summary(hitFitCellMeans)
summary(hitFitCellMeans)$coef

with(subset(miniDat, gene == "theHit"),
     model.matrix(~ devStage))

## nothing below here used in 2013 or 2014
## in general, I am trying to switch from abstract / simulated
## examples to using the photoRec dataset to show as much as possible

## institute color scheme based on original Keynote cartoon
jCols <- c(ref = "gray60", A = "blue", B = "orange", C = "green4")
trellis.par.set(superpose.symbol = list(col = jCols),
                superpose.line = list(col = jCols))
jCex <- 1.3                             # used for sample means

trueGroupMeans <- c(0, 1, 4, 7)
n <- c(5, 4, 7, 5)
nGrp <- length(trueGroupMeans)

set.seed(783)
jDat <- data.frame(grp = rep(paste0("grp", 1:nGrp), n),
                   y = unlist(lapply(1:nGrp, function(j) {
                     return(rnorm(n = n[j], mean = trueGroupMeans[j]))
                   })))

jFudge <- 0.35
stripplot(grp ~ y, jDat,
          type = c("p", "g"), xlab = "",
          groups = grp, jitter.data = TRUE,
          panel = panel.superpose,
          panel.groups = function(x, y, ..., group.number) {
            yo <- group.number
            panel.stripplot(x, y, ...)
            theAvg <- mean(x)
            panel.points(theAvg, y[1], pch = 2, cex = jCex,
                         col = jCols[group.number])
            panel.points(trueGroupMeans[group.number], y[1],
                         pch = 19, cex = jCex,
                         col = jCols[group.number])
          },
          key = list(x = 0.03, y = 0.96,
            text = list(c("truth", "observed")),
            points = list(pch = c(19, 2))))

dev.print(pdf,
          paste0(whereAmI, "figs/fourGrpComp-stripplot.pdf"),
          width = 5, height = 3.5)


stripplot(grp ~ y, jDat,
          type = c("p", "g"), xlab = "",
          groups = grp, jitter.data = TRUE,
          panel = function(x, y, ...) {
            panel.stripplot(x, y, ...)

          })

(theAvgs <- with(jDat,
                 tapply(y, grp, mean)))
##      grp1      grp2      grp3      grp4
## 0.1008892 0.5573339 4.3096060 6.8967497

data.frame(truth = trueGroupMeans,
           obs = round(theAvgs, 4))

##      truth       obs
## grp1     0 0.1008892
## grp2     1 0.5573339
## grp3     4 4.3096060
## grp4     7 6.8967497

summary(lm(y ~ 0 + grp, jDat))

summary(lm(y ~ grp, jDat))

data.frame(truth = trueGroupMeans - trueGroupMeans[1],
           obs = round(theAvgs - theAvgs[1], 4))[-1, ]

trueGroupMeans[2] - trueGroupMeans[1]

kDat <- data.frame(jDat,
                   geneA = factor(rep(c("wt", "mut",
                     "wt", "mut"), n),
                     levels = c("wt", "mut")),
                   geneB = factor(rep(c("wt", "mut"),
                     c(n[1] + n[2], n[3] + n[4])),
                     levels = c("wt", "mut")))

summary(lm(y ~ geneA * geneB, kDat))

theAvgs[4] -
  (theAvgs[1] +
   theAvgs[2] - theAvgs[1] +
   theAvgs[3] - theAvgs[1])

fitGroups <- lm(y ~ 0 + grp, jDat)

fitGpEffects <- lm(y ~ grp, jDat)

fitTwoWay <- lm(y ~ geneA * geneB, kDat)

fitNothing <- lm(y ~ 0, kDat)

fitIntercept <- lm(y ~ 1, kDat)

## F num =
## (RSS(small) - RSS(big)) / (nParams(big) - nParams(small))

## F denom =
## RSS(big) / (n - nParams(big))


names(fitGroups)
names(summary(fitGroups))
names(aov(fitGroups))

sum(resid(fitGroups)^2)

rDat <-
  data.frame(pzation = I(c("groups", "gpEffects", "twoWay")),
             RSSBig = c(sum(resid(fitGroups)^2),
               sum(resid(fitGpEffects)^2),
               sum(resid(fitTwoWay)^2)),
             nParamsBig = sum(n) - c(fitGroups$df.residual,
               fitGpEffects$df.residual,
               fitTwoWay$df.residual),
             RSSSmall = c(sum(kDat$y^2),
               sum((kDat$y - mean(kDat$y))^2),
               sum((kDat$y - mean(kDat$y))^2)),
             nParamsSmall = c(0, 1, 1))

rDat$numDf <- with(rDat, nParamsBig - nParamsSmall)
rDat$denDf <- with(rDat, sum(n) - nParamsBig)

rDat$FStat <-
  with(rDat,
       ( (RSSSmall - RSSBig) / numDf ) /
       (RSSBig / denDf) )
rDat

anova(fitNothing, fitGroups)
summary(fitGroups)

anova(fitIntercept, fitGpEffects)
summary(fitGroups)

anova(fitIntercept, fitTwoWay)
summary(fitTwoWay)
