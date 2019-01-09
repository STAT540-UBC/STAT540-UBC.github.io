## TO DO (?): refactor using ggplot2, plyr

## assume working directory is where this file lives

library(lattice)

## imitate color scheme I keep showing in Keynote two group cartoon
jCols <- c(x = "blue", y = "orange")
trellis.par.set(superpose.symbol = list(col = jCols),
                superpose.line = list(col = jCols))

## load photoRec data and design from the interwebs to delay filepath pain
## use the dput/dget workflow since works better with URLs 
jURL <- "http://www.ugrad.stat.ubc.ca/~stat540/examples/photoRec/data/GSE4051_design_DPUT.txt"
str(prDes <- dget(jURL))
jURL <- "http://www.ugrad.stat.ubc.ca/~stat540/examples/photoRec/data/GSE4051_data.tsv"
str(prDat <- read.table(jURL), max.level = 0)

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

(nWt <- table(prDes$gType)["wt"])
(nKo <- table(prDes$gType)["NrlKO"])

## extract the data just for Irs4 ....

## Irs4 (insulin receptor substrate 4) was selected at random as a
## boring non differentially expressed gene; NrlKO ~= wt

## 1422248_at
##       4608
##       Irs4

prDat[4608, ]

realDat <- data.frame(prDes, gExp = unlist(prDat[4608, ]))

stripplot(gType ~ gExp, realDat, groups = gType, auto.key = TRUE)
densityplot(~ gExp, realDat, groups = gType, auto.key = TRUE)

## 'null-ify' this data, i.e. make sure the null hypothesis of equal
## group means is true ... leave data otherwise alone

(groupAvgs <- with(realDat, tapply(gExp, gType, mean)))
##       wt    NrlKO
## 7.765750 7.739684
groupAvgs <- data.frame(gType = factor(names(groupAvgs),
                        levels = levels(realDat$gType)),
                        groupAvg = groupAvgs)
(globalAvg <- mean(realDat$gExp))
## [1] 7.753051

(theDiff <- groupAvgs$groupAvg["NrlKO"] - groupAvgs$groupAvg["wt"])
##       NrlKO
## -0.02606579

(ttRes <- t.test(gExp ~ gType, realDat))
(welchStat <- -1 * ttRes$stat)          # I'm doing ko - wt

(groupVars <- with(realDat, tapply(gExp, gType, var)))
##         wt      NrlKO
## 0.02403557 0.02332078
(commVar <- sum(groupVars * 1/c(nWt, nKo)))
## [1] 0.002429188


nullDat <- merge(realDat, groupAvgs)
nullDat$gExp <- with(nullDat, gExp - groupAvg + globalAvg)
nullDat$groupAvg <- NULL

## did i really do it?
mean(nullDat$gExp)
globalAvg
with(nullDat, tapply(gExp, gType, mean))
## YES

## make figures of the data-generating distributions

stripplot(gType ~ gExp, nullDat,
          groups = gType, auto.key = TRUE)

densityplot(~ gExp, nullDat, groups = gType, auto.key = TRUE, n = 400)

jFudge <- 0.35
jNudge <- 0.04
jCex <- 1.5
densityplot(~ gExp, nullDat, groups = gType, auto.key = TRUE,
            n = 400, plot.points = FALSE, ref = TRUE,
            panel = panel.superpose,
            panel.groups = function(x, ..., group.number) {
                yo <- group.number
                panel.densityplot(x, ...)
                theAvg <- mean(x)
                panel.points(theAvg, (group.number * 2 - 3) * jNudge,
                             pch = 19, cex = jCex,
                             col = jCols[group.number])
                if(group.number == 1) {
                    jLab <- substitute(paste(mu, " = ", foo),
                                       list(Z = c("Y", "Z")[group.number],
                                            foo = round(theAvg, 2)))
                    panel.text(x = mean(x), y = jFudge, jLab)
                }
            })

dev.print(pdf, "../figs/twoGroupsSimulation/dataGenDensityplot.pdf",
          width = 7, height = 4.5)

## generate the bootstrap data
set.seed(101)
B <- 10000
(nWt <- sum(nullDat$gType == "wt"))     # 20
(nKo <- sum(nullDat$gType == "NrlKO"))  # 19

bootDatWt <-
  matrix(sample(nullDat$gExp[nullDat$gType == "wt"],
                nWt * B, replace = TRUE),
                nrow = nWt, ncol = B)
str(bootDatWt)
##  num [1:20, 1:10000] 7.65 7.85 7.86 7.85 7.72 ...

bootDatKo <-
  matrix(sample(nullDat$gExp[nullDat$gType == "NrlKO"],
                nKo * B, replace = TRUE),
                nrow = nKo, ncol = B)
str(bootDatKo)
## num [1:19, 1:10000] 7.58 7.56 7.78 7.58 7.76 ...

bootDat <- rbind(bootDatWt, bootDatKo)
str(bootDat)
## num [1:39, 1:10000] 7.65 7.85 7.86 7.85 7.72 ...

## let's look at the first 4 datasets

## good example of data reshaping for the purposes of plotting
bootDatToPlot <- stack(as.data.frame(bootDat[ , 1:4]))
str(bootDatToPlot)
## 'data.frame':	156 obs. of  2 variables:
##  $ values: num  7.65 7.85 7.86 7.85 7.72 ...
##  $ ind   : Factor w/ 4 levels "V1","V2","V3",..: 1 1 1 1 1 1 1 1 1 1 ...

bootDatToPlot <-                              # tidy up
  data.frame(gType = nullDat$gType,
             bootIter = unclass(bootDatToPlot$ind),
             gExp = bootDatToPlot$values)
head(bootDatToPlot)

jFudge <- 0.35
jCex <- 3
stripplot(gType ~ gExp | factor(bootIter), bootDatToPlot,
          groups = gType, grid = TRUE, jitter.data = TRUE,
          panel = panel.superpose,
          panel.groups = function(x, y, ..., group.number) {
            yo <- group.number
            panel.stripplot(x, y, ...)
            theAvg <- mean(x)
            panel.points(theAvg, y[1], pch = "|", cex = jCex,
                         col = jCols[group.number])
            jLab <- substitute(paste(bar(Z), " = ", foo),
                               list(Z = c("Y", "Z")[group.number],
                                    foo = round(theAvg, 2)))
            panel.text(x = mean(x), y = y[1] + jFudge,
                       jLab)
          })

dev.print(pdf, "../figs/twoGroupsSimulation/stripplot-firstFewSimDsets.pdf",
          width = 5, height = 5)

computeStats <- function(jObs, jRv) {
  foo <- t.test(jObs ~ jRv)
  return(c(foo$estimate[1] - foo$estimate[2],
           foo$statistic, foo$p.value))
}

system.time(
    bootTestStats <-
            apply(bootDat, 2,
                  computeStats, jRv = nullDat$gType)
            )                           # ~10 secs for JB

str(bootTestStats)
## num [1:3, 1:10000] -0.0492 -1.1866 0.245 -0.0126 -0.2422 ...

bootTestStats <- data.frame(t(bootTestStats))
names(bootTestStats) <- c("smDiff", "tStat", "pVal")

head(bootTestStats)

densityplot(~ smDiff, bootTestStats,
            xlab = expression(bar(Z) - bar(Y)), n = 400,
            panel = function(...) {
              panel.densityplot(...)
              panel.mathdensity(args = list(mean = 0,
                                  sd = sqrt(commVar)),
                                n = 300, col.line = "grey75",
                                lty = "dashed", lwd = 2)
              panel.abline(v = theDiff)
            })

dev.print(pdf, "../figs/twoGroupsSimulation/densityplot-empiDistDiff.pdf",
          width = 5, height = 5)


densityplot(~ tStat, bootTestStats,
            xlab = "Welch's t statistic", n = 400,
            panel = function(...) {
              panel.densityplot(...)
              panel.mathdensity(dmath = dt,
                                args = list(df = nWt + nKo - 2),
                                n = 300, col.line = "grey75",
                                lty = "dashed", lwd = 2)
              panel.abline(v = welchStat)
            })

dev.print(pdf, "../figs/twoGroupsSimulation/densityplot-empiDistTstat.pdf",
          width = 5, height = 5)

mean(abs(bootTestStats$tStat) >= abs(welchStat))
## [1] 0.5942

mean(abs(bootTestStats$smDiff) >= abs(theDiff))
## [1] 0.5818

histogram(~ pVal, bootTestStats,
          xlab = "p-value for Welch's t statistic",
          col = "grey70", border = "grey30",
          endpoints = c(0, 1))

dev.print(pdf, "../figs/twoGroupsSimulation/histogram-empiDistPvals.pdf",
          width = 5, height = 5)

### 2014: nothing below here updated or re-run! originally used to produce
### figures looking at the effect of different variance or sample, holding all
### other things equal, on statistical significance
str(jDat)

lDat <- subset(jDat, n == "big")
str(lDat)

with(lDat,
     by(lDat, sigStat, function(yo) {
       t.test(obs ~ rv, yo)
     }))


jennyStrip <-
  function(which.panel, factor.levels, ...) {
    panel.rect(0, 0, 1, 1, col = "grey90", border = 1)
    foo <- factor.levels[which.panel]
    jLab <- substitute(paste(sigma^2, " ", foo),
                       list(foo = foo))
    panel.text(x = 0.5, y = 0.5, pos = 2, lab = jLab)
  }


stripplot(rv ~ obs | sigStat, lDat,
          type = c("p", "g"),
          groups = rv, layout = c(1, 2),
          as.table = TRUE, strip = jennyStrip,
          par.strip.text = list(lines = 1.25))

dev.print(pdf,
          paste0(whereAmI, "figs/twoGrpComp-stripplot-divideVarBy10.pdf"),
          width = 5, height = 5)


mDat <- subset(jDat, sigStat == "big")
str(mDat)

with(mDat,
     by(mDat, n, function(yo) {
       t.test(obs ~ rv, yo)
     }))

with(mDat,
     by(mDat, n, function(yo) {
       wilcox.test(obs ~ rv, yo)
     }))


with(mDat,
     by(mDat, n, function(yo) {
       ks.test(x = yo$obs[yo$rv == "x"],
               y = yo$obs[yo$rv == "y"])
     }))


n <- with(subset(mDat, rv == "x"),
          table(n))

jennyStrip <-
  function(which.panel, factor.levels, ...) {
    panel.rect(0, 0, 1, 1, col = "grey90", border = 1)
    panel.text(x = 0.5, y = 0.5, pos = 2,
               lab = paste("n = ", n[which.panel]))
  }


stripplot(rv ~ obs | n, mDat,
          type = c("p", "g"),
          groups = rv, layout = c(1, 2),
          as.table = TRUE, strip = jennyStrip,
          par.strip.text = list(lines = 1.25))

dev.print(pdf,
          paste0(whereAmI, "figs/twoGrpComp-stripplot-changeSampleSize.pdf"),
          width = 5, height = 5)

