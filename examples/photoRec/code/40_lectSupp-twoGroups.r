## TO DO (?): refactor using ggplot2, plyr

## assume working directory is where this file lives

library(lattice)
library(latticeExtra)                   # ecdfplot()

## imitate color scheme I keep showing in Keynote two group cartoon
jCols <- c(x = "blue", y = "orange")
trellis.par.set(superpose.symbol = list(col = jCols),
                superpose.line = list(col = jCols))
jCex <- 3                               # used for sample means

## load photoRec data and design from the interwebs to delay filepath pain
## use the dput/dget workflow since works better with URLs 
jURL <- "http://www.ugrad.stat.ubc.ca/~stat540/examples/photoRec/data/GSE4051_design_DPUT.txt"
str(prDes <- dget(jURL))
jURL <- "http://www.ugrad.stat.ubc.ca/~stat540/examples/photoRec/data/GSE4051_data.tsv"
str(prDat <- read.table(jURL), max.level = 0)

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

## Nrl (neural retina leucine zipper gene) is the gene that was
## knocked out in half the mice; obviously should be differentially
## expressed

## Irs4 (insulin receptor substrate 4) was selected as a boring non 
## differentially expressed gene

## Knowledge from the future: From downstream differential expression analysis, 
## I know that the probeset representing Nrl has rank 1 (!!!), i.e. it shows
## most evidence for DE, and the probeset representing Irs4 has rank 21000 
## (that's how I selected it, actually). These ranks refer to the rank variable 
## in the data.frame stored later as deGtype, FYI.

## probeset IDs:
## Nrl = 1450946_at
## Irs4 = 1422248_at

miniDat <- as.vector(t(prDat[c("1422248_at", "1450946_at"), ]))
miniDat <- data.frame(gene = rep(c("Irs4", "Nrl"), each = nrow(prDes)),
                      gExp = miniDat)
miniDat <- data.frame(prDes, miniDat) # ignore the warning about row names
str(miniDat)

###############################################
## making basic "two groups" figures
###############################################

stripplot(gType ~ gExp | gene, miniDat,
          scales = list(x = list(relation = "free")),
          groups = gType, auto.key = TRUE,
          aspect = 0.5, layout = c(2, 1))

dev.print(pdf, file = "../figs/twoGroups/stripplot.pdf",
          width = 7, height = 3.5)

stripplot(gType ~ gExp | gene, miniDat,
          groups = gType, auto.key = TRUE,
          aspect = 0.5, layout = c(2, 1))

dev.print(pdf, "../figs/twoGroups/stripplot-commonScale.pdf",
          width = 7, height = 3.5)

densityplot(~ gExp | gene, miniDat,
            scales = list(x = list(relation = "free")),
            groups = gType, auto.key = TRUE)

dev.print(pdf, "../figs/twoGroups/densityplot.pdf",
          width = 7, height = 4.5)

bwplot(gType ~ gExp | gene, miniDat,
       scales = list(x = list(relation = "free")))

dev.print(pdf, "../figs/twoGroups/boxplot.pdf",
          width = 7, height = 4.5)

bwplot(gType ~ gExp | gene, miniDat,
       scales = list(x = list(relation = "free")),
       panel = function(...) {
         panel.violin(..., col = "grey90")
         panel.stripplot(..., col = "black")
       })

dev.print(pdf, "../figs/twoGroups/violinAndStripplot.pdf",
          width = 7, height = 4.5)

bwplot(gType ~ gExp | gene, miniDat,
       scales = list(x = list(relation = "free")),
       panel = function(..., box.ratio) {
         panel.violin(..., col = "transparent", border = "black",
                      varwidth = FALSE, box.ratio = box.ratio)
         panel.bwplot(..., fill = NULL, box.ratio = .1)
       })

dev.print(pdf, "../figs/twoGroups/violinAndBoxplot.pdf",
          width = 7, height = 4.5)


jFudge <- 0.35
stripplot(gType ~ gExp | gene, miniDat,
          grid = TRUE,
          scales = list(x = list(relation = "free")),
          groups = gType, jitter.data = TRUE,
          aspect = 0.5, layout = c(2, 1),
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

dev.print(pdf, "../figs/twoGroups/stripplot-yBar-zBar.pdf",
          width = 7, height = 3.5)

###############################################
## "two groups" inference: t tests
###############################################

## compute difference of sample averages
(theAvgs <- with(miniDat,
                 tapply(gExp, list(gType, gene), mean)))

(theDiff <- theAvgs["NrlKO", ] - theAvgs["wt", ])

## compute sample variances
(theVars <- with(miniDat,
                 tapply(gExp, list(gType, gene), var)))

## compute estimated variance of zbar - ybar
(nY <- with(miniDat, sum(gType == "wt" & gene == "Nrl")))
(nZ <- with(miniDat, sum(gType == "NrlKO" & gene == "Nrl")))

## assuming unequal true variance
(s2DiffWelch <- colSums(theVars / c(nY, nZ)))

## assuming equal true variance
(s2Pooled <- colSums(theVars * c((nY - 1) / (nY + nZ - 2),
                                 (nZ - 1) / (nY + nZ - 2))))

(s2Diff <- s2Pooled * (1/nY + 1/nZ))

(welchStat <- theDiff / sqrt(s2DiffWelch))

by(miniDat, miniDat$gene, function(theDat) {
  ## to be compatible with the direction of the test, as I've described in
  ## lecture slides, I need to reverse order of levels of gType factor
  theDat$gType <- factor(theDat$gType, rev(levels(theDat$gType)))
  t.test(gExp ~ gType, theDat)
})

(tstStat <- theDiff / sqrt(s2Diff))

by(miniDat, miniDat$gene, function(theDat) {
  theDat$gType <- factor(theDat$gType, rev(levels(theDat$gType)))
  t.test(gExp ~ gType, theDat, var.equal = TRUE)
})

## draw the null dist'n

nullDistn <- data.frame(tVals = seq(-3.5, 3.5, length = 300))
nullDistn$densT <- dt(nullDistn$tVals, df = nY + nZ - 2)
nullDistn$densNorm <- dnorm(nullDistn$tVals)

with(nullDistn,
     matplot(x = nullDistn$tVals,
             y = nullDistn[c("densT", "densNorm")],
             type = "l"))

## depicting only for Irs4!

with(nullDistn,
     plot(densNorm ~ tVals, type = "l"))
abline(v = tstStat["Irs4"], lty = "dashed")

foo <- subset(nullDistn, tVals < tstStat["Irs4"])
polygon(x = c(foo$tVals, rev(foo$tVals)),
        y = c(rep(0, nrow(foo)), rev(foo$densNorm)),
        col = "blue")

polygon(x = -1 * c(foo$tVals, rev(foo$tVals)),
        y = c(rep(0, nrow(foo)), rev(foo$densNorm)),
        col = "blue")

dev.print(pdf, "../figs/twoGroups/pValBasedOnNorm.pdf",
          width = 7, height = 5)


round(pt(-1 * abs(tstStat), df = nY + nZ - 2) * 2, 5)
round(pnorm(-1 * abs(tstStat)) * 2, 5)

###############################################
## "two groups" inference: Wilcoxon
###############################################

(wilcTest <- by(miniDat, miniDat$gene, function(theDat) {
    wilcox.test(gExp ~ gType, theDat)
}))

## Warning messages:
## 1: In wilcox.test.default(x = c(7.867, 7.783, 7.585, 7.4, 7.667, 7.818,  :
##   cannot compute exact p-value with ties
## 2: In wilcox.test.default(x = c(8.761, 9.079, 10.61, 9.796, 11.39,  :
##   cannot compute exact p-value with ties

## hmmm..... there are *ties*? really?
by(miniDat, miniDat$gene, function(theDat) rank(theDat$gExp))
## yep, sure are, for both Irs4 and Nrl .... strange

## let's recreate the Wilcoxon test statistic
(jRankSums <- by(miniDat, miniDat$gene, function(theDat) {
    tapply(rank(theDat$gExp), theDat$gType, sum)
}))

(sampSize <- c(nY, nZ))

## from wilcox.text, we learn this: "The literature is not unanimous
## about the definitions of the Wilcoxon rank sum and Mann-Whitney
## tests.  The two most common definitions correspond to the sum of
## the ranks of the first sample with the minimum value subtracted or
## not: R subtracts and S-PLUS does not, giving a value which is
## larger by m(m+1)/2 for a first sample of size m."

sapply(jRankSums, function(foo) foo - sampSize * (sampSize + 1)/2)

## focus on the first row of that table -- the rank sums for wild type
## with the associated minimum subtracted

sapply(wilcTest, function(foo) foo$stat)

## Compare to output of the "official" test.
## See! we've computed the test statistic for ourselves.

###############################################
## "two groups" inference: Kolmogorov-Smirnov
###############################################

## sadly, ks.test() has no formula interface

(ksTest <- by(miniDat, miniDat$gene, function(theDat) {
    ks.test(x = theDat$gExp[theDat$gType == "wt"],
            y = theDat$gExp[theDat$gType == "NrlKO"])
}))
## more moaning about the ties


## fodder for recreating the KS test stat by hand ... leaving for now
## THIS CODE WILL NOT WORK 'AS IS'!!

## foo <- ecdf(miniDat$gExp)

## x <- rnorm(12)
## Fn <- ecdf(x)
## Fn     # a *function*
## Fn(x)  # returns the percentiles for x
## tt <- seq(-2,2, by = 0.1)
## 12 * Fn(tt) # Fn is a 'simple' function {with values k/12}
## summary(Fn)
## ##--> see below for graphics
## knots(Fn)# the unique data values {12 of them if there were no ties}
## plot(tt, Fn(tt))

ecdfplot(~ gExp | gene, miniDat,
         scales = list(x = list(relation = "free")),
          groups = gType)

dev.print(pdf, "../figs/twoGroups/ecdfplot.pdf",
          width = 7, height = 4.5)

###############################################
## direct support of first little bit of "beyond two groups" lecture
###############################################

## stripplot **for Irs4** by itself
stripplot(gType ~ gExp | gene, miniDat,
          subset = gene == "Irs4",
          groups = gType, auto.key = TRUE,
          aspect = 0.5)

dev.print(pdf, "../figs/twoGroups/stripplot-Irs4.pdf",
          width = 7, height = 3.5)

## stripplot with sample means **for Irs4** by itself
jFudge <- 0.35
stripplot(gType ~ gExp | gene, miniDat,
          subset = gene == "Irs4",
          grid = TRUE, aspect = 0.5,
          groups = gType, jitter.data = TRUE,
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

dev.print(pdf, "../figs/twoGroups/stripplot-yBar-zBar-Irs4.pdf",
          width = 7, height = 3.5)

t.test(gExp ~ gType, miniDat,
       subset = gene == "Irs4", var.equal = TRUE)

summary(aov(gExp ~ gType, miniDat,
            subset = gene == "Irs4"))

summary(lm(gExp ~ gType, miniDat,
           subset = gene == "Irs4"))

## for pasting into slides
theAvgs
7.739684 - 7.765750
t.test(gExp ~ gType, miniDat,
       subset = gene == "Irs4", var.equal = TRUE)$stat
0.5286494 ^ 2

###############################################
## NOTHING BELOW HERE UPDATED FOR CURRENT EXAMPLE!
## FODDER FOR SECTION OR A TAKE HOME PROBLEM?
###############################################

## what happens when we remove the outliers?
str(kDat)

kDat$rankWithinGrp <- unlist(with(kDat,
                                  tapply(obs, rv, rank)))

kDat[which(kDat$rankWithinGrp == 1), ]
##          obs rv   n sigStat nNum sigSqNum rankWithinGrp
## 28 -3.460970  x big     big   30        1             1
## 49 -2.171199  y big     big   30        1             1

kDat$outlier <- 0                       # default, codes for "not an
                                        # outlier"
kDat$outlier[kDat$rankWithinGrp == 1 & kDat$rv == "y"] <- 1 # orange
                                        # outlier is "milder"

kDat$outlier[kDat$rankWithinGrp == 1 & kDat$rv == "x"] <- 2 # blue
                                        # outlier is extreme

table(kDat$rv, kDat$outlier)
##      0  1  2
##   x 29  0  1
##   y 29  1  0

t.test(obs ~ rv, kDat)

t.test(obs ~ rv, kDat,
       subset = outlier < 2)

t.test(obs ~ rv, kDat,
       subset = outlier < 1)

## for first bit of "beyond two groups" lecture, show how two sample t
## test and anova and lm give same results
kDat <- subset(jDat, n == "big" & sigStat == "big")
str(kDat)
kDat <- droplevels(kDat)
str(kDat)

stripplot(rv ~ obs, kDat,
          type = c("p", "g"), xlab = "",
          groups = rv, jitter.data = TRUE)

t.test(obs ~ rv, kDat, var.equal = TRUE)

summary(lm(obs ~ rv, kDat))

summary(aov(obs ~ rv, kDat))

(theAvgs <- with(kDat,
                tapply(obs, rv, mean)))
theAvgs["y"] - theAvgs["x"]

t.test(obs ~ rv, kDat, var.equal = TRUE)$statistic

t.test(obs ~ rv, kDat, var.equal = TRUE)$statistic^2
