library(car)          # recode()
library(lattice)      # stripplot() in stripplotIt()
                      # xyplot() in xyplotIt()

##########################################################
## source in-house functions for post-processing mlm objects
##########################################################

## code written by Rick White
source("80_anova-mlm.r")

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
## prepare data for multivariate regression
##########################################################

## from lm() documentation: "If response is a matrix a linear model is fitted separately by least-squares to each column of the matrix."

## responses must be a matrix, one *column* per probeset
prMat <- t(as.matrix(prDat))

##########################################################
## adding quantitative version of devStage
##########################################################
## recode() is from add-on package 'car'
prDes$age <-
  recode(prDes$devStage,
         "'E16'=-2; 'P2'=2; 'P6'=6; 'P10'=10; '4_weeks'=28",
         as.factor.result = FALSE)
head(prDes)
str(prDes)

##########################################################
## helper functions for finding and featuring genes
##########################################################

source("98_lectSupp-prepareDataFunction.r")
source("99_lectSupp-plotItFunction.r")

##########################################################
## re-fitting some of the models we've worked with
##########################################################

twAnova <- lm(prMat ~ gType * devStage, prDes)
dsOnly <- lm(prMat ~ devStage, prDes)
twAnovaNoInter <- lm(prMat ~ gType + devStage, prDes)
quadBoth <- lm(prMat ~ gType * (age + I(age^2)), prDes)
quadAge <- lm(prMat ~ age + I(age^2), prDes)

twSumm <- summary(twAnova)
dsSumm <- summary(dsOnly)
twNiSumm <- summary(twAnovaNoInter)
qBothSumm <- summary(quadBoth)
qAgeSumm <- summary(quadAge)

##########################################################
## error variance
##########################################################
str(twSumm, max.level = 1)
str(twSumm$Full)
str(twSumm$Full[ , "Res. SE"])

stdDevn <-
    data.frame(twAnova = twSumm$Full[ , "Res. SE"],
               dsOnly = dsSumm$Full[ , "Res. SE"],
               quadBoth = qBothSumm$Full[ , "Res. SE"],
               quadAge = qAgeSumm$Full[ , "Res. SE"])
summary(stdDevn)
sdTall <-
    data.frame(fittedModel =
               factor(rep(names(stdDevn), each = nrow(stdDevn)),
                      levels = names(stdDevn)),
               estResSd = unlist(stdDevn))

splom(stdDevn,
      prepanel.limits = function(x) c(0, 2),
      panel = function(x, y, ... ) {
          panel.smoothScatter(x, y, ...)
          panel.abline(a = 0, b = 1, col = "orange")
      })

dev.print(pdf, "../figs/enMasse/splom-errorStdDevn.pdf",
          width = 7, height = 7)

densityplot(~ estResSd, sdTall,
            group = fittedModel, auto.key = TRUE,
            plot.points = FALSE, n = 300)

dev.print(pdf, "../figs/enMasse/densityplot-errorStdDevn.pdf",
          width = 7, height = 7)

densityplot(~ estResSd, sdTall,
            group = fittedModel, auto.key = TRUE,
            plot.points = FALSE, n = 300, xlim = c(0, 0.3))

dev.print(pdf, "../figs/enMasse/densityplot-errorStdDevn-zoomed.pdf",
          width = 7, height = 7)

sum(highVar <- with(stdDevn,
                    twAnova > quantile(twAnova, 0.95))) # 1498

sum(lowVar <- with(stdDevn,
                    twAnova < quantile(twAnova, 0.05))) # 1498

set.seed(111)
xyplotIt(prepareData(yo <- c(sample(which(highVar), size = 4),
                             sample(which(lowVar), size = 4))),
         layout = c(4, 2))
dev.print(pdf,
          "../figs/enMasse/scatterplot-lowHighVar.pdf",
          width = 7, height = 7)

##########################################################
## interaction between gType and devStage, ANOVA
##########################################################
huh <- twSumm$Full[ , "Pval"]
foo <- anova(twAnova, twAnovaNoInter)[ , "M1:Pval"]
str(foo)

densityplot(~ foo, plot.points = FALSE, n = 300, xlab = "p-value")

dev.print(pdf,
          "../figs/enMasse/densityplot-interactionYesAnova.pdf",
          width = 7, height = 7)

sum(interactionYES <- foo < 0.005 & huh < 0.005) # 1344
sum(interactionNO <- foo > 0.9 & huh < 0.005)    # 155

set.seed(333)
stripplotIt(prepareData(yo <- c(sample(which(interactionYES), size = 4),
                                sample(which(interactionNO), size = 4))),
         layout = c(4, 2))

dev.print(pdf,
          "../figs/enMasse/stripplot-interactionYesNoAnova.pdf",
          width = 7, height = 7)

##########################################################
## does gType matter? ANOVA
##########################################################
huh <- twSumm$Full[ , "Pval"]
foo <- anova(twAnova, dsOnly)[ , "M1:Pval"]
str(foo)

densityplot(~ foo, plot.points = FALSE, n = 300, xlab = "p-value")

dev.print(pdf,
          "../figs/enMasse/densityplot-gTypeMattersAnova.pdf",
          width = 7, height = 7)

sum(gTypeYES <- foo < 0.005 & huh < 0.005) # 2629

set.seed(333)
stripplotIt(prepareData(yo <- sample(which(gTypeYES), size = 8)),
            layout = c(4, 2))

dev.print(pdf,
          "../figs/enMasse/stripplot-gTypeYesAnova.pdf",
          width = 7, height = 7)

##########################################################
## overall significance of model
##########################################################
tw <- twSumm$Full[ , "Pval"]
qb <- qBothSumm$Full[ , "Pval"]

densityplot(~ tw + qb, plot.points = FALSE, n = 300, xlab = "p-value",
            auto.key = TRUE)

dev.print(pdf,
          "../figs/enMasse/densityplot-overallFAnovaVsQuad.pdf",
          width = 7, height = 7)

xyplot(tw ~ qb,
       panel = function(x, y, ...) {
           panel.smoothScatter(x, y, nbin = 200, raster = TRUE, ...)
           panel.abline(a = 0, b = 1, col = "orange")
           panel.smooth(x, y, col = "pink")
       })

dev.print(pdf,
          "../figs/enMasse/stripplot-gTypeYesAnova.pdf",
          width = 7, height = 7)





