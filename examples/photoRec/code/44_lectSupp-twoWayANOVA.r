library(lattice)
library(car)                            # Anova()

## load design data.frame
prDes <- readRDS("../data/GSE4051_design.rds")
str(prDes)

## load gene expression data.frame
prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)
## 'data.frame':  29949 obs. of  39 variables:

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

##########################################################
## source in-house functions for post-processing mlm objects
##########################################################

## code written by Rick White
source("80_anova-mlm.r")

##########################################################
## helper functions for finding and featuring genes
##########################################################

source("98_lectSupp-prepareDataFunction.r")
source("99_lectSupp-stripplotItFunction.r") # not used here!

##########################################################
## do DEA 
########################################################

## re-loading the fitted mlm model and the summary produced by in-house code
## proved too time-consuming; more practical to re-fit the model, although
## violates DRY principles

## gType * devStage
## summary() will hang if you've not sourced our in-house code!
prMat <- t(as.matrix(prDat))
gTypeByDevStageFit <- lm(prMat ~ gType * devStage, prDes)
gTypeByDevStageSumm <- summary(gTypeByDevStageFit)

##########################################################
## create ANOVA-table type of results for these models
##########################################################

gTypeByDevStageAnova <- anova(gTypeByDevStageFit)
gTypeByDevStageAnova                   # first 2 fits, by default
head(gTypeByDevStageAnova)

##########################################################
## marshall p-values to help find various types of 'hits'
##########################################################

## get p-values for the the main effects and interaction
pVals <- gTypeByDevStageAnova[ , c("gType:Pval", "devStage:Pval",
                                   "gType:devStage:Pval")]
colnames(pVals) <- c("gType", "devStage", "gType:devStage")
splom(~ pVals, panel = panel.smoothScatter)
pValsTall <- data.frame(pVals = as.vector(pVals),
                        what = rep(colnames(pVals),
                        each = nrow(pVals)))
densityplot(~ pVals, pValsTall,
            group = what, auto.key = TRUE, plot.points = FALSE, n = 200)

pValCutoff <- 0.15
pHits <- apply(pVals, 2, function(yo) yo <= pValCutoff)
table(gType = pHits[ , "gType"],
      devStage = pHits[ , "devStage"],
      interaction = pHits[ , "gType:devStage"])

##########################################################
## helper functions for finding and featuring genes
##########################################################

## I aimed to get these functions out of lecture support files and to reuse
## them. I would expect to use the function defined in
## 99_lectSupp-stripplotItFunction.r here. But there are some differences, so
## will go ahead and define and use graphIt() here.

graphIt <- function(jDat, yFree = FALSE) {
    thePlot <- stripplot(gExp ~ devStage | gene, jDat,
                         type = c('p', 'a'), grid = TRUE,
                         group = gType, auto.key = TRUE,
                         jitter.data = TRUE)
    if(yFree) {
        thePlot <- update(thePlot,
                          scales = list(y = list(relation = "free")))
    }
    print(thePlot)
}

##########################################################
## finding genes to feature in lecture
##########################################################

## boring
getMe <- which(rowSums(pHits) == 0)
length(getMe) # 4948
set.seed(111)
(boring <- sample(getMe, size = 3))
graphIt(prepareData(boring))
dev.print(pdf,
          "../figs/gType-by-devStage/stripplot-boring.pdf",
          width = 7, height = 7)
print(gTypeByDevStageSumm, show = which(colnames(prMat) == names(boring)[2]))
print(gTypeByDevStageAnova, show = which(colnames(prMat) == names(boring)[2]))

## devStage YES NrlKO NO interaction NO
getMe <- which(pVals[ , "devStage"] < 0.0005 & # demand more here!
               !pHits[ , "gType"] &
               !pHits[ , "gType:devStage"])
length(getMe)                           # 1243
## set.seed(555)
## (devStageOnly <- sample(getMe, size = 12))
## through various random looks, I settled on these 3:
foo <- c("1456219_at", "1455007_s_at", "1445613_at")
(devStageOnly <- which(rownames(prDat) %in% foo))
names(devStageOnly) <- foo
graphIt(prepareData(devStageOnly))
dev.print(pdf,
          "../figs/gType-by-devStage/stripplot-devStage.pdf",
          width = 7, height = 7)
print(gTypeByDevStageSumm,
      show = which(colnames(prMat) == names(devStageOnly)[2]))
print(gTypeByDevStageAnova,
      show = which(colnames(prMat) == names(devStageOnly)[2]))

## devStage NO NrlKO YES interaction NO
getMe <- which(pVals[ , "devStage"] > 0.8 &
               pVals[ , "gType"]  < 0.05 & # demand more here!
               !pHits[ , "gType:devStage"])
length(getMe)                           # 775
set.seed(111)
(gTypeOnly <- sample(getMe, size = 12))
## through various random looks, I settled on these 3:
foo <- c("1452142_at", "1429993_s_at", "1456468_x_at")
(gTypeOnly <- which(rownames(prDat) %in% foo))
names(gTypeOnly) <- foo
graphIt(prepareData(gTypeOnly))
dev.print(pdf,
          "../figs/gType-by-devStage/stripplot-gType.pdf",
          width = 7, height = 7)
print(gTypeByDevStageSumm, show = which(colnames(prMat) == names(gTypeOnly)[2]))
print(gTypeByDevStageAnova, show = which(colnames(prMat) == names(gTypeOnly)[2]))
print(gTypeByDevStageAnova, show = which(colnames(prMat) == names(gTypeOnly)[1]))
print(gTypeByDevStageAnova, show = which(colnames(prMat) == names(gTypeOnly)[3]))

## devStage YES NrlKO YES interaction NO
getMe <- which(pVals[ , "devStage"] < 0.05 &
               pVals[ , "gType"]  < 0.05 &
               pVals[ , "gType:devStage"] > 0.7)
length(getMe)                           # 515
set.seed(111)
(mainEffOnly <- sample(getMe, size = 12))
## through various random looks, I settled on these:
foo <- c("1438815_at", "1438786_a_at", "1448690_at")
## more candidates:
## "1433716_x_at", "1434564_at", , "1427347_s_at", "1435176_a_at"
(mainEffOnly <- which(rownames(prDat) %in% foo))
names(mainEffOnly) <- foo
graphIt(prepareData(mainEffOnly))
dev.print(pdf,
          "../figs/gType-by-devStage/stripplot-mainYESinterNO.pdf",
          width = 7, height = 7)
print(gTypeByDevStageSumm, show = which(colnames(prMat) == names(mainEffOnly)[2]))
print(gTypeByDevStageAnova, show = which(colnames(prMat) %in% names(mainEffOnly)))

## devStage YES NrlKO YES interaction YES
getMe <- which(pVals[ , "devStage"] < 0.05 &
               pVals[ , "gType"]  < 0.05 &
               pVals[ , "gType:devStage"] < 0.05)
length(getMe)                           # 1570
set.seed(111)
(exciting <- sample(getMe, size = 12))
## through various random looks, I settled on these:
foo <- c("1451763_at", "1419655_at", "1429225_at")
## more candidates:
## "1452601_a_at", "1433577_at"
(exciting <- which(rownames(prDat) %in% foo))
names(exciting) <- foo
graphIt(prepareData(exciting))
dev.print(pdf,
          "../figs/gType-by-devStage/stripplot-mainYESinterYES.pdf",
          width = 7, height = 7)
print(gTypeByDevStageSumm, show = which(colnames(prMat) == names(exciting)[2]))
print(gTypeByDevStageAnova, show = which(colnames(prMat) %in% names(exciting)))

##########################################################
## detailed look at a boring gene, to explore lm versus ANOVA conventions in
## reporting; also get into the sum of squares problem for unbalanced designs :(
##########################################################

## http://goanna.cs.rmit.edu.au/~fscholer/anova.php

## https://stat.ethz.ch/pipermail/r-help/2011-April/273863.html

## http://psych.wisc.edu/moore/Rpdf/610-R12_UnbalFactorialBetw.pdf

## referred to multiple times ... server not responding last I checked
## http://prometheus.scp.rochester.edu/zlab/sites/default/files/InteractionsAndTypesOfSS.pdf

(luckyGene <- boring[2])
jDat <- data.frame(gType = prDes$gType,
                   devStage = prDes$devStage,
                   gExp = unlist(prDat[luckyGene, ]))
lmRes <- lm(gExp ~ gType * devStage, jDat)
summary(lmRes)

deviance(lmRes) # residual sum of squares; see bottom row of anova table

anova(lmRes)
##                Df  Sum Sq  Mean Sq F value Pr(>F)
## gType           1 0.02985 0.029848  1.5657 0.2208
## devStage        4 0.10365 0.025914  1.3594 0.2722
## gType:devStage  4 0.07191 0.017977  0.9430 0.4532
## Residuals      29 0.55283 0.019063

lmResAlt <- lm(gExp ~ devStage * gType, jDat) # reverse order of covariates

anova(lmResAlt)
##                Df  Sum Sq  Mean Sq F value Pr(>F)
## devStage        4 0.10328 0.025819  1.3544 0.2739
## gType           1 0.03022 0.030225  1.5855 0.2180
## devStage:gType  4 0.07191 0.017977  0.9430 0.4532
## Residuals      29 0.55283 0.019063

## evidence of the slightly weird and undesirable nature of the Type I
## sums of squares approach ... see how the entries for gType and
## devStage depend on the order in which they appear in the model
## formula?

## Anova() from the car package uses Type II sums of squares
Anova(lmRes)
##                 Sum Sq Df F value Pr(>F)
## gType          0.03022  1  1.5855 0.2180
## devStage       0.10365  4  1.3594 0.2722
## gType:devStage 0.07191  4  0.9430 0.4532
## Residuals      0.55283 29

Anova(lmResAlt)
##                 Sum Sq Df F value Pr(>F)
## devStage       0.10365  4  1.3594 0.2722
## gType          0.03022  1  1.5855 0.2180
## devStage:gType 0.07191  4  0.9430 0.4532
## Residuals      0.55283 29

## see, these Type II sums of squares results are invariant to the
## order in which the factors are given in the model

Anova(lm(gExp ~  gType * devStage, jDat,
         contrasts=list(gType = contr.sum, devStage = contr.sum)),
      type = 3)
##                 Sum Sq Df    F value Pr(>F)
## (Intercept)    2753.78  1 1.4446e+05 <2e-16 ***
## gType             0.03  1 1.6862e+00 0.2043
## devStage          0.10  4 1.3679e+00 0.2693
## gType:devStage    0.07  4 9.4300e-01 0.4532
## Residuals         0.55 29

aovRes <- aov(gExp ~ gType * devStage, jDat)
summary(aovRes)
##                Df Sum Sq Mean Sq F value Pr(>F)
## gType           1 0.0298 0.02985   1.566  0.221
## devStage        4 0.1037 0.02591   1.359  0.272
## gType:devStage  4 0.0719 0.01798   0.943  0.453
## Residuals      29 0.5528 0.01906
## same table as anova(lmRes)

aovResAlt <- aov(gExp ~ devStage * gType, jDat)
summary(aovResAlt)
##                Df Sum Sq Mean Sq F value Pr(>F)
## devStage        4 0.1033 0.02582   1.354  0.274
## gType           1 0.0302 0.03022   1.585  0.218
## devStage:gType  4 0.0719 0.01798   0.943  0.453
## Residuals      29 0.5528 0.01906
## same table as anova(lmResAlt)

anova(lm(gExp ~ gType * devStage, jDat))
anova(lm(gExp ~ devStage * gType, jDat))

Anova(lm(gExp ~ gType * devStage, jDat))
Anova(lm(gExp ~ devStage * gType, jDat))

## building some F stats by hand

lmBig <- lm(gExp ~ gType * devStage, jDat)
deviance(lmBig)

lmSmall <- lm(gExp ~ gType + devStage, jDat)
deviance(lmSmall)

deviance(lmSmall) - deviance(lmBig)
lmBig$rank - lmSmall$rank

(numerator <- (deviance(lmSmall) - deviance(lmBig)) /
 (lmBig$rank - lmSmall$rank))

(denominator <- deviance(lmBig) / lmBig$df.residual)

(Fstat <- numerator / denominator)

pf(Fstat, lmBig$rank - lmSmall$rank, lmBig$df.residual, lower.tail =
   FALSE)

## Type I, also called "sequential" sum of squares:
## SS(A) for factor A.
## SS(B | A) for factor B.
## SS(AB | B, A) for interaction AB.

## This tests the main effect of factor A, followed by the main effect
## of factor B after the main effect of A, followed by the interaction
## effect AB after the main effects.

lmFull <- lm(gExp ~ gType * devStage, jDat)
lmNoInter <- lm(gExp ~ gType + devStage, jDat)
lmGtype <- lm(gExp ~ gType, jDat)
lmInter <- lm(gExp ~ 1, jDat)

myFtest <- function(lmBig, lmSmall) {
    numerator <- (deviance(lmSmall) - deviance(lmBig)) /
        (lmBig$rank - lmSmall$rank)
    denominator <- deviance(lmBig) / lmBig$df.residual
    Fstat <- numerator / denominator
    pVal <- pf(Fstat, lmBig$rank - lmSmall$rank, lmBig$df.residual,
    lower.tail = FALSE)
    return(c(num = numerator, den = denominator, Fstat, pVal))
}

myOtherTest <- function(lmHuge, lmBig, lmSmall) {
    numerator <- (deviance(lmSmall) - deviance(lmBig)) /
        (lmBig$rank - lmSmall$rank)
    denominator <- deviance(lmHuge) / lmHuge$df.residual
    Fstat <- numerator / denominator
    pVal <- pf(Fstat, lmBig$rank - lmSmall$rank, lmHuge$df.residual,
    lower.tail = FALSE)
    return(c(num = numerator, den = denominator, Fstat, pVal))
}

deviance(lmFull)
deviance(lmFull) / lmFull$df.resid

deviance(lmNoInter) - deviance(lmFull)
(deviance(lmNoInter) - deviance(lmFull)) / (lmFull$rank -
                                            lmNoInter$rank)
myFtest(lmFull, lmNoInter)

deviance(lmGtype) - deviance(lmNoInter)
(deviance(lmGtype) - deviance(lmNoInter)) / (lmNoInter$rank -
                                             lmGtype$rank)
myFtest(lmNoInter, lmGtype)
myOtherTest(lmFull, lmNoInter, lmGtype)
myOtherTest(lmFull, lmGtype, lmInter)


anova(lmFull, lmNoInter, lmGtype, lmInter)
anova(lmInter, lmGtype, lmNoInter, lmFull)

anova(lmInter, lmGtype)

