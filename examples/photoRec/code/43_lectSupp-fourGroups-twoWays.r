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

##########################################################
## source in-house functions for post-processing mlm objects
##########################################################

## code written by Rick White
source("80_anova-mlm.r")

##########################################################
## simplying devStage to first and last timepoints
##########################################################
prDes <- 
  droplevels(subset(prDes,
                    subset = devStage %in%
                      levels(devStage)[c(1, nlevels(devStage))]))
str(prDes) # 15 obs. of  4 variables
prDat <- subset(prDat, select = prDes$sidChar)
prMat <- t(as.matrix(prDat))
str(prMat) # num [1:15, 1:29949]

##########################################################
## show an example
##########################################################

## I've already done differential expression analysis for gType * a
## simplified version of devStage, including only E16 and 4_weeks

## I've selected this as our example
(luckyGene <- which(rownames(prDat) == "1455695_at")) # 26861
miniDat <- data.frame(gExp = unlist(prDat[luckyGene, ]))
miniDat <- data.frame(prDes, miniDat)
miniDat$grp <- with(miniDat, interaction(gType, devStage))
str(miniDat)
with(miniDat, table(gType, devStage))
table(miniDat$grp)

stripplot(gExp ~ devStage, miniDat,
          group = gType, auto.key = TRUE,
          grid = TRUE, type = c('p', 'a'),
          jitter.data = TRUE)

dev.print(pdf, "../figs/gType-by-simpleDevStage/stripplot-twoFact.pdf",
          width = 7, height = 4)

stripplot(gExp ~ grp, miniDat,
          grid = TRUE, type = c('p', 'a'),
          jitter.data = TRUE)

dev.print(pdf, "../figs/gType-by-simpleDevStage/stripplot-grp.pdf",
          width = 7, height = 4)

## stuff I use in lecture slides
(theAvgs <- with(miniDat, tapply(gExp, grp, mean)))

with(miniDat, tapply(gExp, list(gType, devStage), mean))

twoFactFit <- lm(gExp ~ gType * devStage, miniDat)
grpFit <- lm(gExp ~ grp, miniDat)

summary(grpFit)
cbind(sampleMeans = theAvgs,
      minuRef = theAvgs - theAvgs["wt.E16"],
      grpFit = coef(grpFit))

summary(twoFactFit)
cbind(sampleMeans = theAvgs,
      minuRef = theAvgs - theAvgs["wt.E16"],
      twoFactFit = coef(twoFactFit))

theAvgs["NrlKO.4_weeks"] -
  (theAvgs["wt.E16"] +
     (theAvgs["NrlKO.E16"] - theAvgs["wt.E16"]) +
     (theAvgs["wt.4_weeks"] - theAvgs["wt.E16"]))

## take-home challenges

## another parametrization to try
summary(lm(gExp ~ gType/devStage, miniDat))

xyplot(fitted(twoFactFit) ~ fitted(grpFit))
xyplot(resid(twoFactFit) ~ resid(grpFit))

## using a contrast matrix
(contMat <- rbind(c(1, 0, 0, 0),
                  c(1, 1, 0, 0),
                  c(1, 0, 1, 0),
                  c(1, 0, 0, 1)))

cbind(matrixResult = as.vector(contMat %*% coef(grpFit)),
      cellMeansFit = coef(lm(gExp ~ 0 + grp, miniDat)))

##########################################################
## do DEA 
########################################################

## re-loading the fitted mlm model and the summary produced by in-house code
## proved too time-consuming; more practical to re-fit the model, although
## violates DRY principles

## gType * simplified devStage
gTypeBySimpleDevStageFit <- lm(prMat ~ gType * devStage, prDes)
gTypeBySimpleDevStageSumm <- summary(gTypeBySimpleDevStageFit)

##########################################################
## marshall p-values to help us find various types of 'hits'
##########################################################

## get p-values for both main effects, interaction, overall F test
pVals <-
  gTypeBySimpleDevStageSumm$Coef[ ,
                                 c("gTypeNrlKO", "devStage4_weeks", "gTypeNrlKO:devStage4_weeks"),
                                 "Pr(>|t|)"]
pVals <- cbind(pVals, gTypeBySimpleDevStageSumm$Full[ , "Pval"])
colnames(pVals) <- c("NrlKO", "4_weeks", "interaction", "overallF")
splom(~ pVals, panel = panel.smoothScatter)
pValsTall <- data.frame(pVals = as.vector(pVals),
                        what = rep(colnames(pVals), each = nrow(pVals)))
densityplot(~ pVals, pValsTall,
            group = what, auto.key = TRUE, plot.points = FALSE, n = 200)

##########################################################
## helper functions for finding and featuring genes
##########################################################

source("98_lectSupp-prepareDataFunction.r")
source("99_lectSupp-stripplotItFunction.r")

##########################################################
## designate "hits"
##########################################################

## NOTE written January 2013 and updated January 2014, for my future self: I got
## better at this task of finding interesting genes, graphing them, and printing
## some model results as I went along. So the next analysis -- DEA by gType *
## devStage -- has a superior approach for using the p-values to find hits. I
## haven't switched over here because I am [1] lazy and [2] want to keep
## the code that produces the figures actually in the lecture slides.

pValCutoff <- 0.15
pHits <- apply(pVals, 2, function(yo) yo <= pValCutoff)
table(NrlKO = pHits[ , "NrlKO"], four_weeks = pHits[ , "4_weeks"],
      interaction = pHits[ , "interaction"],
      overallF = pHits[ , "overallF"])

## boring genes
## not a hit w/r/t either main effect OR the interaction OR overall
getMe <- which(rowSums(pHits) == 0)
length(getMe) # 13315
set.seed(111)
(boring <- sample(getMe, size = 3))
stripplotIt(prepareData(boring))
dev.print(pdf, "../figs/gType-by-simpleDevStage/stripplot-boring.pdf",
          width = 7, height = 7)
print(gTypeBySimpleDevStageSumm,
      show = which(colnames(prMat) == names(boring)[2]))

## devStage YES NrlKO NO interaction NO
getMe <- which(pHits[ , "overallF"] &
                !pHits[ , "interaction"] &
                pHits[ , "4_weeks"] &
                !pHits[ , "NrlKO"])
length(getMe)                           # 3666
set.seed(111)
(devStageOnly <- sample(getMe, size = 3))
stripplotIt(prepareData(devStageOnly))
dev.print(pdf, "../figs/gType-by-simpleDevStage/stripplot-4_weeks.pdf",
          width = 7, height = 7)
print(gTypeBySimpleDevStageSumm,
      show = which(colnames(prMat) == names(devStageOnly)[2]))

## devStage NO NrlKO YES interaction NO
getMe <- which(pHits[ , "overallF"] &
               !pHits[ , "interaction"] &
               !pHits[ , "4_weeks"] &
               pHits[ , "NrlKO"])
length(getMe)                           # 1043
set.seed(111)
(NrlKoOnly <- sample(getMe, size = 3))
stripplotIt(prepareData(NrlKoOnly))
dev.print(pdf, "../figs/gType-by-simpleDevStage/stripplot-NrlKO.pdf",
          width = 7, height = 7)
print(gTypeBySimpleDevStageSumm,
      show = which(colnames(prMat) == names(NrlKoOnly)[2]))

## devStage YES NrlKO YES interaction NO
getMe <- which(pHits[ , "overallF"] &
               !pHits[ , "interaction"] &
               pHits[ , "4_weeks"] &
               pHits[ , "NrlKO"])
length(getMe)                           # 1307
set.seed(222)
(mainEffOnly <- sample(getMe, size = 3))
stripplotIt(prepareData(mainEffOnly))
dev.print(pdf, "../figs/gType-by-simpleDevStage/stripplot-mainYESinterNO.pdf",
          width = 7, height = 7)
print(gTypeBySimpleDevStageSumm,
      show = which(colnames(prMat) == names(mainEffOnly)[2]))

## devStage YES NrlKO YES interaction YES
## the randomly selected ones were not as good as some I had found
## earlier with more stringent methods
(exciting <-
   which(rownames(prDat) %in% c("1434709_at", "1458220_at", "1455695_at")))
stripplotIt(prepareData(exciting))
dev.print(pdf, "../figs/gType-by-simpleDevStage/stripplot-mainYESinterYES.pdf",
          width = 7, height = 7)
## little different here, due to hard-wiring
print(gTypeBySimpleDevStageSumm, show = exciting[2])
