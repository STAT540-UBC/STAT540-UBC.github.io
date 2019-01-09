## does several differential expression analyses

## all use built-in lm() and its ability to do multivariate regression
## no limma, no empirical bayes

## FYI: lm() returns an object of class mlm when doing multivariate regression

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
str(prDes)

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
## DEA by devStage
##########################################################

# see 20_dea-by-devStage.rmd or .html for a more chatty, educational tour
# through this analysis; this is stripped down to bare minimum

## it looks odd to pass the design as the data, but that's where the covariates 
## are to be found; in this context, it is no longer convenient (possible?) to
## keep everything neatly together in a data.frame
devStageFit <- lm(prMat ~ devStage, prDes)

## summary() will hang if you've not sourced our in-house code!
devStageSumm <- summary(devStageFit)
str(devStageSumm, max.level = 1)

## write the fit to file for use elsewhere
dput(devStageFit, "../results/dea-devStage-mlm-DPUT.txt")

## write the summary of the fit
dput(devStageSumm, "../results/dea-devStage-mlm-summary-DPUT.txt")

## write the most useful results to (a much smaller) file
devStageStats <-
  data.frame(F = (devStageSumm$FullModelStats)[ , "Fstat"],
             pVal = (devStageSumm$FullModelStats)[ , "Pval"],
             sigma = (devStageSumm$FullModelStats)[ , "Res. SE"],
             rank = rank((devStageSumm$FullModelStats)[ , "Pval"]))
write.table(devStageStats,
            file = "../results/dea-devStage-mlm-stats.txt", quote = FALSE)

rm("devStageFit", "devStageSumm", "devStageStats")
## related limma code ... parking here for now
# library(limma)
# (dmDevStage <- model.matrix(~ devStage, prDes))
# colnames(dmDevStage) <- levels(prDes$devStage)
# fitDevStage <- lmFit(prDat, dmDevStage)
# ebFitDevStage <- eBayes(fitDevStage)

##########################################################
## DEA gType * simplified version of devStage
##########################################################

# see 21_dea-by-gType-oneDevStage for a more chatty, educational tour through
# this analysis; this is stripped down to bare minimum

# simplying devStage to first and last stages, i.e. E16 and 4_weeks
prDesSimple <-
  droplevels(subset(prDes,
                    subset = devStage %in%
                      levels(devStage)[c(1, nlevels(devStage))]))
str(prDesSimple) # 15 obs. of 4 variables
prMatSimple <- prMat[prDesSimple$sidChar, ]
str(prMatSimple) # num [1:15, 1:29949]

## gType * simplified devStage
gTypeBySimpleDevStageFit <- lm(prMatSimple ~ gType * devStage, prDesSimple)

## summary() will hang if you've not sourced our in-house code!
gTypeBySimpleDevStageSumm <- summary(gTypeBySimpleDevStageFit)
str(gTypeBySimpleDevStageSumm, max.level = 1)

## write the fit to file for use elsewhere
dput(gTypeBySimpleDevStageFit,
     "../results/dea-gType-by-simpleDevStage-mlm-DPUT.txt")

## write the summary of the fit
dput(gTypeBySimpleDevStageSumm,
     "../results/dea-gType-by-simpleDevStage-mlm-summary-DPUT.txt")

## write the most useful results to (a much smaller) file
## ??do I really need/use this anywhere??
est_coefs <- t(coef(gTypeBySimpleDevStageFit))
gTypeBySimpleDevStageStats <-
  data.frame(est_coefs, # est coefs
             gTypeBySimpleDevStageSumm$Coef[ , , "Std. Error"], # est se's
             gTypeBySimpleDevStageSumm$Coef[ , , "t value"], # t stats
             gTypeBySimpleDevStageSumm$Coef[ , , "Pr(>|t|)"], # p vals
             deviance(gTypeBySimpleDevStageFit) /
               df.residual(gTypeBySimpleDevStageFit), # resid var
             gTypeBySimpleDevStageSumm$FullModelStats[ , "Fstat"], # F test
             gTypeBySimpleDevStageSumm$FullModelStats[ , "Pval"])

yo <- c(rep(c("est", "se", "tStat", "pVal"), each = 4))
yo2 <- c(rep(c("int", "NrlKO", "4_weeks", "NrlKO:4_weeks"), 4))
colnames(gTypeBySimpleDevStageStats) <-
  c(paste(yo, yo2, sep = "."), c("sigSq", "Fstat", "Fp"))

head(gTypeBySimpleDevStageStats)
str(gTypeBySimpleDevStageStats)

write.table(gTypeBySimpleDevStageStats,
            file = "../results/dea-gType-by-simpleDevStage-mlm-stats.txt",
            quote = FALSE, row.names = FALSE)

rm("prDesSimple", "prMatSimple",
   "gTypeBySimpleDevStageFit", "gTypeBySimpleDevStageSumm",
   "est_coefs", "yo", "yo2", "gTypeBySimpleDevStageStats")
## related limma code ... parking here for now
# library(limma)
# (dmGtypeStarDevStage <- model.matrix(~ gType * devStage, prDesSimple))
# colnames(dmGtypeStarDevStage) <- levels(prDesSimple$devStage)
# fitGtypeStarDevStage <- lmFit(prDatSimple, dmGtypeStarDevStage)
# ebFitGtypeStarDevStage <- eBayes(fitGtypeStarDevStage)

##########################################################
## DEA gType * devStage
##########################################################
gTypeByDevStageFit <- lm(prMat ~ gType * devStage, prDes)

## write the fit to file for use elsewhere
dput(gTypeByDevStageFit,
     "../results/dea-gType-by-devStage-mlm-DPUT.txt")

##########################################################
## devStage as factor vs. quantitative covariate
##########################################################

# given where and how I use these models, makes no sense to fit here and write
# to file; putting the code here just as a central collection of lm() commands
# we've used

## recode() is from add-on package 'car'
# prDes$age <-
#   recode(prDes$devStage,
#          "'E16'=-2; 'P2'=2; 'P6'=6; 'P10'=10; '4_weeks'=28",
#          as.factor.result = FALSE)

# rFit <- lm(prMat ~ devStage, prDes, subset = gType == "wt")
# rFit <- lm(prMat ~ age, prDes, subset = gType == "wt")
# rFit <- lm(prMat ~ age + I(age^2), prDes, subset = gType == "wt")
