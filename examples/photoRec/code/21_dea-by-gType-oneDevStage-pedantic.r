## This should be converted to R Markdown as I did for DEA w/r/t developmental
## stage in 20_dea-by-devStage.rmd. But this really covers lots of the same
## ground, so doesn't seem mission critical now (January 2014). Just leaving in
## the repo.

library(limma)
library(lattice)

## code Rick developed for working with mlm objects, ie fitted multivariate
## regression
source("80_anova-mlm.r")

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
## simplying devStage to first and last timepoints
##########################################################
keepMe <- with(prDes,
               devStage %in%
               levels(devStage)[c(1, nlevels(devStage))])
prDesSimple <- droplevels(subset(prDes, keepMe))
str(prDesSimple)
prDatSimple <- subset(prDat, select = prDesSimple$sidChar)
str(prDatSimple)

##########################################################
## limma
##########################################################
(dmGtypeStarDevStage <- model.matrix(~ gType * devStage, prDesSimple))

fitGtypeStarDevStage <- lmFit(prDatSimple, dmGtypeStarDevStage)
fitGtypeStarDevStage

## matrix, one row per probeset, four columns
## (Intercept) = wt, E16
## gTypeNrlKO = knockout main effect
## devStage4_weeks = 4_weeks main effect
## gTypeNrlKO:devStage4_weeks = interaction term
head(fitGtypeStarDevStage$coef)
str(fitGtypeStarDevStage$coef)

## square root of diagonal elements of (t(X) %*% X)^{-1}
## same for each probeset, of course, in our case
head(fitGtypeStarDevStage$stdev.unscaled)
sqrt(diag(solve(t(dmGtypeStarDevStage) %*% dmGtypeStarDevStage))) # yep

## estimate residual std dev'n
head(fitGtypeStarDevStage$sigma)

## residual degress of freedom
## n - #params = 15 - 4 = 11
## same for each probeset, of course, in our case
head(fitGtypeStarDevStage$df.residual)

## top of page 62 in 2013 Jan 04 limma User's Guide
## "The ordinary t-statistics can be recovered by ..."
head(fitGtypeStarDevStage$coef/fitGtypeStarDevStage$stdev.unscaled/fitGtypeStarDevStage$sigma)
## one row per probeset, one col for each estimated parameter
##              (Intercept) gTypeNrlKO devStage4_weeks gTypeNrlKO:devStage4_weeks
## 1415670_at     109.16960  0.9672920      -0.2645388                  0.0399375
## 1415671_at      30.15966 -1.4905696      -1.7668220                  2.1175915
## 1415672_at      42.71904 -1.0002450      -1.6205821                  1.4849438
## 1415673_at     116.29623  1.0527907       1.8585579                 -1.6979205
## 1415674_a_at   162.80779  0.7165277      -3.5614565                  0.5166830
## 1415675_at     183.87810 -0.4031307      -2.4245328                  0.9893521

## not populated yet, i.e. not before eBayes()
head(fitGtypeStarDevStage$F)
fitGtypeStarDevStage$s2.prior
fitGtypeStarDevStage$df.prior
head(fitGtypeStarDevStage$var.prior)
head(fitGtypeStarDevStage$t)
## NULL

ebFitGtypeStarDevStage <- eBayes(fitGtypeStarDevStage)

## Yes these are the same
head(cbind(coef(ebFitGtypeStarDevStage),
           coef(fitGtypeStarDevStage)))

## this is the mean of the inverse Chisquare prior for the
## gene-specific variances
ebFitGtypeStarDevStage$s2.prior
## [1] 0.07478654

## this is the degrees of freedom associated with the inverse
## ChiSquare prior for gene-specific variances
ebFitGtypeStarDevStage$df.prior
## [1] 3.367102

##########################################################
## two-way ANOVA "by hand" for one probeset
##########################################################
## working with first probeset
jDat <- data.frame(gType = prDesSimple$gType,
                   devStage = prDesSimple$devStage,
                   gExp = unlist(prDatSimple[1, ]))
lmRes <- lm(gExp ~ gType * devStage, jDat)
summary(lmRes)

(foo <- with(jDat, tapply(gExp, list(gType, devStage), mean)))
foo - foo[1, 1]
coef(lmRes)
fitGtypeStarDevStage$coef[1, ]
## agree agree agree agree!!!!

summary(lmRes)
fitGtypeStarDevStage$sigma[1] * fitGtypeStarDevStage$stdev.unscaled[1, ] # std errors
round((fitGtypeStarDevStage$coef/fitGtypeStarDevStage$stdev.unscaled/
       fitGtypeStarDevStage$sigma)[1, ], 5) # t stats
## agree agree agree agree!!!!

##########################################################
## homegrown mlm code
##########################################################

## responses must be a matrix, one *row* per response
prMatSimple <- t(as.matrix(prDatSimple))
gType <- prDesSimple$gType              # lesser of two evils
devStage <- prDesSimple$devStage        # lesser of two evils
rFit <- lm(prMatSimple ~ gType * devStage)

## summary(rFit) hangs unless you've loaded our homegrown code
## source()ed at top!
rfSumm <- summary(rFit)
rfSumm                                  # by default, info on first 2
                                        # fitted models
print(rfSumm, show = c(2, 4555, 29403)) # show gives more flexibility

str(rfSumm, max.level = 1)
str(rfSumm$Coef)

##########################################################
## extracting interesting bits and verifying agreement, when
## appropriate, between different methods
##########################################################

## estimated coefficients ... manual loop
coefName <- "(Intercept)"
coefName <- "gTypeNrlKO"
coefName <- "devStage4_weeks"
coefName <- "gTypeNrlKO:devStage4_weeks"

foo <-
    cbind(mlm = coef(rFit)[coefName, ], # we want 1 row!
          limmaPre = coef(fitGtypeStarDevStage)[ , coefName],
          limma = coef(ebFitGtypeStarDevStage)[ , coefName])
splom(~ foo, panel = panel.smoothScatter)
max(apply(foo, 1, function(z) max(z) - min(z))) # 0

## estimated standard error ... manual loop
coefName <- "(Intercept)"
coefName <- "gTypeNrlKO"
coefName <- "devStage4_weeks"
coefName <- "gTypeNrlKO:devStage4_weeks"

foo <-
    cbind(mlm = rfSumm$Coef[ , coefName, "Std. Error"],
          limmaPre = fitGtypeStarDevStage$sigma *
             fitGtypeStarDevStage$stdev.unscaled[ , coefName],
          limma = sqrt(ebFitGtypeStarDevStage$s2.post) *
             ebFitGtypeStarDevStage$stdev.unscaled[ , coefName])
splom(~ foo, panel = panel.smoothScatter)
max(apply(foo, 1, function(z) max(z) - min(z)))
## [1] 0.1191471 for intercept
## [1] 0.1820002 for NrlKO
## [1] 0.1684995 for 4_weeks
## [1] 0.2480245 for interaction effect

## t statistics ... manual loop
coefName <- "(Intercept)"
coefName <- "gTypeNrlKO"
coefName <- "devStage4_weeks"
coefName <- "gTypeNrlKO:devStage4_weeks"

foo <-
    cbind(mlm = rfSumm$Coef[ , coefName, "t value"],
          limmaPre = coef(fitGtypeStarDevStage)[ , coefName] /
              (fitGtypeStarDevStage$sigma *
               fitGtypeStarDevStage$stdev.unscaled[ , coefName]),
          limma = ebFitGtypeStarDevStage$t[ , coefName])
splom(~ foo, panel = panel.smoothScatter)
max(apply(foo, 1, function(z) max(z) - min(z)))
## [1] 330.4443 for intercept !!!!
## [1] 4.766263 for NrlKO
## [1] 7.784197 for 4_weeks
## [1] 4.512867 for interaction effect

## p-values ... manual loop
coefName <- "(Intercept)"
coefName <- "gTypeNrlKO"
coefName <- "devStage4_weeks"
coefName <- "gTypeNrlKO:devStage4_weeks"

foo <-
    cbind(mlm = log(rfSumm$Coef[ , coefName, "Pr(>|t|)"]),
          limma = log(ebFitGtypeStarDevStage$p.value[ , coefName]))
splom(~ foo, panel = panel.smoothScatter)

## estimated residual variance
foo <-
    cbind(mlm1 = deviance(rFit)/df.residual(rFit),
          mlm2 = rfSumm$FullModelStats[ , "Res. SE"] ^ 2,
          limmaPre = fitGtypeStarDevStage$sigma ^ 2,
          limma = ebFitGtypeStarDevStage$s2.post)
splom(~ foo, panel = panel.smoothScatter)
max(apply(foo, 1, function(z) max(z) - min(z)))
## [1] 0.8713952

## overall F test ... clearly not comparable
foo <-
    cbind(mlm = rfSumm$FullModelStats[ , "Fstat"],
          limma = topTableF(ebFitGtypeStarDevStage, sort.by = "none",
          adjust.method = "none", n = Inf)[ , "F"])
splom(~ foo, panel = panel.smoothScatter)
max(apply(foo, 1, function(z) max(z) - min(z)))
## [1] 31760.08 !!!
## I think limma is truly comparing to no model at all and I don't
## feel like fixing this right now

## overall F test's p-value
foo <-
    cbind(mlm = rfSumm$FullModelStats[ , "Pval"],
          limma1 = topTableF(ebFitGtypeStarDevStage, sort.by = "none",
                             adjust.method = "none", n = Inf)[ ,
                                                     "P.Value"],
          limma2 = ebFitGtypeStarDevStage$F.p.value)
splom(~ foo, panel = panel.smoothScatter)
max(apply(foo, 1, function(z) max(z) - min(z)))
## [1] 0.9999895
## As above, I think limma is truly comparing to no model at all and I don't
## feel like fixing this right now
max(apply(foo[ , 2:3], 1, function(z) max(z) - min(z)))
## 0

##########################################################
## writing results in case I want to mine this again soon
##########################################################

deGtypeBy4weeks <-
    data.frame(coef(ebFitGtypeStarDevStage)[ , ], # est coefs
               rfSumm$Coef[ , , "Std. Error"], # est se's
               rfSumm$Coef[ , , "t value"], # t stats
               rfSumm$Coef[ , , "Pr(>|t|)"], # p vals
               deviance(rFit)/df.residual(rFit), # resid var
               rfSumm$FullModelStats[ , "Fstat"], # F test
               rfSumm$FullModelStats[ , "Pval"])
peek(deGtypeBy4weeks)
str(deGtypeBy4weeks)

yo <- c(rep(c("est", "se", "tStat", "pVal"), each = 4))
yo2 <- c(rep(c("int", "NrlKO", "4_weeks", "NrlKO:4_weeks"), 4))
colnames(deGtypeBy4weeks) <- c(paste(yo, yo2, sep = "."), c("sigSq", "Fstat", "Fp"))

# write.table(deGtypeBy4weeks,
#             file = "../results/deGtypeBy4weeks.txt",
#             quote = FALSE, row.names = FALSE)
