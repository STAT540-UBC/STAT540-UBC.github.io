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
## finding good exemplars ... line, first try
##########################################################

## I want a gene where linear regression on age will work nicely
rFit <- lm(prMat ~ devStage, prDes, subset = gType == "wt")
rfSumm <- summary(rFit)
str(rfSumm$Coef)
## isolate effect estimates of P2, P6 and P10
str(rfSumm$Coef[ , 2:4, "Estimate"])
litmus <-
    data.frame(maxPval = apply(rfSumm$Coef[ , 2:4, "Pr(>|t|)"], 1, max),
               magnitude = apply(abs(rfSumm$Coef[ , 2:4, "Estimate"]), 1, min),
               sameSign = abs(rowSums(sign(rfSumm$Coef[ , 2:4,
                                                       "Estimate"]))) == 3,
               range = apply(rfSumm$Coef[ , 2:4, "Estimate"], 1,
                             function(z) max(z) - min(z)))
summary(litmus)
sum(getMe <-
      with(litmus,
           maxPval < 0.05 & sameSign &
             magnitude > quantile(magnitude, 0.95)))  # 1116

set.seed(999)
yo <- sample(which(getMe), size = 12)
xyplotIt(prepareData(yo))

## 2014 note: I apparently used various, undocumented methods in 2012 (or 2013?)
## to find these probesets, but they don't come up with the current code. Sigh.
## I think it is a moot point, as none of these are used as examples in the
## lecture slides. This is a dead end, so to speak.
(yo <- c("1424963_at", "1437687_x_at", "1439571_at", "1423506_a_at",
         "1425222_x_at", "1436552_at"))
xyplotIt(prepareData(yo))

##########################################################
## finding good exemplars ... line, second try
##########################################################

## I want a gene where linear regression on age will work nicely
rFit <- lm(prMat ~ age, prDes, subset = gType == "wt")
rfSumm <- summary(rFit)
str(rfSumm$Coef)
## isolate the age param ests
str(rfSumm$Coef[ , "age", "Estimate"])
litmus <- data.frame(rfSumm$Coef[ , "age", ])
names(litmus) <- c("est","se","t","pval")
summary(litmus)
sum(getMe <-
    with(litmus,
         pval < 0.05 &
         est > quantile(est, 0.95))) # 1464

set.seed(999)
yo <- sample(which(getMe), size = 12)
xyplotIt(prepareData(yo), myType = c('p', 'r'))

## interesting
## "1441811_x_at"

## 2014 note: I don't get this probeset again. This IS A probeset used in
## lecture slides.

##########################################################
## finding good exemplars ... parabola
##########################################################

## I want a gene where quadratic regression on age will work nicely
rFit <- lm(prMat ~ age + I(age^2), prDes, subset = gType == "wt")
rfSumm <- summary(rFit)
str(rfSumm$Coef)
## isolate the non-intercept param ests
str(rfSumm$Coef[ , -1, "Estimate"])
litmus <-
    data.frame(quadPval = rfSumm$Coef[ , "I(age^2)", "Pr(>|t|)"],
               quadMagnitude = abs(rfSumm$Coef[ , "I(age^2)", "Estimate"]))
summary(litmus)
sum(getMe <- with(litmus,
                  quadPval < 0.05 &
                  quadMagnitude > quantile(quadMagnitude, 0.9))) # 2714
set.seed(222)
yo <- sample(which(getMe), size = 12)
xyplotIt(prepareData(yo))

## 2014 note: although I get a very different group of probesets than the ones 
## listed below, the probeset featured in the slides -- 1427275_at -- is in the
## list AND in `yo` as defined above. Yay.

## fairly steady decrease in expression
## practically no difference in wt and NrlKO
## "1429688_at"
## "1415840_at"
## "1419740_at"
## "1457682_at"
## "1427275_at"
## "1456341_a_at"
## "1449519_at"
## "1419864_x_at"
## "1455972_x_at"

(yo <- c("1429688_at", "1415840_at", "1419740_at", "1457682_at",
         "1427275_at", "1456341_a_at", "1449519_at", "1419864_x_at",
         "1455972_x_at"))
xyplotIt(prepareData(yo))

##########################################################
## finding good exemplars ... weird
##########################################################

## I want a gene where neither linear nor quadratic regression will
## work well
rFit <- lm(prMat ~ devStage, prDes, subset = gType == "wt")
rfSumm <- summary(rFit)
str(rfSumm$Coef)
## isolate the non-intercept param ests
str(rfSumm$Coef[ , 2:4, "Estimate"])
litmus <-
  data.frame(maxPval = apply(rfSumm$Coef[ , 2:4, "Pr(>|t|)"], 1,
                             max),
             minPval = apply(rfSumm$Coef[ , 2:4, "Pr(>|t|)"], 1, min),
             magnitude = apply(abs(rfSumm$Coef[ , 2:4, "Estimate"]), 1, min),
             signSum = rowSums(sign(rfSumm$Coef[ , 2:4, "Estimate"])),
             range = apply(rfSumm$Coef[ , 2:4, "Estimate"], 1,
                           function(z) max(z) - min(z)))
summary(litmus)
sum(getMe <- with(litmus,
                  minPval < 0.025 &
                  signSum >= -1 &
                  signSum <= 1))        # 413
set.seed(222)
yo <- sample(which(getMe), size = 12)
xyplotIt(prepareData(yo))

## 2014 note: Of the two probesets listed below, I have 1456341_a_at in `yo` as computed above AND this is the probeset featured in the slides. Yay.

## interestingly weird
## "1456341_a_at"
## "1444790_at"

(yo <- c("1456341_a_at", "1444790_at"))
xyplotIt(prepareData(yo))

##########################################################
## fix my exemplars
##########################################################
(luckyGenes <- c("1427275_at", "1456341_a_at", "1441811_x_at"))
jDat <- prepareData(luckyGenes)
xyplotIt(jDat)

##########################################################
## simple linear regression on age for the exemplars
## wild type data only
##########################################################

jDat <- prepareData(luckyGenes)
jDat <- droplevels(subset(jDat, gType == "wt", select = -gType))
str(jDat)

print(jStrip <- stripplot(gExp ~ devStage | gene, jDat,
                    type = c('p', 'a'), grid = TRUE, layout = c(3, 1)))

print(jScatter <- xyplot(gExp ~ age | gene, jDat,
                   type = c('p', 'a'), grid = TRUE, layout = c(3, 1)))

print(jStrip, c(0, 0, 1, 0.5), more = TRUE)
print(jScatter, c(0, 0.5, 1, 1), more = FALSE)

## this reproduces the figure in lecture slides, up to the panel order
dev.print(pdf, "../figs/ageQuantCovariate/exemplars-wtOnly-devStageVsAge.pdf",
          width = 10, height = 7)

linFits <- with(jDat, by(jDat, gene, function(z) lm(gExp ~ age, z)))

## there are easier ways to put simple linear regression lines down
## but they don't generalize to other fits, so I'm doing it the 'hard' way
nPts <- 100
linFitted <- lapply(linFits, function(theFit) {
    newDat <-
        with(prDes,
             data.frame(age = seq(from = min(age),
                                  to = max(age), length = nPts)))
    newDat$fitted <- predict(theFit, newDat)
    newDat
})

xyplot(gExp ~ age | gene, jDat,
       grid = TRUE, layout = c(3, 1),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           with(linFitted[[packet.number()]],
                panel.lines(age, fitted))
       })

## this reproduces the figure in lecture slides, up to the panel order
dev.print(pdf, "../figs/ageQuantCovariate/exemplars-wtOnly-scatterAndLine.pdf",
          width = 10, height = 4)

summary(linFits[["1441811_x_at"]])

##########################################################
## quadratic model for the exemplars
## wild type data only
##########################################################

jDat <- prepareData(luckyGenes)
jDat <- droplevels(subset(jDat, gType == "wt", select = -gType))
str(jDat)

quadFits <- with(jDat, by(jDat, gene,
                          function(z) lm(gExp ~ age + I(age^2), z)))

## there are easier ways to put simple linear regression lines down
## but they don't generalize well, so I'm doing it the 'hard' way
nPts <- 100
quadFitted <- lapply(quadFits, function(theFit) {
    newDat <-
        with(prDes,
             data.frame(age = seq(from = min(age),
                                  to = max(age), length = nPts)))
    newDat$fitted <- predict(theFit, newDat)
    newDat
})

xyplot(gExp ~ age | gene, jDat,
       grid = TRUE, layout = c(3, 1),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           with(quadFitted[[packet.number()]],
                panel.lines(age, fitted))
       })

## this reproduces the figure in lecture slides, up to the panel order
dev.print(pdf,
          "../figs/ageQuantCovariate/exemplars-wtOnly-scatterAndParabola.pdf",
          width = 10, height = 4)

summary(quadFits[["1427275_at"]])
summary(quadFits[["1441811_x_at"]])
summary(quadFits[["1456341_a_at"]])

##########################################################
## compare models for the exemplars
## wild type data only
##########################################################

factFits <- with(jDat, by(jDat, gene, function(z) lm(gExp ~ devStage, z)))

(jGene <- luckyGenes[1]) # "1427275_at"
anova(linFits[[jGene]], quadFits[[jGene]])
# Model 1: gExp ~ age
# Model 2: gExp ~ age + I(age^2)
# Res.Df     RSS Df Sum of Sq      F   Pr(>F)    
# 1     18 17.9021                                 
# 2     17  7.0591  1    10.843 26.113 8.71e-05 ***

AIC(linFits[[jGene]], quadFits[[jGene]], factFits[[jGene]])
#                   df      AIC
# linFits[[jGene]]   3 60.54129
# quadFits[[jGene]]  4 43.92930
# factFits[[jGene]]  6 47.54810

(jGene <- luckyGenes[2]) # "1456341_a_at"
anova(linFits[[jGene]], quadFits[[jGene]])
# Model 1: gExp ~ age
# Model 2: gExp ~ age + I(age^2)
# Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
# 1     18 19.370                              
# 2     17 15.139  1    4.2308 4.7509 0.04363 *

AIC(linFits[[jGene]], quadFits[[jGene]], factFits[[jGene]])
#                   df      AIC
# linFits[[jGene]]   3 62.11743
# quadFits[[jGene]]  4 59.18864
# factFits[[jGene]]  6 48.70210


(jGene <- luckyGenes[3]) # "1441811_x_at"
anova(linFits[[jGene]], quadFits[[jGene]])
# Model 1: gExp ~ age
# Model 2: gExp ~ age + I(age^2)
# Res.Df    RSS Df Sum of Sq      F Pr(>F)  
# 1     18 3.7176                             
# 2     17 3.1120  1   0.60559 3.3081 0.0866 .

AIC(linFits[[jGene]], quadFits[[jGene]], factFits[[jGene]])
#                   df      AIC
# linFits[[jGene]]   3 29.10466
# quadFits[[jGene]]  4 27.54851
# factFits[[jGene]]  6 27.12587

for(jGene in luckyGenes) {
  print(jLinear <-
          xyplot(gExp ~ age, jDat,
                 subset = gene == jGene,
                 grid = TRUE,
                 panel = function(x, y, ...) {
                   panel.xyplot(x, y, ...)
                   with(linFitted[[jGene]],
                        panel.lines(age, fitted))
                 }))  
  
  print(jQuadratic <-
          xyplot(gExp ~ age, jDat,
                 subset = gene == jGene,
                 grid = TRUE,
                 panel = function(x, y, ...) {
                   panel.xyplot(x, y, ...)
                   with(quadFitted[[jGene]],
                        panel.lines(age, fitted))
                 }))
  
  print(jAnova <-
          xyplot(gExp ~ age, jDat,
                 subset = gene == jGene,
                 grid = TRUE, type = c('p','a')))
  
  print(jAnova, c(0, 0, 1, 0.33), more = TRUE)
  print(jQuadratic, c(0, 0.33, 1, 0.67), more = TRUE)
  print(jLinear, c(0, 0.67, 1, 1), more = FALSE)
  
  dev.print(pdf,
            paste0("../figs/ageQuantCovariate/",
                   jGene, "-wtOnly-linearVsQuad.pdf"),
            width = 6, height = 10)
}


##########################################################
## now put knockouts back in, i.e. use gType
## just the linear term
##########################################################

jDat <- prepareData(luckyGenes)
str(jDat) # 117 obs. of  7 variables

xyplot(gExp ~ age | gene, jDat,
       group = gType, grid = TRUE, layout = c(3, 1),
       type = c('p', 'a'))

dev.print(pdf,
          "../figs/ageQuantCovariate/exemplars-scatterAndMeans.pdf",
          width = 10, height = 5)

## let's focus on "1441811_x_at" for linear
jDat <- prepareData("1441811_x_at")
str(jDat)

xyplot(gExp ~ age | gene, jDat,
       group = gType, auto.key = list(x = 0.9, y = 0.9, corner = c(1,1)),
       grid = TRUE,
       type = c('p', 'r'))

dev.print(pdf,
          "../figs/ageQuantCovariate/1441811_x_at-scatterAndLines.pdf",
          width = 7, height = 5)

jFit <- lm(gExp ~ gType * age, jDat)
summary(jFit)

jFitAlt <- lm(gExp ~ gType/age - 1, jDat)
summary(jFitAlt)

coef(jFit)
coef(jFitAlt)

(contMat <- rbind(c(1, 0, 0, 0),
                  c(1, 1, 0, 0),
                  c(0, 0, 1, 0),
                  c(0, 0, 1, 1)))

cbind(coefDefault = coef(jFit),
      coefAlt = coef(jFitAlt),
      matrixResult = as.vector(contMat %*% coef(jFit)))

anova(lm(gExp ~ age, jDat), jFit)

##########################################################
## now put knockouts back in, i.e. use gType
## include quadratic term
##########################################################

## let's focus on "1427275_at" for quadratic
jDat <- prepareData("1427275_at")
str(jDat)

jFit <- lm(gExp ~ gType * (age + I(age^2)), jDat)
summary(jFit)

nPts <- 100
quadFitted <-
    with(jDat,
         data.frame(age = seq(from = min(age), to = max(age),
                    length = nPts),
                    gType = factor(rep(levels(gType), each = nPts),
                    levels(gType))))
quadFitted$fitted <- predict(jFit, quadFitted)

xyplot(gExp ~ age | gene, jDat,
       group = gType, grid = TRUE,
       panel = panel.superpose,
       panel.group = function(x, y, ..., group.number) {
           panel.xyplot(x, y, ...)
           with(quadFitted,
                panel.lines(age[gType == levels(gType)[group.number]],
                            fitted[gType == levels(gType)[group.number]],
                            col =
                            trellis.par.get("superpose.line")$col[group.number]))
       })

dev.print(pdf,
          "../figs/ageQuantCovariate/1427275_at-scatterAndParabolas.pdf",
          width = 7, height = 5)

anova(lm(gExp ~ age + I(age^2), jDat),
      lm(gExp ~ gType * (age + I(age^2)), jDat))
