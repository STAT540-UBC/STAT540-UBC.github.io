##########################################################
## code below here came from 45_lectSupp-ageQuantCovariate; when re-visited Jan 2014, I found it described as "rough, not currently used but I'm afraid to delete"
##########################################################

rFactFit <- lm(prMat ~ gType * devStage)
rQuadFullFit <- lm(prMat ~ gType * (age + I(age^2)))
rQuadRestrFit <- lm(prMat ~ age + I(age^2))

rffSumm <- summary(rFactFit)

rqfSumm <- summary(rQuadFullFit)
rqrSumm <- summary(rQuadRestrFit)

##########################################################
## testing for interaction
##########################################################

rfAnova <- anova(rQuadFullFit, rQuadRestrFit)
plot(rfAnova)                           # things whizz by ....
rfAnova                                 # first 2 fits, by default

peek(rfAnova)

##########################################################
## make it easy to find various types of 'hits'
##########################################################

## get p-values for the the main effects and interaction
pVals <- cbind(interaction = rfAnova[ , "M1:Pval"],
               rqfSumm$Coef[ , c("gTypeNrlKO:I(age^2)",
                                 "gTypeNrlKO:age",
                                 "gTypeNrlKO",
                                 "I(age^2)",
                                 "age"), "Pr(>|t|)"])
splom(~ pVals, panel = panel.smoothScatter)
pValsTall <- data.frame(pVals = as.vector(pVals),
                        what = rep(colnames(pVals),
                                   each = nrow(pVals)))
densityplot(~ pVals, pValsTall,
            group = what, auto.key = TRUE,
            plot.points = FALSE)


##########################################################
## finding some genes to feature in lecture
##########################################################

## boring
getMe <- which(rowSums(pVals < 0.05) == 0)
length(getMe)                           # 11244
set.seed(111)
(boring <- sample(getMe, size = 3))
graphIt(extractIt(boring))
print(rqfSumm, show = which(colnames(prMat) %in% names(boring)))
print(rqfSumm, show = which(colnames(prMat) == names(boring)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(boring)))

#dev.print(pdf,
#          paste0(whereLectureLives, "figs/gTypeByDevStage/stripplot-boring.pdf"),
#          width = 7, height = 7)

## age YES NrlKO NO interaction NO
getMe <- which(pVals[ , "interaction"] > 0.2 &
                 (pVals[ , "age"] < 0.05 | pVals[ , "I(age^2)"] < 0.05) &
                 (pVals[ , "gTypeNrlKO"] > 0.15 &
                    pVals[ , "gTypeNrlKO:age"] > 0.15 &
                    pVals[ , "gTypeNrlKO:I(age^2)"] > 0.15))
length(getMe)                           # 3550
set.seed(111)
(ageOnly <- sample(getMe, size = 3))
graphIt(extractIt(ageOnly))
print(rqfSumm, show = which(colnames(prMat) %in% names(ageOnly)))
print(rqfSumm, show = which(colnames(prMat) == names(ageOnly)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(ageOnly)))

## age NO NrlKO YES interaction NO
getMe <- which(pVals[ , "interaction"] > 0.2 &
                 pVals[ , "age"] > 0.20 & pVals[ , "I(age^2)"] > 0.2 &
                 pVals[ , "gTypeNrlKO"] < 0.05)
length(getMe)                           # 28
set.seed(111)
(koOnly <- sample(getMe, size = 3))
graphIt(extractIt(koOnly))
print(rqfSumm, show = which(colnames(prMat) %in% names(koOnly)))
print(rqfSumm, show = which(colnames(prMat) == names(koOnly)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(koOnly)))

## age YES NrlKO YES interaction NO
getMe <- which(pVals[ , "interaction"] > 0.2 &
                 (pVals[ , "age"] < 0.05 | pVals[ , "I(age^2)"] < 0.05) &
                 pVals[ , "gTypeNrlKO"] < 0.05)
length(getMe)                           # 90
set.seed(111)
(mainEffOnly <- sample(getMe, size = 3))
graphIt(extractIt(mainEffOnly))
print(rqfSumm, show = which(colnames(prMat) %in% names(mainEffOnly)))
print(rqfSumm, show = which(colnames(prMat) == names(mainEffOnly)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(mainEffOnly)))

## devStage YES NrlKO YES interaction YES
getMe <- which(pVals[ , "interaction"] < 0.05 &
                 (pVals[ , "age"] < 0.05 | pVals[ , "I(age^2)"] < 0.05) &
                 pVals[ , "gTypeNrlKO"] < 0.05)
length(getMe)                           # 1771
set.seed(111)
(exciting <- sample(getMe, size = 3))
graphIt(extractIt(exciting))
print(rqfSumm, show = which(colnames(prMat) %in% names(exciting)))
print(rqfSumm, show = which(colnames(prMat) == names(exciting)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(exciting)))

##########################################################
## writing results in case I want to mine this again soon
##########################################################

## not edited for this particular diff exp anal yet!

## deGtypeByDevStage <-
##     data.frame(coef(ebFitGtypeStarDevStage), # est coefs
##                rfSumm$Coef[ , , "Std. Error"], # est se's
##                rfSumm$Coef[ , , "t value"], # t stats
##                rfSumm$Coef[ , , "Pr(>|t|)"], # p vals
##                deviance(rFit)/df.residual(rFit), # resid var
##                rfSumm$FullModelStats[ , "Fstat"], # overall F test
##                rfSumm$FullModelStats[ , "Pval"],  # overall F pval
##                rfAnova[ , c("gType:Pval", "devStage:Pval", # anova
##                             "gType:devStage:Pval")]) # table stuff
## str(deGtypeByDevStage)
## peek(deGtypeByDevStage)

## yo <- c(rep(c("est", "se", "tStat", "pVal"), each = 10))
## yo2 <- c(rep(c("int", "NrlKO", "4_weeks", "NrlKO:4_weeks"), 10))
## colnames(deGtypeByDevStage) <- c(paste(yo, yo2, sep = "."),
##                                  c("sigSq", "overallFstat",
##                                    "overallFp",
##                                    "gTypeMainEffPval",
##                                    "devStageMainEffPval",
##                                    "interactionPval"))
## peek(deGtypeByDevStage)

## ## good news: writes a plain text version to file, easy to read
## ## by both humans and machines (e.g. Excel)
## ## bad news: order of factor levels not captured
## write.table(deGtypeByDevStage,
##             file = file.path(whereAmI,
##             "rmd/data/photoRec/deGtypeDevStage.txt"),
##             quote = FALSE, row.names = FALSE)

## ## good news: writes an easy-to-reload-in-R version that
## ## captures order of factor levels
## ## bad news: the file is binary and totally specific to R
## save(deGtypeByDevStage,
##      file = file.path(whereAmI,
##      "rmd/data/photoRec/deGtypeByDevStage.robj"))





## factFitted <-
##     with(prDes,
##          data.frame(gType = factor(rep(levels(gType), each = nlevels(devStage)),
##                     levels = levels(gType)),
##                     devStage = factor(levels(devStage),
##     levels(devStage))))
## factFitted$fitted <- predict(factFit, factFitted)

## jStrip <-
##     stripplot(gExp ~ devStage, jDat, group = gType, grid = TRUE,
##               panel = panel.superpose,
##               panel.groups = function(x, y , ..., group.number) {
##                   panel.stripplot(x, y, ...)
##                   panel.lines(seq_len(nlevels(prDes$devStage)),
##                               factFitted$fitted[factFitted$gType ==
##                                                 levels(factFitted$gType)[group.number]],
##                               col = trellis.par.get("superpose.line")$col[group.number])
##               })



## code below here came out of 48_lectSupp-multiple-testing.r

##########################################################
## trying to understand the utility, if any, of a 'sum to zero' parametrization
##########################################################
miniDes <- with(prDes, expand.grid(levels(gType), levels(devStage)))
names(miniDes) <- c("gType", "devStage")

## "the usual"
model.matrix(~ gType * devStage, miniDes)

## "sum to zero"
jDesMat2 <- model.matrix(~ gType * devStage, prDes,
                         contrasts.arg = list(gType = "contr.sum",
                                              devStage = "contr.sum"))
system.time(jFit2 <- lmFit(prDat, jDesMat2))
ebFit2 <- eBayes(jFit2)
colnames(coef(ebFit2))                # parameters in the linear model
grep("gType1", colnames(coef(ebFit2)))  # isolate the genotype terms
hits2 <- topTable(ebFit2, coef = grep("gType1",
                                      colnames(coef(ebFit2))),
                  n = Inf)
tail(hits2)

colnames(coef(ebFit))                 # parameters in the linear model
grep("gType", colnames(coef(ebFit)))  # isolate the genotype terms
hits3 <- topTable(ebFit,
                  coef = grep("gType", colnames(coef(ebFit))),
                  n = Inf)

head(hits2)
head(hits3)

stripplotIt(hitDat <- prepareData(tail(rownames(hits2), 6)),
            scales = list(rot = c(45, 0)))






## below was removed from 48_lectSupp-multiple-testing.r

##########################################################
## trying to understand limma manual around page 43 (different ways to
## do 2x2 factorial
##########################################################

limDat <- data.frame(strain = factor(c("wt", "wt", "mut", "mut", "mut"),
                                     levels = c("wt", "mut")),
                     tx = factor(c("untx", "stim", "untx", "stim", "stim"),
                                 levels = c("untx", "stim")))

model.matrix(~ strain + strain:tx, limDat)
##   (Intercept) strainmut strainwt:txstim strainmut:txstim
## 1           1         0               0                0
## 2           1         0               1                0
## 3           1         1               0                0
## 4           1         1               0                1
## 5           1         1               0                1

model.matrix(~ strain/tx, limDat)
##   (Intercept) strainmut strainwt:txstim strainmut:txstim
## 1           1         0               0                0
## 2           1         0               1                0
## 3           1         1               0                0
## 4           1         1               0                1
## 5           1         1               0                1

## Yep, the same result.
