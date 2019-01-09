library(limma)
library(hexbin)

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
## helper functions for finding and featuring genes
##########################################################

source("98_lectSupp-prepareDataFunction.r")
source("99_lectSupp-plotItFunction.r")

##########################################################
## limma
##########################################################
jDesMat <- model.matrix(~ gType * devStage, prDes)
jFit <- lmFit(prDat, jDesMat)
ebFit <- eBayes(jFit)
str(ebFit, max.level = 3)
summary(ebFit)

topTable(ebFit, coef = "gTypeNrlKO")

colnames(coef(ebFit))                 # parameters in the linear model
grep("gType", colnames(coef(ebFit)))  # isolate the genotype terms
(gTypeHits <- topTable(ebFit,
                       coef = grep("gType", colnames(coef(ebFit)))))

stripplotIt(prepareData(head(rownames(gTypeHits), 3)),
            scales = list(rot = c(45, 0)))

## figure in lecture slides is a masked version of a figure made in
## 47_lectSupp-limma.r

##########################################################
## experimenting with contrasts.fit and decideTests for the multiple
## testing lecture
##########################################################

## get my head straight with just 1 probeset
stripplotIt(tinyDat <- prepareData("1448904_at"))

dev.print(pdf,
          "../figs/limma/stripplot-1448904_at.pdf",
          width = 7, height = 5)

tinyDat$grp <- with(tinyDat, interaction(gType, devStage))
summary(lm(gExp ~ grp + 0, tinyDat))

## set up design matrix to get "cell means" fit for devStage, within gType
prDes$grp <- with(prDes, interaction(gType, devStage))
kDesMat <- model.matrix(~ grp + 0, prDes)

## fitting two-way ANOVA to all probesets at once
kFit <- lmFit(prDat, kDesMat)
colnames(coef(kFit))                  # parameters in the linear model
                                      # = columns of the design matrix
(cont.matrix <- makeContrasts(
    wt.P10vsE16 = grpwt.P10 - grpwt.E16,
    NrlKO.P10vsE16 = grpNrlKO.P10 - grpNrlKO.E16,
    levels = kDesMat))
kFitCont <- contrasts.fit(kFit, cont.matrix)
kebFit <- eBayes(kFitCont)

cutoff <- 1e-4

## use topTable to get separate hits for each of the contrasts
wtHits <-
    topTable(kebFit, coef = "wt.P10vsE16", p.value = cutoff, n = Inf)
NrlKOHits <-
    topTable(kebFit, coef = "NrlKO.P10vsE16", p.value = cutoff, n = Inf)

## same topTable -- same order -- but don't enforce a cutoff
wtHitsAll <- topTable(kebFit, coef = "wt.P10vsE16", n = Inf)
NrlKOHitsAll <- topTable(kebFit, coef = "NrlKO.P10vsE16", n = Inf)

## same topTable -- original row order -- but don't enforce a cutoff
wtAllUnsrtd <- topTable(kebFit, coef = "wt.P10vsE16", n = Inf, sort = "none")
NrlKOAllUnsrtd <- topTable(kebFit, coef = "NrlKO.P10vsE16", n = Inf,
                           sort = "none")

nrow(wtHitsAll)                         # yeah 29949 as expected
## make sure order is same
head(wtHits, 4)
head(wtHitsAll, 4)

nrow(NrlKOHitsAll)                      # yeah 29949 as expected
## make sure order is same
head(NrlKOHits, 4)
head(NrlKOHitsAll, 4)

## how many "hits" when we consider the contrasts separately?
nrow(wtHits)                            # 349
nrow(NrlKOHits)                         # 196

nrow(wtHits) + nrow(NrlKOHits)          # 545 total if no overlap
sepHits <- union(rownames(wtHits), rownames(NrlKOHits))
length(sepHits)                         # 482 actual total

## what are the associated p-value cutoffs?
wtHits[nrow(wtHits), ] # 1.162691e-06
NrlKOHits[nrow(NrlKOHits), ] # 6.409428e-07

## show the top hit for wt and NrlKO
stripplotIt(tinyDat <- prepareData(c("1425171_at", "1425530_a_at")))

dev.print(pdf,
          "../figs/limma/stripplot-topHitsE16vsP10ByBtype_at.pdf",
          width = 7, height = 5)

## use decideTests, optionally with method = "global"
kRes <- decideTests(kebFit, p.value = cutoff, method = "separate")
kResGlobal <- decideTests(kebFit, p.value = cutoff, method = "global")

str(kRes)
str(kResGlobal)

## decideTests returns a matrix, one row per probeset, one column per contrast
## elements are -1 = statistically significant and negative
##               0 = not statistically significant
##               1 = statistically significant and positive

summary(kRes)
##    wt.P10vsE16 NrlKO.P10vsE16
## -1         256             52
## 0        29600          29753
## 1           93            144

summary(kResGlobal)
##    wt.P10vsE16 NrlKO.P10vsE16
## -1         246             61
## 0        29615          29735
## 1           88            153

## verifying that the separate hits obtained via two topTable calls
## and via one decideTests call are the same
identical(sort(rownames(wtHits)),
          sort(rownames(kRes)[which(kRes[ , "wt.P10vsE16"] != 0)]))
identical(sort(rownames(NrlKOHits)),
          sort(rownames(kRes)[which(kRes[ , "NrlKO.P10vsE16"] != 0)]))

## how many "global" hits are there?
length(kHitsGlobal <- which(rowSums(abs(kResGlobal)) > 0)) # 485

head(kResGlobal[kHitsGlobal, ])

## how many "hits" when we consider the contrasts separately?
length(which(rowSums(abs(kRes)) > 0)) # 482

## I would have expected the "global" hit list to be *shorter* than
## the separate but apparently my intuition is faulty
## https://stat.ethz.ch/pipermail/bioconductor/2008-July/023593.html
## http://thread.gmane.org/gmane.science.biology.informatics.conductor/17847

## how do these "global" hist relate to the "separate" hits?

## first we've got to parse the global hits more finely and track if
## they are hits w.r.t. wt or NrlKO
length(globHitsWt <-
       names(which(abs(kResGlobal[ , "wt.P10vsE16"]) == 1))) # 334
length(globHitsNrlKO <-
       names(which(abs(kResGlobal[ , "NrlKO.P10vsE16"]) == 1))) # 214

## overlap between separate and global results for wt
sum(abs(kResGlobal[ ,"wt.P10vsE16"]))
vennCounts(cbind(separate = wtAllUnsrtd$adj.P.Val < cutoff,
                 global = kResGlobal[ ,"wt.P10vsE16"]),
           include = "both")
#   separate global Counts
# 1        0      0  29600
# 2        0      1      0
# 3        1      0     15
# 4        1      1    334

## so in this case, the separate p-value adjustment DOES lead to more hits

## overlap between separate and global results for NrlKO
sum(abs(kResGlobal[ ,"NrlKO.P10vsE16"]))
vennCounts(cbind(separate = NrlKOAllUnsrtd$adj.P.Val < cutoff,
                 global = kResGlobal[ ,"NrlKO.P10vsE16"]),
           include = "both")
#   separate global Counts
# 1        0      0  29735
# 2        0      1     18
# 3        1      0      0
# 4        1      1    196

## but in this case, the separate p-value adjustment leads to fewer hits

## painstaking checks to make sure I understand what's going on

## the "wt.P10vsE16" contrast
wt <- data.frame(ID = I(rownames(wtHitsAll)),
                 adj.P.Val = wtHitsAll$adj.P.Val)
wt$sepHit <- wt$adj.P.Val < cutoff
wt$globHit <- wt$ID %in% globHitsWt

## checking that sepHit and globHit is TRUEs followed by FALSEs
xyplot(sepHit + globHit ~ seq_len(nrow(wt)), wt,
       outer = TRUE, xlim = c(0, 500), type = "l")
with(wt,
     c(max(which(sepHit)), min(which(!sepHit)),
       max(which(globHit)), min(which(!globHit))))
## [1] 349 350 334 335

## I want to see the global hit status for NrlKOHits
ko <- data.frame(ID = I(rownames(NrlKOHitsAll)),
                 adj.P.Val = NrlKOHitsAll$adj.P.Val)
ko$sepHit <- ko$adj.P.Val < cutoff
ko$globHit <- ko$ID %in% globHitsNrlKO

## checking that sepHit and globHit is TRUEs followed by FALSEs
xyplot(sepHit + globHit ~ seq_len(nrow(ko)), ko,
       outer = TRUE, xlim = c(0, 500), type = "l")
with(ko,
     c(max(which(sepHit)), min(which(!sepHit)),
       max(which(globHit)), min(which(!globHit))))
## [1] 196 197 214 215

## combine all hit info together so I can choose some genes to
## spotlight
names(wt) <- paste0("wt.", names(wt))
names(ko) <- paste0("ko.", names(ko))
jRes <- merge(wt, ko, by.x = "wt.ID", by.y = "ko.ID")
str(jRes) # 29949 obs. of  7 variables
with(jRes,
     table(wt.sepHit, wt.globHit, ko.sepHit, ko.globHit))

## picking probes to feature
set.seed(123)
(boring <- with(jRes,
                sample(which(!wt.sepHit & !wt.globHit &
                             !ko.sepHit & !ko.globHit),
                       size = 1)))      # 8620

(yesYesYes <- with(jRes,
                   sample(which(wt.sepHit & wt.globHit &
                                ko.sepHit & ko.globHit),
                          size = 1)))      # 22171

(kodisc <- with(jRes,
                sample(which(!wt.sepHit &
                               !ko.sepHit & ko.globHit),
                       size = 1)))      # 15969

(wtdisc <- with(jRes,
                sample(which(wt.sepHit & !wt.globHit &
                !ko.globHit & !ko.globHit),
                       size = 1)))      # 28080

(getMe <- rownames(prDat)[c(boring, yesYesYes, NrlKOdisc, wtdisc)])
## "1427795_s_at" "1448904_at"   "1438294_at"   "1457558_at" 
hitDat <- prepareData(getMe)

stripplot(gExp ~ devStage | gene, hitDat,
          group = gType, subset = devStage %in% c("E16", "P10"),
          jitter.data = TRUE,
          auto.key = TRUE, type = c('p', 'a'), grid = TRUE,
          scales = list(x = list(rot = c(45, 0)),
          y = list(relation = "free")))

dev.print(pdf,
          "../figs/limma/stripplot-fourExamplesE16vsP10ByBtype_at.pdf",
          width = 10, height = 6)

## I can't say I learned very much from that at all. :(
