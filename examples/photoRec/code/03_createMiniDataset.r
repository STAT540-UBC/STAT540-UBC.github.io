## working directory assumed to be the parent of the code directory where this
## file lies; at the time of writing, that means the 'photoRec' directory

## design data.frame
prDes <- readRDS("data/GSE4051_design.rds")
str(prDes)
# 'data.frame':  39 obs. of  4 variables:
# $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
# $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
# $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
# $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...

## gene expression data.frame
prDat <- read.table("data/GSE4051_data.tsv")
str(prDat, max.level = 0)
## 'data.frame':	29949 obs. of  39 variables:
##  $ Sample_20: num  7.24 9.48 10.01 8.36 8.59 ...
##  $ Sample_21: num  7.41 10.02 10.04 8.37 8.62 ...
## ...
##  $ Sample_2 : num  7.35 9.66 9.91 8.4 8.37 ...
##  $ Sample_9 : num  7.32 9.8 9.85 8.4 8.46 ...

## pick a small number of probesets to bring over into a mini-dataset

## set.seed() is how to make something random repeatable
set.seed(5782)
(getMe <- sample(1:nrow(prDat), size = 3))
## [1]  2668 20534 17554

## I really wanted to demonstrate the use of subset() but the subset
## argument expects a logical, which is awkward for me here ...
jDat <- prDat[getMe, ]
str(jDat)
## 'data.frame':	3 obs. of  39 variables:
##  $ Sample_20: num  10.22 7.46 7.37
##  $ Sample_21: num  10.02 6.89 7.18

## transpose this data so that expression data for each probeset /
## gene will ultimately be a variable in our data.frame

## renaming the probesets from the Pokemon move list to make them more
## catchy!
## http://pokemondb.net/move/all

jDat <- data.frame(t(jDat))
names(jDat) <- c("crabHammer", "eggBomb", "poisonFang")
str(jDat)
## 'data.frame':	39 obs. of  3 variables:
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...

## reconfirm that the row order of prDes and the column order of prDat
## and, therefore, jDat are SAME

## visual check .... looks good
cbind(prDes$sidChar, rownames(jDat))

## automated check .... YES
identical(prDes$sidChar, rownames(jDat))
## [1] TRUE

## merge the info on samples from prDes and on gene expression from
## jDat

kDat <- data.frame(prDes, jDat)
str(kDat)
# 'data.frame':  39 obs. of  7 variables:
# $ sidChar   : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
# $ sidNum    : num  20 21 22 23 16 17 6 24 25 26 ...
# $ devStage  : Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
# $ gType     : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
# $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
# $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
# $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...

## let's write this minidataset to file for future use

## good news: writes a plain text version of prDes to file, easy read
## by both humans and machines (e.g. Excel)
## bad news: all our work on factor levels is not captured
write.table(kDat, file = "data/GSE4051_MINI.tsv", quote = FALSE,
            sep = "\t", row.names = FALSE)

## good news: writes a plain text representation of prDes that *does*
## capture our hard-won factor levels
## bad news: the file itself is ugly and R specific
dput(kDat, file = "data/GSE4051_MINI_DPUT.txt")

## good news: writes an easy-to-reload-in-R version of prDes that
## captures our hard-won factor levels
## bad news: the file is binary and totally specific to R
saveRDS(kDat, file = "data/GSE4051_MINI.rds")
