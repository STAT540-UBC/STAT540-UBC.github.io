## working directory assumed to be the parent of the code directory where this
## file lies; at the time of writing, that means the 'photoRec' directory

prDat <- read.table("data/GSE4051_data_RAW.txt",
                    sep = "\t", header = T, row.names = 1)

str(prDat)
## 'data.frame':	29949 obs. of  39 variables:
##  $ Sample_35: num  7.15 9.22 10.06 8.35 8.45 ...
##  $ Sample_32: num  7.54 9.53 9.92 8.78 8.57 ...
## ...
##  $ Sample_3 : num  7.16 9.55 9.84 8.33 8.5 ...
##  $ Sample_14: num  7.09 9.56 9.88 8.57 8.59 ...

## bring the design data.frame in

## showing various options, given all the different ways we preserved
## it after cleaning

## write.table / read.table
## protect the two sample id variables from conversion to factor
prDes <- read.table("data/GSE4051_design.tsv", as.is = 1:2, header = TRUE)

str(prDes)
# 'data.frame':  39 obs. of  4 variables:
# $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
# $ sidNum  : int  20 21 22 23 16 17 6 24 25 26 ...
# $ devStage: Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
# $ gType   : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...

## PROBLEM: OUR FACTOR LEVELS ARE NOT IN THE ORDER WE WANTED!
## one solution:
prDes$devStage <- factor(prDes$devStage,
                         levels = c("E16", "P2", "P6", "P10", "4_weeks"))
prDes$gType <- factor(prDes$gType, levels = c("wt", "NrlKO"))
str(prDes)
# 'data.frame':  39 obs. of  4 variables:
# $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
# $ sidNum  : int  20 21 22 23 16 17 6 24 25 26 ...
# $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
# $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
## PROBLEM SOLVED.

## another solution is to just import from an R-native format
## dput / dget
prDes <- dget("data/GSE4051_design_DPUT.txt")
str(prDes)
# 'data.frame':  39 obs. of  4 variables:
# $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
# $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
# $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
# $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
## NOTICE the factor levels set during cleaning are preserved.

## another R-native format
## saveRDS / readRDS
prDes <- readRDS("data/GSE4051_design.rds")
str(prDes)
# 'data.frame':  39 obs. of  4 variables:
# $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
# $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
# $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
# $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
## NOTICE the factor levels set during cleaning are preserved.

## I want to rearrange the columns / variables of the gene expression
## data to match the order of the rows in the design data.frame

## see they are not the same!
cbind(design = prDes$sidChar, data = names(prDat))
## DANGER! DANGER!

## find out where the rownames of jDat appear in those of prDes
(alpha <- match(names(prDat), prDes$sidChar))

## if I order variables of prDat with alpha, will that work?
## see if these two columns are "harmonized"

## visual check   ..... NO!
cbind(names(prDat)[alpha], prDes$sidChar)

## automated check
identical(names(prDat)[alpha], rownames(prDes))
## [1] FALSE

## let's try it the other way around
beta <- match(prDes$sidChar, names(prDat))

## if I order jDat with beta, will that work?

## visual check   ..... YES!
cbind(names(prDat)[beta], prDes$sidChar)

## automated check
identical(names(prDat)[beta], prDes$sidChar)
## [1] TRUE

prDat <- prDat[, beta]

head(prDat)

## I consider this data.frame clean and ready for downstream analyses
## .... time to save our work

## there is no compelling reason to use any strategy other than
## write.table here (e.g. no factor levels to preserve)
write.table(prDat, file = "data/GSE4051_data.tsv", sep = "\t", quote = FALSE)
