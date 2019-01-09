library(car)                            # recode() used below but I
                                        # also present alternatives,
                                        # i.e. the car package is not
                                        # mission critical

## working directory assumed to be the parent of the code directory where this
## file lies; at the time of writing, that means the 'photoRec' directory

prDes <- read.table("data/GSE4051_design_RAW.txt",
                    sep = "\t", header = T, row.names = 1)

str(prDes)
## 'data.frame':	39 obs. of  2 variables:
##  $ developmentStage   : Factor w/ 5 levels "4_weeks","E16",..: 1 1 1 1 1 1 1 1..
##  $ genotypeOrVariation: Factor w/ 2 levels "Nrl_deficient",..: 1 1 1 1 2 2 2 2..

summary(prDes)
 ## developmentStage    genotypeOrVariation
 ## 4_weeks:8        Nrl_deficient:19
 ## E16    :7        wild_type    :20
 ## P10    :8
 ## P2     :8
 ## P6     :8

## give variables shorter names
names(prDes) <- c("devStage", "gType")

## convert rownames = sample labels to proper variable
prDes$sidChar <- as.character(rownames(prDes))
rownames(prDes) <- NULL

## both factors need their levels set more rationally!

## for gType, the first level -- the reference level, by convention --
## should be wild type
(foo <- levels(prDes$gType))
## [1] "Nrl_deficient" "wild_type"

## there are various ways to change the order of the levels
## just showing some here

prDes$gType <- factor(as.character(prDes$gType),
                      levels = foo[2:1])

summary(prDes$gType)
##     wild_type Nrl_deficient
##            20            19

## the recode() function in the `car` package is handy for changing
## the factor levels
## here I'm using it to switch to shorter level names and enforce
## level order

## we could also shorten the factor levels
## the recode() function is in the car package loaded at top of script
prDes$gType <-
    recode(prDes$gType, "'wild_type'='wt'; 'Nrl_deficient'='NrlKO'",
           levels = c('wt', 'NrlKO'))

summary(prDes$gType)
##    wt NrlKO
##    20    19

## TO DO: PRESENT A RECODE/CAR FREE WAY TO DO THIS

## the developmental stage factor has even more levels and they are in
## alphabetical order, which makes no sense

## " ... 5 different developmental stages: day 16 of embryonic
## development (E16), postnatal days 2,6 and 10 (P2, P6 and p10) as
## well as 4_weeks"

levels(prDes$devStage)
## [1] "4_weeks" "E16"     "P10"     "P2"      "P6"
## currently in nonsensical alphabetical order :(

## re-make it into a factor and specify the levels in logical,
## developmental order
prDes$devStage <-
    factor(prDes$devStage,
           levels = c("E16", "P2", "P6", "P10", "4_weeks"))

## TO DO: could expand here with some cautionary tales about how wrong
## things will go if one includes a level above that was not already
## present

summary(prDes$devStage)
##     E16      P2      P6     P10 4_weeks
##       7       8       8       8       8

## I rarely use rownames
## would prefer to have sample appear as a proper variable in the
## data.frame

jExtract <- function(charVec, split, whichElem = 1) {
  return(sapply(strsplit(charVec, split = split),
                function(x) return(x[whichElem])))
}

prDes$sidNum <- as.numeric(jExtract(prDes$sidChar, split = "_", whichElem = 2))

## reorder the variables
prDes <- prDes[c("sidChar", "sidNum", "devStage", "gType")]

## reorder the rows
prDes <- prDes[with(prDes, order(devStage, gType)), ]

str(prDes)
# 'data.frame':  39 obs. of  4 variables:
# $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
# $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
# $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
# $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...

## JB private function
## peek(prDes)
##           sample devStage gType
## Sample_20     20      E16    wt
## Sample_26     26       P2    wt
## Sample_5       5       P2 NrlKO
## Sample_1       1       P6 NrlKO
## Sample_10     10       P6 NrlKO
## Sample_15     15      P10 NrlKO
## Sample_36     36  4_weeks    wt

## I consider this data.frame clean and ready for downstream analyses
## .... time to save our work

## good news: writes a plain text version of prDes to file, easy to read
## by both humans and machines (e.g. Excel)
## bad news: all our work on factor levels is not captured
write.table(prDes, file = "data/GSE4051_design.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

## good news: writes a plain text representation of prDes that *does*
## capture our hard-won factor levels
## bad news: the file itself is ugly and R specific
dput(prDes, file = "data/GSE4051_design_DPUT.txt")

## good news: writes an easy-to-reload-in-R version of prDes that
## captures our hard-won factor levels
## bad news: the file is binary and totally specific to R
saveRDS(prDes, file = "data/GSE4051_design.rds")

