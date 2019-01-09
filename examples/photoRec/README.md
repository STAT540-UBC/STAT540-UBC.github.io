The `photoRec` dataset
=================================================================

The aim of the study was to "generate gene expression profiles of purified photoreceptors at distinct developmental stages and from different genetic backgrounds". The experimental units were mice and the microarray platform was Affymetrix mouse genomic expression array 430 2.0.

For more information on this study, please refer to the 2006 publication:

http://www.ncbi.nlm.nih.gov/pubmed/16505381

The data is also directly accessible from GEO:

http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4051

TO DO: can we track down provenance of the files we start with? And/or can we start over with data from a link above and preserve the entire process?

There are two main sources of information in the `data` directory.

1. The gene expression data itself. The `GSE4051_data_RAW.txt` file contains expression values of 29949 probes from photoreceptor cells in 39 mice samples. A "cleaned" version, where the columns = variables have been rearranged rationally, is given in `GSE4051_data.tsv`. See the script `code/02_cleanData.R` for details.

2. The metadata file `GSE4051_design_RAW.txt` describes the experimental condition for each sample. Gene expression was studied at 5 different developmental stages: day 16 of embryonic development (E16), postnatal days 2,6 and 10 (P2, P6 and p10) as well as 4 weeks (`4_weeks`). Each of these 5 experimental conditions was studied in wild type mice and Nrl knockout mice. A "cleaned" version, with variables renamed and row order rationalized w/r/t developmental stage, is given in `GSE4051_design.tsv`. The same info is preserved in two R-specific formats in files `GSE4051_design_DPUT.txt` and `GSE4051_design.rds`. These files have an advantage over the plain text delimited file because, upon import, the levels of the developmental stage factor will be in chronological order (not alphabetical). See the script `code/01_cleanDesign.R` for details.

There are several derived datasets, created from processing the above.

`GSE4051_MINI.txt` holds data for 3 randomly selected probesets, renamed for fun, transposed into convenient column variables and stored together in a data.frame with the experimental condition information for each sample. See the script `code/03_createMiniDataset.R` for how it was created.

TO DO: discuss differential expression analysis

TO DO: use diff exp analysis to pick 3 more interesting probesets for a new mini dataset

### Reading the raw data and design

WARNING: It is your responsibility to make sure the working directory is set to where these files live or to edit paths accordingly below!

__Raw__ data and design:


```r
prDat <- read.table("data/GSE4051_data_RAW.txt", sep = "\t", header = T, row.names = 1)
str(prDat, max.level = 0)
```

```
## 'data.frame':	29949 obs. of  39 variables:
```

```r
prDes <- read.table("data/GSE4051_design_RAW.txt", sep = "\t", header = T, row.names = 1)
str(prDes)
```

```
## 'data.frame':	39 obs. of  2 variables:
##  $ developmentStage   : Factor w/ 5 levels "4_weeks","E16",..: 1 1 1 1 1 1 1 1 2 2 ...
##  $ genotypeOrVariation: Factor w/ 2 levels "Nrl_deficient",..: 1 1 1 1 2 2 2 2 1 1 ...
```


Do columns of the data match the rows of the design?


```r
head(cbind(names(prDat), rownames(prDes)))
```

```
##      [,1]        [,2]       
## [1,] "Sample_35" "Sample_11"
## [2,] "Sample_32" "Sample_12"
## [3,] "Sample_34" "Sample_2" 
## [4,] "Sample_33" "Sample_9" 
## [5,] "Sample_28" "Sample_36"
## [6,] "Sample_29" "Sample_37"
```


Heck, no! That would be too easy and transparent. Watch out for such things!


### Reading the cleaned data and design

__Cleaned__ data and design (saved in various formats):


```r
prDat <- read.table("data/GSE4051_data.tsv")
str(prDat, max.level = 0)
```

```
## 'data.frame':	29949 obs. of  39 variables:
```

```r
prDes <- read.table("data/GSE4051_design.tsv", header = TRUE, as.is = 1)
str(prDes)
```

```
## 'data.frame':	39 obs. of  4 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType   : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
```

```r
head(cbind(names(prDat), prDes$sidChar))
```

```
##      [,1]        [,2]       
## [1,] "Sample_20" "Sample_20"
## [2,] "Sample_21" "Sample_21"
## [3,] "Sample_22" "Sample_22"
## [4,] "Sample_23" "Sample_23"
## [5,] "Sample_16" "Sample_16"
## [6,] "Sample_17" "Sample_17"
```

```r
identical(names(prDat), prDes$sidChar)
```

```
## [1] TRUE
```


In the above case, note that the factor levels for `devStage` and `gType` may not be as you want. Wild type ('wt') is not the reference level and the developmental stages are not in chronological order. Set explicitly or import the design from a format that preserves factor levels, i.e. from a cleaned design written to file via `dput()` or `saveRDS()`.

### Reading the design with sane factor levels


```r
prDes <- dget("data/GSE4051_design_DPUT.txt")
str(prDes)
```

```
## 'data.frame':	39 obs. of  4 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
```

```r
prDes <- readRDS("data/GSE4051_design.rds")
str(prDes)
```

```
## 'data.frame':	39 obs. of  4 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
```


### Loading and getting to know the mini dataset:

The usual problem with factor level order occurs here:


```r
kDat <- read.delim("data/GSE4051_MINI.tsv")
str(kDat)
```

```
## 'data.frame':	39 obs. of  7 variables:
##  $ sidChar   : Factor w/ 39 levels "Sample_1","Sample_10",..: 13 14 15 16 8 9 36 17 18 19 ...
##  $ sidNum    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```


Importing from these R-specific formats preserves factor levels:


```r
kDat <- dget("data/GSE4051_MINI_DPUT.txt")
str(kDat)
```

```
## 'data.frame':	39 obs. of  7 variables:
##  $ sidChar   : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum    : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType     : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```

```r
kDat <- readRDS("data/GSE4051_MINI.rds")
str(kDat)
```

```
## 'data.frame':	39 obs. of  7 variables:
##  $ sidChar   : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum    : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType     : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```

