Git(Hub) Intro
========================================================
A quick introduction to Git(Hub) is available [here](githubPres.pdf). We are going to cover this topic further in the next seminar. You can also visit Prof Jenny Bryan's [STAT545A webpage](https://stat545-ubc.github.io/git01_git-install.html) for further details on installation. 

Explore a small gene expression dataset
========================================================



### Optional getting started advice

*Ignore if you don't need this bit of support.*

This is one in a series of tutorials in which we work with the `photoRec` dataset. Now is the time to make sure you are working in an appropriate directory on your computer, probably through the use of an [RStudio project](https://stat545-ubc.github.io/block002_hello-r-workspace-wd-project.html). To ensure a clean slate, you may wish to clean out your workspace and restart R (both available from the RStudio Session menu, among other methods). Confirm that the new R process has the desired working directory, for example, with the `getwd()` command or by glancing at the top of RStudio's Console pane.

Open a new R script (in RStudio, File > New > R Script). Develop and run your code from there (recommended) or periodicially copy "good" commands from the history. In due course, save this script with a name ending in `.r` or `.R`, __containing no spaces or other funny stuff__, and evoking "seminar01" and "data exploration".

See [R basics, workspace and working directory, RStudio projects](https://stat545-ubc.github.io/block002_hello-r-workspace-wd-project.html) if you need some help developing your basic R mojo.

### Obtain the data file `GSE4051_MINI.txt`

We will work with a small data excerpt. Please save this file locally, for example, in the directory associated with your RStudio Project:

  * <https://github.com/STAT540-UBC/STAT540-UBC.github.io/blob/master/examples/photoRec/data/GSE4051_MINI.tsv>
  
### The `photoRec` data

See the [`photoRec` README](https://github.com/STAT540-UBC/STAT540-UBC.github.io/tree/master/examples/photoRec).

Gene expression data from photo receptor cells in mice. Various developmental stages, two genotypes.

The excerpt we work with here contains gene expression data for 3 genes. Their arcane names have been replaced with 3 randomly selected Pokemon attacks: crabHammer, eggBomb, and poisonFang.

### Load `GSE4051_MINI.txt`



We load the mini gene expression dataset into a `data.frame`, which is the preferred way to handle tabular data in R.

If you have not already done so, now is a good time to work through [Basic care and feeding of data in R](https://stat545-ubc.github.io/block006_care-feeding-data.html). The exercises below will test whether your basic understanding of `data.frames` is developed enough to explore the mini gene expression data we work with here.

If you are new to R, you'll want to work through [R objects (beyond data.frames) and indexing](https://stat545-ubc.github.io/block004_basic-r-objects.html) at some in the near future as well.

Read the data from a local file in the current working directory and use `str()` to get a basic overview of the resulting `prDat` object:


```r
prDat <- read.table("GSE4051_MINI.txt", header = TRUE, row.names = 1)
```

Does the `read.table()` command fail for you? If so, have you taken control of your working directory? Have you placed `GSE4051_MINI.txt` in that directory or edited the `read.table()` command to reflect its location? Mastering the working directory and filepaths is a surprisingly important part of the learning curve. Submit commands like these to do some troubleshooting:


```r
getwd()      # is working directory what you think it should be?
list.files() # do see GSE4051_MINI.txt sitting there?!?
setwd()      # use with moderation
```

Upon successful import, `str()` the object as a sanity check:

```r
str(prDat)
```

```
## 'data.frame':	39 obs. of  6 variables:
##  $ sample    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```

Exercise: Why do we use the `header =` and `row.names =` arguments above upon import? Submit the command *without* these arguments and note any difference in the result. Form the habit of reading error messages *carefully* and working the problem. Mastering the arguments of `read.table()` and friends is time well spent.




### Basic exploration of `prDat`

Use R to answer these questions. You could create a script with these questions as comments, then insert actual code to get the answers. You could use RStudio's "Compile Notebook" command to create a nice little report!

How many rows are there? Hint: `nrow()`, `dim()`.

How many columns or variables are there? Hint: `ncol()`, `length()`, `dim()`.

Inspect the first few observations or the last few or a random sample. Hint: `head()`, `tail()`, `x[i, j]` combined with `sample()`.

What does row correspond to -- different genes or different mice?

What are the variable names? Hint: `names()`, `dimnames()`.

What "flavor" is each variable, i.e. numeric, character, factor? Hint: `str()`.

For `sample`, do a sanity check that each integer between 1 and the number of rows in the dataset occurs exactly once. Hint: `a:b`,  `seq()`, `seq_len()`, `sort()`, `table()`, `==`, `all()`, `all.equal()`, `identical()`.

For each factor variable, what are the levels? Hint: `levels()`, `str()`.

How many observations do we have for each level of `devStage`? For `gType`? Hint: `summary()`, `table()`.

Perform a cross-tabulation of `devStage` and `gType`. Hint: `table()`.

If you had to take a wild guess, what do you think the *intended* experimental design was? What actually happened in real life?

For each quantitative variable, what are the extremes? How about average or median? Hint: `min()`, `max()`, `range()`, `summary()`, `fivenum()`, `mean()`, `median()`, `quantile()`.



### Indexing and subsetting `prDat`

You'll definitely need the knowledge in [R objects (beyond data.frames) and indexing](https://stat545-ubc.github.io/block004_basic-r-objects.html) now, so go work through it if you need to.

Create a new data.frame called `weeDat` only containing observations for which expression of poisonFang is above 7.5.

For how many observations poisonFang > 7.5? How do they break down by genotype and developmental stage? 



Print the observations with row names "Sample_16" and "Sample_38" to screen, showing only the 3 gene expression variables.

Which samples have expression of eggBomb less than the 0.10 quantile?



### Bonus

Try solving some of the previous questions using `dplyr` package. 

### Answers!

Code for the exercises is shown in the source document that makes this page, but is hidden here. Go read [the R Markdown source]("https://github.com/STAT540-UBC/STAT540-UBC.github.io/blob/master/sm01b_gitIntro-basic-data-exploration.rmd") if you want to see it.
