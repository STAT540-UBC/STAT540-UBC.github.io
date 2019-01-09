Play with probability distributions and simulations
========================================================

This seminar was originally written by Dr.Jenny Bryan and modified by Alice Zhu.

There are 9 exercises embedded for you to practice. It is recommended for you to attempt at least half of them to familiarize yourself with probability-related functions





### Optional getting started advice

*Ignore if you don't need this bit of support.*

This is a self-contained activity, so you won't need to download and import any data. But you will still want to make a deliberate choice about where to do your work. Pick a directory on your computer, perhaps within an existing directory and/or RStudio project for STAT540, where you will save the script you develop and that will be the working directory for these activities.

To ensure a clean slate, it is **highly recommended** that you clean out your workspace and restart R (both available from the RStudio Session menu, among other methods). You can use `ls()` to confirm you have a clean workspace and `getwd()` to confirm your working directory. You can get the same info from the RStudio IDE by looking at the Environment tab (by default, in the same pane as the History) and at the "title bar" of the Console.

Open a new R script (in RStudio, File > New > R Script). Develop and run your code from there (recommended) or periodicially copy "good" commands from the history. Before long, save this script in your chosen directory with a name ending in `.r` or `.R`, __containing no spaces or other funny stuff__, and evoking "seminar02" and "probability".

See [R basics, workspace and working directory, RStudio projects](https://stat545-ubc.github.io/block002_hello-r-workspace-wd-project.html) if you need some help developing your basic R mojo.

### Prerequisites

These activities will utilize stuff you learned here:

  * [Basics of R/RStudio, workspaces, and projects](https://stat545-ubc.github.io/block002_hello-r-workspace-wd-project.html), borrowed from STAT 545A.
  * [Basic care and feeding of data in R](https://stat545-ubc.github.io/block006_care-feeding-data.html), borrowed from STAT 545A.
  * [R objects (beyond data.frames) and indexing](https://stat545-ubc.github.io/block004_basic-r-objects.html), borrowed from STAT 545A.

### Generate some data

R has built-in functions for generating data from many continuous and discrete distributions. Here's how we generate data from the standard normal distribution:


```r
rnorm(n = 10)
```

```
##  [1] -0.57701 -1.60965  0.69499  2.41824 -0.36879  0.88706  0.83328
##  [8]  0.56332  0.05497 -1.66002
```


The "random data generating" function have a common scheme: the first letter is `r` for *random*, then comes something that evokes the family of the distribution, e.g. `norm` for the normal distribution.

#### Sidebar: calling functions

Call up the [help on `rnorm()`](http://stat.ethz.ch/R-manual/R-patched/library/stats/html/Normal.html) and focus on its usage and arguments. Here's the usage:

```
rnorm(n, mean = 0, sd = 1)
```

This tells us that we must specify `n`, which is the number of observations, since it has no default value. There are two other arguments specifying the mean, `mean`, and the standard deviation, `sd`. We see that each of these has a default value. So, if the user does not specify the mean and the variance, she will get observations from the standard normal. In general, if you are happy with the default value of an argument, just omit it from your function call.

Recall how R resolves function arguments. You can always specify in `name = value` form, as we did with `rnorm(n = 10)`. But if you do not, _R attempts to resolve by position_. So a readable, equivalent command is this:


```r
rnorm(10)
```

```
##  [1] -1.1363 -1.4743  1.7301  1.1048  1.2447 -0.9635 -1.0839  1.6628
##  [9] -1.7273  0.6833
```


For functions you call often and/or where the first 1 or 2 arguments are "obvious", it's ok to rely on this "resolve by position" functionality.  Otherwise, it's best to use `name = value`.

Exercises: Generate some data from a different normal distribution, i.e. where the mean is not 0 and/or the standard deviation is not 1. Make the new parameter values extreme, so the data will be obviously different. Experiment with giving the arguments in different order, accepting defaults, using "resolve by position" and `name = value`.

### Repeatably random

It can be incredibly helpful to make your "random" work *repeatable*. This may sound like an oxymoron, but it's not. R, like almost any other system you might use, is actually generating *pseudo-random* numbers and, if you exert control over the generator's initial state, you can obtain the same (pseudo-)random numbers all over again. We use the function `set.seed()` for this:


```r
set.seed(540)
rnorm(10)
```

```
##  [1] -1.6145  0.7340  2.0947  0.5009  1.7355 -2.0584 -1.0339  1.4914
##  [9] -0.5068 -0.7001
```


So, in practice, at the top of a script that has some random work in it, you might want to use `set.seed()` near the very top. As for the seed itself, it doesn't matter what it is. Just pick a positive integer.

Main contexts where I am glad to have made things repeatably random:

  * In teaching, when I am using simulated data. By setting seed, I can make sure my random data sets are the same and it's easy to keep various tables, figures, etc. "in sync" without constantly remaking all of them.
  * In data analysis, when a complicated pipeline has some random components. My collaborators would kill me if our "hit list" of genes changed subtly every time I revisited part of the analysis. By setting seed, I can lock things down.
  * In simulation studies, when I am studying the performance of some novel procedure. When things go really wrong, i.e. something that should be "impossible" happens, by setting seed, it is easy for me to get the offending random dataset in front of my eyeballs and begin the debugging process.
  
### Generate lots of data, the R way

If we want to explore some of the probability facts and theorems discussed in lecture, we need to generate many samples. The natural object to store these in is a *matrix*. Unlike data analysis, where the `data.frame` reigns supreme, simulation studies often utilize matrices. Why? Because simulated data is so homogeneous! In our current example, it's all numeric. In contrast, real datasets often have variables of different flavors.

Now is a good time to run through [R objects (beyond data.frames) and indexing](https://stat545-ubc.github.io/block004_basic-r-objects.html), if you have not already done so.

Here is the best, most R native way to generate many samples of a common size:


```r
n <- 10
B <- 4

x <- matrix(rnorm(n * B), nrow = n)

str(x)
```

```
##  num [1:10, 1:4] -0.00874 0.1083 0.806 -1.5045 -0.41763 ...
```

```r
head(x)
```

```
##           [,1]    [,2]    [,3]     [,4]
## [1,] -0.008744 -0.8413 -1.2752  1.13066
## [2,]  0.108303 -1.0581  0.3885 -1.46676
## [3,]  0.806000 -0.2552 -0.1950  0.07621
## [4,] -1.504498  0.3213  1.4259 -0.11602
## [5,] -0.417628  0.7739 -0.6323 -0.73436
## [6,]  0.480787 -0.1932  0.6684  0.55222
```


We used `rnorm()` to generate all the data we needed -- enough for 4 samples, each of size 10. We provided that directly as the first argument to `matrix()`, which was resolved by position to be the `data` for populating the matrix. Then we specified that the matrix should have `n` rows, implying that each column of the matrix represents one sample.
  
### Generate lots of data, the awkward ways

Many people new to R would generate the samples a different way. I'll show some of these approaches and discuss why they're not ideal.

#### The `for` loop

Programmers from other languages might turn to a `for` loop:


```r
x <- matrix(0, nrow = n, ncol = B)
for(j in 1:B) {
  x[ , j] <- rnorm(n)
}
```


We still get what we need but with more code, less readable code, and, in larger examples, this uses more RAM and compute time. In real life, *programmer time and psychic energy* is your main limiting factor, so it is good to develop a style that leads to less code and more elegant code. A secondary consideration is efficient use of computing resources. Because of R's insistence that it be able to "back out" of `for` loops if things go wrong, they can consume alot of memory, but that really only matters when things get really large. You will gradually learn the R native ways for iterating; we don't use `for` loops very often, although we have them.

#### Gluing stuff together

People who haven't previously programmed at all are prone to "copy-and-paste" parties like this:


```r
sample1 <- rnorm(n)
sample2 <- rnorm(n)
sample3 <- rnorm(n)
sample4 <- rnorm(n)
x <- cbind(sample1, sample2, sample3, sample4)
```


We still get what we need but this does not scale well and it's easy to make mistakes. I made 4 typos just typing that code. Repetitive code like this is a *sure sign* you're attacking from the wrong angle. But the lack of scalability is the main problem here. What happens when you scale up to `B = 1000` samples???

#### Gluing stuff together inside a `for` loop

With this hybrid of the previous two approaches, we have reached the absolute nadir. This is called "growing an object" and, in R, it's usually a big no-no. It works, but it's memory inefficient and opaque.


```r
x <- rnorm(n)
for(j in 1:(B - 1)) {
  x <- cbind(x, rnorm(n))
}
```


We create the first sample and then, inside the `for` loop, keep appending new samples as columns in the matrix. If you find yourself using constructor functions, such as `c()`, `cbind()`, or `rbind()`, recursively like this, you should take a step back and consider a different attack.

### Compute sample means and explore

Let's re-create our samples, with the preferred code. I also apply row and column names. This is somewhat artificial in a simulation study, but I do it here to demonstrate how nicely informative row and column names propagate through downstream computations.


```r
x <- matrix(rnorm(n * B), nrow = n)
rownames(x) <- sprintf("obs%02d", 1:n)
colnames(x) <- sprintf("samp%02d", 1:B)
x
```

```
##        samp01   samp02   samp03    samp04
## obs01 -0.6625  3.42818 -1.31936 -0.100429
## obs02  0.7743 -0.52823  0.70759  0.001006
## obs03 -1.5737  0.47803 -0.42372 -0.552930
## obs04 -0.2075  1.34923 -0.75977 -0.209284
## obs05 -0.6732 -0.06336  0.12161  0.292004
## obs06 -0.4987 -0.23567 -1.33881  0.240990
## obs07  0.3075 -0.67465 -0.06423  0.952912
## obs08  0.2241  0.47812  1.54898  0.643902
## obs09 -0.5773  0.30075  0.51883 -1.265136
## obs10  1.1534  0.57059  2.36283 -0.779013
```

```r
dim(x)
```

```
## [1] 10  4
```


Let's take the mean of each sample. Here we show how to operate on one sample and then two approaches for doing the same for all samples (look Ma, no top-level `for` loops!).


```r
mean(x[ , 2]) # sample mean for 2nd sample
```

```
## [1] 0.5103
```

```r
colMeans(x)
```

```
##  samp01  samp02  samp03  samp04 
## -0.1734  0.5103  0.1354 -0.0776
```

```r
apply(x, 2, mean)
```

```
##  samp01  samp02  samp03  samp04 
## -0.1734  0.5103  0.1354 -0.0776
```


Computing sums and averages for the rows and columns of a matrix is such a common task that there are special purpose built-in functions for this, such as `colMeans()`. They are super-fast, but you will only notice this on truly large datasets. If you want to apply an arbitrary function to rows or columns of a matrix, you can use the built-in function `apply()`. The first argument is the matrix, the second specifies which dimension to retain (1 means row dimension, 2 means column dimension), and the third argument `FUN` specifies the function. Above, we just request the built-in function `mean()`.

Exercise: Recall the claim that the expected value of the sample mean is the true mean. Compute the average of the 4 sample means we have. Is it (sort of) close the true mean? Feel free to change `n` or `B` at any point.




Exercise: Recall the Weak Law of Large Numbers said that, as the sample size gets bigger, the distribution of the sample means gets more concentrated around the true mean. Recall also that the variance of the sample mean is equal to the true data-generating variance divided by the sample size `n`. Explore these probability facts empirically. Don't go crazy -- just pick a few different sample sizes, compute sample means, and explore the variability of the sample means as a function of sample size. Don't stress out about your code, i.e. some repetition is OK. Hint: `var()` and `sd()` will compute variance and standard deviation; `IQR()` and `mad()` will compute other measures of spread, namely, the inter-quartile range and the median absolute deviation.




Here's a table that might be your end result (once you've had a try, you can read the hidden source [here](https://github.com/jennybc/stat540_2014/blob/master/seminars/seminar02_playing-with-probability.rmd))


```
##        sampSize trueSEM   obsSEM sampMeanIQR sampMeanMad
## n10          10 0.31623 0.310619     0.42673     0.31442
## n100        100 0.10000 0.101024     0.13172     0.09806
## n1000      1000 0.03162 0.031598     0.04244     0.03157
## n10000    10000 0.01000 0.009825     0.01414     0.01048
```


Exercise: Repeat the above for a different distribution. Call up the [help on distributions](http://stat.ethz.ch/R-manual/R-patched/library/stats/html/Distributions.html) for an overview of other readily available distributions.

### Compare probabilities and observed relative frequencies

As you've probably already discovered, R's built-in functions can give you probabilities from many distributions. In the normal case, the function `pnorm()` provides the cumulative distribution function (CDF). Go read [the documention on it](http://stat.ethz.ch/R-manual/R-patched/library/stats/html/Normal.html). I'll wait.

Exercise: Generate a reasonably large sample from some normal distribution (it need not be standard normal!). Pick a threshhold. What is the CDF at that threshhold, i.e. what's the true probability of seeing an observation less than or equal to the threshhold? Use your large sample to compute the observed proportion of observations that are less than threshhold. Are the two numbers sort of close? Hint: If `x` is a numeric vector, then `mean(x <= threshhold)` computes the proportion of values less than or equal to `threshhold`.

Exercise: Do the same for a different distribution.

Exercise: Do the same for a variety of sample sizes. Do the two numbers tend to be closer for larger samples?

Exercise: Instead of focusing on values less than the threshhold, focus on values *greater than* the threshhold.

Exercise: Instead of focusing on tail probabilities, focus on the probability of the observed values falling in an interval.

### Disclaimer

Computing the average for each of many samples is an example of what we call *data aggregation*. The built-in `apply` family of functions are one way to approach data aggregation. But I feel that the add-on package `plyr` is substantilly easier to use in the long-run and we will discuss that soon. Similarly, the add-on package `reshape2` is handy for reshaping datasets, which is a closely related and important issue. I deliberately use built-in functions here, to not overwhelm novices, but I will recommend other approaches for general long-term use.

### Explore the distribution of sample means and the CLT

We don't really have the tools yet! Soon we will go over how to create figures, which is really the best way to explore a distribution. Then we can revisit some of the work above and empirically verify, e.g., the Central Limit Theorem. If you already have some graphing prowess, go ahead and try. Here's some code I present without explanation (yet) if you're eager to play around.


```r
library(lattice)

## theoretical vs. empirical distriubution for a single sample
## demo of the Central Limit Theorem
n <- 35
x  <- rnorm(n)

densityplot(~ x)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-141.png) 

```r
densityplot(~x, n = 200, ylim = dnorm(0) * c(-0.1, 1.15),
            panel = function(x, ...) {
              panel.densityplot(x, ...)
              panel.mathdensity(n = 200, col.line = "grey74")
            })
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-142.png) 



Here, let's take a look at Lattice which is convenient for visualizing high-dimensional multivariate data. In this function call for densityplot, try decrease the value of `n=200` (the 2nd argument), what do you observe? The blue curve becomes less smooth, why? This is because `n` defines the number at which the kernel density will be evaluated. `ylim` defines the range of the y axis. The `panel` defines how the two curves overlay: firstly, the empirical distribution, secondly, the theoretical distribution. Now replace `n <- 35` with a larger number, such as 500, what do you observe comparing the theoretical and empirical distribution?




```r
## empirical distribution of sample means for various sample sizes
B <- 1000
n <- round(10^(seq(from = 1, to = 2.5, length = 4)), 0)
names(n) <- paste0("n", n)
getSampleMeans <- function(n, B) colMeans(matrix(rnorm(n * B), nrow = n))
x <- data.frame(sapply(n, getSampleMeans, B))
```


Note that `getSampleMeans` is a user-specified function which takes 2 arguments as input: `n` and `B`. In `sapply(n, getSampleMeans, B)`, the 4 numbers of `n` are iterated one by one, and each number is passed to `getSampleMeans`. When `getSampleMeans` requires more than 1 argument, it will look at the arguments provided at the 3rd position, and so on. Applying `getSampleMeans` on each number returns a vector of size `B`. `sapply` binds the 4 returned vectors together, resulting in a matrix with `B` rows and 4 columns.





```r
## using the "extended formula interface" in lattice
jFormula <- as.formula(paste("~", paste(names(n), sep = "", collapse = " + ")))
## building the formula programmatically is slicker than the alternative, which
## is hard wiring to "~ n10 + n32 + n100 + n316", which is not a crime
densityplot(jFormula, x, xlab = "sample means",
            auto.key = list(x = 0.9, y = 0.9, corner = c(1, 1),
                            reverse.rows = TRUE))
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-161.png) 

```r

## keeping the data "tidy", i.e. tall and skinny, for a happier life
xTallSkinny <- stack(x)
names(xTallSkinny) <- c("x","n")
xTallSkinny$n <- factor(xTallSkinny$n, levels = colnames(x))
densityplot(~ x, xTallSkinny, xlab = "sample means", groups = n,
            auto.key = list(x = 0.9, y = 0.9, corner = c(1, 1),
                            reverse.rows = TRUE))
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-162.png) 

Note that the density plots generated by these 2 alternative methods are actually the same.

