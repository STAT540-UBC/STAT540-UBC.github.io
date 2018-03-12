# Introduction to dplyr

Introduction to data aggregation - `dplyr` version
======================================================================

This content is borrowed from the following two tutorials of Dr.Jenny Bryan's STAT545 course:

  * [Introduction to dplyr](https://github.com/STAT545-UBC/STAT545-UBC.github.io/blob/master/block009_dplyr-intro.md)
  * [dplyr functions for a single dataset](https://github.com/STAT545-UBC/STAT545-UBC.github.io/blob/master/block010_dplyr-end-single-table.md)





### Intro

`dplyr` is a new package for data manipulation. It is built to be fast, highly expressive, and open-minded about how your data is stored. It is developed by Hadley Wickham and Romain Francois.

`dplyr`'s roots are in an earlier, still-very-useful package called [`plyr`](http://plyr.had.co.nz), which implements the "split-apply-combine" strategy for data analysis. Where `plyr` covers a diverse set of inputs and outputs (e.g., arrays, data.frames, lists), `dplyr` has a laser-like focus on data.frames and related structures.

Have no idea what I'm talking about? Not sure if you care? If you use these base R functions: `subset()`, `apply()`, `[sl]apply()`, `tapply()`, `aggregate()`, `split()`, `do.call()`, then you should keep reading.

#### Load `dplyr`


```r
## install if you do not already have

## from CRAN:
## install.packages('dplyr')

## from GitHub using devtools (which you also might need to install!):
## devtools::install_github("hadley/lazyeval")
## devtools::install_github("hadley/dplyr")
suppressPackageStartupMessages(library(dplyr))
```

#### Load the Gapminder data

An excerpt of the Gapminder data which we work with alot.


```r
gd_url <- "http://tiny.cc/gapminder"
gdf <- read.delim(file = gd_url)
str(gdf)
```

```
## 'data.frame':	1704 obs. of  6 variables:
##  $ country  : Factor w/ 142 levels "Afghanistan",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ year     : int  1952 1957 1962 1967 1972 1977 1982 1987 1992 1997 ...
##  $ pop      : num  8425333 9240934 10267083 11537966 13079460 ...
##  $ continent: Factor w/ 5 levels "Africa","Americas",..: 3 3 3 3 3 3 3 3 3 3 ...
##  $ lifeExp  : num  28.8 30.3 32 34 36.1 ...
##  $ gdpPercap: num  779 821 853 836 740 ...
```

```r
head(gdf)
```

```
##       country year      pop continent lifeExp gdpPercap
## 1 Afghanistan 1952  8425333      Asia  28.801  779.4453
## 2 Afghanistan 1957  9240934      Asia  30.332  820.8530
## 3 Afghanistan 1962 10267083      Asia  31.997  853.1007
## 4 Afghanistan 1967 11537966      Asia  34.020  836.1971
## 5 Afghanistan 1972 13079460      Asia  36.088  739.9811
## 6 Afghanistan 1977 14880372      Asia  38.438  786.1134
```

### Meet `tbl_df`, an upgrade to `data.frame`


```r
gtbl <- tbl_df(gdf)
gtbl
```

```
## Source: local data frame [1,704 x 6]
## 
##        country year      pop continent lifeExp gdpPercap
## 1  Afghanistan 1952  8425333      Asia  28.801  779.4453
## 2  Afghanistan 1957  9240934      Asia  30.332  820.8530
## 3  Afghanistan 1962 10267083      Asia  31.997  853.1007
## 4  Afghanistan 1967 11537966      Asia  34.020  836.1971
## 5  Afghanistan 1972 13079460      Asia  36.088  739.9811
## 6  Afghanistan 1977 14880372      Asia  38.438  786.1134
## 7  Afghanistan 1982 12881816      Asia  39.854  978.0114
## 8  Afghanistan 1987 13867957      Asia  40.822  852.3959
## 9  Afghanistan 1992 16317921      Asia  41.674  649.3414
## 10 Afghanistan 1997 22227415      Asia  41.763  635.3414
## ..         ...  ...      ...       ...     ...       ...
```

```r
glimpse(gtbl)
```

```
## Observations: 1704
## Variables:
## $ country   (fctr) Afghanistan, Afghanistan, Afghanistan, Afghanistan,...
## $ year      (int) 1952, 1957, 1962, 1967, 1972, 1977, 1982, 1987, 1992...
## $ pop       (dbl) 8425333, 9240934, 10267083, 11537966, 13079460, 1488...
## $ continent (fctr) Asia, Asia, Asia, Asia, Asia, Asia, Asia, Asia, Asi...
## $ lifeExp   (dbl) 28.801, 30.332, 31.997, 34.020, 36.088, 38.438, 39.8...
## $ gdpPercap (dbl) 779.4453, 820.8530, 853.1007, 836.1971, 739.9811, 78...
```

A `tbl_df` is basically an improved data.frame, for which `dplyr` provides nice methods for high-level inspection. Specifically, these methods do something sensible for datasets with many observations and/or variables. You do __NOT__ need to turn your data.frames into `tbl_df`s to use `plyr`. I do so here for demonstration purposes only.

### Think before you create excerpts of your data ...

If you feel the urge to store a little snippet of your data:


```r
(snippet <- subset(gdf, country == "Canada"))
```

```
##     country year      pop continent lifeExp gdpPercap
## 241  Canada 1952 14785584  Americas  68.750  11367.16
## 242  Canada 1957 17010154  Americas  69.960  12489.95
## 243  Canada 1962 18985849  Americas  71.300  13462.49
## 244  Canada 1967 20819767  Americas  72.130  16076.59
## 245  Canada 1972 22284500  Americas  72.880  18970.57
## 246  Canada 1977 23796400  Americas  74.210  22090.88
## 247  Canada 1982 25201900  Americas  75.760  22898.79
## 248  Canada 1987 26549700  Americas  76.860  26626.52
## 249  Canada 1992 28523502  Americas  77.950  26342.88
## 250  Canada 1997 30305843  Americas  78.610  28954.93
## 251  Canada 2002 31902268  Americas  79.770  33328.97
## 252  Canada 2007 33390141  Americas  80.653  36319.24
```

Stop and ask yourself ...

> Do I want to create mini datasets for each level of some factor (or unique combination of several factors) ... in order to compute or graph something?  

If YES, __use proper data aggregation techniques__ or facetting in `ggplot2` plots or conditioning in `lattice` -- __don’t subset the data__. Or, more realistic, only subset the data as a temporary measure while you develop your elegant code for computing on or visualizing these data subsets.

If NO, then maybe you really do need to store a copy of a subset of the data. But seriously consider whether you can achieve your goals by simply using the `subset =` argument of, e.g., the `lm()` function, to limit computation to your excerpt of choice. Lots of functions offer a `subset =` argument!

Copies and excerpts of your data clutter your workspace, invite mistakes, and sow general confusion. Avoid whenever possible.

Reality can also lie somewhere in between. You will find the workflows presented below can help you accomplish your goals with minimal creation of temporary, intermediate objects.

### Use `filter()` to subset data row-wise.

`filter()` takes logical expressions and returns the rows for which all are `TRUE`.


```r
filter(gtbl, lifeExp < 29)
```

```
## Source: local data frame [2 x 6]
## 
##       country year     pop continent lifeExp gdpPercap
## 1 Afghanistan 1952 8425333      Asia  28.801  779.4453
## 2      Rwanda 1992 7290203    Africa  23.599  737.0686
```

```r
filter(gtbl, country == "Rwanda")
```

```
## Source: local data frame [12 x 6]
## 
##    country year     pop continent lifeExp gdpPercap
## 1   Rwanda 1952 2534927    Africa  40.000  493.3239
## 2   Rwanda 1957 2822082    Africa  41.500  540.2894
## 3   Rwanda 1962 3051242    Africa  43.000  597.4731
## 4   Rwanda 1967 3451079    Africa  44.100  510.9637
## 5   Rwanda 1972 3992121    Africa  44.600  590.5807
## 6   Rwanda 1977 4657072    Africa  45.000  670.0806
## 7   Rwanda 1982 5507565    Africa  46.218  881.5706
## 8   Rwanda 1987 6349365    Africa  44.020  847.9912
## 9   Rwanda 1992 7290203    Africa  23.599  737.0686
## 10  Rwanda 1997 7212583    Africa  36.087  589.9445
## 11  Rwanda 2002 7852401    Africa  43.413  785.6538
## 12  Rwanda 2007 8860588    Africa  46.242  863.0885
```

```r
filter(gtbl, country %in% c("Rwanda", "Afghanistan"))
```

```
## Source: local data frame [24 x 6]
## 
##        country year      pop continent lifeExp gdpPercap
## 1  Afghanistan 1952  8425333      Asia  28.801  779.4453
## 2  Afghanistan 1957  9240934      Asia  30.332  820.8530
## 3  Afghanistan 1962 10267083      Asia  31.997  853.1007
## 4  Afghanistan 1967 11537966      Asia  34.020  836.1971
## 5  Afghanistan 1972 13079460      Asia  36.088  739.9811
## 6  Afghanistan 1977 14880372      Asia  38.438  786.1134
## 7  Afghanistan 1982 12881816      Asia  39.854  978.0114
## 8  Afghanistan 1987 13867957      Asia  40.822  852.3959
## 9  Afghanistan 1992 16317921      Asia  41.674  649.3414
## 10 Afghanistan 1997 22227415      Asia  41.763  635.3414
## ..         ...  ...      ...       ...     ...       ...
```

Compare with some base R code to accomplish the same things

```r
gdf[gdf$lifeExp < 29, ] ## repeat `gdf`, [i, j] indexing is distracting
subset(gdf, country == "Rwanda") ## almost same as filter ... but wait ...
```

### Meet the new pipe operator

Before we go any further, we should exploit the new pipe operator that `dplyr` imports from the [`magrittr`](https://github.com/smbache/magrittr) package. This is going to change your data analytical life. You no longer need to enact multi-operation commands by nesting them inside each other, like so many [Russian nesting dolls](http://blogue.us/wp-content/uploads/2009/07/Unknown-21.jpeg). This new syntax leads to code that is much easier to write and to read.

Here's what it looks like: `%>%`. The RStudio keyboard shortcut: Ctrl + Shift + M (Windows), Cmd + Shift + M (Mac), according to [this tweet](https://twitter.com/rstudiotips/status/514094879316525058).

<!-- alt-shift-. works for me but I'm not running bleeding edge RStudio -->

Let's demo then I'll explain:


```r
gdf %>% head
```

```
##       country year      pop continent lifeExp gdpPercap
## 1 Afghanistan 1952  8425333      Asia  28.801  779.4453
## 2 Afghanistan 1957  9240934      Asia  30.332  820.8530
## 3 Afghanistan 1962 10267083      Asia  31.997  853.1007
## 4 Afghanistan 1967 11537966      Asia  34.020  836.1971
## 5 Afghanistan 1972 13079460      Asia  36.088  739.9811
## 6 Afghanistan 1977 14880372      Asia  38.438  786.1134
```

This is equivalent to `head(gdf)`. This pipe operator takes the thing on the left-hand-side and __pipes__ it into the function call on the right-hand-side -- literally, drops it in as the first argument.

Never fear, you can still specify other arguments to this function! To see the first 3 rows of Gapminder, we could say `head(gdf, 3)` or this:


```r
gdf %>% head(3)
```

```
##       country year      pop continent lifeExp gdpPercap
## 1 Afghanistan 1952  8425333      Asia  28.801  779.4453
## 2 Afghanistan 1957  9240934      Asia  30.332  820.8530
## 3 Afghanistan 1962 10267083      Asia  31.997  853.1007
```

I've advised you to think "gets" whenever you see the assignment operator, `<-`. Similary, you should think "then" whenever you see the pipe operator, `%>%`.

You are probably not impressed yet, but the magic will soon happen.

### Use `select()` to subset the data on variables or columns.

Back to `dplyr` ...

Use `select()` to subset the data on variables or columns. Here's a conventional call:


```r
select(gtbl, year, lifeExp) ## tbl_df prevents TMI from printing
```

```
## Source: local data frame [1,704 x 2]
## 
##    year lifeExp
## 1  1952  28.801
## 2  1957  30.332
## 3  1962  31.997
## 4  1967  34.020
## 5  1972  36.088
## 6  1977  38.438
## 7  1982  39.854
## 8  1987  40.822
## 9  1992  41.674
## 10 1997  41.763
## ..  ...     ...
```

And here's similar operation, but written with the pipe operator and piped through `head`:

```r
gtbl %>%
  select(year, lifeExp) %>%
  head(4)
```

```
## Source: local data frame [4 x 2]
## 
##   year lifeExp
## 1 1952  28.801
## 2 1957  30.332
## 3 1962  31.997
## 4 1967  34.020
```

Think: "Take `gtbl`, then select the variables year and lifeExp, then show the first 4 rows."

### Revel in the convenience

Here's the data for Cambodia, but only certain variables:


```r
gtbl %>%
  filter(country == "Cambodia") %>%
  select(year, lifeExp)
```

```
## Source: local data frame [12 x 2]
## 
##    year lifeExp
## 1  1952  39.417
## 2  1957  41.366
## 3  1962  43.415
## 4  1967  45.415
## 5  1972  40.317
## 6  1977  31.220
## 7  1982  50.957
## 8  1987  53.914
## 9  1992  55.803
## 10 1997  56.534
## 11 2002  56.752
## 12 2007  59.723
```

and what a typical base R call would look like:


```r
gdf[gdf$country == "Cambodia", c("year", "lifeExp")]
```

```
##     year lifeExp
## 217 1952  39.417
## 218 1957  41.366
## 219 1962  43.415
## 220 1967  45.415
## 221 1972  40.317
## 222 1977  31.220
## 223 1982  50.957
## 224 1987  53.914
## 225 1992  55.803
## 226 1997  56.534
## 227 2002  56.752
## 228 2007  59.723
```

or, possibly?, a nicer look using base R's `subset()` function:


```r
subset(gdf, country == "Cambodia", select = c(year, lifeExp))
```

```
##     year lifeExp
## 217 1952  39.417
## 218 1957  41.366
## 219 1962  43.415
## 220 1967  45.415
## 221 1972  40.317
## 222 1977  31.220
## 223 1982  50.957
## 224 1987  53.914
## 225 1992  55.803
## 226 1997  56.534
## 227 2002  56.752
## 228 2007  59.723
```

### Pause to reflect

We've barely scratched the surface of `dplyr` but I want to point out key principles you may start to appreciate. If you're new to R or "programing with data", feel free skip this section and [move on](block010_dplyr-end-single-table.html).

`dplyr`'s verbs, such as `filter()` and `select()`, are what's called [pure functions](http://en.wikipedia.org/wiki/Pure_function). To quote from Wickham's [Advanced R Programming book](http://adv-r.had.co.nz/Functions.html):

> The functions that are the easiest to understand and reason about are pure functions: functions that always map the same input to the same output and have no other impact on the workspace. In other words, pure functions have no side effects: they don’t affect the state of the world in any way apart from the value they return.

In fact, these verbs are a special case of pure functions: they take the same flavor of object as input and output. Namely, a data.frame or one of the other data receptacles `dplyr` supports. And finally, the data is __always__ the very first argument of the verb functions.

This set of deliberate design choices, together with the new pipe operator, produces a highly effective, low friction [domain-specific language](http://adv-r.had.co.nz/dsl.html) for data analysis.

Go to the next block, [`dplyr` functions for a single dataset](block010_dplyr-end-single-table.html), for more `dplyr`!



### Where were we?

So far, we used two very important verbs and an operator:

  * `filter()` for subsetting data row-wise
  * `select()` for subsetting data variable- or column-wise
  * the pipe operator `%>%`, which feeds the LHS as the first argument to the expression on the RHS
  
Here we explore other `dplyr` functions, especially more verbs, for working with a single dataset.

#### Load `dplyr` and the Gapminder data

We use an excerpt of the Gapminder data and store it as a `tbl_df` object, basically an enhanced data.frame. I'll use the pipe operator even here, to demonstrate its utility outside of `dplyr`.


```r
suppressPackageStartupMessages(library(dplyr))
gd_url <- "http://tiny.cc/gapminder"
gtbl <- gd_url %>% read.delim %>% tbl_df
gtbl %>% glimpse
```

```
## Observations: 1704
## Variables:
## $ country   (fctr) Afghanistan, Afghanistan, Afghanistan, Afghanistan,...
## $ year      (int) 1952, 1957, 1962, 1967, 1972, 1977, 1982, 1987, 1992...
## $ pop       (dbl) 8425333, 9240934, 10267083, 11537966, 13079460, 1488...
## $ continent (fctr) Asia, Asia, Asia, Asia, Asia, Asia, Asia, Asia, Asi...
## $ lifeExp   (dbl) 28.801, 30.332, 31.997, 34.020, 36.088, 38.438, 39.8...
## $ gdpPercap (dbl) 779.4453, 820.8530, 853.1007, 836.1971, 739.9811, 78...
```

### Use `mutate()` to add new variables

Imagine we wanted to recover each country's GDP. After all, the Gapminder data has a variable for population and GDP per capita. Let's multiply them together.


```r
gtbl <- gtbl %>%
  mutate(gdp = pop * gdpPercap)
gtbl %>% glimpse
```

```
## Observations: 1704
## Variables:
## $ country   (fctr) Afghanistan, Afghanistan, Afghanistan, Afghanistan,...
## $ year      (int) 1952, 1957, 1962, 1967, 1972, 1977, 1982, 1987, 1992...
## $ pop       (dbl) 8425333, 9240934, 10267083, 11537966, 13079460, 1488...
## $ continent (fctr) Asia, Asia, Asia, Asia, Asia, Asia, Asia, Asia, Asi...
## $ lifeExp   (dbl) 28.801, 30.332, 31.997, 34.020, 36.088, 38.438, 39.8...
## $ gdpPercap (dbl) 779.4453, 820.8530, 853.1007, 836.1971, 739.9811, 78...
## $ gdp       (dbl) 6567086330, 7585448670, 8758855797, 9648014150, 9678...
```

Hmmmm ... those GDP numbers are almost uselessly large and abstract. Consider the [advice of Randall Munroe of xkcd](http://fivethirtyeight.com/datalab/xkcd-randall-munroe-qanda-what-if/): "One thing that bothers me is large numbers presented without context... 'If I added a zero to this number, would the sentence containing it mean something different to me?' If the answer is 'no,' maybe the number has no business being in the sentence in the first place." Maybe it would be more meaningful to consumers of my tables and figures if I reported GDP per capita, *relative to some benchmark country*. Since Canada is my adopted home, I'll go with that.


```r
just_canada <- gtbl %>% filter(country == "Canada")
gtbl <- gtbl %>%
  mutate(canada = just_canada$gdpPercap[match(year, just_canada$year)],
         gdpPercapRel = gdpPercap / canada)
gtbl %>%
  select(country, year, gdpPercap, canada, gdpPercapRel)
```

```
## Source: local data frame [1,704 x 5]
## 
##        country year gdpPercap   canada gdpPercapRel
## 1  Afghanistan 1952  779.4453 11367.16   0.06856992
## 2  Afghanistan 1957  820.8530 12489.95   0.06572108
## 3  Afghanistan 1962  853.1007 13462.49   0.06336874
## 4  Afghanistan 1967  836.1971 16076.59   0.05201335
## 5  Afghanistan 1972  739.9811 18970.57   0.03900679
## 6  Afghanistan 1977  786.1134 22090.88   0.03558542
## 7  Afghanistan 1982  978.0114 22898.79   0.04271018
## 8  Afghanistan 1987  852.3959 26626.52   0.03201305
## 9  Afghanistan 1992  649.3414 26342.88   0.02464959
## 10 Afghanistan 1997  635.3414 28954.93   0.02194243
## ..         ...  ...       ...      ...          ...
```

```r
gtbl %>%
  select(gdpPercapRel) %>%
  summary
```

```
##   gdpPercapRel     
##  Min.   :0.007236  
##  1st Qu.:0.061648  
##  Median :0.171521  
##  Mean   :0.326659  
##  3rd Qu.:0.446564  
##  Max.   :9.534690
```

Note that, `mutate()` builds new variables sequentially so you can reference earlier ones (like `canada`) when defining later ones (like `gdpPercapRel`). (I got a little off topic here using `match()` to do table look up, but [you can figure that out](http://www.rdocumentation.org/packages/base/functions/match).)

The relative GDP per capita numbers are, in general, well below 1. We see that most of the countries covered by this dataset have substantially lower GDP per capita, relative to Canada, across the entire time period.

### Use `arrange()` to row-order data in a principled way

Imagine you wanted this data ordered by year then country, as opposed to by country then year.


```r
gtbl %>%
  arrange(year, country)
```

```
## Source: local data frame [1,704 x 9]
## 
##        country year      pop continent lifeExp  gdpPercap          gdp
## 1  Afghanistan 1952  8425333      Asia  28.801   779.4453   6567086330
## 2      Albania 1952  1282697    Europe  55.230  1601.0561   2053669902
## 3      Algeria 1952  9279525    Africa  43.077  2449.0082  22725632678
## 4       Angola 1952  4232095    Africa  30.015  3520.6103  14899557133
## 5    Argentina 1952 17876956  Americas  62.485  5911.3151 105676319105
## 6    Australia 1952  8691212   Oceania  69.120 10039.5956  87256254102
## 7      Austria 1952  6927772    Europe  66.800  6137.0765  42516266683
## 8      Bahrain 1952   120447      Asia  50.939  9867.0848   1188460759
## 9   Bangladesh 1952 46886859      Asia  37.484   684.2442  32082059995
## 10     Belgium 1952  8730405    Europe  68.000  8343.1051  72838686716
## ..         ...  ...      ...       ...     ...        ...          ...
## Variables not shown: canada (dbl), gdpPercapRel (dbl)
```

Or maybe you want just the data from 2007, sorted on life expectancy?


```r
gtbl %>%
  filter(year == 2007) %>%
  arrange(lifeExp)
```

```
## Source: local data frame [142 x 9]
## 
##                     country year      pop continent lifeExp gdpPercap
## 1                 Swaziland 2007  1133066    Africa  39.613 4513.4806
## 2                Mozambique 2007 19951656    Africa  42.082  823.6856
## 3                    Zambia 2007 11746035    Africa  42.384 1271.2116
## 4              Sierra Leone 2007  6144562    Africa  42.568  862.5408
## 5                   Lesotho 2007  2012649    Africa  42.592 1569.3314
## 6                    Angola 2007 12420476    Africa  42.731 4797.2313
## 7                  Zimbabwe 2007 12311143    Africa  43.487  469.7093
## 8               Afghanistan 2007 31889923      Asia  43.828  974.5803
## 9  Central African Republic 2007  4369038    Africa  44.741  706.0165
## 10                  Liberia 2007  3193942    Africa  45.678  414.5073
## ..                      ...  ...      ...       ...     ...       ...
## Variables not shown: gdp (dbl), canada (dbl), gdpPercapRel (dbl)
```

Oh, you'd like to sort on life expectancy in __desc__ending order? Then use `desc()`.


```r
gtbl %>%
  filter(year == 2007) %>%
  arrange(desc(lifeExp))
```

```
## Source: local data frame [142 x 9]
## 
##             country year       pop continent lifeExp gdpPercap
## 1             Japan 2007 127467972      Asia  82.603  31656.07
## 2  Hong Kong, China 2007   6980412      Asia  82.208  39724.98
## 3           Iceland 2007    301931    Europe  81.757  36180.79
## 4       Switzerland 2007   7554661    Europe  81.701  37506.42
## 5         Australia 2007  20434176   Oceania  81.235  34435.37
## 6             Spain 2007  40448191    Europe  80.941  28821.06
## 7            Sweden 2007   9031088    Europe  80.884  33859.75
## 8            Israel 2007   6426679      Asia  80.745  25523.28
## 9            France 2007  61083916    Europe  80.657  30470.02
## 10           Canada 2007  33390141  Americas  80.653  36319.24
## ..              ...  ...       ...       ...     ...       ...
## Variables not shown: gdp (dbl), canada (dbl), gdpPercapRel (dbl)
```

I advise that your analyses NEVER rely on rows or variables being in a specific order. But it's still true that human beings write the code and the interactive development process can be much nicer if you reorder the rows of your data as you go along. Also, once you are preparing tables for human eyeballs, it is imperative that you step up and take control of row order.

### Use `rename()` to rename variables

*NOTE: I am using the development version of `dplyr` which will soon become the official release 0.3. If `rename()` does not work for you, try `rename_vars()`, which is what this function is called in version 0.2 on CRAN. You could also use `plyr::rename()`, but then you have to be careful to always load `plyr` before `dplyr`.*

I am in the awkward life stage of switching from [`camelCase`](http://en.wikipedia.org/wiki/CamelCase) to [`snake_case`](http://en.wikipedia.org/wiki/Snake_case), so I am vexed by the variable names I chose when I cleaned this data years ago. Let's rename some variables!


```r
gtbl %>%
  rename(life_exp = lifeExp, gdp_percap = gdpPercap,
         gdp_percap_rel = gdpPercapRel)
```

```
## Source: local data frame [1,704 x 9]
## 
##        country year      pop continent life_exp gdp_percap         gdp
## 1  Afghanistan 1952  8425333      Asia   28.801   779.4453  6567086330
## 2  Afghanistan 1957  9240934      Asia   30.332   820.8530  7585448670
## 3  Afghanistan 1962 10267083      Asia   31.997   853.1007  8758855797
## 4  Afghanistan 1967 11537966      Asia   34.020   836.1971  9648014150
## 5  Afghanistan 1972 13079460      Asia   36.088   739.9811  9678553274
## 6  Afghanistan 1977 14880372      Asia   38.438   786.1134 11697659231
## 7  Afghanistan 1982 12881816      Asia   39.854   978.0114 12598563401
## 8  Afghanistan 1987 13867957      Asia   40.822   852.3959 11820990309
## 9  Afghanistan 1992 16317921      Asia   41.674   649.3414 10595901589
## 10 Afghanistan 1997 22227415      Asia   41.763   635.3414 14121995875
## ..         ...  ...      ...       ...      ...        ...         ...
## Variables not shown: canada (dbl), gdp_percap_rel (dbl)
```

I did NOT assign the post-rename object back to `gtbl` because that would make the chunks in this tutorial harder to copy/paste and run out of order. In real life, I would probably assign this back to `gtbl`, in a data preparation script, and proceed with the new variable names.

### `group_by()` is a mighty weapon

I have found friends and family love to ask seemingly innocuous questions like, "which country experienced the sharpest 5-year drop in life expectancy?". In fact, that is a totally natural question to ask. But if you are using a language that doesn't know about data, it's an incredibly annoying question to answer.

`dplyr` offers powerful tools to solve this class of problem.

  * `group_by()` adds extra structure to your dataset -- grouping information -- which lays the groundwork for computations within the groups.
  * `summarize()` takes a dataset with $n$ observations, computes requested summaries, and returns a dataset with 1 observation.
  * window functions take a dataset with $n$ observations and return a dataset with $n$ observations.
  
Combined with the verbs you already know, these new tools allow you to solve an extremely diverse set of problems with relative ease.

#### Counting things up

Let's start with simple counting.  How many observations do we have per continent?


```r
gtbl %>%
  group_by(continent) %>%
  summarize(n_obs = n())
```

```
## Source: local data frame [5 x 2]
## 
##   continent n_obs
## 1    Africa   624
## 2  Americas   300
## 3      Asia   396
## 4    Europe   360
## 5   Oceania    24
```

The `tally()` function is a convenience function for this sort of thing.


```r
gtbl %>%
  group_by(continent) %>%
  tally
```

```
## Source: local data frame [5 x 2]
## 
##   continent   n
## 1    Africa 624
## 2  Americas 300
## 3      Asia 396
## 4    Europe 360
## 5   Oceania  24
```

What if we wanted to add the number of unique countries for each continent?


```r
gtbl %>%
  group_by(continent) %>%
  summarize(n_obs = n(), n_countries = n_distinct(country))
```

```
## Source: local data frame [5 x 3]
## 
##   continent n_obs n_countries
## 1    Africa   624          52
## 2  Americas   300          25
## 3      Asia   396          33
## 4    Europe   360          30
## 5   Oceania    24           2
```

#### General summarization

The functions you'll apply within `summarize()` include classical statistical summaries, like `mean()`, `median()`, `sd()`, and `IQR`. Remember they are functions that take $n$ inputs and distill them down into 1 output.

Although this may be statistically ill-advised, let's compute the average life expectancy by continent.


```r
gtbl %>%
  group_by(continent) %>%
  summarize(avg_lifeExp = mean(lifeExp))
```

```
## Source: local data frame [5 x 2]
## 
##   continent avg_lifeExp
## 1    Africa    48.86533
## 2  Americas    64.65874
## 3      Asia    60.06490
## 4    Europe    71.90369
## 5   Oceania    74.32621
```

`summarize_each()` applies the same summary function(s) to multiple variables. Let's compute average and median life expectancy and GDP per capita by continent by year ... but only for 1952 and 2007.

*NOTE: you won't have `summarize_each()` if you're using `dplyr` version 0.2. Just wait for it.*


```r
gtbl %>%
  filter(year %in% c(1952, 2007)) %>%
  group_by(continent, year) %>%
  summarise_each(funs(mean, median), lifeExp, gdpPercap)
```

```
## Source: local data frame [10 x 6]
## Groups: continent
## 
##    continent year lifeExp_mean gdpPercap_mean lifeExp_median
## 1     Africa 1952     39.13550       1252.572        38.8330
## 2     Africa 2007     54.80604       3089.033        52.9265
## 3   Americas 1952     53.27984       4079.063        54.7450
## 4   Americas 2007     73.60812      11003.032        72.8990
## 5       Asia 1952     46.31439       5195.484        44.8690
## 6       Asia 2007     70.72848      12473.027        72.3960
## 7     Europe 1952     64.40850       5661.057        65.9000
## 8     Europe 2007     77.64860      25054.482        78.6085
## 9    Oceania 1952     69.25500      10298.086        69.2550
## 10   Oceania 2007     80.71950      29810.188        80.7195
## Variables not shown: gdpPercap_median (dbl)
```

Let's focus just on Asia. What are the minimum and maximum life expectancies seen by year?

```r
gtbl %>%
  filter(continent == "Asia") %>%
  group_by(year) %>%
  summarize(min_lifeExp = min(lifeExp), max_lifeExp = max(lifeExp))
```

```
## Source: local data frame [12 x 3]
## 
##    year min_lifeExp max_lifeExp
## 1  1952      28.801      65.390
## 2  1957      30.332      67.840
## 3  1962      31.997      69.390
## 4  1967      34.020      71.430
## 5  1972      36.088      73.420
## 6  1977      31.220      75.380
## 7  1982      39.854      77.110
## 8  1987      40.822      78.670
## 9  1992      41.674      79.360
## 10 1997      41.763      80.690
## 11 2002      42.129      82.000
## 12 2007      43.828      82.603
```

Of course it would be much more interesting to see *which* country contributed these extreme observations. Is the minimum (maximum) always coming from the same country? That's where window functions come in.

#### Window functions

Recall that window functions take $n$ inputs and give back $n$ outputs. Here we use window functions based on ranks and offsets.

Let's revisit the worst and best life expectancies in Asia over time, but retaining info about *which* country contributes these extreme values.


```r
gtbl %>%
  filter(continent == "Asia") %>%
  select(year, country, lifeExp) %>%
  arrange(year) %>%
  group_by(year) %>%
  filter(min_rank(desc(lifeExp)) < 2 | min_rank(lifeExp) < 2)
```

```
## Source: local data frame [24 x 3]
## Groups: year
## 
##    year     country lifeExp
## 1  1952 Afghanistan  28.801
## 2  1952      Israel  65.390
## 3  1957 Afghanistan  30.332
## 4  1957      Israel  67.840
## 5  1962 Afghanistan  31.997
## 6  1962      Israel  69.390
## 7  1967 Afghanistan  34.020
## 8  1967       Japan  71.430
## 9  1972 Afghanistan  36.088
## 10 1972       Japan  73.420
## ..  ...         ...     ...
```

We see that (min = Agfhanistan, max = Japan) is the most frequent result, but Cambodia and Israel pop up at least once each as the min or max, respectively. That table should make you impatient for our upcoming work on tidying and reshaping data! Wouldn't it be nice to have one row per year?

How did that actually work? First, I store and view the result including everything but the last `filter()` statement. All of these operations are familiar.


```r
asia <- gtbl %>%
  filter(continent == "Asia") %>%
  select(year, country, lifeExp) %>%
  arrange(year) %>%
  group_by(year)
asia
```

```
## Source: local data frame [396 x 3]
## Groups: year
## 
##    year          country lifeExp
## 1  1952      Afghanistan  28.801
## 2  1952          Bahrain  50.939
## 3  1952       Bangladesh  37.484
## 4  1952         Cambodia  39.417
## 5  1952            China  44.000
## 6  1952 Hong Kong, China  60.960
## 7  1952            India  37.373
## 8  1952        Indonesia  37.468
## 9  1952             Iran  44.869
## 10 1952             Iraq  45.320
## ..  ...              ...     ...
```

Now we apply a window function -- `min_rank()`. Since `asia` is grouped by year, `min_rank()` operates within mini-datasets, each for a specific year. Applied to the variable `lifeExp`, `min_rank()` returns the rank of each country's observed life expectancy. FYI, the `min` part just specifies how ties are broken. Here is an explicit peek at these within-year life expectancy ranks, in both the (default) ascending and descending order.


```r
asia %>%
  mutate(le_rank = min_rank(lifeExp),
         le_desc_rank = min_rank(desc(lifeExp)))
```

```
## Source: local data frame [396 x 5]
## Groups: year
## 
##    year          country lifeExp le_rank le_desc_rank
## 1  1952      Afghanistan  28.801       1           33
## 2  1952          Bahrain  50.939      25            9
## 3  1952       Bangladesh  37.484       7           27
## 4  1952         Cambodia  39.417       9           25
## 5  1952            China  44.000      16           18
## 6  1952 Hong Kong, China  60.960      31            3
## 7  1952            India  37.373       5           29
## 8  1952        Indonesia  37.468       6           28
## 9  1952             Iran  44.869      17           17
## 10 1952             Iraq  45.320      18           16
## ..  ...              ...     ...     ...          ...
```

You can understand the original `filter()` statement now:


```r
filter(min_rank(desc(lifeExp)) < 2 | min_rank(lifeExp) < 2)
```

These two sets of ranks are formed, within year group, and `filter()` retains rows with rank less than 2, which means ... the row with rank = 1. Since we do for ascending and descending ranks, we get both the min and the max.

If we had wanted just the min OR the max, an alternative approach using `top_n()` would have worked.

```r
gtbl %>%
  filter(continent == "Asia") %>%
  select(year, country, lifeExp) %>%
  arrange(year) %>%
  group_by(year) %>%
  #top_n(1)               ## gets the min
  top_n(1, desc(lifeExp)) ## gets the max
```

```
## Source: local data frame [12 x 3]
## Groups: year
## 
##    year     country lifeExp
## 1  1952 Afghanistan  28.801
## 2  1957 Afghanistan  30.332
## 3  1962 Afghanistan  31.997
## 4  1967 Afghanistan  34.020
## 5  1972 Afghanistan  36.088
## 6  1977    Cambodia  31.220
## 7  1982 Afghanistan  39.854
## 8  1987 Afghanistan  40.822
## 9  1992 Afghanistan  41.674
## 10 1997 Afghanistan  41.763
## 11 2002 Afghanistan  42.129
## 12 2007 Afghanistan  43.828
```

#### Grand Finale

So let's answer that "simple" question: which country experienced the sharpest 5-year drop in life expectancy? Recall that this excerpt of the Gapminder data only has data every five years, e.g. for 1952, 1957, etc. So this really means looking at life expectancy changes between adjacent timepoints.

At this point, that's just too easy, so let's do it by continent while we're at it.


```r
gtbl %>%
  group_by(continent, country) %>%
  select(country, year, continent, lifeExp) %>%
  mutate(le_delta = lifeExp - lag(lifeExp)) %>%
  summarize(worst_le_delta = min(le_delta, na.rm = TRUE)) %>%
  filter(min_rank(worst_le_delta) < 2) %>%
  arrange(worst_le_delta)
```

```
## Source: local data frame [5 x 3]
## Groups: continent
## 
##   continent     country worst_le_delta
## 1    Africa      Rwanda        -20.421
## 2  Americas El Salvador         -1.511
## 3      Asia    Cambodia         -9.097
## 4    Europe  Montenegro         -1.464
## 5   Oceania   Australia          0.170
```

Ponder that for a while. The subject matter and the code. Mostly you're seeing what genocide looks like in dry statistics on average life expectancy.

Break the code into pieces, starting at the top, and inspect the intermediate results. That's certainly how I was able to *write* such a thing. These commands do not [leap fully formed out of anyone's forehead](http://tinyurl.com/athenaforehead) -- they are built up gradually, with lots of errors and refinements along the way. I'm not even sure it's a great idea to do so much manipulation in one fell swoop. Is the statement above really hard for you to read? If yes, then by all means break it into pieces and make some intermediate objects. Your code should be easy to write and read when you're done.

In later tutorials, we'll explore more of `dplyr`, such as operations based on two datasets.


### Take-home Exercise

Can you modify the code above to answer this questions: which country gains the most growth in GDP in a 5-year interval? Between which 2 years?




### Resources

`dplyr` official stuff

  * package home [on CRAN](http://cran.r-project.org/web/packages/dplyr/index.html)
    - note there are several vignettes, with the [introduction](http://cran.r-project.org/web/packages/dplyr/vignettes/introduction.html) being the most relevant right now
    - the [one on window functions](http://cran.rstudio.com/web/packages/dplyr/vignettes/window-functions.html) will also be interesting to you now
  * development home [on GitHub](https://github.com/hadley/dplyr)
  * [tutorial HW delivered](https://www.dropbox.com/sh/i8qnluwmuieicxc/AAAgt9tIKoIm7WZKIyK25lh6a) (note this links to a DropBox folder) at useR! 2014 conference

Blog post [Hands-on dplyr tutorial for faster data manipulation in R](http://www.dataschool.io/dplyr-tutorial-for-faster-data-manipulation-in-r/) by Data School, that includes a link to an R Markdown document and links to videos

[Cheatsheet](bit001_dplyr-cheatsheet.html) I made for `dplyr` join functions (not relevant yet but soon)

### Resources on data aggregation via 'plyr'

`plyr` is the older iteration of `dplyr`. There is a nice tutorial from STAT545 for [data aggregation using `plyr`](http://stat545-ubc.github.io/block013_plyr-ddply.html) and some helpful lecture slides [Intro to data aggregation](http://www.slideshare.net/jenniferbryan5811/cm009-data-aggregation).









