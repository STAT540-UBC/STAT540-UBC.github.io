STAT 540 - Assignment 1 - Quality Control
================

The dataset used for this assignment has been published by Scheffer et al. in 2015. The raw RNA-Seq reads have been submitted to GEO under the series ID [GSE60019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60019). Be sure to review the [paper](https://www.ncbi.nlm.nih.gov/pubmed/25904789) to gain some familiarity with the study before you start.

The transcriptomic data and the samples metadata can be downloaded here: [Samples Metadata](data/gse60019_expression_matrix.RDS), and [Expression Matrix](data/gse60019_experiment_design.RDS).

The raw reads have been mapped and processed into gene expression values using an established RNA-Seq pipeline at the [Pavlidis Lab](http://pavlab.msl.ubc.ca/) during data curation in [Gemma](https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=9818). Expression values are given in Counts per Million (CPM). Quantile normalization has been done as part of the data processing pipeline.

### Question 1: Data Inspection and Basic Manipulation

#### Q1.1 Importing the data and getting familiar with it (1 POINT)

-   Read the datasets into R-Studio.
-   How many genes are there?
-   How many samples are there?
-   How many factors are there? How may levels per factor? List them out.

#### Q1.2 Data manipulation (2 POINTS)

The time\_point factor levels actually refer to points on a continous axis; it is not strictly a categorical variable. In order to make graphing easier, it is helpful to restore the numeric property of this variable.

-   Create a new column in the samples metadata tibble. Call it "age" and populate it with the appropriate numeric values. Hint: Assume that the mouse gestation length is 18 days (ie. P0 = 18).

#### Q1.3 Single gene graphing (2 POINTS)

-   Find the expression profile for the gene **Vegfa**. Make a scatterplot with age on the x-axis and expression value in CPM on the y-axis. Color the data points by cell\_type. Add in regression lines separately for each cell\_type.
-   Is there sign of interaction between cell\_type and age for **Vegfa**? Explain.

### Question 2: Assessing overall data quality

#### Q2.1 Overall distributions (2 POINTS)

-   The expression values are currently in CPM. Log2 transform them so that the distribution is more evenly spread out and can be examined more easily.
-   Examine the distribution of gene expressions across all samples using 1. box plots and 2. overlapping density plots.
    -   For the box plots, samples should be on the x-axis and expression should be on the y-axis.
    -   For the overlapping density plots, expression should be on the x-axis and density should be on the y-axis. Lines should be colored by sample (i.e. one line per sample).
    -   Hint: There are a number of data manipulation steps required. Look at the melt() function in reshape2.
-   Which two samples stand out as different, in terms of the distribution of expression values, compared to the rest?

#### Q2.2 How do the samples correlate with one another? (2 POINTS)

-   Examine the correlation between samples using one or more heatmaps. Again, use the log2 transformed expression values. Display cell\_type, organism\_part, age, and batch for each sample in the heatmap. Hint: Consider using pheatmap() with annotations.
-   Among the factors cell\_type, organism\_part, age, and batch, which one seems to be most strongly correlated with clusters in gene expression data? Hint: Consider using 'cluster\_rows=TRUE' in pheatmap().
-   There is a sample whose expression values correlate with the samples of the different cell\_type just as well as with the samples of the same cell\_type. Identify this sample by its ID.

### Question 3: Using PCA to dig deeper

#### Q3.1 Perform PCA to summarize the samples. (2 POINTS)

-   Perform PCA. Hint: Use svd() and remember to scale and center the expression data for **all genes** by using scale() and t().
-   Make a bar graph showing the amount of variance explained by each PC.

#### Q3.2 Confirm your suspicion. Is cell\_type the biggest contributor to variation? (2 POINTS)

-   Which PCs seem to be associated with the cell\_type variable? Hint: Explore using scatterplots.
-   Quantitatively assess the association of the cell\_type variable with all PCs up to PC10. Hint: Fit a linear model and look at the coefficient of determination (R-Squared); make a bar graph where x is PCs and y is R-Squared.

#### Q3.3 Characterizing batch effects. (2 POINTS)

-   Quantitatively assess the association of the batch variable with all PCs up to PC10. Hint: Fit a linear model and look at the coefficient of determination (R-Squared); make a bar graph where x is PCs and y is R-Squared.
-   How much of the variation captured by PC2 is due to batch effects? How does that compare with PC1?
