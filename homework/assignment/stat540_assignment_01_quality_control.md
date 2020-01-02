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
-   How many factors are in our experimental design? How may levels per factor? List out the levels for each factor.

#### Q1.2 Data manipulation (2 POINTS)

The levels of the factor `time_point` actually refer to points on a continous axis. In other words, it doesn't have to be interpreted as strictly categorical variable. In order to make graphing easier, it will be helpful to convert this variable to a numeric representation.

-   Create a new column in the samples metadata tibble. Call it "age" and populate it with the appropriate numeric values. Hint: Assume that the mouse gestation length is 18 days (ie. P0 = 18).

#### Q1.3 Single gene graphing (2 POINTS)

-   Find the expression profile for the gene **Vegfa**. Make a scatterplot with age on the x-axis and expression value in CPM on the y-axis. Color the data points by cell\_type. Add in a regression line for each cell type.

-   Is there sign of interaction between cell\_type and age for **Vegfa**? Explain using what you observed in your graph from the previous question.

### Question 2: Assessing overall data quality

#### Q2.1 Overall distributions (2 POINTS)

-   The expression values are currently in CPM. Log2 transform them so that the distribution is more evenly spread out and can be examined more easily.
-   Examine the distribution of gene expression across all samples using 1. box plots and 2. overlapping density plots.
    -   For the box plots, samples should be on the x-axis and expression should be on the y-axis.
    -   For the overlapping density plots, expression should be on the x-axis and density should be on the y-axis. Lines should be colored by sample (i.e. one line per sample).
    -   Hint: There are a number of data manipulation steps required. Look at the melt() function in reshape2.
-   Which two samples stand out as different, in terms of the distribution of expression values, compared to the rest?

#### Q2.2 How do the samples correlate with one another? (2 POINTS)

-   Examine the correlation **between samples** using one or more heatmaps (i.e. samples should be on the x axis and the y axis, and the values in the heatmap should be correlations). Again, use the log2 transformed expression values. Display cell\_type, organism\_part, age, and batch for each sample in the heatmap. Hint: Consider using pheatmap() with annotations and cor to correlate gene expression between each pair of samples.
-   Among the factors cell\_type, organism\_part, age, and batch, which one seems to be most strongly correlated with clusters in gene expression data? Hint: Consider using 'cluster\_rows=TRUE' in pheatmap().
-   There is a sample whose expression values correlate with the samples of the different cell\_type just as well as with the samples of the same cell\_type. Identify this sample by its ID.

### Question 4: Conducting differential expression analysis
-------------------------------------------------------

#### 4.1 Remove lowly expressed genes (2 POINTS)

-   Remove lowly expressed genes by retaining genes that have CPM &gt; 1 in at least as many samples as the *smallest group size* (i.e use table() to identify the number of samples belonging to each treatment group. The *smallest group size* is the smallest number in that table).
-   How many genes are there after filtering?

#### 4.2 Construct linear model (2 POINTS)

-   Use limma-voom to fit a linear model with cell type, organism part, age and the interaction between age and cell type as covariates (hint: use lmFit, voom and eBayes). Before you do this, reformat the data frame so that gene IDs are row names, and not a column (limma requires the dataset in this format).


#### 4.3: Interpret model (3 POINTS)

-   For the gene Eva1a, what is the numeric value of the coeffcient of the age term? What does it mean?
-   Write down an equation describing the linear model you specified in 4.2. Hint: your equation should have a term for each column of the design matrix.
-   When you create a linear model, what is the underlying assumptions you are making about the distribution of the data (this question is multiple choice)?
    1.  All expression values are normally distributed
    2.  The residuals of the fitted model are normally distributed
    
### Question 5: Evaluating the results
----------------------------------

#### 5.1: Quantifying the number of genes differentially expressed (3 POINTS)

-   Using the linear model defined above, determine the number of genes differentially expressed by cell type at an FDR (use adjust.method = "fdr" in topTable()) less than 0.05.
-   Although an FDR cutoff of 0.05 was used, many of the identified genes have smaller FDRs. By taking an average of the FDR across the set of differentially expressed genes, determine the number of genes that we expect to be false discoveries on average.
-   Use decideTests() to quantify the number of genes that increase, decrease or don't change by cell type, organism part and age. Which variable is associated with the largest number of differentially expressed genes?

#### 5.2: Interpret the interaction term (2 POINTS)

-   Explain what you are modeling with this interaction term. For a particular gene, what does a signifcant interaction term mean?
-   For how many probes is the interaction effect significant (FDR less than 0.05)?

#### **Bonus Question** (2 POINTS)

-   Compare your results to those obtained by Scheffer et al (2015). Discuss any discrepancies. List at least three explanations for these discrepancies.
