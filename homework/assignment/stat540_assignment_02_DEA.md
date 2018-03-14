STAT 540 - Assignment 2 - Differential Expression Analysis
================

Question 4: Conducting differential expression analysis
-------------------------------------------------------

**4.1**: Remove lowly expressed genes (2 POINT)

-   Remove lowly expressed genes by retaining genes that have CPM &gt; 1 in at least as many samples as the *smallest group size* (i.e use table() to identify the number of samples belonging to each treatment group. The *smallest group size* is the smallest number in that table).
-   How many genes are there after filtering?

**4.2**: Construct linear model (2 POINTS)

-   Use limma-voom to fit a linear model with cell type, organism part, age and the interaction between age and cell type as covariates (hint: use lmFit, voom and eBayes). Before you do this, reformat the data frame so that gene IDs are row names, and not a column (limma requires the dataset in this format).

**4.3**: Interpret model (3 POINTS)

-   For the gene Eva1a, what is the numeric value of the coeffcient of the age term? What does it mean?
-   Write down an equation describing the linear model you specified in 4.2. Hint: your equation should have a term for each column of the design matrix.
-   When you create a linear model, what is the underlying assumptions you are making about the distribution of the data (this question is multiple choice)?
    1.  All expression values are normally distributed
    2.  The residuals of the fitted model are normally distributed

Question 5: Evaluating the results
----------------------------------

**5.1**: Quantifying the number of genes differentially expressed (3 POINTS)

-   Using the linear model defined above, determine the number of genes differentially expressed by cell type at an FDR (use adjust.method = "fdr" in topTable()) less than 0.05.
-   Although an FDR cutoff of 0.05 was used, many of the identified genes have smaller FDRs. By taking an average of the FDR across the set of differentially expressed genes, determine the number of genes that we expect to be false discoveries on average.
-   Use decideTests() to quantify the number of genes that increase, decrease or don't change by cell type, organism part and age. Which variable is associated with the largest number of differentially expressed genes?

**5.2**: Effect of age - visualize top differentially expressed genes (2 POINTS)

-   Take the top 50 genes differentially expressed by age and create a heatmap of their expression levels in logCPM. Sort the genes by p-values and group the samples by time point.

**5.3**: Interpret the interaction term (2 POINTS)

-   Explain what you are modeling with this interaction term. For a particular gene, what does a signifcant interaction term mean?
-   For how many probes is the interaction effect significant (FDR less than 0.05)?

**5.4**: Plot three genes where the interaction does matter (1 POINT)

-   Plot the top three genes with the most significant interaction term. Make a scatterplot with log CPM on the y-axis and age on the x-axis. The shape of each point should correspond to the organism part, and the cell type should correspond to the colour. Note: some conditions have multiple samples, so it is best to plot the mean expression of each treatment group.

**Bonus Question** (2 POINTS)

-   Compare your results to those obtained by Scheffer et al (2015). Discuss any discrepancies. List at least three explanations for these discrepancies.
