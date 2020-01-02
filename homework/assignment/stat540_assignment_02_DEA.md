STAT 540 - Assignment 2 - Differential Expression Analysis
================

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
