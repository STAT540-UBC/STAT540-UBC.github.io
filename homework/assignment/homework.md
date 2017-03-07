## THE assignment (fun!) 

In this assignment, the questions are structured to test your knowledge on some statistical concepts and R techniques, such as data visualization, data wrangling, etc. If you find yourself referring to the seminars and perhaps googling things a lot, you're on the right track. Asking the TAs for help during the seminars is also a good option, too! 

Also check out [this list of tips](http://stat540-ubc.github.io/subpages/assignment_tips.html).

## Evaluation 

The assignment has a total of 30 points and represents 30% of your overall course mark. Overall presentation and mechanics is worth 5 points, with full points awarded if it's exceptional. Other 25 points are spread among the specific questions, as detailed below. We will try to give partial credit when/if we can figure out what you're trying to do. 

Each student must submit their own work: no group efforts. It's fine to talk to fellow students and have discussion on the Discussion repository. If someone is really helpful, mention them in your homework to give them some credit. It's crossing the line to to copy and paste any code or wording. 

## Submission 

Please refer to the [submission instructions](http://stat540-ubc.github.io/subpages/assignments.html#stat-540-homework-submission-instructions) on the website. Questions 1-5 are due on **March 1** and the second half is due on **March 22**. Late submission will cost 5 points per day. 

## About R session 

It might not matter but just in case: this assignment assumes that you have 

- R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree" 
- Bioconductor 3.2 
- limma 3.26.7
- edgeR 3.12.0

Note that `limma`, `edgeR` need to be download via bioconductor (eg `biocLite("limma")`; `biocLite("edgeR")`), the rest are installed by runnnig `install.packages()` (eg `install.packages("gridExtra")`, etc)

Also another tip: if you're using RStudio, you can use `tab` to autocomplete any filepath, variable names, function names, etc. 

---

> For questions 1-5, you will be analyzing a publicly-available gene expression study of normal human bronchial epithelial (NHBE) cells, run on the Affymetrix GeneChip Human Genome U133 Plus 2.0 Array. Please briefly read about the study: [Time course of NHBE cells exposed to whole cigarette smoke (full flavor)](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10718). We will use "probes" and "genes" interchangeably, and also "agent" and "treatment".

---



### Question 1: Data inspection 

#### 1.1 Download and inspect the data (1 point)

The transcriptome data and the meta-data for the design details can be downloaded here : [data](homework_data/NHBE_transcriptome_data.txt.gz?raw=true), [meta-data](homework_data/NHBE_design.txt?raw=true). The data has been preprocessed and it is on a log2 scale. One sample with NA values on all probes is removed from both the data and the meta-data, so overall it’s a little different than the processed data provided via GEO.

- Please load and inspect the data and the meta-data. Fiddle with factor levels and do any other light cleaning you deem necessary.
- How many genes are in our data? How many samples? 
- Check the breakdown for treatment and time. Do you think it's a good experiment design?
	- Hint: `table()` and `addmargins()` will help.

#### 1.2 Basic data manipulation (1 point)

Note that the current variable time in this dataset is a factor (time column of design.txt). Create a new variable which is a numeric version of time (basically a new column in the metadata's dataframe). It should represent a numeric count of hours. You will need this new quantitative variable later in the analysis.
 
#### 1.3 Basic graphing (1 point)

Create a plot showing the gene expression data for a random probe and the averages for all possible combinations of agent and time. 

Hint: Choose a random probe first. Use position, panels or facets, and color to convey the treatment, time of the samples. Report the average expression of the selected probe for all possible combinations of treatment and time. Treat time like a factor in this question.

---



### Question 2: Assessing data quality 

After loading the data, the second must-do thing is to sanity check the data.

#### 2.1 Examine the sample-to-sample correlations in a heatmap (2 points)

- Create the following sample-sample correlation matrices. 
	- Order the samples by `time`; within each time group, sort on the other factor, i.e., `agent` (treatment). 
	- Order the samples by `agent`; and within each agent group, sort by `time`.  
- Interpret your results. What does the sample correlation matrix tell us about the overall impact of time and agent?

#### 2.2 Assess the presence of outlier samples (1 point)

- Comment on the potential presence of outliers (if any) based on the sample correlation matrices created in Q2a.
- Try to go beyond merely eyeballing, and make a quantitative statement about presence of an outlier: e.g., quantify for each sample, whether it “sticks out” compared to the other samples.  
- If any sample does “stick out” from the previous step, examine it in the context of its experimental group (e.g., given a sample that is least similar to all other samples: does it correlate with samples treated with same agent better than other samples treated with a different agent?)

---



### Question 3: Differential expression with respect to treatment

#### 3.1 Linear model (1 point)

- Fit a linear model, modeling expression level of each probe using treatment as a single covariate.  (Hint: Lecture 10 on limma)

- Write out in English and as an equation the model you are using to assess differential expression. In the context of that model, what statistical test are you performing?

#### 3.2 Look at the hits (1 point)

- How many hits (probes) are associated with treatment at unadjusted p-value 1e-3? How may are associated with treatment at FDR 0.05?
- Take the top 50 probes as your “hits” and create a heatmap of their expression levels. Sort the hits by p-values and the samples by treatment.
- What is the (estimated) false discovery rate of this "hits" list? How many of these hits do we expect to be false discoveries?

---



### Question 4: Differential expression with respect to time 

For this question, time is treated as a quantitative covariate (unit: hours). 

#### 4.1 Linear model (1 point)

- You know the drill! Fit a linear model for assessing the effect of time on gene expression
- How many hits are associated with time at unadjusted p-value 1e-3? At FDR 0.05?

---



### Question 5: Differential expression analysis with a full model

Like question 4, treat time as a quantitative variable in this question. We're using both treatment and time as covariates in the full model. 

#### 5.1 Quantify the number of hits for treatment (2 point)

- For how many probes is treatment a significant factor at the unadjusted p-value 1e-3, and at FDR 0.05 level?

- Is this number different from what you reported in 3.2? Why? Quantify the proportion of overlapping probes among your hits, when using the unadjusted p-value threshold of 1e-3.

- Plot the distributions of all the p-values for treatment when using both models, i.e., one from the model in Q3 and one from the full model in this question. Compare and comment on the similarity/differences in the shape of the distributions.

#### 5.2 Test the null hypothesis (1 point)

Null hypothesis: there's no significant interaction between time and treatment.

- Explain in English what you are modeling with this interaction term (what does it represent?).
- For how many probes is the interaction effect significant at the unadjusted p-value 1e-3, and at FDR 0.05 level?

#### 5.3 Plot a few probes where the interaction does and does not matter (2 point)

Plot the expression levels for each sample on the y-axis and time on the x-axis. Color the points based on treatment group. Include the fitted regression lines in your plots.

---


> For questions 6 - 8, you will work with data originally published in [A comprehensive comparison of RNA-Seq-based transcriptome analysis from reads to differential gene expression and cross-comparison with microarrays: a case study in _Saccharomyces cerevisiae_.](http://dx.doi.org/10.1093%2Fnar%2Fgks804). The authors used two different platforms -- microarrays and RNA-Seq -- to obtain gene expression data for yeast grown in two conditions: batch medium and chemostat. They then compared the results of differential expression analysis (DEA) across platforms. We will work with their data obtained from the NCBI SRA and GEO repositories. Because of some differences in pre-processing your results will differ from theirs.

---



### Question 6: Microarray analysis 

The six samples in this study, 3 replicates for each of the two conditions, were analysed on the Affymetrix Yeast Genome Array 2.0 platform. We have already downloaded the raw CEL files from GEO and normalised them. The normalized data is saved in the file [GSE37599-data.tsv](homework_data/GSE37599-data.tsv?raw=true)

#### 6.1 Data loading and QC (2 points)

The labels on two of the samples have been swapped, that is one of the batch samples has been labelled as chemostat and vice-versa. 

- First, load the normalized data. What are dimensions of the dataset? In addition to reporting number of rows and columns, make it clear what rows and columns represent and how you're interpreting column names.

- Produce the plots described below and explain how they allow you to identify the swapped samples.
	- (High volume) scatter plot matrix
	- A heatmap of the first 100 genes
	- Compute the Pearson correlation of the samples and plot the results using a heatmap.
	- Scatterplot the six data samples with respect to the first two principal components and label the samples.

#### 6.2 Microarray DEA (2 points)

- Fix the label swap identified in the previous question. 

- Now use this data to do a differential expression analysis with `limma`. Package these results in a data frame with six columns:
	- probe.id - The array probe id.
	- gene.id - The id of the gene which the probe overlaps (see below).
	- p.value - The raw p-value for the probe.
	- q.value - The BH corrected p-value, aka the q-value.
	- log.fc - The log fold change which is the column called "logFC" in the limma results table.
	- test.stat - The test statistics which for limma is the moderated t statistic. This is the column called "t" in the limma results table.

>  You can retrieve the gene IDs using the `yeast2.db` package from Bioconductor:

```r
source("http://bioconductor.org/biocLite.R")    # Installing biocLite package
biocLite() 										# Load the package 
biocLite("yeast2.db") 							# Download the package using biocLite
library("yeast2.db") 							# Load the yeast2.db package 
mget("1772391_at", yeast2ORF)  					# See what gene is that probe for 
gene.ids <- unlist(mget(probe.ids, yeast2ORF))  # See all the genes for the list of probes
```

#### 6.3 Microarray DEA continue (1 point)

- Remove any rows with probes which don't map to genes. You'll be able to find these because they will have `NA` as their gene id. 
- How many probes did we start with and how many remain after removing probes without gene ids?
- Illustrate the differential expression between the batch and the chemostat samples for the top hit (i.e., probe with the lowest p- or q-value).
- How many probes are identified as differentially expressed at a false discovery rate (FDR) of 1e-5 (note: this is a FDR cutoff used in the original paper)?
- Save this data.frame for later with `write.table()`

> When using write.table to save the data, you need to pass the arguments `row.names = TRUE, col.names = NA` to make sure the headers are written correctly.

---



### Question 7: RNA-seq analysis 

We have aligned the RNA-Seq library using the [Stampy](http://www.well.ox.ac.uk/project-stampy) aligner and generated count data. This particular RNA-Seq dataset was sequenced to quite a high depth. The data file is available as [stampy.deep.counts.tsv](homework_data/stampy.deep.counts.tsv?raw=true). In this question you will use this data to do a differential expression analysis.

#### 7.1 Load RNA Count Data and Sanity Check (1 point)

Load the count data using `read.table`; you will need to pass the arguments `header=TRUE` and `row.names=1`. 

- What are dimensions of the dataset? In addition to reporting number of rows and columns, make it clear what rows and columns represent. What is the difference between the rows of this dataset versus rows of the array data in question 1a?
- Do a sanity check to make sure there is no sample swap by plotting a heatmap of the sample correlations.


#### 7.2 DEA of deep sequencing data (1 point)

Now you will use `edgeR` to identify differentially expressed genes between the batch medium vs. chemostat conditions.

- Recall that `edgeR` needs to estimate the dispersion parameter in the negative binomial model using an empirical Bayes method. Estimate the dispersion parameters using `estimateGLMCommonDisp`, `estimateGLMTrendedDisp` and `estimateGLMTagwiseDisp`. 

- Use the glm functionality of `edgeR`, i.e. use the `glmFit` function, to identify differentially expressed genes between conditions. 

- How many genes are differentially expressed between conditions at a false discovery rate (FDR) of 1e-5? 

- Package these results in a data.frame called 'edger.deep.results' with five columns:
	- gene.id - The id of the gene which reads were aligned to.
	- p.value - The raw p-value for the gene.
	- q.value - The BH corrected p-value, aka the q-value.
	- log.fc - The log fold change which is the column called "logFC" in the `edgeR` results table.
	- test.stat - The test statistic, which for `edgeR` is a likelihood ratio. This is the column called "LR" in the `edgeR` results table.

- Save this data.frame for later with `write.table()`

#### 7.3 DEA of low sequencing data (1 point)

Sequencing depth is a major factor in RNA-Seq analysis; as rule of thumb, the deeper you sequence the more power you have. To investigate this claim, we have created a low depth dataset by randomly selecting 1% of the reads and discarding the remaining 99%. This low depth count data was also aligned with Stampy.

- Load the low depth count data using `read.table`; you will need to pass the arguments `header=TRUE` and `row.names=1`. The data file is available as [stampy.low.counts.tsv](homework_data/stampy.low.counts.tsv?raw=true)

- Repeat the same steps as in 7.2 except package your results in a data.frame called "edger.low.results". Answer the same questions. 


#### 7.4 Deep vs low sequencing (1 point)

Create a Venn diagram showing all genes identified as differentially expressed (at FDR of 1e-5) in the two previous RNA-Seq analyses and answer the following questions based on the diagram: 

- How many genes were identified by edgeR in both low and deep count data?
- How many genes were identified in all the analyses?
- Comment on the effect of sequencing depth on the DEA results. 

---



### Question 8: Compare DEA results from RNA-Seq and arrays

In question 6, you performed a DEA of array data using limma. In question 7, you used edgeR to perform DEA of RNA-Seq data. In this question, you will compare the results of RNA-Seq DEA with those of the array DEA.

> Remember that you've packaged your results and saved them using write.table(). If the data.frames containing the results of these analyses are no longer in your workspace, you can load them using read.table().


#### 8.1 Plots of interesting and boring genes (2 point)

Plot the gene expression (i.e., from array data) and the logged counts (i.e., from deep RNA-Seq data) of: two interesting genes identified as DE by both analyses; one DE gene identified as DE only in the array analysis; one DE gene only in the edgeR analysis; one boring gene in both analyses (i.e., 5 genes total measured with 2 platforms)

Hint: Write a plotting function and reuse it 

---


### Question 9 (OPTIONAL): Try another package 

In bioinformatics, there are so many tools and methods out there to analyze the same data set for the same purpose. So in this optional question, you would use `voom+limma` to identify differentially expressed genes between the batch medium vs. chemostat conditions. 

#### 9.1 DEA with `voom+limma` on deep seq data 
- `voom` normalizes the counts before it converts counts to log2-cpm. Use `calcNormFactors` to normalize counts.
- Use `voom` to convert count data into logged CPM data and then use `limma` to identify differentially expressed genes between conditions. Package these results in a data.frame called 'voom.limma.deep.results' with five columns:
  - gene.id - The id of the gene which reads were aligned to.
  - p.value - The raw p-value for the gene.
  - q.value - The BH corrected p-value, aka the q-value.
  - log.fc - The log fold change which is the column called "logFC" in the `edgeR` results table.
  - test.stat - The test statistic, which is the column called "t".
- How many genes are differentially expressed between conditions at a false discovery rate (FDR) of 1e-5?

#### 9.2 DEA with `voom+limma` on low seq data 

Repeat 9.1 with the low seq data and store your result in `voom.limma.low.results` 

#### 9.3 Compare your results with edgeR

Now that we have the results of the differential expression analysis performed by `voom+limma` and `edgeR` methods on both low and deep count data, we are going to compare and illustrate the results.

Create a Venn diagram showing all genes identified as differentially expressed (at FDR of 1e-5) in the four previous RNA-Seq analyses. If the comparison of 4 sets gets very confusing, you can also create different pairs of Venn diagrams of interest.
- How many genes were identified by `voom+limma` in both low and deep count data? 
- How many genes were identified by `edgeR` and `voom+limma` in deep count data? 
- What can you say about these two different methods? 

**Congrats! You're done :clap: Don't forget to open an issue to submit the homework. See [submission instructions](http://stat540-ubc.github.io/subpages/assignments.html#stat-540-homework-submission-instructions)** 

