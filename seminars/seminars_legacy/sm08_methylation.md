DNA methylation analysis with Illumina's Infinium HumanMethylation450K array
========================================================

Contributors: Gloria Li

> This seminar is inspired by a STAT 540 project in 2013: *Analysis of Gene Expression Omnibus Leukaemia Data from the Illumina HumanMethylation450 Array* by Alice Zhu, Rachel Edgar, Shaun Jackman and Nick Fishbane. See their project website [here](https://sites.google.com/site/stat540diffmethleuk/).   

The Illumina HumanMethylation450 Array is a microarray based high throughput platform for methylation profiling on 450,000 pre-selected probes for CpGs across the human genome. A great number of datasets have been made publicly available through Gene Expression Omnibus (GEO). Each dataset in GEO has a unique GSE ID that we can use to retrieve it from GEO. There are several Bioconductor packages that are very useful for obtaining data from GEO and microarray analysis. If you haven't done so yet, install and load these packages.   


```r
source('http://bioconductor.org/biocLite.R')
biocLite('GEOquery')
biocLite('wateRmelon')
biocLite("IlluminaHumanMethylation450k.db")
```

```r
library(GEOquery)
library(wateRmelon)
library(IlluminaHumanMethylation450k.db)
```

## Explore 450k methylation array data
In this seminar, we are going to perform differential methylation analysis for Acute lymphoblastic leukemia (ALL) cases and healthy B cells as control group. The datasets we are going to use are: GSE39141 which contains 29 ALL bone marrow samples and 4 healthy B cell samples, and to make up for the small sample size for healthy cells,  we will also use GSE42865 which contains 9 healthy B cell samples and 7 samples with other conditions.    

First, let's retrieve our datasets from GEO with `getGEO` from `GEOquery` package. Warning: this may take several minutes! So to avoid re-downloading in the future, save the data once you get it into a good shape.    

As you can see, the datasets we got are large lists including Beta values as well as metadata. We can extract this information separately. Then we can isolate a subset of GSE42865 containing only the healthy cells.    


```r
if(file.exists("methyl_ALL.Rdata")){ # if previously downloaded
  load("methyl_ALL.Rdata")
} else { # if downloading for the first time
  GSE39141 <- getGEO('GSE39141')
  show(GSE39141) ## 33 samples (29 ALL and 4 healthy B cells)
  GSE42865 <- getGEO('GSE42865') # took ~2 mins for JB
  show(GSE42865) ## 16 samples (9 healthy cells B cells and 7 other cells)
  
  # Extract expression matrices (turn into data frames at once) 
  ALL.dat <- as.data.frame(exprs(GSE39141[[1]]))
  CTRL.dat <- as.data.frame(exprs(GSE42865[[1]]))
  
  # Obtain the meta-data for the samples and rename them perhaps?
  ALL.meta <- pData(phenoData(GSE39141[[1]]))
  CTRL.meta <- pData(phenoData(GSE42865[[1]]))
  
  # create some labels
  ALL.meta$Group<- c(rep('ALL', 29), rep('HBC', 4)) 
  ## ALL: Case; HBC: Healthy B Cells
  
  # Subset both meta-data and data for control (healthy) donors
  CTRL.meta <- droplevels(subset(CTRL.meta,
                                 grepl("Healthy donor", characteristics_ch1.1)))
  CTRL.dat <- subset(CTRL.dat, select = as.character(CTRL.meta$geo_accession)) 
  
  # Rename variables
  names(ALL.dat) <- paste(ALL.meta$Group,
                          gsub("GSM", "", names(ALL.dat)), sep = '_')
  names(CTRL.dat) <- paste('HBC', gsub("GSM", "", names(CTRL.dat)), sep = '_')
  
  # save the data to avoid future re-downloading
  save(ALL.dat, CTRL.dat, ALL.meta, CTRL.meta, file = "methyl_ALL.Rdata")
}
```

Now, we can do some exploratory analysis of the data, for examples, looking at distribution of Beta values in each sample or each probe. Here is a density plot of average Beta values for probes in the two datasets. See if you get similar results. 
> Hint: you can use `rowMeans` function to calculate the average Beta value for each probe. Don't forget to use `na.rm = T` to exclude NAs. 

![](sm08_methylation_files/figure-html/exploratory-1.png) 

## Normalization
You can see the distribution of Beta values are not the same between different experiments. That is why it is important to normalize the data before proceed. The Bioconductor package `wateRmelon` offers 15 methods for 450K methylation array normalization. If you are interested in comparing these methods, check out [this paper](http://www.biomedcentral.com/1471-2164/14/293). For the sake of convenience, we will just use the simplest method, `betaqn`, which performs quantile normalization on the beta values.        


```r
# combine data from two experiments into one matrix, each column represents beta
# values of one sample
beta.matrix <- as.matrix(cbind(ALL.dat, CTRL.dat))
str(beta.matrix, max.level = 0)
```

```
##  num [1:485577, 1:42] 0.512 0.911 0.857 0.149 0.729 ...
##  - attr(*, "dimnames")=List of 2
```

```r
# quantile normalization
system.time(beta.norm <- betaqn(beta.matrix))
```

```
##    user  system elapsed 
##  36.186   0.228  36.385
```
![](sm08_methylation_files/figure-html/postNorm-1.png) 

## M values
As you can see, after normalization, beta values from the different experiments have more similar distributions. However, this distribution is on the interval [0, 1], with two modes near the two extremes. This type of distribution is not very compatible with the typical assumptions of linear models. A common solution is to apply a logit transformation on the data to convert it to a continuous variable that spans $(-\infty, \infty)$, i.e. compute the M value. The `beta2m()` function in the `wateRmelon` package does this conversion. 
$$M = log({Beta\over 1-Beta})$$


```r
M.norm <- beta2m(beta.norm)
```

## CpG Islands
Now we can go ahead with differential methylation analysis -- cases vs. controls or leukemia vs. healthy. You can perform this analysis on each probe, but in this seminar, we will aggregate the probes into CpG Islands (CGIs), i.e. the CpG dense regions of the genome, and then detect differentially methylated CGIs. The biological function of CGIs is better studied, so the interpretation of our results will be easier if we focus on CGIs. Conveniently, the Bioconductor package `IlluminaHumanMethylation450k.db` provides all sorts of annotation information for the 450K methylation array including the association of probes to CGIs. We will extract this information and use the mean of M values for all probes in a specific CGI to represent its methylation level. Check if you get the same boxplot for CGI M values.    


```r
# Extracting probe ID to CpG islands association
cginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
names(cginame) <- c('Probe_ID', 'cginame')
rownames(cginame) <- cginame$Probe_ID
length(levels(factor(cginame$cginame)))   # No. of CGIs
```

```
## [1] 27176
```

```r
# restrict probes to those within CGIs
beta.inCGI <- beta.norm[cginame$Probe_ID, ]
M.inCGI <- M.norm[cginame$Probe_ID, ]
nrow(M.inCGI)  # No. of probes within CGIs
```

```
## [1] 309465
```

```r
# aggregate probes to CGIs
beta.CGI <- aggregate(beta.inCGI, by = list(cginame$cginame), mean, na.rm = T)
rownames(beta.CGI) <- beta.CGI[, "Group.1"]
beta.CGI <- subset(beta.CGI, select = - Group.1)
str(beta.CGI, max.level = 0)
```

```
## 'data.frame':	27176 obs. of  42 variables:
```

```r
M.CGI <- aggregate(M.inCGI, by = list(cginame$cginame), mean, na.rm = T)
rownames(M.CGI) <- M.CGI[, "Group.1"]
M.CGI <- subset(M.CGI, select = - Group.1)
str(M.CGI, max.level = 0)
```

```
## 'data.frame':	27176 obs. of  42 variables:
```
![](sm08_methylation_files/figure-html/M.CGI.boxplot-1.png) 

## Differential methylation analysis with limma
Next, we can use a linear model to identify differentially methylated CGIs with `limma`. You are already familiar with this part.    


```r
library(limma)
design <-
  data.frame(Group = relevel(factor(gsub("_[0-9]+", "", colnames(M.CGI))),
                             ref = "HBC"), row.names = colnames(M.CGI))
str(design)
```

```
## 'data.frame':	42 obs. of  1 variable:
##  $ Group: Factor w/ 2 levels "HBC","ALL": 2 2 2 2 2 2 2 2 2 2 ...
```

```r
(DesMat <- model.matrix(~ Group, design))
```

```
##             (Intercept) GroupALL
## ALL_956761            1        1
## ALL_956762            1        1
## ALL_956763            1        1
## ALL_956764            1        1
## ALL_956765            1        1
## ALL_956766            1        1
## ALL_956767            1        1
## ALL_956768            1        1
## ALL_956769            1        1
## ALL_956770            1        1
## ALL_956771            1        1
## ALL_956772            1        1
## ALL_956773            1        1
## ALL_956774            1        1
## ALL_956775            1        1
## ALL_956776            1        1
## ALL_956777            1        1
## ALL_956778            1        1
## ALL_956779            1        1
## ALL_956780            1        1
## ALL_956781            1        1
## ALL_956782            1        1
## ALL_956783            1        1
## ALL_956784            1        1
## ALL_956785            1        1
## ALL_956786            1        1
## ALL_956787            1        1
## ALL_956788            1        1
## ALL_956789            1        1
## HBC_956790            1        0
## HBC_956791            1        0
## HBC_956792            1        0
## HBC_956793            1        0
## HBC_1052420           1        0
## HBC_1052421           1        0
## HBC_1052422           1        0
## HBC_1052423           1        0
## HBC_1052424           1        0
## HBC_1052425           1        0
## HBC_1052426           1        0
## HBC_1052427           1        0
## HBC_1052428           1        0
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$Group
## [1] "contr.treatment"
```

```r
DMRfit <- lmFit(M.CGI, DesMat)
DMRfitEb <- eBayes(DMRfit)
cutoff <- 0.01
DMR <- topTable(DMRfitEb, coef = 'GroupALL', number = Inf, p.value = cutoff)
head(DMR)   # top hits 
```

```
##                               logFC   AveExpr         t      P.Value
## chr19:49828412-49828668   1.3084022  2.283294  12.05206 2.798199e-15
## chr4:156680095-156681386 -1.1115033 -2.521386 -10.94730 6.244596e-14
## chr20:11871374-11872207  -1.2638610 -2.284080 -10.48924 2.370313e-13
## chr19:33726654-33726946   0.9428988  1.886933  10.39039 3.172261e-13
## chr18:77166704-77167043   0.8103505  3.198909  10.36406 3.428991e-13
## chr18:46447718-46448083  -0.8990291  2.034496 -10.31377 3.979679e-13
##                             adj.P.Val        B
## chr19:49828412-49828668  7.604386e-11 24.30792
## chr4:156680095-156681386 8.485157e-10 21.38964
## chr20:11871374-11872207  1.802529e-09 20.12796
## chr19:33726654-33726946  1.802529e-09 19.85172
## chr18:77166704-77167043  1.802529e-09 19.77792
## chr18:46447718-46448083  1.802529e-09 19.63664
```

So using a cutoff of FDR = 0.01, we identified 4115 CGIs that are differentially methylated between ALL and control group. Now we can make some plots to check these hits. 

First, heatmap of beta values of top 100 hits.    


```r
library(gplots)
DMR100 <- topTable(DMRfitEb, coef = 'GroupALL', number = 100)
DMR.CGI <- t(as.matrix(subset(beta.CGI,
                              rownames(beta.CGI) %in% rownames(DMR100))))
str(DMR.CGI, max.level = 0)
```

```
##  num [1:42, 1:100] 0.751 0.735 0.72 0.702 0.737 ...
##  - attr(*, "dimnames")=List of 2
```

```r
col <- c(rep("darkgoldenrod1", times = nrow(DMR.CGI))) 
col[grepl("HBC", rownames(DMR.CGI))] <- "forestgreen"
op <- par(mai = rep(0.5, 4))
heatmap.2(DMR.CGI, col = redblue(256), RowSideColors = col,
          density.info = "none", trace = "none", Rowv = TRUE, Colv = TRUE,
          labCol = FALSE, labRow = FALSE, dendrogram="row",
          margins = c(1, 5))
legend("topright", c("ALL", "HBC"),
       col = c("darkgoldenrod1", "forestgreen"), pch = 15)
```

![](sm08_methylation_files/figure-html/heatmap-1.png) 

```r
par(op)
```

Next, stripplot of beta values of probes within top 5 CGI hits.   


```r
DMR5 <- topTable(DMRfitEb, coef = 'GroupALL', number = 5)
beta.DMR5probe <-
  beta.inCGI[cginame[rownames(beta.inCGI),]$cginame %in% rownames(DMR5),]
beta.DMR5probe.tall <-
  melt(beta.DMR5probe, value.name = 'M', varnames = c('Probe_ID', 'Sample'))
beta.DMR5probe.tall$Group <-
  factor(gsub("_[0-9]+", "", beta.DMR5probe.tall$Sample))
beta.DMR5probe.tall$CGI <-
  factor(cginame[as.character(beta.DMR5probe.tall$Probe_ID),]$cginame)
(beta.DMR5.stripplot <-
   ggplot(data = beta.DMR5probe.tall, aes(x = Group, y = M, color = Group)) + 
   geom_point(position = position_jitter(width = 0.05), na.rm = T) + 
   stat_summary(fun.y = mean, aes(group = 1), geom = "line", color = "black") + 
   facet_grid(. ~ CGI) + 
   ggtitle("Probe beta values within top 5 DM CGIs") + 
   xlab("Group") + 
   ylab("beta") + 
   theme_bw())
```

Finally, plot location of differential methylated probes along each chromosome.   


```r
# get the length of chromosome 1-22 and X
chrlen <-
  unlist(as.list(IlluminaHumanMethylation450kCHRLENGTHS)[c(as.character(1:22),
                                                           "X")])   
chrlen <- data.frame(chr = factor(names(chrlen)), length = chrlen)
chr <- IlluminaHumanMethylation450kCHR        # get the chromosome of each probe
# get the probe identifiers that are mapped to chromosome
chr <- unlist(as.list(chr[mappedkeys(chr)]))
# get chromosome coordinate of each probe
coord <- IlluminaHumanMethylation450kCPGCOORDINATE   
# get the probe identifiers that are mapped to coordinate
coord <- unlist(as.list(coord[mappedkeys(coord)]))      
coord <- data.frame(chr = chr[intersect(names(chr), names(coord))],
                    coord = coord[intersect(names(chr), names(coord))])
# coordinates of probes in DM CGIs
coordDMRprobe <-
  droplevels(na.omit(coord[cginame[cginame$cginame %in%
                                     rownames(DMR),]$Probe_ID,])) 
(coord.plot <- ggplot(data = coordDMRprobe) + 
   geom_linerange(aes(factor(chr, levels = c("X", as.character(22:1))),
                      ymin = 0, ymax = length), data = chrlen, alpha = 0.5) + 
   geom_point(aes(x = factor(chr,
                             levels = c("X", as.character(22:1))), y = coord),
              position = position_jitter(width = 0.03), na.rm = T) + 
   ggtitle("DMR positions on chromosomes") + 
   ylab("Position of DMRs") +
   xlab("chr") +
   coord_flip() + 
   theme_bw())
```

![](sm08_methylation_files/figure-html/coord-1.png) 

## Interpretation and functional enrichment analysis
The next step in the analysis pipeline would be to interpret our result by associating these differentially methylated regions (DMRs) with biological features. DNA methylation in different genomic regions, i.e. promoter regions, enhancers, gene body etc., has different impact on gene transcription. A common practice in DMR analysis is to separate DMRs associated with different types of genomic regions and study them separately. Our key interest is to see how these DMRs affect gene expression level and consequently biological function. For that purpose, after we separate DMRs into different genomic regions, we will associate them with genes and then biological functions. There are many ways to do this. The `IlluminaHumanMethylation450k.db` package itself contains some annotation information, including associating probes with different types of genomic regions and gene IDs. You can obtain this information similarly to how we extracted chromosome coordinates above. Here is a full list of all objects available in this package. Their names are self-explanatory, but if want more detailed description, check the [package manual document](http://www.bioconductor.org/packages/release/data/annotation/manuals/IlluminaHumanMethylation450k.db/man/IlluminaHumanMethylation450k.db.pdf). 

```r
ls("package:IlluminaHumanMethylation450k.db")
```

```
##  [1] "IlluminaHumanMethylation450k"                  
##  [2] "IlluminaHumanMethylation450kACCNUM"            
##  [3] "IlluminaHumanMethylation450kALIAS2PROBE"       
##  [4] "IlluminaHumanMethylation450kBLAME"             
##  [5] "IlluminaHumanMethylation450kBUILD"             
##  [6] "IlluminaHumanMethylation450kCHR"               
##  [7] "IlluminaHumanMethylation450kCHR36"             
##  [8] "IlluminaHumanMethylation450kCHR37"             
##  [9] "IlluminaHumanMethylation450kCHRLENGTHS"        
## [10] "IlluminaHumanMethylation450kCHRLOC"            
## [11] "IlluminaHumanMethylation450kCHRLOCEND"         
## [12] "IlluminaHumanMethylation450kCOLORCHANNEL"      
## [13] "IlluminaHumanMethylation450kCPG36"             
## [14] "IlluminaHumanMethylation450kCPG37"             
## [15] "IlluminaHumanMethylation450kCPGCOORDINATE"     
## [16] "IlluminaHumanMethylation450kCPGILOCATION"      
## [17] "IlluminaHumanMethylation450kCPGINAME"          
## [18] "IlluminaHumanMethylation450kCPGIRELATION"      
## [19] "IlluminaHumanMethylation450kCPGS"              
## [20] "IlluminaHumanMethylation450k_dbconn"           
## [21] "IlluminaHumanMethylation450k_dbfile"           
## [22] "IlluminaHumanMethylation450k_dbInfo"           
## [23] "IlluminaHumanMethylation450k_dbschema"         
## [24] "IlluminaHumanMethylation450kDESIGN"            
## [25] "IlluminaHumanMethylation450kDHS"               
## [26] "IlluminaHumanMethylation450kDMR"               
## [27] "IlluminaHumanMethylation450kENHANCER"          
## [28] "IlluminaHumanMethylation450kENSEMBL"           
## [29] "IlluminaHumanMethylation450kENSEMBL2PROBE"     
## [30] "IlluminaHumanMethylation450kENTREZID"          
## [31] "IlluminaHumanMethylation450kENZYME"            
## [32] "IlluminaHumanMethylation450kENZYME2PROBE"      
## [33] "IlluminaHumanMethylation450kFANTOM"            
## [34] "IlluminaHumanMethylation450kGENEBODY"          
## [35] "IlluminaHumanMethylation450kGENENAME"          
## [36] "IlluminaHumanMethylation450k_get27k"           
## [37] "IlluminaHumanMethylation450k_getControls"      
## [38] "IlluminaHumanMethylation450k_getProbeOrdering" 
## [39] "IlluminaHumanMethylation450k_getProbes"        
## [40] "IlluminaHumanMethylation450kGO"                
## [41] "IlluminaHumanMethylation450kGO2ALLPROBES"      
## [42] "IlluminaHumanMethylation450kGO2PROBE"          
## [43] "IlluminaHumanMethylation450kISCPGISLAND"       
## [44] "IlluminaHumanMethylation450kMAP"               
## [45] "IlluminaHumanMethylation450kMAPCOUNTS"         
## [46] "IlluminaHumanMethylation450kMETHYL27"          
## [47] "IlluminaHumanMethylation450kNUID"              
## [48] "IlluminaHumanMethylation450kOMIM"              
## [49] "IlluminaHumanMethylation450kORGANISM"          
## [50] "IlluminaHumanMethylation450kORGPKG"            
## [51] "IlluminaHumanMethylation450kPATH"              
## [52] "IlluminaHumanMethylation450kPATH2PROBE"        
## [53] "IlluminaHumanMethylation450kPFAM"              
## [54] "IlluminaHumanMethylation450kPMID"              
## [55] "IlluminaHumanMethylation450kPMID2PROBE"        
## [56] "IlluminaHumanMethylation450kPROBELOCATION"     
## [57] "IlluminaHumanMethylation450kPROSITE"           
## [58] "IlluminaHumanMethylation450kRANDOM"            
## [59] "IlluminaHumanMethylation450kREFSEQ"            
## [60] "IlluminaHumanMethylation450kREGULATORYGROUP"   
## [61] "IlluminaHumanMethylation450kREGULATORYLOCATION"
## [62] "IlluminaHumanMethylation450kREGULATORYNAME"    
## [63] "IlluminaHumanMethylation450kSVNID"             
## [64] "IlluminaHumanMethylation450kSYMBOL"            
## [65] "IlluminaHumanMethylation450kUNIGENE"           
## [66] "IlluminaHumanMethylation450kUNIPROT"
```

There are a lot of resources that provide association between genes and biological function. Gene Ontology (GO) is a popular database that provides a controlled vocabulary of terms for describing gene product characteristics and functions including molecular function, biological process and cellular component. Also, there are databases like KEGG and PANTHER that provide association of genes with pathways for pathway analysis, and InterPro and Pfam etc. that provide protein domain analysis. You might have noticed from the list above that `IlluminaHumanMethylation450k.db` also provides some of these association analyses like GO, but you would have more control and flexibility if you use other tools more specialized for functional enrichment analysis. There are two online tools that I found very helpful, both extremely easy to use. [DAVID](http://david.abcc.ncifcrf.gov/tools.jsp) takes a list of genes and outputs annotation and enrichment analysis from all sorts of databases including all the ones we discussed above. The results are available for downloading in tab delimited text format. Another more recent tool is [GREAT](http://bejerano.stanford.edu/great/public/html/). It was originally developed to associate regulatory regions with genes and biological functions and the required input is a list of coordinates. It also performs functional enrichment analysis from different databases and it can make some fancy plots for output.       

The interpretation of differential methylation analysis results is a very flexible and most often context-specific process. Functional enrichment analysis is just one popular way to approach this but there is really no *gold-standard pipeline* we can follow. Although we won't be able to go into it in this seminar, you should understand that identifying DMRs is really the easy part and digging up the biology behind it and make a story of it would be something you need to spend a lot of time considering.       

## Takehome assignment 

Using plots generated above (or via) the above analysis, describe how and which linear modelling assumptions are imperfectly satisfied. If assumptions are violated, how does this affect the usage of Limma?
