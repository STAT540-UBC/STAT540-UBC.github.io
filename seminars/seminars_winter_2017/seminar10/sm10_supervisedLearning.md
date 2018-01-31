
Supervised learning, classification, cross validation, variable selection
=========================================================================

Contributors: Gabriela Cohen Freue, W. Evan Durno, Jasleen Grewal

Learning objectives
-------------------

By the end of this tutorial, you should be able to
- Filter gene expression data to remove uninformative features.
- Further subselect a set of variables using a fixed information criteria, and know when it may not be required/possible to do so.
- Understand the concept of held-out test set, training set, and validation set.
- Select a suitable classifier (based on size of training dataset) and train a model with the input training data.
- Get predicted labels for a validation set, from the trained classifier model.
- Understand the utility of cross validation in selecting an optimal model
- Be able to identify at what part of the supervised learning process CV should be used. - Explain how cross validation is different from bootstrapping.
- Assess the best selected model's performance on the held out test set.
- Recognize why we can't go back and reselect models after assessment on the held out test set.
- Be able to distinguish between the following scenarios for using supervised machine learning, and enumerate metrics for assessing performance in each case:
- 2 class classification problem
- Multi-clas classification problem
- Be able to define accuracy, precision, sensitivity, specificity, F1-score, Kappa score.
- Replicate analysis in either MLInterface, GLMnet, CMA, or Caret (take home exercise, optional).

Introduction
============

In this Seminar we go over packages and codes for performing supervised learning and evaluation in R. In supervised learning, one is given a training data set with known response and covariates and our goal is to predict the response in a new test set for which we only know the covariates. A supervised learning process can be decomposed into the following steps:

*Step 1*: Data preprocessing. Before getting started with feature selection and training, we must make sure we have adjusted the data to account for any batch effects or biases from outliers.

*Step 2*: Select Features. Before training a model, in many applications, it is usually important to perform a pre-filtering step in which one retains only the most informative features (e.g., genes) as candidate "biomarkers". The amount of features retained to select and train a model is up to the analyst and the methods used in the next steps. For example, some methods may be unfeasible or too slow to run with a large number of features.

*Step 3*: Select and train a classifier. Once the set of candidate markers have been selected, the next step is to select and train a model to predict the labels of a test data. We will also tune the parameters for each classifier using cross-validation.

*Step 4*: Test. Finally, a model is chosen and used to predict labels of a test data.

R packages
----------

There are many packages for performing supervised learning in R, each of which may implement one or more algorithms. There have also been at least two major efforts to unify these libraries under a common framework to make them easier to use: `MLInterfaces` and `CMA`. Although these may be useful and save you a lot of time in your analysis, it is important that you understand what these packages are doing and what they are *not* doing. Thus, I will not use these packages in this Seminar but I encourage you to reproduce the analysis using at least one of them! (I recommend `CMA`).

Install the following packages from Bioconductor: `CMA` and `GEOquery`, and from CRAN: `ROCR`, `car`, `e1071` (for SVM), and `glmnet` along with ther dependencies.

``` r
source('http://bioconductor.org/biocLite.R')
biocLite('GEOquery')
biocLite('CMA')

install.packages('ROCR')

install.packages(c('e1071','glmnet','caret','mlbench','gbm'))
```

``` r
library(MASS)
library(reshape)
library(car)
library(limma)
library(e1071)
library(glmnet)
library(ROCR)
library(CMA)
library(lattice)
library(class)
library(RCurl)
options('download.file.method'='curl')
library(GEOquery)
```

Data Set
--------

This seminar is based on a dataset that comes from a paper by [Smeets et al. 2010](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23177), who studied Affymetrix expression profiles from primary breast tumors. Smeets group was interested in whether tumors which had spread to lymph nodes (LN positive, generally a bad sign) have different gene expression profiles than LN negative tumors. If so, a gene expression signature can be use to predict tumor class.

Their data set contains 24236 genes on 116 samples. The status of the lymph node is known for each sample, with 59 LN positive and 57 LN negative. Samples were divided into two parts: 96 samples (48 LN positive and 48 LN negative) were used as a "training" set and 20 samples (11 LN positive and 9 LN negative) were used as a "test" set. There is also a quantitative measure, "LnRatio", the fraction of affected lymph nodes, presumably reflecting "how bad" the LnStatus is. Thus, we can use this dataset to illustrate classification and regularization methods! This seminar will focus on the first task, i.e., classification. In the past, Paul selected this dataset to illustrate the challenges of supervised learning tasks!

In the paper, the authors trained a support vector machine classifier to distinguish between LN positive and LN negative tumors (i.e., classification), and evaluated the results using ROC curves. After some optimization, they got an area under the ROC curve (AUC) of 0.66 on the training set and 0.65 on the test data. This is better than chance, but still not very convincing results about the relevance of the derived molecular signature (random chance would give an AUC of 0.5; perfect classification would give an AUC of 1.0).

### Data Preparation

First, let's retrieve our datasets from GEO with `getGEO` from `GEOquery` package. Warning: this may take several minutes! So to avoid re-downloading in the future, save the data once you get it into a good shape.

``` r
if(file.exists("class_LNstatus.Rdata")) { # if previously downloaded
  load("class_LNstatus.Rdata")
} else { # if downloading for the first time
  # takes a several mins!; returns a list
  datgeo <- getGEO('GSE23177', GSEMatrix = TRUE) 
  dat <- datgeo[[1]]   #Note that dat is an ExpressionSets
  
  str(pData(dat), max.level = 0)
  
  # extract only those variables of interest 
  pData(dat) <-
    subset(pData(dat),
           select = c("characteristics_ch1.2",
                      "characteristics_ch1.3","characteristics_ch1"))
  names(pData(dat))<-c("LnStatus", "LnRatio", "Set")

  #Note: LNRatio will not be used in this Seminar. However, you can use it to try some of the regularization techniques learned in class
  
  # split the ExpressionSet into training and test sets. 
  train.es <- dat[, dat$Set == "patient type: training set"]
  test.es <- dat[ , dat$Set != "patient type: training set"]

  #Re-label factor
  pData(train.es)$LnStatus <-
      recode(pData(train.es)$LnStatus, "levels(pData(train.es)$LnStatus)[1]='neg'; else='pos'", levels = c('neg', 'pos'))

  pData(test.es)$LnStatus <-
      recode(pData(test.es)$LnStatus, "levels(pData(test.es)$LnStatus)[1]='neg'; else='pos'",
             levels = c('neg', 'pos'))

  # create data matrices with expression values (probesets in rows). Some of the functions we will use do not take ExpressionSets as objects
  trainDat <- exprs(train.es)
  testDat <- exprs(test.es)

  # Redefine the quantitative variable LnRatio to make it a numeric variable.
  ntrain <- dim(pData(train.es))[1]
  ntest <- dim(pData(test.es))[1]
  
  pData(train.es)$LnRatio <- as.numeric(unlist(strsplit(as.vector(unlist(pData(train.es)$LnRatio)), ":", fixed = TRUE))[(1:ntrain)*2])
  pData(test.es)$LnRatio <- as.numeric(unlist(strsplit(as.vector(unlist(pData(test.es)$LnRatio)), ":", fixed = TRUE))[(1:ntest)*2])

  # save the data to avoid future re-downloading
  save(dat,trainDat,testDat,train.es,test.es, file = "class_LNstatus.Rdata")
}
```

Now, we can do some exploratory analysis of the data before trying some classification methods.

``` r
# undestand your data for classification
table(pData(train.es)$LnStatus)
```

    ## 
    ## neg pos 
    ##  48  48

``` r
table(pData(test.es)$LnStatus)
```

    ## 
    ## neg pos 
    ##   9  11

``` r
# understand the continuous response
tapply(pData(train.es)$LnRatio,pData(train.es)$LnStatus,summary)
```

    ## $neg
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0       0       0       0       0       0 
    ## 
    ## $pos
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0400  0.0700  0.1100  0.1935  0.2275  0.9600

``` r
tapply(pData(test.es)$LnRatio,pData(test.es)$LnStatus,summary)
```

    ## $neg
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0       0       0       0       0       0 
    ## 
    ## $pos
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0500  0.1800  0.5000  0.4573  0.6100  0.9400

``` r
# look at the expression of 3 randomly picked genes in both training and test sets
set.seed(1234)
(getMe <- sample(1:nrow(train.es), size = 3)) ## [1]   2756 15082 14766
```

    ## [1]  2756 15082 14766

``` r
# training data
trDat <- trainDat[getMe, ]
str(trDat)
```

    ##  num [1:3, 1:96] 7.19 9.17 8.38 7.13 9.38 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:3] "201513_at" "223139_s_at" "222698_s_at"
    ##   ..$ : chr [1:96] "GSM570518" "GSM570519" "GSM570520" "GSM570521" ...

``` r
trDat <- data.frame(LnStatus=pData(train.es)$LnStatus, Set=rep("train",nrow(pData(train.es))),t(trDat))
str(trDat)
```

    ## 'data.frame':    96 obs. of  5 variables:
    ##  $ LnStatus    : Factor w/ 2 levels "neg","pos": 1 1 2 1 1 2 1 1 2 1 ...
    ##  $ Set         : Factor w/ 1 level "train": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ X201513_at  : num  7.19 7.13 7.39 6.86 6.96 ...
    ##  $ X223139_s_at: num  9.17 9.38 9.03 9.55 9.5 ...
    ##  $ X222698_s_at: num  8.38 8.24 7.23 7.87 8.45 ...

``` r
plotDat.train <- melt(trDat, id=c("LnStatus","Set"),variable_name="gene")
colnames(plotDat.train)[colnames(plotDat.train)=="value"]="gExp"

# test data
tDat <- testDat[getMe, ]
str(tDat)
```

    ##  num [1:3, 1:20] 6.05 9.15 7.55 6.87 8.95 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:3] "201513_at" "223139_s_at" "222698_s_at"
    ##   ..$ : chr [1:20] "GSM570498" "GSM570499" "GSM570500" "GSM570501" ...

``` r
tDat <- data.frame(LnStatus=pData(test.es)$LnStatus,Set=rep("test",nrow(pData(test.es))), t(tDat))
str(tDat)
```

    ## 'data.frame':    20 obs. of  5 variables:
    ##  $ LnStatus    : Factor w/ 2 levels "neg","pos": 1 1 1 1 1 1 1 1 1 2 ...
    ##  $ Set         : Factor w/ 1 level "test": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ X201513_at  : num  6.05 6.87 6.71 8 6.54 ...
    ##  $ X223139_s_at: num  9.15 8.95 9.09 9.81 9.2 ...
    ##  $ X222698_s_at: num  7.55 8.34 8.32 7.33 8.14 ...

``` r
plotDat.test <- melt(tDat, id=c("LnStatus","Set"),variable_name="gene")
colnames(plotDat.test)[colnames(plotDat.test)=="value"]="gExp"

plotDat <- rbind(plotDat.train,plotDat.test)

# plot 3 randomly picked genes in both training and test sets
stripplot(gExp ~ LnStatus | gene+Set, plotDat,
                       grid = TRUE,
                       group = LnStatus, auto.key = TRUE,
                       jitter.data = TRUE)
```

![](sm10_supervisedLearning_files/figure-markdown_github/unnamed-chunk-3-1.png)

Classification
==============

The prediction of a discrete response is usually refer to as *classification*. A response taking values over a finite set of labels is essentially the same thing as a factor. We will use the dataset from Smeets et al. to find the *best-trained* classifier and use it to predict the `LnStatus` of the 20 samples in the test set, i.e., classify those as "lymph node positive" or "negative".

Data preprocessing
------------------

We should check to ensure there are no missing values in our data.

``` r
sum(is.na(trainDat))
```

    ## [1] 0

``` r
sum(is.na(testDat))
```

    ## [1] 0

Here we see there are no missing values in our dataset, so we don't have to worry about that.

Other pre-processing operations that can be done are:
- centering, scaling, normalizing
- imputing missing data
- transforming individual features (like boolean measurements)

When you are working with expression data, you may also need to use some normalization methods to ensure your variables are comparable across all your samples. For example, when working with count data, it is advised to log transformed and quantile-normalize your data, so that your samples have similar distributions.

Feature and Model Selection
---------------------------

We will now identify the best set of features that we will use to train the model using a cross-validation. Thus, I will divide the training set into 6 folds (the authors used 10 folds). We also want the proportion of positive and negative examples in each split to be approximately the same as for the full data set (i.e., stratified 6-fold CV with 8 positive and 8 negative samples within each fold). For each round of cross-validation, we use one fold as the test data and the rest of the data as training to select features and train different classifier.

### Cross-validation

Although it makes sense to proceed as described above, many methods are available and many constants within methods need to be selected in these steps. Thus, a *cross-validation* is usually required to *evaluate* how well different trained models work and select the *best* model to proceed. Note that although you may only want to select among different choices available in Step 2, the cross-validation needs to start in Step 1. Why? The results of the cross-validation will be over-optimistic and biased if the samples in the test sets of the cross-validation (i.e., left-out folds) were used to *select* the most promising features in Step 1!! For example, if the performance of a complex model is (artificially) good, you may not penalize regression coefficients enough in Step 2, and may yield to a poor performance in Step 3.

In many studies, in the absence of a test set, cross-validation is used to estimate performance. In those cases, a *nested cross-validation* is required! The inner cross-validation will be used to select features and tune parameters, the outer cross-validation will be used to a selected test model.

In this seminar, we will work with a dataset that has both a training and a test set. Thus, we will not do a nested cross-validation. However, keep it in mind for your project or future work!

#### Making cross validation splits

This is not the only way to create splits of the training data to run a cross-validation. Note that if the samples can not be evenly divided into the nfolds you specified, then you need to complete the matrices below with NAs and call for entries different from NA at those folds.

``` r
nfold <- 6

tabTrain <- table(train.es$LnStatus)

indlist <- sapply(names(tabTrain), function(z) which(train.es$LnStatus == z), simplify = FALSE)

set.seed(1234)

#Each row contains 8 pos and 8 negative samples. 

fold.pos <- matrix(sample(indlist[["pos"]]),nrow=nfold)
fold.neg <- matrix(sample(indlist[["neg"]]),nrow=nfold)
```

*Note*: with `CMA` you can use the command `GenerateLearningsets` to split the training data into folds. However, it does not show you how the data was split. Thus, you either use CMA for all or you write your own script.

``` r
splits <- GenerateLearningsets(y = train.es$LnStatus, method="CV", fold=6, strat= TRUE)
```

### Loop for feature selection and modeling

To illustrate how to select a model, I will use the top-50 genes selected by `limma` (within each fold). Note that this number is very arbitrary and other options may make more sense like using a p-value threshold or testing different options with this CV. For this example, I'm using only the top-50 genes as methods like LDA and Logit can not be run on more features than samples. However, other methods like kNN or SVM will do well with more features.

In this example, I will compare 7 different models: kNN for k={1,5,10,15}, LDA, Logit, SVM. Feel free to add other methods to the list!

``` r
#Define here the constants that you will not evaluate. For example, I will use the top-50 limma genes

ngenes <- 50
nmethod <- 7 #number of methods you plan to compare. 

#Define here an output objects to store results
pr.err <- matrix(-1, nfold,nmethod, dimnames=list(paste0("Fold",1:nfold),c("1NN","5NN","10NN", "15NN","LDA","Logit","SVM")))

for(i in 1:nfold){

  #Test Fold for the i-th step
  testdat.fold<-trainDat[,c(fold.pos[i,],fold.neg[i,])]
  #I will create a factor of classes for the test set of the i_th fold
  testclass.fold<-train.es$LnStatus[c(fold.pos[i,],fold.neg[i,])]
  
    
  #The rest of the samples are the training set for the i-th step
  traindat.fold<-trainDat[,-c(fold.pos[i,],fold.neg[i,])]
  trainclass.fold<-train.es$LnStatus[-c(fold.pos[i,],fold.neg[i,])]

  #Step 1: feature selection (do you remember limma?). 

  # Note that a different set of genes will be selected for each fold! you can then compare how consistent these sets were.

  limma.dat<-as.data.frame(traindat.fold)
  desMat <- model.matrix(~ trainclass.fold, limma.dat) #design matrix
  trainFit <- lmFit(limma.dat, desMat)
  eBtrainFit <- eBayes(trainFit)
  
  # top-50 limma genes
  top.fold <- topTable(eBtrainFit, coef = which(colnames(coef(trainFit)) != "(Intercept)"),
                       n = ngenes,sort.by="P")
  
  #Retain the top-50 limma genes from the train and test sets
  traindat.fold <- traindat.fold[rownames(top.fold),]
  testdat.fold <-  testdat.fold[rownames(top.fold),]

  
  #STEP 2: select a classifier
  #Set a counter for the method tested
  l <- 0

  #kNN classifiers
  for(kk in c(1,5,10,15)) {
    #every time you get inside this loop, the l counter gets redefined (i.e., 1, 2, etc for         method 1, method 2, etc)
    l <- l+1

    #knn needs samples in rows
    yhat.knn <- knn(train=t(traindat.fold), test=t(testdat.fold), cl=trainclass.fold,
                    k = kk)
    #Store the prediction error for each kk within this fold
    pr.err[i,l]<- mean(testclass.fold != yhat.knn)
                          } #end of kNN loop

  #LDA method. Note that you can change the prior parameter to reflect a different proportion of case and control samples. The default is to use the class proportions from the training set.
  
  m.lda <- lda(x=t(traindat.fold), group=trainclass.fold, prior=c(.5, .5))
  yhat.lda <- predict(m.lda, newdata=t(testdat.fold))$class
  pr.err[i,"LDA"] <-mean(testclass.fold != yhat.lda)
   
  #Logit
  glm.dat <- data.frame(t(traindat.fold), group=trainclass.fold)
  
  # 50 factors still will cause optimization warnings  
  # Try without warning suppression to see 
  # To further reduce parameters, regularized regression can be used 
  # To use regularized regression uncomment lines followed by "uncomment for regularized regression" 
  suppressWarnings( m.log <- glm(group ~ ., data=glm.dat,family=binomial) ) 
  
  # uncomment for regularized regression 
  # m.log <- glmnet( t(traindat.fold) , trainclass.fold ,family="binomial") 

  pr.log <- predict(m.log,newdata=data.frame(t(testdat.fold)),type="response")
  
  # uncomment for regularized regression 
  # pr.log <- predict(m.log,newdata=data.frame(t(testdat.fold)),type="response",newx=t(testdat.fold)) 
  
  pr.cl <- rep(0,length(testclass.fold))
  pr.cl[pr.log > 1/2] <- "pos"
  pr.cl[pr.log <= 1/2] <- "neg"

  pr.cl <- factor(pr.cl)
  pr.err[i,"Logit"] <- mean( pr.cl != testclass.fold )

  #SVM
  m.svm <- svm(x=t(traindat.fold), y=trainclass.fold, cost=1, type="C-classification", 
               kernel="linear")
  pr.svm <- predict(m.svm,newdata=t(testdat.fold)) 
   
  pr.err[i,"SVM"] <- mean( pr.svm != testclass.fold )
  } #end of CV loop
```

### Error Rates

Now you can get the average prediction error for all methods. Note that the prediction errors are high! not too much hope for the real test run!

``` r
cv.err <- colMeans(pr.err)

# mean - 1 sd (sd of the 6 error rates)
ls <- cv.err - apply(pr.err, 2, sd)

# mean + 1 sd (sd of the 6 error rates)
us <- cv.err + apply(pr.err, 2, sd)

# plot the results
plot(1:nmethod, cv.err, ylim=c(0, 1), xlim=c(1, (nmethod+.5)),type='n', 
axes=FALSE, xlab='Classifier', ylab='Error rate',main="6-fold CV Error")

for(j in 1:ncol(pr.err)) 
   points(jitter(rep(j, 6), factor=2), jitter(pr.err[,j]), cex=0.8, pch='X', col='gray')

for(i in 1:nmethod)
   lines(c(i, i), c(ls[i], us[i]), lwd=2, col='gray')
points(1:nmethod, ls, pch=19, col='red')
points(1:nmethod, us, pch=19, col='green')
points(1:nmethod, cv.err, pch=19, cex=1.5, col='black')
axis(2, ylab='Error rate')
axis(1, 1:nmethod, colnames(pr.err))

box()
```

![](sm10_supervisedLearning_files/figure-markdown_github/unnamed-chunk-8-1.png)

### Results of the CV

According to these results, LDA and 10NN may be the better classifier to try in the test data. However, remember that this CV results depend on the first split of the data we did. Thus, we need to repeat this CV

**Exercise 1**: perform 100 runs of this CV before selecting a model to test! Add at least one rule to select the list of models, e.g., use genes with a p-val threshold &lt; cutoff.

**Exercise 2**: Use AUC as a criteria to select a model based on the training data! Tip: extract the predicted probabilities from each method and use the roc function in ROCR.

Testing the selected model
--------------------------

Now that we decided on which method we are going to use to classify samples in the test set, we need to train the model using the *FULL* training set and then classify samples of the test set. I will use the 10NN model.

``` r
yhat.knn <- knn(train=t(trainDat), test=t(testDat), cl=train.es$LnStatus,
                     k = 10)
#Store the prediction error for each kk within this fold
pr.errTest<- mean(test.es$LnStatus != yhat.knn)
pr.errTest
```

    ## [1] 0.45

What does the prediction error mean?
In this instance, we have evaluated how often the prediction matched the actual lymph node status, against the total number of cases. This is the **accuracy** metric.

Not good! In real practice, you should not keep trying until we get a good result! In fact, you must use cross-validation on the training dataset to evaluate different parameters and classifiers, and only evaluate the generalizability of the *best* model on the test set.
However, in this seminar, I encourage you to try different options **as an exercise** and to see how much the results can change.

CMA
---

Many steps of the CV defined above can be easily done with CMA. For example, Step 1 in the loop above can also be done using 'CMA' with the function 'GeneSelection', which selects the most informative features (e.g., gene) to build a classifier within each of the splits generated by 'GenerateLearningsets'. Some learning algorithms do better if you only give them "useful" features.

``` r
featureScores<-GeneSelection(X=t(trainDat), y=train.es$LnStatus, learningsets=splits, method="limma")
```

    ## GeneSelection: iteration 1 
    ## GeneSelection: iteration 2 
    ## GeneSelection: iteration 3 
    ## GeneSelection: iteration 4 
    ## GeneSelection: iteration 5 
    ## GeneSelection: iteration 6

``` r
#Compare list of selected genes using:
toplist(featureScores)
```

    ## top  10  genes for iteration  1 
    ##  
    ##    index importance
    ## 1   9265   25.37579
    ## 2   6936   24.85400
    ## 3  21592   24.55555
    ## 4   1702   23.99424
    ## 5  21919   23.98643
    ## 6  19932   21.29348
    ## 7   6938   20.60153
    ## 8  22524   20.18380
    ## 9  17847   18.77700
    ## 10  6937   18.74738

``` r
#We can aggregate the results across the 6 splits.

seliter<-numeric()
for(i in 1:nfold) seliter<-c(seliter, toplist(featureScores, iter=i, top = 10, show=FALSE)$index)
(sort(table(seliter), dec=T)) # summarize
```

    ## seliter
    ##  1702  9265  6936 21919   808  6938 18958 19932 20571 21592 23206 23567 
    ##     6     5     3     3     2     2     2     2     2     2     2     2 
    ##    18   377   767  2690  3386  6937  7182  7183  8447 10254 10292 13581 
    ##     1     1     1     1     1     1     1     1     1     1     1     1 
    ## 13620 13802 15997 16094 16165 17847 18668 19152 19265 19402 19526 19577 
    ##     1     1     1     1     1     1     1     1     1     1     1     1 
    ## 21533 22524 22943 
    ##     1     1     1

``` r
# Choose the 20 probes which are chosen most commonly in the 6 splits
bestprobes<-as.numeric(names(sort(table(seliter), dec=T)))[1:20]

# examine the annotations. I just selected a few columns from the fData of the eSet.
(fData(dat)[bestprobes,c("Gene Symbol", "Gene Title", "ENTREZ_GENE_ID", "Representative Public ID" )])
```

    ##                                         Gene Symbol
    ## 1569472_s_at                        TTC3 /// TTC3P1
    ## 212384_at    ATP6V1G2-DDX39B /// DDX39B /// SNORD84
    ## 208661_s_at                         TTC3 /// TTC3P1
    ## 237746_at                                          
    ## 1556088_at                                    RPAIN
    ## 208663_s_at                         TTC3 /// TTC3P1
    ## 228510_at                                     ATAT1
    ## 230609_at                                    CLINT1
    ## 232740_at                                MCM3AP-AS1
    ## 236196_at                                    ZNF326
    ## 242562_at                                          
    ## 243751_at                                      CHD2
    ## 1552283_s_at                   ZDHHC11 /// ZDHHC11B
    ## 1554182_at                        TRIM73 /// TRIM74
    ## 1555920_at                                     CBX3
    ## 201440_at                                     DDX23
    ## 202182_at                                     KAT2A
    ## 208662_s_at                         TTC3 /// TTC3P1
    ## 208920_at                                       SRI
    ## 208921_s_at                                     SRI
    ##                                                                                                                                  Gene Title
    ## 1569472_s_at                                           tetratricopeptide repeat domain 3 /// tetratricopeptide repeat domain 3 pseudogene 1
    ## 212384_at    ATP6V1G2-DDX39B readthrough (NMD candidate) /// DEAD (Asp-Glu-Ala-Asp) box polypeptide 39B /// small nucleolar RNA, C/D box 84
    ## 208661_s_at                                            tetratricopeptide repeat domain 3 /// tetratricopeptide repeat domain 3 pseudogene 1
    ## 237746_at                                                                                                                                  
    ## 1556088_at                                                                                                          RPA interacting protein
    ## 208663_s_at                                            tetratricopeptide repeat domain 3 /// tetratricopeptide repeat domain 3 pseudogene 1
    ## 228510_at                                                                                                 alpha tubulin acetyltransferase 1
    ## 230609_at                                                                                                             clathrin interactor 1
    ## 232740_at                                                                                                            MCM3AP antisense RNA 1
    ## 236196_at                                                                                                           zinc finger protein 326
    ## 242562_at                                                                                                                                  
    ## 243751_at                                                                                       chromodomain helicase DNA binding protein 2
    ## 1552283_s_at                                                 zinc finger, DHHC-type containing 11 /// zinc finger, DHHC-type containing 11B
    ## 1554182_at                                                                tripartite motif containing 73 /// tripartite motif containing 74
    ## 1555920_at                                                                                                              chromobox homolog 3
    ## 201440_at                                                                                         DEAD (Asp-Glu-Ala-Asp) box polypeptide 23
    ## 202182_at                                                                                                    K(lysine) acetyltransferase 2A
    ## 208662_s_at                                            tetratricopeptide repeat domain 3 /// tetratricopeptide repeat domain 3 pseudogene 1
    ## 208920_at                                                                                                                            sorcin
    ## 208921_s_at                                                                                                                          sorcin
    ##                             ENTREZ_GENE_ID Representative Public ID
    ## 1569472_s_at               7267 /// 286495                 BC026260
    ## 212384_at    7919 /// 692199 /// 100532737                 AI282485
    ## 208661_s_at                7267 /// 286495                 AW510696
    ## 237746_at                                                  AI168187
    ## 1556088_at                           84268                 AK098491
    ## 208663_s_at                7267 /// 286495                 AI652848
    ## 228510_at                            79969                 AL566825
    ## 230609_at                             9685                 BF510429
    ## 232740_at                           114044                 BC002458
    ## 236196_at                           284695                 BF939032
    ## 242562_at                                                  AW772288
    ## 243751_at                             1106                 AA709148
    ## 1552283_s_at              79844 /// 653082                NM_024786
    ## 1554182_at               375593 /// 378108                 BC033871
    ## 1555920_at                           11335                 BU683892
    ## 201440_at                             9416                NM_004818
    ## 202182_at                             2648                NM_021078
    ## 208662_s_at                7267 /// 286495                 AI885338
    ## 208920_at                             6717                 AV752215
    ## 208921_s_at                           6717                   L12387

This looks promising since I get TTC3 and at least a couple of other genes that show up on Table 3 of the paper.

Similarly, you can use CMA to train and test a classifier within each CV fold (learningsets). However, there are things you can not do within CMA or that CMA is not doing right. For example, CMA can not do a full nested cross-validation. Additionally, it is not trivial to train the selected in the full dataset and then test it in the test set. CMA is more designed for CV. Thus, it is good to know how to do this things by hand as well.

Paul solved this problem in the following way: he made a `learningsets` object that has just one "split" defined by the samples in the training set.

``` r
m<-matrix(which(dat$Set=="patient type: training set"), 1)

full.learningset<-new("learningsets", learnmatrix=m, method="my own", ntrain=96, iter=1)

fullFeatureScores<-GeneSelection(X=t(exprs(dat)), learningsets= full.learningset, y=dat$LnStatus, method="t.test")
```

    ## GeneSelection: iteration 1

``` r
testclassif<-classification(X=t(exprs(dat)), y=dat$LnStatus, learningsets= full.learningset, genesel=fullFeatureScores, nbgene = 100, classifier =pknnCMA, k=5)
```

    ## iteration 1

``` r
#Evaluation:
tres<-testclassif[[1]]
ftable(tres)
```

    ## number of missclassifications:  11 
    ## missclassification rate:  0.55 
    ## sensitivity: 0.545 
    ## specificity: 0.333 
    ##     predicted
    ## true 0 1
    ##    0 3 6
    ##    1 5 6

``` r
roc(tres)
```

![](sm10_supervisedLearning_files/figure-markdown_github/unnamed-chunk-11-1.png)

Note: his optimized classifier did terribly as well.

Multiclass learning
-------------------

You won't always have a binary learning problem, where you are classifying a sample into 1 of 2 classes. Sometimes we might want to train a classifier a classifier with more than two classes.
Here we will use the *caret* and *mlbench* packages for classification on a multi-class problem.

``` r
library(caret)
```

    ## Loading required package: ggplot2

    ## 
    ## Attaching package: 'caret'

    ## The following objects are masked from 'package:CMA':
    ## 
    ##     best, rfe

``` r
library(mlbench)
```

We will be using the Soybean dataset, where our prediction is for the different problems associated with soybean crops. Our dataset has 683 samples, and 35 features being measured for each sample.
Our categories are 19.

``` r
cv_folds <- trainControl(method="cv", number=5)
data(Soybean)
unique(Soybean$Class)
```

    ##  [1] diaporthe-stem-canker       charcoal-rot               
    ##  [3] rhizoctonia-root-rot        phytophthora-rot           
    ##  [5] brown-stem-rot              powdery-mildew             
    ##  [7] downy-mildew                brown-spot                 
    ##  [9] bacterial-blight            bacterial-pustule          
    ## [11] purple-seed-stain           anthracnose                
    ## [13] phyllosticta-leaf-spot      alternarialeaf-spot        
    ## [15] frog-eye-leaf-spot          diaporthe-pod-&-stem-blight
    ## [17] cyst-nematode               2-4-d-injury               
    ## [19] herbicide-injury           
    ## 19 Levels: 2-4-d-injury alternarialeaf-spot ... rhizoctonia-root-rot

``` r
summary(Soybean)
```

    ##                  Class          date     plant.stand  precip      temp    
    ##  brown-spot         : 92   5      :149   0   :354    0   : 74   0   : 80  
    ##  alternarialeaf-spot: 91   4      :131   1   :293    1   :112   1   :374  
    ##  frog-eye-leaf-spot : 91   3      :118   NA's: 36    2   :459   2   :199  
    ##  phytophthora-rot   : 88   2      : 93               NA's: 38   NA's: 30  
    ##  anthracnose        : 44   6      : 90                                    
    ##  brown-stem-rot     : 44   (Other):101                                    
    ##  (Other)            :233   NA's   :  1                                    
    ##    hail     crop.hist  area.dam    sever     seed.tmt     germ    
    ##  0   :435   0   : 65   0   :123   0   :195   0   :305   0   :165  
    ##  1   :127   1   :165   1   :227   1   :322   1   :222   1   :213  
    ##  NA's:121   2   :219   2   :145   2   : 45   2   : 35   2   :193  
    ##             3   :218   3   :187   NA's:121   NA's:121   NA's:112  
    ##             NA's: 16   NA's:  1                                   
    ##                                                                   
    ##                                                                   
    ##  plant.growth leaves  leaf.halo  leaf.marg  leaf.size  leaf.shread
    ##  0   :441     0: 77   0   :221   0   :357   0   : 51   0   :487   
    ##  1   :226     1:606   1   : 36   1   : 21   1   :327   1   : 96   
    ##  NA's: 16             2   :342   2   :221   2   :221   NA's:100   
    ##                       NA's: 84   NA's: 84   NA's: 84              
    ##                                                                   
    ##                                                                   
    ##                                                                   
    ##  leaf.malf  leaf.mild    stem     lodging    stem.cankers canker.lesion
    ##  0   :554   0   :535   0   :296   0   :520   0   :379     0   :320     
    ##  1   : 45   1   : 20   1   :371   1   : 42   1   : 39     1   : 83     
    ##  NA's: 84   2   : 20   NA's: 16   NA's:121   2   : 36     2   :177     
    ##             NA's:108                         3   :191     3   : 65     
    ##                                              NA's: 38     NA's: 38     
    ##                                                                        
    ##                                                                        
    ##  fruiting.bodies ext.decay  mycelium   int.discolor sclerotia  fruit.pods
    ##  0   :473        0   :497   0   :639   0   :581     0   :625   0   :407  
    ##  1   :104        1   :135   1   :  6   1   : 44     1   : 20   1   :130  
    ##  NA's:106        2   : 13   NA's: 38   2   : 20     NA's: 38   2   : 14  
    ##                  NA's: 38              NA's: 38                3   : 48  
    ##                                                                NA's: 84  
    ##                                                                          
    ##                                                                          
    ##  fruit.spots   seed     mold.growth seed.discolor seed.size  shriveling
    ##  0   :345    0   :476   0   :524    0   :513      0   :532   0   :539  
    ##  1   : 75    1   :115   1   : 67    1   : 64      1   : 59   1   : 38  
    ##  2   : 57    NA's: 92   NA's: 92    NA's:106      NA's: 92   NA's:106  
    ##  4   :100                                                              
    ##  NA's:106                                                              
    ##                                                                        
    ##                                                                        
    ##   roots    
    ##  0   :551  
    ##  1   : 86  
    ##  2   : 15  
    ##  NA's: 31  
    ##            
    ##            
    ## 

Let us pre-process our data.

``` r
#Remove rows (samples) with missing values  
soybean_x = Soybean[rowSums(is.na(Soybean)) == 0,]
#Then remove columns (features) with missing values  
soybean_x = soybean_x[,colSums(is.na(soybean_x)) == 0]
dim(soybean_x)
```

    ## [1] 562  36

We are left with 562 samples and 35 attributes. The first column, `Class`, describes the categories. We will refactor this column since we have removed certain columns (and possibly some of the 19 classes)

``` r
soybean_x$Class = (as.factor(as.character(soybean_x$Class)))
unique(soybean_x$Class)
```

    ##  [1] diaporthe-stem-canker  charcoal-rot           rhizoctonia-root-rot  
    ##  [4] phytophthora-rot       brown-stem-rot         powdery-mildew        
    ##  [7] downy-mildew           brown-spot             bacterial-blight      
    ## [10] bacterial-pustule      purple-seed-stain      anthracnose           
    ## [13] phyllosticta-leaf-spot alternarialeaf-spot    frog-eye-leaf-spot    
    ## 15 Levels: alternarialeaf-spot anthracnose ... rhizoctonia-root-rot

Now we have 15 classes!

In this instance, we don't have an external test set for evaluation, so we will be assessing the performance on a held-out test set. First, we create the held-out test set.
*Note* that we are holding out this test set from our training data, and then performing data-splitting for validation within our training subset. In practise, like we did earlier, you would want to loop this multiple times with holding out a test set and training on the remainder dataset (using cross validation or bootstrapping to estimate your model accuracy on the test set).

``` r
trainIndex <- createDataPartition(soybean_x$Class, p = .8,  list = FALSE,  times = 1)
soyTrain <- soybean_x[ trainIndex,]
soyTest  <- soybean_x[-trainIndex,]
```

**Cross validation results**
We set up our cross-validation folds. Note we can also choose the option 'cv' or 'LOOCV' instead of 'repeatedcv'.
With 'repeatedcv', we don't have to manually set up the multiple loops we did when we were using CMA.

``` r
#Prepare resamling method for Cross-Validation folds  
set.seed(7)
cv_control = trainControl(method="repeatedcv", number=5, repeats=10, classProbs=FALSE) #summaryFunction=mnLogLoss)
modelSvm_cv <- train(Class~., data=soyTrain, method="svmRadial", metric="accuracy", trControl=cv_control, verbose=FALSE)
```

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in train.default(x, y, weights = w, ...): The metric "accuracy" was
    ## not in the result set. Accuracy will be used instead.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

``` r
# display results
print(modelSvm_cv)
```

    ## Support Vector Machines with Radial Basis Function Kernel 
    ## 
    ## 452 samples
    ##  35 predictors
    ##  15 classes: 'alternarialeaf-spot', 'anthracnose', 'bacterial-blight', 'bacterial-pustule', 'brown-spot', 'brown-stem-rot', 'charcoal-rot', 'diaporthe-stem-canker', 'downy-mildew', 'frog-eye-leaf-spot', 'phyllosticta-leaf-spot', 'phytophthora-rot', 'powdery-mildew', 'purple-seed-stain', 'rhizoctonia-root-rot' 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (5 fold, repeated 10 times) 
    ## Summary of sample sizes: 359, 362, 363, 359, 365, 361, ... 
    ## Resampling results across tuning parameters:
    ## 
    ##   C     Accuracy   Kappa    
    ##   0.25  0.7188972  0.6787205
    ##   0.50  0.8878138  0.8740985
    ##   1.00  0.9222808  0.9130760
    ## 
    ## Tuning parameter 'sigma' was held constant at a value of 0.07
    ## Accuracy was used to select the optimal model using the largest value.
    ## The final values used for the model were sigma = 0.07 and C = 1.

``` r
y_pred_cv = predict(modelSvm_cv, soyTest[,2:36])
```

When fitting this model with SVM, you might get the following request: '1 package is needed for this model and is not installed. (kernlab). Would you like to try to install it now?'. Enter *yes*.

We can also use bootstrap re-sampling instead of cross-validation folds. This means we take random samples from the dataset (with re-selection), against which to evaluate our model. This gives us an idea about hte variance of the model itself.

Bootstrapping is different from cross validation in that the latter splits the entire dataset into folds *without re-selection*.It is a robust method to estimate the accuracy of our model.

**Bootstrapped results**

``` r
#Prepare resamling method for Bootstrapping
set.seed(7)
cv_control = trainControl(method="boot", number=100, classProbs=FALSE) 
modelSvm_boot <- train(Class~., data=soyTrain, method="svmRadial", metric="accuracy", trControl=cv_control, verbose=FALSE)
```

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in train.default(x, y, weights = w, ...): The metric "accuracy" was
    ## not in the result set. Accuracy will be used instead.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

``` r
# display results
print(modelSvm_boot)
```

    ## Support Vector Machines with Radial Basis Function Kernel 
    ## 
    ## 452 samples
    ##  35 predictors
    ##  15 classes: 'alternarialeaf-spot', 'anthracnose', 'bacterial-blight', 'bacterial-pustule', 'brown-spot', 'brown-stem-rot', 'charcoal-rot', 'diaporthe-stem-canker', 'downy-mildew', 'frog-eye-leaf-spot', 'phyllosticta-leaf-spot', 'phytophthora-rot', 'powdery-mildew', 'purple-seed-stain', 'rhizoctonia-root-rot' 
    ## 
    ## No pre-processing
    ## Resampling: Bootstrapped (100 reps) 
    ## Summary of sample sizes: 452, 452, 452, 452, 452, 452, ... 
    ## Resampling results across tuning parameters:
    ## 
    ##   C     Accuracy   Kappa    
    ##   0.25  0.7498099  0.7156718
    ##   0.50  0.8844782  0.8700662
    ##   1.00  0.9134670  0.9029037
    ## 
    ## Tuning parameter 'sigma' was held constant at a value of 0.07
    ## Accuracy was used to select the optimal model using the largest value.
    ## The final values used for the model were sigma = 0.07 and C = 1.

``` r
y_pred_boot = predict(modelSvm_boot, soyTest[,2:36])
```

*Fitting other models*
We can consider other algorithmic approaches for classifiers too, as follows:

``` r
# train the GLM model
set.seed(7)
modelNN <- train(Class~., data=soyTrain, method="nnet", metric="accuracy", trControl=cv_control, verbose=FALSE)
```

    ## Warning in train.default(x, y, weights = w, ...): The metric "accuracy" was
    ## not in the result set. Accuracy will be used instead.

    ## # weights:  95
    ## initial  value 1327.679959 
    ## iter  10 value 982.170176
    ## iter  20 value 835.519629
    ## iter  30 value 768.761550
    ## iter  40 value 718.525648
    ## iter  50 value 658.795486
    ## iter  60 value 618.464227
    ## iter  70 value 566.702345
    ## iter  80 value 520.418568
    ## iter  90 value 491.830936
    ## iter 100 value 473.711406
    ## final  value 473.711406 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1334.809255 
    ## iter  10 value 663.320864
    ## iter  20 value 522.444709
    ## iter  30 value 420.854689
    ## iter  40 value 334.757615
    ## iter  50 value 283.638942
    ## iter  60 value 249.792281
    ## iter  70 value 220.469904
    ## iter  80 value 182.764805
    ## iter  90 value 154.097431
    ## iter 100 value 134.623459
    ## final  value 134.623459 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1316.527953 
    ## iter  10 value 418.881758
    ## iter  20 value 177.733835
    ## iter  30 value 105.115709
    ## iter  40 value 62.877600
    ## iter  50 value 46.118487
    ## iter  60 value 39.631082
    ## iter  70 value 33.230157
    ## iter  80 value 31.069831
    ## iter  90 value 29.970791
    ## iter 100 value 29.160790
    ## final  value 29.160790 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1300.318660 
    ## iter  10 value 1067.844566
    ## iter  20 value 926.120180
    ## iter  30 value 863.917746
    ## iter  40 value 839.546972
    ## iter  50 value 825.072397
    ## iter  60 value 805.025842
    ## iter  70 value 791.896759
    ## iter  80 value 783.632939
    ## iter  90 value 778.002468
    ## iter 100 value 772.079332
    ## final  value 772.079332 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1255.782989 
    ## iter  10 value 715.691028
    ## iter  20 value 515.256191
    ## iter  30 value 400.705835
    ## iter  40 value 352.479210
    ## iter  50 value 328.681994
    ## iter  60 value 318.976919
    ## iter  70 value 312.053669
    ## iter  80 value 307.972090
    ## iter  90 value 307.091034
    ## iter 100 value 306.838678
    ## final  value 306.838678 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1312.576894 
    ## iter  10 value 616.825793
    ## iter  20 value 398.321508
    ## iter  30 value 289.973161
    ## iter  40 value 212.034690
    ## iter  50 value 178.409651
    ## iter  60 value 169.698716
    ## iter  70 value 166.062355
    ## iter  80 value 164.842670
    ## iter  90 value 164.241196
    ## iter 100 value 163.860382
    ## final  value 163.860382 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1212.804742 
    ## iter  10 value 830.813091
    ## iter  20 value 764.372190
    ## iter  30 value 719.285409
    ## iter  40 value 688.112947
    ## iter  50 value 664.547004
    ## iter  60 value 626.278073
    ## iter  70 value 615.118777
    ## iter  80 value 613.921584
    ## iter  90 value 613.613106
    ## iter 100 value 612.934510
    ## final  value 612.934510 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1288.070337 
    ## iter  10 value 789.697799
    ## iter  20 value 396.692136
    ## iter  30 value 250.837186
    ## iter  40 value 204.255441
    ## iter  50 value 163.554489
    ## iter  60 value 136.873093
    ## iter  70 value 118.542178
    ## iter  80 value 108.416804
    ## iter  90 value 100.688495
    ## iter 100 value 96.049496
    ## final  value 96.049496 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1328.551719 
    ## iter  10 value 497.189185
    ## iter  20 value 174.185419
    ## iter  30 value 83.395431
    ## iter  40 value 48.270890
    ## iter  50 value 29.904182
    ## iter  60 value 24.297297
    ## iter  70 value 19.840660
    ## iter  80 value 19.118021
    ## iter  90 value 18.156610
    ## iter 100 value 16.966931
    ## final  value 16.966931 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1263.291487 
    ## iter  10 value 847.668026
    ## iter  20 value 786.831827
    ## iter  30 value 745.166592
    ## iter  40 value 711.986561
    ## iter  50 value 680.204903
    ## iter  60 value 654.092570
    ## iter  70 value 639.163048
    ## iter  80 value 630.494030
    ## iter  90 value 622.390337
    ## iter 100 value 618.352152
    ## final  value 618.352152 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1191.611780 
    ## iter  10 value 596.975465
    ## iter  20 value 433.327739
    ## iter  30 value 312.544838
    ## iter  40 value 260.225134
    ## iter  50 value 223.481073
    ## iter  60 value 190.526529
    ## iter  70 value 176.370058
    ## iter  80 value 157.794177
    ## iter  90 value 147.431180
    ## iter 100 value 142.520948
    ## final  value 142.520948 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1305.360318 
    ## iter  10 value 395.592260
    ## iter  20 value 101.223822
    ## iter  30 value 51.680663
    ## iter  40 value 32.786067
    ## iter  50 value 19.852612
    ## iter  60 value 13.154035
    ## iter  70 value 6.960582
    ## iter  80 value 0.924599
    ## iter  90 value 0.269172
    ## iter 100 value 0.087041
    ## final  value 0.087041 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1262.407231 
    ## iter  10 value 1032.478034
    ## iter  20 value 983.837643
    ## iter  30 value 926.137416
    ## iter  40 value 850.737629
    ## iter  50 value 824.613227
    ## iter  60 value 804.948050
    ## iter  70 value 794.460607
    ## iter  80 value 786.517587
    ## iter  90 value 783.196186
    ## iter 100 value 782.093245
    ## final  value 782.093245 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1329.911296 
    ## iter  10 value 1016.037124
    ## iter  20 value 919.390411
    ## iter  30 value 724.659622
    ## iter  40 value 631.711891
    ## iter  50 value 582.474778
    ## iter  60 value 451.558512
    ## iter  70 value 370.457094
    ## iter  80 value 340.974421
    ## iter  90 value 329.328601
    ## iter 100 value 316.042009
    ## final  value 316.042009 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1219.679657 
    ## iter  10 value 541.717514
    ## iter  20 value 348.517309
    ## iter  30 value 264.514317
    ## iter  40 value 210.405441
    ## iter  50 value 186.196268
    ## iter  60 value 177.781774
    ## iter  70 value 174.767747
    ## iter  80 value 165.618157
    ## iter  90 value 159.595340
    ## iter 100 value 157.493659
    ## final  value 157.493659 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1260.460538 
    ## iter  10 value 861.348392
    ## iter  20 value 777.520832
    ## iter  30 value 719.832821
    ## iter  40 value 689.365430
    ## iter  50 value 663.635345
    ## iter  60 value 653.310635
    ## iter  70 value 648.393241
    ## iter  80 value 646.038758
    ## iter  90 value 645.019345
    ## iter 100 value 644.734730
    ## final  value 644.734730 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1249.158623 
    ## iter  10 value 645.663209
    ## iter  20 value 452.327396
    ## iter  30 value 370.114179
    ## iter  40 value 315.721643
    ## iter  50 value 283.356213
    ## iter  60 value 260.852477
    ## iter  70 value 246.814034
    ## iter  80 value 229.972047
    ## iter  90 value 224.587279
    ## iter 100 value 220.840197
    ## final  value 220.840197 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1366.805777 
    ## iter  10 value 463.074935
    ## iter  20 value 132.275766
    ## iter  30 value 58.484326
    ## iter  40 value 28.604772
    ## iter  50 value 19.588712
    ## iter  60 value 14.054791
    ## iter  70 value 12.428132
    ## iter  80 value 11.112760
    ## iter  90 value 10.572892
    ## iter 100 value 7.440712
    ## final  value 7.440712 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1236.639234 
    ## iter  10 value 802.791370
    ## iter  20 value 762.789606
    ## iter  30 value 727.669749
    ## iter  40 value 691.369582
    ## iter  50 value 663.145678
    ## iter  60 value 645.470418
    ## iter  70 value 639.882739
    ## iter  80 value 637.199354
    ## iter  90 value 633.050013
    ## iter 100 value 631.940124
    ## final  value 631.940124 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1283.739095 
    ## iter  10 value 625.822873
    ## iter  20 value 405.622056
    ## iter  30 value 326.288676
    ## iter  40 value 290.540275
    ## iter  50 value 265.423635
    ## iter  60 value 243.957180
    ## iter  70 value 226.554470
    ## iter  80 value 219.118705
    ## iter  90 value 211.117922
    ## iter 100 value 207.409092
    ## final  value 207.409092 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1306.316783 
    ## iter  10 value 331.854165
    ## iter  20 value 96.869965
    ## iter  30 value 65.585765
    ## iter  40 value 47.238068
    ## iter  50 value 34.013989
    ## iter  60 value 25.986470
    ## iter  70 value 21.487820
    ## iter  80 value 18.024792
    ## iter  90 value 16.457645
    ## iter 100 value 16.387794
    ## final  value 16.387794 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1253.235990 
    ## iter  10 value 949.985023
    ## iter  20 value 827.136915
    ## iter  30 value 781.693759
    ## iter  40 value 763.677167
    ## iter  50 value 742.530542
    ## iter  60 value 724.343287
    ## iter  70 value 706.508732
    ## iter  80 value 702.575079
    ## iter  90 value 700.024135
    ## iter 100 value 699.511897
    ## final  value 699.511897 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1314.841247 
    ## iter  10 value 1020.422239
    ## iter  20 value 870.267863
    ## iter  30 value 655.194179
    ## iter  40 value 449.760192
    ## iter  50 value 388.197077
    ## iter  60 value 343.102264
    ## iter  70 value 322.802172
    ## iter  80 value 303.040564
    ## iter  90 value 290.529464
    ## iter 100 value 286.620258
    ## final  value 286.620258 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1231.162065 
    ## iter  10 value 560.609123
    ## iter  20 value 309.402006
    ## iter  30 value 200.001249
    ## iter  40 value 159.821838
    ## iter  50 value 148.037383
    ## iter  60 value 143.761504
    ## iter  70 value 142.283477
    ## iter  80 value 141.335979
    ## iter  90 value 140.173344
    ## iter 100 value 139.859314
    ## final  value 139.859314 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1234.303389 
    ## iter  10 value 1055.665354
    ## iter  20 value 1052.704202
    ## iter  30 value 1052.634006
    ## final  value 1052.633920 
    ## converged
    ## # weights:  255
    ## initial  value 1248.684566 
    ## iter  10 value 647.201030
    ## iter  20 value 355.485685
    ## iter  30 value 280.144450
    ## iter  40 value 245.578504
    ## iter  50 value 221.574496
    ## iter  60 value 206.468784
    ## iter  70 value 201.860565
    ## iter  80 value 200.945351
    ## iter  90 value 199.632950
    ## iter 100 value 195.453577
    ## final  value 195.453577 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1308.111449 
    ## iter  10 value 456.145713
    ## iter  20 value 94.045316
    ## iter  30 value 38.588658
    ## iter  40 value 11.221074
    ## iter  50 value 7.372763
    ## iter  60 value 7.199632
    ## iter  70 value 7.048880
    ## iter  80 value 6.864288
    ## iter  90 value 6.727944
    ## iter 100 value 6.596727
    ## final  value 6.596727 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1190.669929 
    ## iter  10 value 907.915265
    ## iter  20 value 798.137418
    ## iter  30 value 761.100753
    ## iter  40 value 729.544778
    ## iter  50 value 716.594639
    ## iter  60 value 707.416716
    ## iter  70 value 704.152707
    ## iter  80 value 702.580190
    ## iter  90 value 701.964906
    ## iter 100 value 701.632831
    ## final  value 701.632831 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1268.748618 
    ## iter  10 value 671.401463
    ## iter  20 value 533.271697
    ## iter  30 value 447.949396
    ## iter  40 value 388.816303
    ## iter  50 value 340.285114
    ## iter  60 value 299.906339
    ## iter  70 value 270.509673
    ## iter  80 value 240.711789
    ## iter  90 value 231.699898
    ## iter 100 value 229.701598
    ## final  value 229.701598 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1346.959343 
    ## iter  10 value 547.261903
    ## iter  20 value 229.588426
    ## iter  30 value 106.863052
    ## iter  40 value 71.730285
    ## iter  50 value 49.876037
    ## iter  60 value 34.625225
    ## iter  70 value 29.777272
    ## iter  80 value 25.119958
    ## iter  90 value 20.605841
    ## iter 100 value 16.575136
    ## final  value 16.575136 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1269.331072 
    ## iter  10 value 1032.052327
    ## iter  20 value 899.615768
    ## iter  30 value 814.603844
    ## iter  40 value 792.476973
    ## iter  50 value 763.384731
    ## iter  60 value 749.059481
    ## iter  70 value 744.153804
    ## iter  80 value 743.168995
    ## iter  90 value 739.013729
    ## iter 100 value 737.075786
    ## final  value 737.075786 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1376.308866 
    ## iter  10 value 824.940393
    ## iter  20 value 582.871633
    ## iter  30 value 475.020566
    ## iter  40 value 419.890207
    ## iter  50 value 377.197907
    ## iter  60 value 350.455174
    ## iter  70 value 317.985483
    ## iter  80 value 306.138271
    ## iter  90 value 301.050520
    ## iter 100 value 298.987817
    ## final  value 298.987817 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1217.519275 
    ## iter  10 value 526.326510
    ## iter  20 value 313.203656
    ## iter  30 value 205.205057
    ## iter  40 value 166.296259
    ## iter  50 value 159.897892
    ## iter  60 value 153.647500
    ## iter  70 value 151.042828
    ## iter  80 value 150.161776
    ## iter  90 value 150.014199
    ## iter 100 value 149.946876
    ## final  value 149.946876 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1335.859568 
    ## iter  10 value 839.493747
    ## iter  20 value 763.157700
    ## iter  30 value 725.351244
    ## iter  40 value 689.967888
    ## iter  50 value 670.435158
    ## iter  60 value 668.019971
    ## iter  70 value 665.216158
    ## iter  80 value 663.741682
    ## iter  90 value 662.969714
    ## iter 100 value 658.240187
    ## final  value 658.240187 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1279.905030 
    ## iter  10 value 642.945403
    ## iter  20 value 407.784634
    ## iter  30 value 317.673406
    ## iter  40 value 265.649427
    ## iter  50 value 235.854937
    ## iter  60 value 205.031947
    ## iter  70 value 195.934994
    ## iter  80 value 186.926609
    ## iter  90 value 182.266149
    ## iter 100 value 180.455042
    ## final  value 180.455042 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1259.803411 
    ## iter  10 value 424.437378
    ## iter  20 value 155.168291
    ## iter  30 value 76.170932
    ## iter  40 value 51.535656
    ## iter  50 value 38.400628
    ## iter  60 value 32.470886
    ## iter  70 value 25.345286
    ## iter  80 value 24.606918
    ## iter  90 value 21.385528
    ## iter 100 value 20.317741
    ## final  value 20.317741 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1236.500974 
    ## iter  10 value 986.644300
    ## iter  20 value 856.282589
    ## iter  30 value 783.721988
    ## iter  40 value 734.960356
    ## iter  50 value 697.480653
    ## iter  60 value 681.881728
    ## iter  70 value 673.745374
    ## iter  80 value 665.576610
    ## iter  90 value 663.661774
    ## iter 100 value 662.074336
    ## final  value 662.074336 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1350.492393 
    ## iter  10 value 908.312869
    ## iter  20 value 701.143446
    ## iter  30 value 591.568532
    ## iter  40 value 466.376051
    ## iter  50 value 399.256982
    ## iter  60 value 362.582394
    ## iter  70 value 336.636880
    ## iter  80 value 320.237388
    ## iter  90 value 305.971404
    ## iter 100 value 292.581170
    ## final  value 292.581170 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1333.959327 
    ## iter  10 value 424.362198
    ## iter  20 value 130.603804
    ## iter  30 value 63.115754
    ## iter  40 value 36.400232
    ## iter  50 value 20.446185
    ## iter  60 value 17.365135
    ## iter  70 value 17.192545
    ## iter  80 value 17.091824
    ## iter  90 value 16.306367
    ## iter 100 value 16.303680
    ## final  value 16.303680 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1181.326310 
    ## iter  10 value 931.629793
    ## iter  20 value 828.341109
    ## iter  30 value 794.214679
    ## iter  40 value 771.552979
    ## iter  50 value 760.326583
    ## iter  60 value 754.667100
    ## iter  70 value 752.776104
    ## iter  80 value 752.530398
    ## iter  90 value 752.520020
    ## final  value 752.519633 
    ## converged
    ## # weights:  255
    ## initial  value 1197.131570 
    ## iter  10 value 531.181325
    ## iter  20 value 408.664682
    ## iter  30 value 358.156658
    ## iter  40 value 329.604304
    ## iter  50 value 316.260416
    ## iter  60 value 310.533389
    ## iter  70 value 305.950948
    ## iter  80 value 300.212128
    ## iter  90 value 296.144577
    ## iter 100 value 294.376757
    ## final  value 294.376757 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1275.477806 
    ## iter  10 value 554.598995
    ## iter  20 value 295.999461
    ## iter  30 value 197.384797
    ## iter  40 value 167.009874
    ## iter  50 value 159.224415
    ## iter  60 value 157.912412
    ## iter  70 value 156.976411
    ## iter  80 value 156.223725
    ## iter  90 value 155.485686
    ## iter 100 value 154.268325
    ## final  value 154.268325 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1296.344084 
    ## iter  10 value 932.609052
    ## iter  20 value 878.679572
    ## iter  30 value 837.915172
    ## iter  40 value 818.865956
    ## iter  50 value 808.419673
    ## iter  60 value 790.342308
    ## iter  70 value 788.740930
    ## iter  80 value 788.068687
    ## iter  90 value 787.692784
    ## iter 100 value 787.606132
    ## final  value 787.606132 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1302.507870 
    ## iter  10 value 834.255004
    ## iter  20 value 570.124882
    ## iter  30 value 444.540714
    ## iter  40 value 379.012389
    ## iter  50 value 311.576735
    ## iter  60 value 240.257059
    ## iter  70 value 220.129820
    ## iter  80 value 205.512781
    ## iter  90 value 200.647169
    ## iter 100 value 192.625341
    ## final  value 192.625341 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1389.664334 
    ## iter  10 value 373.174324
    ## iter  20 value 146.221332
    ## iter  30 value 74.379195
    ## iter  40 value 22.429567
    ## iter  50 value 8.244566
    ## iter  60 value 3.146601
    ## iter  70 value 2.825518
    ## iter  80 value 2.508665
    ## iter  90 value 2.325092
    ## iter 100 value 2.147266
    ## final  value 2.147266 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1233.533223 
    ## iter  10 value 1000.060218
    ## iter  20 value 805.542781
    ## iter  30 value 752.496540
    ## iter  40 value 722.598933
    ## iter  50 value 703.886872
    ## iter  60 value 683.827733
    ## iter  70 value 670.939503
    ## iter  80 value 662.632350
    ## iter  90 value 657.838938
    ## iter 100 value 651.994408
    ## final  value 651.994408 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1228.929573 
    ## iter  10 value 589.427272
    ## iter  20 value 333.576619
    ## iter  30 value 267.941123
    ## iter  40 value 230.167285
    ## iter  50 value 201.413993
    ## iter  60 value 179.524225
    ## iter  70 value 168.969150
    ## iter  80 value 163.366116
    ## iter  90 value 157.303634
    ## iter 100 value 154.564093
    ## final  value 154.564093 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1349.296558 
    ## iter  10 value 534.578537
    ## iter  20 value 122.308357
    ## iter  30 value 67.548319
    ## iter  40 value 29.343072
    ## iter  50 value 22.280611
    ## iter  60 value 17.319088
    ## iter  70 value 16.608734
    ## iter  80 value 13.650062
    ## iter  90 value 13.286442
    ## iter 100 value 13.238839
    ## final  value 13.238839 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1323.319475 
    ## iter  10 value 1128.435057
    ## iter  20 value 954.298635
    ## iter  30 value 911.719853
    ## iter  40 value 858.878597
    ## iter  50 value 822.299187
    ## iter  60 value 807.702398
    ## iter  70 value 803.614093
    ## iter  80 value 802.962121
    ## iter  90 value 801.651245
    ## iter 100 value 801.082823
    ## final  value 801.082823 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1332.759953 
    ## iter  10 value 841.531616
    ## iter  20 value 575.863479
    ## iter  30 value 441.501336
    ## iter  40 value 381.835204
    ## iter  50 value 344.128647
    ## iter  60 value 328.439238
    ## iter  70 value 322.191601
    ## iter  80 value 319.014518
    ## iter  90 value 315.003701
    ## iter 100 value 312.909339
    ## final  value 312.909339 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1296.300900 
    ## iter  10 value 471.196523
    ## iter  20 value 243.105889
    ## iter  30 value 174.224992
    ## iter  40 value 161.729088
    ## iter  50 value 159.601206
    ## iter  60 value 157.954106
    ## iter  70 value 157.532126
    ## iter  80 value 155.867448
    ## iter  90 value 153.081609
    ## iter 100 value 151.979211
    ## final  value 151.979211 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1240.099356 
    ## iter  10 value 922.586782
    ## iter  20 value 802.299725
    ## iter  30 value 755.167884
    ## iter  40 value 716.892600
    ## iter  50 value 695.286177
    ## iter  60 value 664.369787
    ## iter  70 value 636.159832
    ## iter  80 value 619.381006
    ## iter  90 value 614.056181
    ## iter 100 value 611.150470
    ## final  value 611.150470 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1377.362843 
    ## iter  10 value 774.946702
    ## iter  20 value 608.640668
    ## iter  30 value 513.431617
    ## iter  40 value 419.726038
    ## iter  50 value 360.712507
    ## iter  60 value 323.793557
    ## iter  70 value 299.561025
    ## iter  80 value 286.250765
    ## iter  90 value 278.896480
    ## iter 100 value 275.551437
    ## final  value 275.551437 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1288.002368 
    ## iter  10 value 410.103952
    ## iter  20 value 121.034067
    ## iter  30 value 56.488019
    ## iter  40 value 45.878188
    ## iter  50 value 41.084232
    ## iter  60 value 34.609388
    ## iter  70 value 28.690482
    ## iter  80 value 25.957356
    ## iter  90 value 24.692993
    ## iter 100 value 24.319782
    ## final  value 24.319782 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1334.965169 
    ## iter  10 value 1097.969048
    ## iter  20 value 1088.661459
    ## iter  30 value 928.107706
    ## iter  40 value 826.702489
    ## iter  50 value 778.523161
    ## iter  60 value 754.238535
    ## iter  70 value 731.771170
    ## iter  80 value 718.662385
    ## iter  90 value 711.325419
    ## iter 100 value 703.951740
    ## final  value 703.951740 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1244.471961 
    ## iter  10 value 490.572382
    ## iter  20 value 298.379002
    ## iter  30 value 240.054224
    ## iter  40 value 208.661692
    ## iter  50 value 176.900239
    ## iter  60 value 153.169695
    ## iter  70 value 135.617219
    ## iter  80 value 123.540253
    ## iter  90 value 118.747981
    ## iter 100 value 116.485748
    ## final  value 116.485748 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1254.127173 
    ## iter  10 value 290.030948
    ## iter  20 value 80.163149
    ## iter  30 value 31.769237
    ## iter  40 value 19.951629
    ## iter  50 value 12.132559
    ## iter  60 value 5.418627
    ## iter  70 value 5.353429
    ## iter  80 value 5.339809
    ## iter  90 value 5.337468
    ## iter 100 value 5.337397
    ## final  value 5.337397 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1242.596441 
    ## iter  10 value 936.017212
    ## iter  20 value 838.146806
    ## iter  30 value 790.808603
    ## iter  40 value 775.213661
    ## iter  50 value 763.093598
    ## iter  60 value 753.334820
    ## iter  70 value 749.516830
    ## iter  80 value 748.414943
    ## iter  90 value 744.763528
    ## iter 100 value 743.291587
    ## final  value 743.291587 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1366.216346 
    ## iter  10 value 734.934164
    ## iter  20 value 605.902165
    ## iter  30 value 553.847965
    ## iter  40 value 527.577812
    ## iter  50 value 507.865814
    ## iter  60 value 491.324350
    ## iter  70 value 426.928966
    ## iter  80 value 364.523422
    ## iter  90 value 340.503693
    ## iter 100 value 333.026396
    ## final  value 333.026396 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1270.385318 
    ## iter  10 value 436.014332
    ## iter  20 value 231.875356
    ## iter  30 value 174.546944
    ## iter  40 value 155.962368
    ## iter  50 value 152.963080
    ## iter  60 value 152.071321
    ## iter  70 value 151.865710
    ## iter  80 value 151.707314
    ## iter  90 value 151.643340
    ## iter 100 value 151.633445
    ## final  value 151.633445 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1225.338045 
    ## iter  10 value 1097.428022
    ## iter  20 value 1092.809893
    ## iter  30 value 1065.099086
    ## iter  40 value 908.014321
    ## iter  50 value 813.052611
    ## iter  60 value 719.819299
    ## iter  70 value 672.705009
    ## iter  80 value 634.255808
    ## iter  90 value 582.324536
    ## iter 100 value 554.410939
    ## final  value 554.410939 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1401.602259 
    ## iter  10 value 983.033605
    ## iter  20 value 849.515329
    ## iter  30 value 648.629958
    ## iter  40 value 554.012020
    ## iter  50 value 503.220394
    ## iter  60 value 477.369041
    ## iter  70 value 454.603966
    ## iter  80 value 435.380278
    ## iter  90 value 406.036969
    ## iter 100 value 389.271346
    ## final  value 389.271346 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1256.455318 
    ## iter  10 value 380.893178
    ## iter  20 value 138.169872
    ## iter  30 value 46.625508
    ## iter  40 value 31.708432
    ## iter  50 value 20.233449
    ## iter  60 value 8.789157
    ## iter  70 value 8.063852
    ## iter  80 value 7.677477
    ## iter  90 value 7.411745
    ## iter 100 value 7.251520
    ## final  value 7.251520 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1286.661649 
    ## iter  10 value 1062.625250
    ## iter  20 value 851.363733
    ## iter  30 value 801.275263
    ## iter  40 value 774.444192
    ## iter  50 value 763.290672
    ## iter  60 value 756.721814
    ## iter  70 value 747.267512
    ## iter  80 value 735.942198
    ## iter  90 value 729.829926
    ## iter 100 value 725.486175
    ## final  value 725.486175 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1279.530398 
    ## iter  10 value 677.810684
    ## iter  20 value 387.030312
    ## iter  30 value 281.010602
    ## iter  40 value 233.771857
    ## iter  50 value 202.083702
    ## iter  60 value 180.695672
    ## iter  70 value 161.693436
    ## iter  80 value 150.001947
    ## iter  90 value 142.281380
    ## iter 100 value 120.152921
    ## final  value 120.152921 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1322.855719 
    ## iter  10 value 373.447750
    ## iter  20 value 118.627105
    ## iter  30 value 58.259893
    ## iter  40 value 38.902613
    ## iter  50 value 23.468251
    ## iter  60 value 18.204370
    ## iter  70 value 17.020523
    ## iter  80 value 16.605071
    ## iter  90 value 16.020421
    ## iter 100 value 15.128250
    ## final  value 15.128250 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1237.146488 
    ## iter  10 value 907.449742
    ## iter  20 value 841.111171
    ## iter  30 value 797.024674
    ## iter  40 value 780.616886
    ## iter  50 value 773.613652
    ## iter  60 value 769.679598
    ## iter  70 value 769.514706
    ## iter  80 value 769.505435
    ## iter  90 value 769.503916
    ## final  value 769.503900 
    ## converged
    ## # weights:  255
    ## initial  value 1264.024935 
    ## iter  10 value 612.924498
    ## iter  20 value 435.260588
    ## iter  30 value 377.112688
    ## iter  40 value 341.946769
    ## iter  50 value 329.536415
    ## iter  60 value 322.064223
    ## iter  70 value 314.943525
    ## iter  80 value 309.164091
    ## iter  90 value 302.662363
    ## iter 100 value 302.188994
    ## final  value 302.188994 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1396.631430 
    ## iter  10 value 536.007428
    ## iter  20 value 303.996332
    ## iter  30 value 217.302972
    ## iter  40 value 175.422007
    ## iter  50 value 160.331499
    ## iter  60 value 150.877133
    ## iter  70 value 145.868672
    ## iter  80 value 144.117508
    ## iter  90 value 143.761948
    ## iter 100 value 143.643892
    ## final  value 143.643892 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1294.124168 
    ## iter  10 value 1118.913836
    ## iter  20 value 962.171796
    ## iter  30 value 915.728785
    ## iter  40 value 853.366419
    ## iter  50 value 807.115106
    ## iter  60 value 715.118237
    ## iter  70 value 632.448675
    ## iter  80 value 580.111537
    ## iter  90 value 517.584193
    ## iter 100 value 446.827727
    ## final  value 446.827727 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1325.005942 
    ## iter  10 value 493.194140
    ## iter  20 value 307.191300
    ## iter  30 value 219.087952
    ## iter  40 value 190.428910
    ## iter  50 value 169.483767
    ## iter  60 value 150.736268
    ## iter  70 value 139.922111
    ## iter  80 value 133.515960
    ## iter  90 value 130.691770
    ## iter 100 value 123.523450
    ## final  value 123.523450 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1243.755530 
    ## iter  10 value 371.428128
    ## iter  20 value 101.755945
    ## iter  30 value 21.006798
    ## iter  40 value 5.640131
    ## iter  50 value 1.980081
    ## iter  60 value 1.627796
    ## iter  70 value 1.460802
    ## iter  80 value 1.359755
    ## iter  90 value 1.285108
    ## iter 100 value 1.235259
    ## final  value 1.235259 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1246.164414 
    ## iter  10 value 1068.631440
    ## iter  20 value 1002.045343
    ## iter  30 value 912.107908
    ## iter  40 value 810.641781
    ## iter  50 value 764.972824
    ## iter  60 value 746.439652
    ## iter  70 value 729.755871
    ## iter  80 value 717.541092
    ## iter  90 value 704.621091
    ## iter 100 value 694.621657
    ## final  value 694.621657 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1263.928369 
    ## iter  10 value 658.816862
    ## iter  20 value 452.255458
    ## iter  30 value 395.306392
    ## iter  40 value 350.182956
    ## iter  50 value 318.340767
    ## iter  60 value 294.723369
    ## iter  70 value 262.705816
    ## iter  80 value 235.003522
    ## iter  90 value 214.081988
    ## iter 100 value 192.493084
    ## final  value 192.493084 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1345.179108 
    ## iter  10 value 534.751032
    ## iter  20 value 218.775084
    ## iter  30 value 133.796982
    ## iter  40 value 43.887993
    ## iter  50 value 19.487193
    ## iter  60 value 7.461190
    ## iter  70 value 5.685370
    ## iter  80 value 5.393481
    ## iter  90 value 5.375502
    ## iter 100 value 5.372306
    ## final  value 5.372306 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1211.707365 
    ## iter  10 value 920.230919
    ## iter  20 value 852.826509
    ## iter  30 value 828.937097
    ## iter  40 value 819.052437
    ## iter  50 value 804.167043
    ## iter  60 value 798.245481
    ## iter  70 value 797.279768
    ## iter  80 value 795.946689
    ## iter  90 value 794.611216
    ## iter 100 value 790.312668
    ## final  value 790.312668 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1302.626811 
    ## iter  10 value 799.927901
    ## iter  20 value 707.877258
    ## iter  30 value 648.521913
    ## iter  40 value 613.393416
    ## iter  50 value 558.560690
    ## iter  60 value 433.882894
    ## iter  70 value 342.349685
    ## iter  80 value 323.923344
    ## iter  90 value 312.894721
    ## iter 100 value 308.332809
    ## final  value 308.332809 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1394.403349 
    ## iter  10 value 547.034949
    ## iter  20 value 350.446569
    ## iter  30 value 250.518905
    ## iter  40 value 201.494646
    ## iter  50 value 177.283034
    ## iter  60 value 167.701145
    ## iter  70 value 164.306738
    ## iter  80 value 163.112593
    ## iter  90 value 162.449229
    ## iter 100 value 162.089874
    ## final  value 162.089874 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1324.133153 
    ## iter  10 value 1055.458987
    ## iter  20 value 891.859670
    ## iter  30 value 817.011722
    ## iter  40 value 777.228582
    ## iter  50 value 752.486795
    ## iter  60 value 746.847665
    ## iter  70 value 745.563870
    ## iter  80 value 742.528506
    ## iter  90 value 737.577698
    ## iter 100 value 736.667800
    ## final  value 736.667800 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1351.436778 
    ## iter  10 value 573.081077
    ## iter  20 value 403.090303
    ## iter  30 value 340.897926
    ## iter  40 value 266.511372
    ## iter  50 value 208.858270
    ## iter  60 value 138.738041
    ## iter  70 value 104.926342
    ## iter  80 value 82.397345
    ## iter  90 value 64.477332
    ## iter 100 value 55.874196
    ## final  value 55.874196 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1412.594616 
    ## iter  10 value 461.779675
    ## iter  20 value 97.981111
    ## iter  30 value 33.818450
    ## iter  40 value 17.151177
    ## iter  50 value 9.665330
    ## iter  60 value 8.927919
    ## iter  70 value 8.582584
    ## iter  80 value 8.202644
    ## iter  90 value 7.610701
    ## iter 100 value 6.659918
    ## final  value 6.659918 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1306.903518 
    ## iter  10 value 1105.721329
    ## iter  20 value 1005.508106
    ## iter  30 value 974.828399
    ## iter  40 value 902.852167
    ## iter  50 value 833.809416
    ## iter  60 value 794.148730
    ## iter  70 value 767.069111
    ## iter  80 value 748.134167
    ## iter  90 value 729.074706
    ## iter 100 value 705.301435
    ## final  value 705.301435 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1191.926349 
    ## iter  10 value 634.899245
    ## iter  20 value 391.081102
    ## iter  30 value 298.890434
    ## iter  40 value 237.972291
    ## iter  50 value 189.917774
    ## iter  60 value 159.340042
    ## iter  70 value 138.251606
    ## iter  80 value 130.619554
    ## iter  90 value 124.120073
    ## iter 100 value 109.035750
    ## final  value 109.035750 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1285.799474 
    ## iter  10 value 427.734260
    ## iter  20 value 172.209972
    ## iter  30 value 98.016576
    ## iter  40 value 66.527811
    ## iter  50 value 50.710112
    ## iter  60 value 43.121924
    ## iter  70 value 38.021165
    ## iter  80 value 31.946994
    ## iter  90 value 28.829037
    ## iter 100 value 24.887766
    ## final  value 24.887766 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1253.464681 
    ## iter  10 value 1073.306429
    ## iter  20 value 948.461751
    ## iter  30 value 842.265862
    ## iter  40 value 810.258640
    ## iter  50 value 801.112341
    ## iter  60 value 791.979080
    ## iter  70 value 785.494920
    ## iter  80 value 770.457497
    ## iter  90 value 768.192647
    ## iter 100 value 767.614695
    ## final  value 767.614695 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1262.909903 
    ## iter  10 value 707.733946
    ## iter  20 value 533.178809
    ## iter  30 value 438.993796
    ## iter  40 value 391.506198
    ## iter  50 value 375.868104
    ## iter  60 value 366.534216
    ## iter  70 value 350.317129
    ## iter  80 value 338.606600
    ## iter  90 value 329.201343
    ## iter 100 value 307.840896
    ## final  value 307.840896 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1309.689087 
    ## iter  10 value 587.443814
    ## iter  20 value 337.415134
    ## iter  30 value 264.096733
    ## iter  40 value 195.878113
    ## iter  50 value 168.979732
    ## iter  60 value 162.166204
    ## iter  70 value 159.786187
    ## iter  80 value 159.034880
    ## iter  90 value 158.889438
    ## iter 100 value 158.807707
    ## final  value 158.807707 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1338.482061 
    ## iter  10 value 1123.221473
    ## iter  20 value 897.503375
    ## iter  30 value 793.795034
    ## iter  40 value 772.341802
    ## iter  50 value 753.779100
    ## iter  60 value 750.130324
    ## iter  70 value 749.831618
    ## iter  80 value 748.395864
    ## iter  90 value 748.029236
    ## iter 100 value 747.976784
    ## final  value 747.976784 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1448.953325 
    ## iter  10 value 654.381873
    ## iter  20 value 408.841103
    ## iter  30 value 325.896660
    ## iter  40 value 273.917455
    ## iter  50 value 247.296894
    ## iter  60 value 229.865065
    ## iter  70 value 219.821151
    ## iter  80 value 217.450942
    ## iter  90 value 215.133537
    ## iter 100 value 213.071336
    ## final  value 213.071336 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1299.308790 
    ## iter  10 value 576.101416
    ## iter  20 value 173.018774
    ## iter  30 value 60.305114
    ## iter  40 value 24.800176
    ## iter  50 value 17.892163
    ## iter  60 value 13.567045
    ## iter  70 value 11.832200
    ## iter  80 value 11.230668
    ## iter  90 value 4.035860
    ## iter 100 value 2.384514
    ## final  value 2.384514 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1262.480414 
    ## iter  10 value 942.506757
    ## iter  20 value 820.917944
    ## iter  30 value 800.950274
    ## iter  40 value 785.622240
    ## iter  50 value 774.571066
    ## iter  60 value 770.831370
    ## iter  70 value 770.623555
    ## iter  80 value 770.204302
    ## iter  90 value 769.971906
    ## iter 100 value 769.652233
    ## final  value 769.652233 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1315.182140 
    ## iter  10 value 582.613313
    ## iter  20 value 363.698993
    ## iter  30 value 247.031324
    ## iter  40 value 182.916440
    ## iter  50 value 136.444140
    ## iter  60 value 112.973402
    ## iter  70 value 97.901007
    ## iter  80 value 90.556185
    ## iter  90 value 84.264817
    ## iter 100 value 80.069944
    ## final  value 80.069944 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1315.238963 
    ## iter  10 value 463.795739
    ## iter  20 value 147.371542
    ## iter  30 value 71.953882
    ## iter  40 value 43.361148
    ## iter  50 value 29.974025
    ## iter  60 value 24.916636
    ## iter  70 value 22.066554
    ## iter  80 value 18.402958
    ## iter  90 value 16.647148
    ## iter 100 value 13.825343
    ## final  value 13.825343 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1223.143183 
    ## iter  10 value 889.056120
    ## iter  20 value 830.185765
    ## iter  30 value 806.675112
    ## iter  40 value 797.620254
    ## iter  50 value 790.359184
    ## iter  60 value 788.518956
    ## iter  70 value 779.676537
    ## iter  80 value 779.295552
    ## iter  90 value 779.292200
    ## iter  90 value 779.292195
    ## iter  90 value 779.292195
    ## final  value 779.292195 
    ## converged
    ## # weights:  255
    ## initial  value 1328.936392 
    ## iter  10 value 1028.229725
    ## iter  20 value 677.672025
    ## iter  30 value 526.434197
    ## iter  40 value 408.971223
    ## iter  50 value 367.065357
    ## iter  60 value 342.624492
    ## iter  70 value 328.663965
    ## iter  80 value 319.803963
    ## iter  90 value 316.477980
    ## iter 100 value 315.375330
    ## final  value 315.375330 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1312.751820 
    ## iter  10 value 683.659464
    ## iter  20 value 355.251631
    ## iter  30 value 232.135962
    ## iter  40 value 193.156438
    ## iter  50 value 168.297867
    ## iter  60 value 160.309374
    ## iter  70 value 157.884913
    ## iter  80 value 156.442718
    ## iter  90 value 155.845151
    ## iter 100 value 155.444869
    ## final  value 155.444869 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1238.281951 
    ## iter  10 value 1009.186420
    ## iter  20 value 895.186282
    ## iter  30 value 837.641145
    ## iter  40 value 805.980649
    ## iter  50 value 770.365162
    ## iter  60 value 755.476482
    ## iter  70 value 747.275854
    ## iter  80 value 739.367441
    ## iter  90 value 737.016829
    ## iter 100 value 730.562995
    ## final  value 730.562995 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1328.702956 
    ## iter  10 value 638.857801
    ## iter  20 value 338.319071
    ## iter  30 value 229.641643
    ## iter  40 value 164.789480
    ## iter  50 value 135.475469
    ## iter  60 value 116.455940
    ## iter  70 value 92.371695
    ## iter  80 value 84.245244
    ## iter  90 value 82.331959
    ## iter 100 value 80.926267
    ## final  value 80.926267 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1299.812818 
    ## iter  10 value 591.517684
    ## iter  20 value 418.776404
    ## iter  30 value 335.464975
    ## iter  40 value 241.137487
    ## iter  50 value 184.001757
    ## iter  60 value 157.874957
    ## iter  70 value 145.555926
    ## iter  80 value 137.585691
    ## iter  90 value 131.067905
    ## iter 100 value 126.702714
    ## final  value 126.702714 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1240.703902 
    ## iter  10 value 831.064307
    ## iter  20 value 741.968333
    ## iter  30 value 671.335777
    ## iter  40 value 628.080679
    ## iter  50 value 599.532803
    ## iter  60 value 576.095830
    ## iter  70 value 567.681464
    ## iter  80 value 560.349763
    ## iter  90 value 558.028309
    ## iter 100 value 554.418162
    ## final  value 554.418162 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1334.114155 
    ## iter  10 value 589.854494
    ## iter  20 value 353.324064
    ## iter  30 value 276.636186
    ## iter  40 value 227.557996
    ## iter  50 value 198.485161
    ## iter  60 value 177.783694
    ## iter  70 value 166.666391
    ## iter  80 value 159.818731
    ## iter  90 value 158.442710
    ## iter 100 value 157.165341
    ## final  value 157.165341 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1283.348600 
    ## iter  10 value 576.632060
    ## iter  20 value 165.524481
    ## iter  30 value 52.589508
    ## iter  40 value 11.514134
    ## iter  50 value 1.427805
    ## iter  60 value 0.103947
    ## iter  70 value 0.015276
    ## iter  80 value 0.004775
    ## iter  90 value 0.001642
    ## iter 100 value 0.000527
    ## final  value 0.000527 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1259.145948 
    ## iter  10 value 987.406168
    ## iter  20 value 858.486771
    ## iter  30 value 821.250141
    ## iter  40 value 798.513806
    ## iter  50 value 784.793583
    ## iter  60 value 771.714819
    ## iter  70 value 767.067391
    ## iter  80 value 766.306763
    ## iter  90 value 765.799134
    ## iter 100 value 765.443260
    ## final  value 765.443260 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1275.181953 
    ## iter  10 value 773.778460
    ## iter  20 value 693.693802
    ## iter  30 value 623.344620
    ## iter  40 value 585.177550
    ## iter  50 value 558.186542
    ## iter  60 value 536.193051
    ## iter  70 value 517.565637
    ## iter  80 value 428.444248
    ## iter  90 value 359.338347
    ## iter 100 value 331.999605
    ## final  value 331.999605 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1467.114373 
    ## iter  10 value 637.738895
    ## iter  20 value 370.968003
    ## iter  30 value 232.823556
    ## iter  40 value 182.737634
    ## iter  50 value 163.744932
    ## iter  60 value 156.286933
    ## iter  70 value 151.164236
    ## iter  80 value 149.865851
    ## iter  90 value 149.408490
    ## iter 100 value 149.035613
    ## final  value 149.035613 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1284.746414 
    ## iter  10 value 937.784370
    ## iter  20 value 868.734742
    ## iter  30 value 765.992382
    ## iter  40 value 694.508764
    ## iter  50 value 650.756709
    ## iter  60 value 607.513834
    ## iter  70 value 574.188329
    ## iter  80 value 533.641243
    ## iter  90 value 517.572909
    ## iter 100 value 510.686472
    ## final  value 510.686472 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1375.171912 
    ## iter  10 value 854.916423
    ## iter  20 value 627.276142
    ## iter  30 value 573.764525
    ## iter  40 value 536.630575
    ## iter  50 value 500.967312
    ## iter  60 value 470.189198
    ## iter  70 value 454.978516
    ## iter  80 value 445.662469
    ## iter  90 value 441.766879
    ## iter 100 value 434.724906
    ## final  value 434.724906 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1251.575515 
    ## iter  10 value 400.271618
    ## iter  20 value 128.567544
    ## iter  30 value 54.513358
    ## iter  40 value 33.219418
    ## iter  50 value 23.258599
    ## iter  60 value 15.875232
    ## iter  70 value 13.596759
    ## iter  80 value 11.377785
    ## iter  90 value 9.624069
    ## iter 100 value 8.822294
    ## final  value 8.822294 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1255.991412 
    ## iter  10 value 1123.149961
    ## iter  20 value 1122.678514
    ## iter  30 value 865.642527
    ## iter  40 value 800.845881
    ## iter  50 value 778.363561
    ## iter  60 value 737.757644
    ## iter  70 value 722.462938
    ## iter  80 value 709.312473
    ## iter  90 value 707.191254
    ## iter 100 value 705.687503
    ## final  value 705.687503 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1278.711210 
    ## iter  10 value 555.473891
    ## iter  20 value 355.365978
    ## iter  30 value 304.398379
    ## iter  40 value 271.863491
    ## iter  50 value 229.539818
    ## iter  60 value 204.394532
    ## iter  70 value 184.252090
    ## iter  80 value 177.153280
    ## iter  90 value 174.605268
    ## iter 100 value 173.901678
    ## final  value 173.901678 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1330.390581 
    ## iter  10 value 432.858830
    ## iter  20 value 113.273753
    ## iter  30 value 41.765959
    ## iter  40 value 17.635326
    ## iter  50 value 4.912656
    ## iter  60 value 1.783717
    ## iter  70 value 1.437332
    ## iter  80 value 1.393788
    ## iter  90 value 1.389629
    ## iter 100 value 1.387347
    ## final  value 1.387347 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1247.046331 
    ## iter  10 value 971.189511
    ## iter  20 value 892.410695
    ## iter  30 value 864.466328
    ## iter  40 value 839.584997
    ## iter  50 value 815.276425
    ## iter  60 value 785.125604
    ## iter  70 value 777.007156
    ## iter  80 value 775.067466
    ## iter  90 value 774.663916
    ## iter 100 value 773.361919
    ## final  value 773.361919 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1228.159060 
    ## iter  10 value 674.191860
    ## iter  20 value 478.381479
    ## iter  30 value 391.480231
    ## iter  40 value 366.591481
    ## iter  50 value 354.477369
    ## iter  60 value 347.443803
    ## iter  70 value 343.302103
    ## iter  80 value 337.266360
    ## iter  90 value 334.217835
    ## iter 100 value 331.973407
    ## final  value 331.973407 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1320.128078 
    ## iter  10 value 588.082781
    ## iter  20 value 334.064386
    ## iter  30 value 223.937843
    ## iter  40 value 189.364934
    ## iter  50 value 171.577627
    ## iter  60 value 163.823289
    ## iter  70 value 160.792554
    ## iter  80 value 159.012233
    ## iter  90 value 158.415799
    ## iter 100 value 157.978488
    ## final  value 157.978488 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1302.604606 
    ## iter  10 value 1125.547006
    ## iter  20 value 992.067123
    ## iter  30 value 911.901479
    ## iter  40 value 853.058521
    ## iter  50 value 805.114409
    ## iter  60 value 786.454948
    ## iter  70 value 783.020631
    ## iter  80 value 781.859905
    ## iter  90 value 780.706760
    ## iter 100 value 780.028187
    ## final  value 780.028187 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1174.832398 
    ## iter  10 value 533.797642
    ## iter  20 value 391.388972
    ## iter  30 value 335.852530
    ## iter  40 value 285.694728
    ## iter  50 value 264.146821
    ## iter  60 value 245.727230
    ## iter  70 value 229.338630
    ## iter  80 value 218.864823
    ## iter  90 value 215.159542
    ## iter 100 value 213.076562
    ## final  value 213.076562 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1313.530895 
    ## iter  10 value 674.472915
    ## iter  20 value 215.140501
    ## iter  30 value 138.744017
    ## iter  40 value 89.666253
    ## iter  50 value 49.260313
    ## iter  60 value 29.830343
    ## iter  70 value 23.386042
    ## iter  80 value 20.660737
    ## iter  90 value 15.448821
    ## iter 100 value 10.914067
    ## final  value 10.914067 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1355.042304 
    ## iter  10 value 1122.335268
    ## iter  20 value 1120.274585
    ## iter  30 value 1106.945601
    ## iter  40 value 967.276967
    ## iter  50 value 912.385339
    ## iter  60 value 897.541969
    ## iter  70 value 889.553694
    ## iter  80 value 886.648048
    ## iter  90 value 882.381543
    ## iter 100 value 874.936677
    ## final  value 874.936677 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1245.668427 
    ## iter  10 value 699.313296
    ## iter  20 value 392.525908
    ## iter  30 value 339.491964
    ## iter  40 value 297.530666
    ## iter  50 value 265.507240
    ## iter  60 value 251.826053
    ## iter  70 value 242.808242
    ## iter  80 value 237.603413
    ## iter  90 value 235.947670
    ## iter 100 value 235.078820
    ## final  value 235.078820 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1315.843197 
    ## iter  10 value 482.774231
    ## iter  20 value 230.566738
    ## iter  30 value 130.697306
    ## iter  40 value 95.077141
    ## iter  50 value 70.847271
    ## iter  60 value 47.901365
    ## iter  70 value 25.119874
    ## iter  80 value 16.365820
    ## iter  90 value 13.779805
    ## iter 100 value 12.856989
    ## final  value 12.856989 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1244.706318 
    ## iter  10 value 1016.700593
    ## iter  20 value 879.104963
    ## iter  30 value 833.597636
    ## iter  40 value 805.566073
    ## iter  50 value 795.496934
    ## iter  60 value 791.676100
    ## iter  70 value 789.737913
    ## iter  80 value 789.248634
    ## iter  90 value 787.295538
    ## iter 100 value 781.812992
    ## final  value 781.812992 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1300.251600 
    ## iter  10 value 719.541515
    ## iter  20 value 494.413662
    ## iter  30 value 403.707250
    ## iter  40 value 357.260268
    ## iter  50 value 332.261234
    ## iter  60 value 316.954275
    ## iter  70 value 311.640707
    ## iter  80 value 305.804684
    ## iter  90 value 300.264304
    ## iter 100 value 293.443839
    ## final  value 293.443839 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1294.958819 
    ## iter  10 value 451.003619
    ## iter  20 value 278.767869
    ## iter  30 value 200.489048
    ## iter  40 value 165.746904
    ## iter  50 value 151.250003
    ## iter  60 value 144.849045
    ## iter  70 value 143.751765
    ## iter  80 value 143.483746
    ## iter  90 value 143.389155
    ## iter 100 value 143.356170
    ## final  value 143.356170 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1249.771589 
    ## iter  10 value 928.623540
    ## iter  20 value 770.889021
    ## iter  30 value 702.221365
    ## iter  40 value 652.419096
    ## iter  50 value 603.261808
    ## iter  60 value 568.666773
    ## iter  70 value 558.111365
    ## iter  80 value 552.868142
    ## iter  90 value 546.571634
    ## iter 100 value 538.452394
    ## final  value 538.452394 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1284.610493 
    ## iter  10 value 563.679298
    ## iter  20 value 362.986537
    ## iter  30 value 310.120682
    ## iter  40 value 279.736158
    ## iter  50 value 245.360097
    ## iter  60 value 225.816579
    ## iter  70 value 221.671118
    ## iter  80 value 221.042912
    ## iter  90 value 220.529660
    ## iter 100 value 220.036068
    ## final  value 220.036068 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1252.137887 
    ## iter  10 value 446.347957
    ## iter  20 value 134.729950
    ## iter  30 value 68.843411
    ## iter  40 value 47.290460
    ## iter  50 value 10.949471
    ## iter  60 value 6.122285
    ## iter  70 value 5.278065
    ## iter  80 value 4.611819
    ## iter  90 value 4.097148
    ## iter 100 value 3.770556
    ## final  value 3.770556 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1299.412900 
    ## iter  10 value 896.147349
    ## iter  20 value 795.660881
    ## iter  30 value 761.364206
    ## iter  40 value 714.307990
    ## iter  50 value 688.947213
    ## iter  60 value 679.802100
    ## iter  70 value 675.180668
    ## iter  80 value 669.792776
    ## iter  90 value 665.955742
    ## iter 100 value 656.412404
    ## final  value 656.412404 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1208.719528 
    ## iter  10 value 707.477184
    ## iter  20 value 391.044289
    ## iter  30 value 306.591025
    ## iter  40 value 256.792377
    ## iter  50 value 229.548925
    ## iter  60 value 216.517875
    ## iter  70 value 209.222725
    ## iter  80 value 205.654603
    ## iter  90 value 202.880845
    ## iter 100 value 201.458021
    ## final  value 201.458021 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1274.816293 
    ## iter  10 value 388.173576
    ## iter  20 value 93.008143
    ## iter  30 value 45.225196
    ## iter  40 value 26.110124
    ## iter  50 value 20.906405
    ## iter  60 value 20.745428
    ## iter  70 value 20.733536
    ## iter  80 value 20.731185
    ## iter  90 value 20.730729
    ## iter 100 value 20.730442
    ## final  value 20.730442 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1261.341962 
    ## iter  10 value 1106.458589
    ## iter  20 value 1102.519070
    ## iter  30 value 892.617465
    ## iter  40 value 832.584580
    ## iter  50 value 804.494282
    ## iter  60 value 783.871104
    ## iter  70 value 772.123366
    ## iter  80 value 760.751587
    ## iter  90 value 757.545806
    ## iter 100 value 756.973852
    ## final  value 756.973852 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1380.216408 
    ## iter  10 value 661.866603
    ## iter  20 value 503.457084
    ## iter  30 value 433.764590
    ## iter  40 value 398.879069
    ## iter  50 value 374.541713
    ## iter  60 value 335.410532
    ## iter  70 value 323.727164
    ## iter  80 value 314.488281
    ## iter  90 value 307.913492
    ## iter 100 value 305.033469
    ## final  value 305.033469 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1270.071948 
    ## iter  10 value 618.976160
    ## iter  20 value 308.860824
    ## iter  30 value 219.481764
    ## iter  40 value 196.231911
    ## iter  50 value 182.738508
    ## iter  60 value 169.255390
    ## iter  70 value 163.820260
    ## iter  80 value 158.099213
    ## iter  90 value 156.283537
    ## iter 100 value 155.629575
    ## final  value 155.629575 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1280.208967 
    ## iter  10 value 1001.347959
    ## iter  20 value 793.098729
    ## iter  30 value 733.018412
    ## iter  40 value 690.818522
    ## iter  50 value 664.653631
    ## iter  60 value 659.654224
    ## iter  70 value 655.798125
    ## iter  80 value 652.290907
    ## iter  90 value 650.346401
    ## iter 100 value 649.107750
    ## final  value 649.107750 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1263.350440 
    ## iter  10 value 961.090442
    ## iter  20 value 844.938097
    ## iter  30 value 715.583502
    ## iter  40 value 603.563878
    ## iter  50 value 572.106524
    ## iter  60 value 544.137752
    ## iter  70 value 527.798794
    ## iter  80 value 521.323439
    ## iter  90 value 513.536843
    ## iter 100 value 512.065956
    ## final  value 512.065956 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1312.537317 
    ## iter  10 value 558.852373
    ## iter  20 value 242.784390
    ## iter  30 value 180.909429
    ## iter  40 value 115.210793
    ## iter  50 value 90.244610
    ## iter  60 value 63.004435
    ## iter  70 value 53.867002
    ## iter  80 value 49.523927
    ## iter  90 value 48.365621
    ## iter 100 value 47.621030
    ## final  value 47.621030 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1251.627306 
    ## iter  10 value 966.491302
    ## iter  20 value 936.694253
    ## iter  30 value 908.007453
    ## iter  40 value 871.829398
    ## iter  50 value 858.405037
    ## iter  60 value 848.394626
    ## iter  70 value 845.562158
    ## iter  80 value 830.820765
    ## iter  90 value 786.869134
    ## iter 100 value 735.616622
    ## final  value 735.616622 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1224.825749 
    ## iter  10 value 620.309497
    ## iter  20 value 308.112004
    ## iter  30 value 229.796890
    ## iter  40 value 187.171061
    ## iter  50 value 149.990667
    ## iter  60 value 117.143640
    ## iter  70 value 95.857028
    ## iter  80 value 78.050209
    ## iter  90 value 67.075233
    ## iter 100 value 57.230072
    ## final  value 57.230072 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1241.522789 
    ## iter  10 value 438.881057
    ## iter  20 value 97.459202
    ## iter  30 value 35.272731
    ## iter  40 value 14.706026
    ## iter  50 value 5.728480
    ## iter  60 value 2.052470
    ## iter  70 value 1.923156
    ## iter  80 value 1.914707
    ## iter  90 value 1.911081
    ## iter 100 value 1.910159
    ## final  value 1.910159 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1286.496552 
    ## iter  10 value 1065.739713
    ## iter  20 value 922.320581
    ## iter  30 value 841.265336
    ## iter  40 value 815.498919
    ## iter  50 value 800.368914
    ## iter  60 value 787.811712
    ## iter  70 value 780.875891
    ## iter  80 value 776.429746
    ## iter  90 value 776.091576
    ## iter 100 value 775.984409
    ## final  value 775.984409 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1202.657621 
    ## iter  10 value 682.917814
    ## iter  20 value 527.504660
    ## iter  30 value 436.786188
    ## iter  40 value 395.879840
    ## iter  50 value 351.802963
    ## iter  60 value 335.733214
    ## iter  70 value 327.958090
    ## iter  80 value 323.241536
    ## iter  90 value 319.836615
    ## iter 100 value 318.101957
    ## final  value 318.101957 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1303.654345 
    ## iter  10 value 486.306103
    ## iter  20 value 282.977055
    ## iter  30 value 200.366608
    ## iter  40 value 165.514827
    ## iter  50 value 158.648754
    ## iter  60 value 154.911143
    ## iter  70 value 153.589645
    ## iter  80 value 153.095619
    ## iter  90 value 152.566788
    ## iter 100 value 152.478297
    ## final  value 152.478297 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1311.749486 
    ## iter  10 value 870.097628
    ## iter  20 value 753.521559
    ## iter  30 value 697.059947
    ## iter  40 value 655.945574
    ## iter  50 value 611.342253
    ## iter  60 value 573.959014
    ## iter  70 value 545.598427
    ## iter  80 value 538.585388
    ## iter  90 value 535.170688
    ## iter 100 value 534.238146
    ## final  value 534.238146 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1291.035671 
    ## iter  10 value 735.690985
    ## iter  20 value 562.507428
    ## iter  30 value 485.283433
    ## iter  40 value 438.824165
    ## iter  50 value 398.546669
    ## iter  60 value 360.919764
    ## iter  70 value 339.806704
    ## iter  80 value 322.724239
    ## iter  90 value 317.759957
    ## iter 100 value 312.262168
    ## final  value 312.262168 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1345.544618 
    ## iter  10 value 460.110755
    ## iter  20 value 156.732059
    ## iter  30 value 76.124957
    ## iter  40 value 50.365624
    ## iter  50 value 34.659609
    ## iter  60 value 28.665376
    ## iter  70 value 27.618650
    ## iter  80 value 27.143183
    ## iter  90 value 24.858853
    ## iter 100 value 21.518767
    ## final  value 21.518767 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1374.147154 
    ## iter  10 value 890.089706
    ## iter  20 value 803.785215
    ## iter  30 value 756.510655
    ## iter  40 value 706.494056
    ## iter  50 value 662.344799
    ## iter  60 value 641.271912
    ## iter  70 value 630.584157
    ## iter  80 value 624.106575
    ## iter  90 value 612.990970
    ## iter 100 value 609.671834
    ## final  value 609.671834 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1289.128463 
    ## iter  10 value 615.117575
    ## iter  20 value 417.460405
    ## iter  30 value 371.473008
    ## iter  40 value 324.284570
    ## iter  50 value 298.477117
    ## iter  60 value 283.365832
    ## iter  70 value 275.831928
    ## iter  80 value 272.359884
    ## iter  90 value 270.476063
    ## iter 100 value 270.113978
    ## final  value 270.113978 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1337.669202 
    ## iter  10 value 484.592950
    ## iter  20 value 267.535601
    ## iter  30 value 141.139054
    ## iter  40 value 81.324944
    ## iter  50 value 60.116345
    ## iter  60 value 48.298075
    ## iter  70 value 42.316884
    ## iter  80 value 39.195014
    ## iter  90 value 37.365281
    ## iter 100 value 36.490493
    ## final  value 36.490493 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1318.229409 
    ## iter  10 value 882.910054
    ## iter  20 value 827.751985
    ## iter  30 value 800.699252
    ## iter  40 value 788.291705
    ## iter  50 value 780.889964
    ## iter  60 value 769.598582
    ## iter  70 value 764.441208
    ## iter  80 value 762.087743
    ## iter  90 value 761.080882
    ## iter 100 value 760.895687
    ## final  value 760.895687 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1279.746254 
    ## iter  10 value 634.192997
    ## iter  20 value 469.616052
    ## iter  30 value 425.465028
    ## iter  40 value 393.877852
    ## iter  50 value 358.782631
    ## iter  60 value 339.419693
    ## iter  70 value 327.166888
    ## iter  80 value 321.178694
    ## iter  90 value 318.875583
    ## iter 100 value 318.005424
    ## final  value 318.005424 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1473.234182 
    ## iter  10 value 602.732199
    ## iter  20 value 326.289391
    ## iter  30 value 235.595788
    ## iter  40 value 202.700792
    ## iter  50 value 174.861972
    ## iter  60 value 161.524442
    ## iter  70 value 155.718112
    ## iter  80 value 152.606053
    ## iter  90 value 150.422091
    ## iter 100 value 148.990048
    ## final  value 148.990048 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1391.434886 
    ## iter  10 value 943.196166
    ## iter  20 value 835.221927
    ## iter  30 value 819.485240
    ## iter  40 value 792.280744
    ## iter  50 value 767.622851
    ## iter  60 value 760.027843
    ## iter  70 value 758.959310
    ## iter  80 value 757.930395
    ## iter  90 value 757.803307
    ## iter 100 value 757.555363
    ## final  value 757.555363 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1274.590873 
    ## iter  10 value 618.693132
    ## iter  20 value 410.686578
    ## iter  30 value 326.531830
    ## iter  40 value 285.320039
    ## iter  50 value 239.985103
    ## iter  60 value 211.845007
    ## iter  70 value 192.196093
    ## iter  80 value 174.931027
    ## iter  90 value 169.730369
    ## iter 100 value 164.291159
    ## final  value 164.291159 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1247.086860 
    ## iter  10 value 407.325059
    ## iter  20 value 149.484676
    ## iter  30 value 67.053022
    ## iter  40 value 30.393375
    ## iter  50 value 20.676437
    ## iter  60 value 12.317860
    ## iter  70 value 7.223599
    ## iter  80 value 6.578581
    ## iter  90 value 6.314844
    ## iter 100 value 6.175931
    ## final  value 6.175931 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1220.888326 
    ## iter  10 value 915.602793
    ## iter  20 value 839.975631
    ## iter  30 value 828.582211
    ## iter  40 value 818.793657
    ## iter  50 value 777.127752
    ## iter  60 value 764.273535
    ## iter  70 value 753.472463
    ## iter  80 value 696.954748
    ## iter  90 value 655.187194
    ## iter 100 value 633.624271
    ## final  value 633.624271 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1364.697384 
    ## iter  10 value 760.016731
    ## iter  20 value 556.728758
    ## iter  30 value 411.790977
    ## iter  40 value 338.137625
    ## iter  50 value 273.239063
    ## iter  60 value 245.805810
    ## iter  70 value 217.248769
    ## iter  80 value 198.765300
    ## iter  90 value 182.823084
    ## iter 100 value 164.509457
    ## final  value 164.509457 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1241.661530 
    ## iter  10 value 519.558506
    ## iter  20 value 171.966098
    ## iter  30 value 105.940980
    ## iter  40 value 70.408609
    ## iter  50 value 39.958211
    ## iter  60 value 25.190543
    ## iter  70 value 15.451420
    ## iter  80 value 9.547582
    ## iter  90 value 8.481313
    ## iter 100 value 8.207274
    ## final  value 8.207274 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1224.255293 
    ## iter  10 value 991.690052
    ## iter  20 value 889.557768
    ## iter  30 value 838.922986
    ## iter  40 value 798.338342
    ## iter  50 value 763.475306
    ## iter  60 value 757.752937
    ## iter  70 value 752.607951
    ## iter  80 value 749.691217
    ## iter  90 value 749.465249
    ## iter 100 value 749.452225
    ## final  value 749.452225 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1261.458796 
    ## iter  10 value 612.089587
    ## iter  20 value 444.431715
    ## iter  30 value 355.242824
    ## iter  40 value 332.727472
    ## iter  50 value 323.772519
    ## iter  60 value 318.849750
    ## iter  70 value 312.761678
    ## iter  80 value 303.384638
    ## iter  90 value 298.130006
    ## iter 100 value 296.604167
    ## final  value 296.604167 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1277.913041 
    ## iter  10 value 556.072163
    ## iter  20 value 299.939393
    ## iter  30 value 211.026613
    ## iter  40 value 172.066540
    ## iter  50 value 163.053690
    ## iter  60 value 160.992040
    ## iter  70 value 157.063468
    ## iter  80 value 152.187330
    ## iter  90 value 151.329456
    ## iter 100 value 151.060722
    ## final  value 151.060722 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1231.313993 
    ## iter  10 value 881.629745
    ## iter  20 value 758.082746
    ## iter  30 value 688.341088
    ## iter  40 value 628.379031
    ## iter  50 value 589.755881
    ## iter  60 value 565.979742
    ## iter  70 value 547.052636
    ## iter  80 value 537.005649
    ## iter  90 value 533.431168
    ## iter 100 value 527.280862
    ## final  value 527.280862 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1275.953367 
    ## iter  10 value 630.111028
    ## iter  20 value 300.862872
    ## iter  30 value 229.461623
    ## iter  40 value 176.442254
    ## iter  50 value 144.676316
    ## iter  60 value 125.717930
    ## iter  70 value 120.351262
    ## iter  80 value 118.657141
    ## iter  90 value 117.124902
    ## iter 100 value 116.450619
    ## final  value 116.450619 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1206.337805 
    ## iter  10 value 451.536103
    ## iter  20 value 196.938903
    ## iter  30 value 95.268570
    ## iter  40 value 52.113127
    ## iter  50 value 32.242528
    ## iter  60 value 14.949459
    ## iter  70 value 6.306749
    ## iter  80 value 5.683824
    ## iter  90 value 5.459072
    ## iter 100 value 5.219202
    ## final  value 5.219202 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1159.406598 
    ## iter  10 value 904.122828
    ## iter  20 value 805.446994
    ## iter  30 value 789.254993
    ## iter  40 value 782.316071
    ## iter  50 value 773.563078
    ## iter  60 value 756.521184
    ## iter  70 value 752.351137
    ## iter  80 value 749.130384
    ## iter  90 value 747.217679
    ## iter 100 value 743.976055
    ## final  value 743.976055 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1215.847806 
    ## iter  10 value 562.122061
    ## iter  20 value 407.523937
    ## iter  30 value 338.350617
    ## iter  40 value 275.509814
    ## iter  50 value 238.203595
    ## iter  60 value 210.054493
    ## iter  70 value 191.682347
    ## iter  80 value 178.119729
    ## iter  90 value 163.735660
    ## iter 100 value 161.475415
    ## final  value 161.475415 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1330.102828 
    ## iter  10 value 448.722766
    ## iter  20 value 175.627377
    ## iter  30 value 129.133122
    ## iter  40 value 110.267416
    ## iter  50 value 93.435905
    ## iter  60 value 56.322007
    ## iter  70 value 47.830033
    ## iter  80 value 42.711403
    ## iter  90 value 38.858609
    ## iter 100 value 35.489582
    ## final  value 35.489582 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1285.040374 
    ## iter  10 value 885.342493
    ## iter  20 value 814.256674
    ## iter  30 value 786.075963
    ## iter  40 value 762.556578
    ## iter  50 value 756.820597
    ## iter  60 value 753.665824
    ## iter  70 value 752.204435
    ## iter  80 value 751.945198
    ## final  value 751.939137 
    ## converged
    ## # weights:  255
    ## initial  value 1341.251867 
    ## iter  10 value 790.005793
    ## iter  20 value 600.930670
    ## iter  30 value 474.792775
    ## iter  40 value 408.219800
    ## iter  50 value 363.919650
    ## iter  60 value 347.762734
    ## iter  70 value 335.410253
    ## iter  80 value 327.823815
    ## iter  90 value 322.009235
    ## iter 100 value 319.643476
    ## final  value 319.643476 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1295.303804 
    ## iter  10 value 559.925579
    ## iter  20 value 319.211515
    ## iter  30 value 225.652303
    ## iter  40 value 186.569047
    ## iter  50 value 167.253543
    ## iter  60 value 162.132737
    ## iter  70 value 160.870634
    ## iter  80 value 160.607832
    ## iter  90 value 160.551514
    ## iter 100 value 160.523530
    ## final  value 160.523530 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1263.000964 
    ## iter  10 value 856.280769
    ## iter  20 value 765.563777
    ## iter  30 value 727.442260
    ## iter  40 value 698.386645
    ## iter  50 value 686.795588
    ## iter  60 value 686.014828
    ## iter  70 value 682.224347
    ## iter  80 value 679.682162
    ## iter  90 value 678.280708
    ## iter 100 value 676.976555
    ## final  value 676.976555 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1349.761828 
    ## iter  10 value 644.738071
    ## iter  20 value 447.055145
    ## iter  30 value 376.183640
    ## iter  40 value 326.758793
    ## iter  50 value 250.852925
    ## iter  60 value 215.818185
    ## iter  70 value 190.726157
    ## iter  80 value 160.389131
    ## iter  90 value 148.423598
    ## iter 100 value 144.205452
    ## final  value 144.205452 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1296.615341 
    ## iter  10 value 463.721685
    ## iter  20 value 168.482923
    ## iter  30 value 98.835622
    ## iter  40 value 46.061713
    ## iter  50 value 26.102810
    ## iter  60 value 17.003527
    ## iter  70 value 13.651141
    ## iter  80 value 9.629303
    ## iter  90 value 8.270559
    ## iter 100 value 7.819010
    ## final  value 7.819010 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1243.619236 
    ## iter  10 value 853.528270
    ## iter  20 value 764.178728
    ## iter  30 value 713.092754
    ## iter  40 value 675.566500
    ## iter  50 value 656.461990
    ## iter  60 value 647.160125
    ## iter  70 value 641.006799
    ## iter  80 value 639.288882
    ## iter  90 value 638.684963
    ## iter 100 value 638.332069
    ## final  value 638.332069 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1276.892600 
    ## iter  10 value 635.429728
    ## iter  20 value 397.082804
    ## iter  30 value 340.492268
    ## iter  40 value 272.166101
    ## iter  50 value 239.140605
    ## iter  60 value 220.439495
    ## iter  70 value 206.041679
    ## iter  80 value 195.830414
    ## iter  90 value 194.127137
    ## iter 100 value 193.974551
    ## final  value 193.974551 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1310.339115 
    ## iter  10 value 481.455901
    ## iter  20 value 254.763511
    ## iter  30 value 158.015522
    ## iter  40 value 117.562495
    ## iter  50 value 95.825185
    ## iter  60 value 76.681564
    ## iter  70 value 52.955483
    ## iter  80 value 47.965677
    ## iter  90 value 43.254038
    ## iter 100 value 41.393714
    ## final  value 41.393714 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1286.816819 
    ## iter  10 value 1096.941373
    ## iter  20 value 1039.859582
    ## iter  30 value 841.316742
    ## iter  40 value 780.560269
    ## iter  50 value 760.030844
    ## iter  60 value 742.912683
    ## iter  70 value 727.549565
    ## iter  80 value 723.701479
    ## iter  90 value 723.461570
    ## iter 100 value 723.448661
    ## final  value 723.448661 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1249.972033 
    ## iter  10 value 622.051354
    ## iter  20 value 458.434732
    ## iter  30 value 380.206418
    ## iter  40 value 342.178974
    ## iter  50 value 328.437796
    ## iter  60 value 317.982262
    ## iter  70 value 312.314739
    ## iter  80 value 308.948669
    ## iter  90 value 307.528937
    ## iter 100 value 302.958992
    ## final  value 302.958992 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1354.247193 
    ## iter  10 value 497.587274
    ## iter  20 value 274.008924
    ## iter  30 value 212.427980
    ## iter  40 value 173.313203
    ## iter  50 value 164.730665
    ## iter  60 value 160.197371
    ## iter  70 value 156.846139
    ## iter  80 value 156.099593
    ## iter  90 value 155.930902
    ## iter 100 value 155.852379
    ## final  value 155.852379 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1246.797029 
    ## iter  10 value 809.155457
    ## iter  20 value 740.336249
    ## iter  30 value 699.512963
    ## iter  40 value 667.223409
    ## iter  50 value 643.523134
    ## iter  60 value 629.612808
    ## iter  70 value 627.033311
    ## iter  80 value 626.716367
    ## iter  90 value 626.524200
    ## iter 100 value 623.967481
    ## final  value 623.967481 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1329.922715 
    ## iter  10 value 1086.946147
    ## iter  20 value 792.791240
    ## iter  30 value 723.863466
    ## iter  40 value 664.110290
    ## iter  50 value 606.450797
    ## iter  60 value 416.046589
    ## iter  70 value 278.892031
    ## iter  80 value 243.703614
    ## iter  90 value 217.325933
    ## iter 100 value 199.774349
    ## final  value 199.774349 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1386.801455 
    ## iter  10 value 534.247982
    ## iter  20 value 152.083342
    ## iter  30 value 75.720287
    ## iter  40 value 45.852759
    ## iter  50 value 33.362729
    ## iter  60 value 22.833169
    ## iter  70 value 19.276777
    ## iter  80 value 18.781189
    ## iter  90 value 17.516550
    ## iter 100 value 16.663826
    ## final  value 16.663826 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1301.678653 
    ## iter  10 value 1111.115839
    ## iter  20 value 1110.430996
    ## iter  30 value 1110.404236
    ## final  value 1110.404082 
    ## converged
    ## # weights:  255
    ## initial  value 1327.750636 
    ## iter  10 value 831.957206
    ## iter  20 value 542.053129
    ## iter  30 value 389.897746
    ## iter  40 value 344.876119
    ## iter  50 value 308.442516
    ## iter  60 value 277.363108
    ## iter  70 value 250.722618
    ## iter  80 value 235.413367
    ## iter  90 value 231.595662
    ## iter 100 value 222.445989
    ## final  value 222.445989 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1323.117966 
    ## iter  10 value 544.832708
    ## iter  20 value 242.420548
    ## iter  30 value 122.278909
    ## iter  40 value 79.486930
    ## iter  50 value 68.059089
    ## iter  60 value 60.462847
    ## iter  70 value 54.025407
    ## iter  80 value 39.564087
    ## iter  90 value 28.394153
    ## iter 100 value 27.557765
    ## final  value 27.557765 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1354.417500 
    ## iter  10 value 1051.954913
    ## iter  20 value 880.790870
    ## iter  30 value 841.087337
    ## iter  40 value 821.194436
    ## iter  50 value 791.963445
    ## iter  60 value 779.794875
    ## iter  70 value 774.645526
    ## iter  80 value 771.947682
    ## iter  90 value 771.196167
    ## iter 100 value 766.687378
    ## final  value 766.687378 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1277.882066 
    ## iter  10 value 783.059468
    ## iter  20 value 522.104168
    ## iter  30 value 433.373165
    ## iter  40 value 391.029221
    ## iter  50 value 361.267748
    ## iter  60 value 331.787646
    ## iter  70 value 320.732908
    ## iter  80 value 312.414975
    ## iter  90 value 308.442833
    ## iter 100 value 305.456364
    ## final  value 305.456364 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1296.870462 
    ## iter  10 value 570.749142
    ## iter  20 value 310.771671
    ## iter  30 value 238.181341
    ## iter  40 value 178.827452
    ## iter  50 value 154.796612
    ## iter  60 value 151.095968
    ## iter  70 value 147.477334
    ## iter  80 value 145.605198
    ## iter  90 value 144.291655
    ## iter 100 value 143.797275
    ## final  value 143.797275 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1286.432834 
    ## iter  10 value 1113.501525
    ## iter  20 value 1110.420354
    ## iter  30 value 1110.384346
    ## iter  40 value 1007.061733
    ## iter  50 value 911.048918
    ## iter  60 value 878.687524
    ## iter  70 value 866.201946
    ## iter  80 value 857.825849
    ## iter  90 value 852.105285
    ## iter 100 value 850.796099
    ## final  value 850.796099 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1292.119037 
    ## iter  10 value 625.320066
    ## iter  20 value 408.052884
    ## iter  30 value 337.925129
    ## iter  40 value 273.649149
    ## iter  50 value 221.962224
    ## iter  60 value 202.445257
    ## iter  70 value 190.325144
    ## iter  80 value 166.151181
    ## iter  90 value 134.785521
    ## iter 100 value 118.984723
    ## final  value 118.984723 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1233.791463 
    ## iter  10 value 310.385950
    ## iter  20 value 75.794553
    ## iter  30 value 34.399436
    ## iter  40 value 12.334201
    ## iter  50 value 3.302099
    ## iter  60 value 1.960303
    ## iter  70 value 1.801614
    ## iter  80 value 1.657778
    ## iter  90 value 1.537982
    ## iter 100 value 1.448299
    ## final  value 1.448299 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1235.749434 
    ## iter  10 value 1100.181133
    ## iter  20 value 955.657153
    ## iter  30 value 873.659448
    ## iter  40 value 827.583193
    ## iter  50 value 789.876511
    ## iter  60 value 749.112827
    ## iter  70 value 720.673479
    ## iter  80 value 699.743975
    ## iter  90 value 687.296202
    ## iter 100 value 680.643680
    ## final  value 680.643680 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1246.862838 
    ## iter  10 value 608.230409
    ## iter  20 value 339.147264
    ## iter  30 value 267.896584
    ## iter  40 value 203.234365
    ## iter  50 value 138.391679
    ## iter  60 value 110.476136
    ## iter  70 value 92.304353
    ## iter  80 value 82.357750
    ## iter  90 value 79.296678
    ## iter 100 value 78.760717
    ## final  value 78.760717 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1329.068151 
    ## iter  10 value 899.562559
    ## iter  20 value 447.377708
    ## iter  30 value 291.831493
    ## iter  40 value 194.769875
    ## iter  50 value 128.013290
    ## iter  60 value 80.828160
    ## iter  70 value 53.473260
    ## iter  80 value 40.832649
    ## iter  90 value 35.041553
    ## iter 100 value 31.266820
    ## final  value 31.266820 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1271.621832 
    ## iter  10 value 942.352832
    ## iter  20 value 835.247360
    ## iter  30 value 811.512699
    ## iter  40 value 784.906612
    ## iter  50 value 770.424627
    ## iter  60 value 764.348313
    ## iter  70 value 762.572246
    ## iter  80 value 761.777826
    ## iter  90 value 759.900633
    ## iter 100 value 758.617040
    ## final  value 758.617040 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1293.883761 
    ## iter  10 value 783.871223
    ## iter  20 value 634.469365
    ## iter  30 value 583.860188
    ## iter  40 value 553.197646
    ## iter  50 value 527.703146
    ## iter  60 value 496.305353
    ## iter  70 value 405.549608
    ## iter  80 value 353.818282
    ## iter  90 value 330.943652
    ## iter 100 value 315.348191
    ## final  value 315.348191 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1263.661500 
    ## iter  10 value 543.977468
    ## iter  20 value 299.703165
    ## iter  30 value 193.449178
    ## iter  40 value 172.920132
    ## iter  50 value 157.538035
    ## iter  60 value 150.471103
    ## iter  70 value 147.597499
    ## iter  80 value 145.707974
    ## iter  90 value 142.719856
    ## iter 100 value 141.937030
    ## final  value 141.937030 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1221.464718 
    ## iter  10 value 948.499846
    ## iter  20 value 921.786191
    ## iter  30 value 888.460504
    ## iter  40 value 860.125164
    ## iter  50 value 822.679146
    ## iter  60 value 759.777560
    ## iter  70 value 720.168343
    ## iter  80 value 702.812459
    ## iter  90 value 682.752535
    ## iter 100 value 673.881864
    ## final  value 673.881864 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1243.607571 
    ## iter  10 value 505.985811
    ## iter  20 value 278.810829
    ## iter  30 value 211.968441
    ## iter  40 value 175.392741
    ## iter  50 value 150.209357
    ## iter  60 value 120.593717
    ## iter  70 value 105.264767
    ## iter  80 value 92.517280
    ## iter  90 value 86.993737
    ## iter 100 value 85.185908
    ## final  value 85.185908 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1395.383320 
    ## iter  10 value 926.172730
    ## iter  20 value 538.366679
    ## iter  30 value 442.554004
    ## iter  40 value 355.525110
    ## iter  50 value 305.917327
    ## iter  60 value 269.453004
    ## iter  70 value 235.983034
    ## iter  80 value 198.982944
    ## iter  90 value 172.474963
    ## iter 100 value 160.283052
    ## final  value 160.283052 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1264.895526 
    ## iter  10 value 984.077722
    ## iter  20 value 867.843677
    ## iter  30 value 822.990339
    ## iter  40 value 777.103813
    ## iter  50 value 724.064287
    ## iter  60 value 699.818002
    ## iter  70 value 690.642088
    ## iter  80 value 679.819121
    ## iter  90 value 671.392085
    ## iter 100 value 669.737551
    ## final  value 669.737551 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1333.485716 
    ## iter  10 value 991.051915
    ## iter  20 value 873.657124
    ## iter  30 value 602.997159
    ## iter  40 value 421.294280
    ## iter  50 value 366.889865
    ## iter  60 value 334.256065
    ## iter  70 value 310.127647
    ## iter  80 value 299.242712
    ## iter  90 value 294.116196
    ## iter 100 value 292.775371
    ## final  value 292.775371 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1358.705282 
    ## iter  10 value 387.827721
    ## iter  20 value 121.012825
    ## iter  30 value 60.774561
    ## iter  40 value 33.989932
    ## iter  50 value 25.556881
    ## iter  60 value 23.641389
    ## iter  70 value 23.517402
    ## iter  80 value 23.507903
    ## iter  90 value 23.506134
    ## iter 100 value 23.505222
    ## final  value 23.505222 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1297.025560 
    ## iter  10 value 1135.376034
    ## iter  20 value 1123.867178
    ## iter  30 value 895.884373
    ## iter  40 value 832.952758
    ## iter  50 value 805.706899
    ## iter  60 value 784.016641
    ## iter  70 value 774.688004
    ## iter  80 value 767.531321
    ## iter  90 value 764.384749
    ## iter 100 value 763.803970
    ## final  value 763.803970 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1318.841967 
    ## iter  10 value 652.046433
    ## iter  20 value 466.007786
    ## iter  30 value 390.613224
    ## iter  40 value 352.132668
    ## iter  50 value 334.715678
    ## iter  60 value 325.595398
    ## iter  70 value 320.369306
    ## iter  80 value 318.290969
    ## iter  90 value 317.233096
    ## iter 100 value 316.418781
    ## final  value 316.418781 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1205.621047 
    ## iter  10 value 503.216234
    ## iter  20 value 307.099265
    ## iter  30 value 211.191003
    ## iter  40 value 165.675739
    ## iter  50 value 153.578851
    ## iter  60 value 150.726298
    ## iter  70 value 147.978334
    ## iter  80 value 146.652912
    ## iter  90 value 145.744132
    ## iter 100 value 145.551112
    ## final  value 145.551112 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1206.248883 
    ## iter  10 value 925.104295
    ## iter  20 value 805.858317
    ## iter  30 value 776.534440
    ## iter  40 value 755.220513
    ## iter  50 value 730.565186
    ## iter  60 value 703.469517
    ## iter  70 value 689.301645
    ## iter  80 value 672.663378
    ## iter  90 value 632.853101
    ## iter 100 value 610.096690
    ## final  value 610.096690 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1214.277723 
    ## iter  10 value 714.165508
    ## iter  20 value 583.377158
    ## iter  30 value 523.111976
    ## iter  40 value 469.908476
    ## iter  50 value 423.777850
    ## iter  60 value 373.400432
    ## iter  70 value 331.422049
    ## iter  80 value 292.238529
    ## iter  90 value 269.996848
    ## iter 100 value 262.633222
    ## final  value 262.633222 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1298.194975 
    ## iter  10 value 453.762434
    ## iter  20 value 143.359087
    ## iter  30 value 39.709599
    ## iter  40 value 11.538724
    ## iter  50 value 6.572053
    ## iter  60 value 4.789245
    ## iter  70 value 4.292869
    ## iter  80 value 4.050785
    ## iter  90 value 3.776137
    ## iter 100 value 3.601243
    ## final  value 3.601243 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1188.327505 
    ## iter  10 value 869.079184
    ## iter  20 value 786.707912
    ## iter  30 value 740.697811
    ## iter  40 value 683.424574
    ## iter  50 value 647.575717
    ## iter  60 value 618.417966
    ## iter  70 value 593.129517
    ## iter  80 value 570.628639
    ## iter  90 value 558.173019
    ## iter 100 value 553.673093
    ## final  value 553.673093 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1283.754256 
    ## iter  10 value 686.022973
    ## iter  20 value 407.868400
    ## iter  30 value 325.539310
    ## iter  40 value 271.660699
    ## iter  50 value 239.048974
    ## iter  60 value 217.326711
    ## iter  70 value 196.584207
    ## iter  80 value 181.976324
    ## iter  90 value 175.503818
    ## iter 100 value 174.499895
    ## final  value 174.499895 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1262.057225 
    ## iter  10 value 457.535653
    ## iter  20 value 220.688637
    ## iter  30 value 149.769424
    ## iter  40 value 119.828835
    ## iter  50 value 100.740770
    ## iter  60 value 83.955571
    ## iter  70 value 76.241373
    ## iter  80 value 67.158700
    ## iter  90 value 64.460312
    ## iter 100 value 64.259189
    ## final  value 64.259189 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1205.979804 
    ## iter  10 value 998.607230
    ## iter  20 value 865.365740
    ## iter  30 value 819.469053
    ## iter  40 value 789.650085
    ## iter  50 value 780.465191
    ## iter  60 value 772.613264
    ## iter  70 value 768.891490
    ## iter  80 value 762.049520
    ## iter  90 value 760.199933
    ## iter 100 value 760.121687
    ## final  value 760.121687 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1239.353592 
    ## iter  10 value 663.408582
    ## iter  20 value 463.639508
    ## iter  30 value 396.662920
    ## iter  40 value 333.151837
    ## iter  50 value 311.516872
    ## iter  60 value 301.202634
    ## iter  70 value 296.931137
    ## iter  80 value 294.703036
    ## iter  90 value 294.164474
    ## iter 100 value 293.921369
    ## final  value 293.921369 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1341.807248 
    ## iter  10 value 429.430646
    ## iter  20 value 248.332683
    ## iter  30 value 179.889810
    ## iter  40 value 159.912200
    ## iter  50 value 155.141832
    ## iter  60 value 150.713331
    ## iter  70 value 147.924598
    ## iter  80 value 146.874189
    ## iter  90 value 146.586026
    ## iter 100 value 146.386855
    ## final  value 146.386855 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1250.998820 
    ## iter  10 value 1028.141801
    ## iter  20 value 956.460783
    ## iter  30 value 919.848134
    ## iter  40 value 907.190620
    ## iter  50 value 877.644730
    ## iter  60 value 861.391565
    ## iter  70 value 853.319893
    ## iter  80 value 846.777761
    ## iter  90 value 841.283991
    ## iter 100 value 832.909083
    ## final  value 832.909083 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1278.930421 
    ## iter  10 value 703.813190
    ## iter  20 value 563.885787
    ## iter  30 value 508.651707
    ## iter  40 value 444.398101
    ## iter  50 value 390.730956
    ## iter  60 value 361.842111
    ## iter  70 value 341.463591
    ## iter  80 value 335.400483
    ## iter  90 value 322.099584
    ## iter 100 value 307.548857
    ## final  value 307.548857 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1375.369247 
    ## iter  10 value 355.856735
    ## iter  20 value 83.610670
    ## iter  30 value 25.806923
    ## iter  40 value 4.067471
    ## iter  50 value 1.833343
    ## iter  60 value 1.692069
    ## iter  70 value 1.557506
    ## iter  80 value 1.466717
    ## iter  90 value 1.379969
    ## iter 100 value 1.319491
    ## final  value 1.319491 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1315.354558 
    ## iter  10 value 1062.627530
    ## iter  20 value 1002.778482
    ## iter  30 value 970.796750
    ## iter  40 value 952.643069
    ## iter  50 value 934.032203
    ## iter  60 value 923.920533
    ## iter  70 value 919.259418
    ## iter  80 value 914.266998
    ## iter  90 value 912.465753
    ## iter 100 value 911.797519
    ## final  value 911.797519 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1229.635973 
    ## iter  10 value 600.426716
    ## iter  20 value 373.609285
    ## iter  30 value 287.868112
    ## iter  40 value 260.972892
    ## iter  50 value 240.299473
    ## iter  60 value 205.755779
    ## iter  70 value 192.794635
    ## iter  80 value 173.094427
    ## iter  90 value 139.001204
    ## iter 100 value 121.053656
    ## final  value 121.053656 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1527.395403 
    ## iter  10 value 1107.843271
    ## iter  20 value 1025.078804
    ## iter  30 value 924.412723
    ## iter  40 value 656.206097
    ## iter  50 value 593.927256
    ## iter  60 value 555.225531
    ## iter  70 value 527.123766
    ## iter  80 value 513.511128
    ## iter  90 value 488.124395
    ## iter 100 value 401.240624
    ## final  value 401.240624 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1213.851334 
    ## iter  10 value 906.088951
    ## iter  20 value 831.579916
    ## iter  30 value 797.802758
    ## iter  40 value 778.674522
    ## iter  50 value 771.353626
    ## iter  60 value 768.374401
    ## iter  70 value 765.558953
    ## iter  80 value 749.570191
    ## iter  90 value 739.464263
    ## iter 100 value 737.029030
    ## final  value 737.029030 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1231.685395 
    ## iter  10 value 680.753371
    ## iter  20 value 469.041110
    ## iter  30 value 375.492395
    ## iter  40 value 349.430033
    ## iter  50 value 338.114252
    ## iter  60 value 330.418680
    ## iter  70 value 322.290361
    ## iter  80 value 316.119796
    ## iter  90 value 311.592770
    ## iter 100 value 307.933944
    ## final  value 307.933944 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1329.971881 
    ## iter  10 value 694.831726
    ## iter  20 value 344.294038
    ## iter  30 value 204.364101
    ## iter  40 value 168.128551
    ## iter  50 value 160.020625
    ## iter  60 value 156.943734
    ## iter  70 value 154.396200
    ## iter  80 value 153.324778
    ## iter  90 value 152.708295
    ## iter 100 value 152.567888
    ## final  value 152.567888 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1267.037019 
    ## iter  10 value 884.662713
    ## iter  20 value 800.447750
    ## iter  30 value 762.984848
    ## iter  40 value 726.701968
    ## iter  50 value 697.625332
    ## iter  60 value 692.531541
    ## iter  70 value 691.352306
    ## iter  80 value 690.103867
    ## iter  90 value 689.896030
    ## iter 100 value 688.554588
    ## final  value 688.554588 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1217.924741 
    ## iter  10 value 623.817211
    ## iter  20 value 395.516499
    ## iter  30 value 303.271071
    ## iter  40 value 268.932567
    ## iter  50 value 241.060748
    ## iter  60 value 213.970121
    ## iter  70 value 201.842892
    ## iter  80 value 194.335831
    ## iter  90 value 189.197340
    ## iter 100 value 185.203137
    ## final  value 185.203137 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1215.223823 
    ## iter  10 value 401.755065
    ## iter  20 value 117.980328
    ## iter  30 value 49.365952
    ## iter  40 value 15.530029
    ## iter  50 value 6.800260
    ## iter  60 value 4.226320
    ## iter  70 value 3.492125
    ## iter  80 value 3.323281
    ## iter  90 value 3.173475
    ## iter 100 value 3.056852
    ## final  value 3.056852 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1262.349225 
    ## iter  10 value 933.140522
    ## iter  20 value 828.917845
    ## iter  30 value 763.399510
    ## iter  40 value 719.076231
    ## iter  50 value 669.883168
    ## iter  60 value 640.225741
    ## iter  70 value 626.481898
    ## iter  80 value 626.097165
    ## iter  90 value 626.085530
    ## iter 100 value 626.083172
    ## final  value 626.083172 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1203.621774 
    ## iter  10 value 566.453897
    ## iter  20 value 373.388223
    ## iter  30 value 329.682720
    ## iter  40 value 285.022052
    ## iter  50 value 257.474456
    ## iter  60 value 234.609126
    ## iter  70 value 224.135761
    ## iter  80 value 220.364238
    ## iter  90 value 216.931851
    ## iter 100 value 216.513279
    ## final  value 216.513279 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1229.013196 
    ## iter  10 value 628.566074
    ## iter  20 value 226.303384
    ## iter  30 value 97.212201
    ## iter  40 value 60.078531
    ## iter  50 value 41.756175
    ## iter  60 value 27.865216
    ## iter  70 value 20.573396
    ## iter  80 value 15.530776
    ## iter  90 value 13.127932
    ## iter 100 value 9.895187
    ## final  value 9.895187 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1253.892074 
    ## iter  10 value 1095.492309
    ## iter  20 value 905.058372
    ## iter  30 value 856.770813
    ## iter  40 value 826.567324
    ## iter  50 value 813.802115
    ## iter  60 value 802.102614
    ## iter  70 value 794.096523
    ## iter  80 value 789.191316
    ## iter  90 value 787.852513
    ## iter 100 value 786.820696
    ## final  value 786.820696 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1223.181947 
    ## iter  10 value 780.817080
    ## iter  20 value 556.658716
    ## iter  30 value 460.699248
    ## iter  40 value 397.238154
    ## iter  50 value 372.384370
    ## iter  60 value 352.267426
    ## iter  70 value 342.870472
    ## iter  80 value 336.290105
    ## iter  90 value 327.397346
    ## iter 100 value 320.200832
    ## final  value 320.200832 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1457.967333 
    ## iter  10 value 1048.879346
    ## iter  20 value 774.476829
    ## iter  30 value 651.311492
    ## iter  40 value 595.561030
    ## iter  50 value 465.121853
    ## iter  60 value 321.657515
    ## iter  70 value 211.862720
    ## iter  80 value 189.727576
    ## iter  90 value 173.715980
    ## iter 100 value 168.589670
    ## final  value 168.589670 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1220.307366 
    ## iter  10 value 876.802459
    ## iter  20 value 782.999545
    ## iter  30 value 743.303302
    ## iter  40 value 705.693952
    ## iter  50 value 681.017071
    ## iter  60 value 669.500773
    ## iter  70 value 658.165600
    ## iter  80 value 623.349036
    ## iter  90 value 606.150166
    ## iter 100 value 573.309224
    ## final  value 573.309224 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1297.701347 
    ## iter  10 value 683.559978
    ## iter  20 value 478.466887
    ## iter  30 value 417.929975
    ## iter  40 value 356.665095
    ## iter  50 value 305.900990
    ## iter  60 value 284.349798
    ## iter  70 value 269.446874
    ## iter  80 value 260.869434
    ## iter  90 value 250.804009
    ## iter 100 value 232.766984
    ## final  value 232.766984 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1270.242820 
    ## iter  10 value 415.070301
    ## iter  20 value 191.167502
    ## iter  30 value 117.066532
    ## iter  40 value 82.501643
    ## iter  50 value 53.116420
    ## iter  60 value 36.685115
    ## iter  70 value 27.939111
    ## iter  80 value 26.127025
    ## iter  90 value 24.966647
    ## iter 100 value 23.860198
    ## final  value 23.860198 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1265.496299 
    ## iter  10 value 1087.793018
    ## iter  20 value 915.656108
    ## iter  30 value 804.963803
    ## iter  40 value 751.969047
    ## iter  50 value 698.738104
    ## iter  60 value 651.443621
    ## iter  70 value 607.317942
    ## iter  80 value 575.258404
    ## iter  90 value 541.328772
    ## iter 100 value 522.942371
    ## final  value 522.942371 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1382.019621 
    ## iter  10 value 939.120064
    ## iter  20 value 817.524303
    ## iter  30 value 700.485399
    ## iter  40 value 585.118984
    ## iter  50 value 523.691792
    ## iter  60 value 478.314728
    ## iter  70 value 437.308905
    ## iter  80 value 394.618895
    ## iter  90 value 375.472424
    ## iter 100 value 370.600676
    ## final  value 370.600676 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1241.776954 
    ## iter  10 value 542.314942
    ## iter  20 value 136.605977
    ## iter  30 value 75.472300
    ## iter  40 value 48.985261
    ## iter  50 value 29.193811
    ## iter  60 value 20.615071
    ## iter  70 value 15.437392
    ## iter  80 value 14.941974
    ## iter  90 value 12.927914
    ## iter 100 value 12.011960
    ## final  value 12.011960 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1278.414039 
    ## iter  10 value 1001.595651
    ## iter  20 value 926.687586
    ## iter  30 value 881.992654
    ## iter  40 value 825.182303
    ## iter  50 value 801.416119
    ## iter  60 value 791.893933
    ## iter  70 value 789.109885
    ## iter  80 value 788.674468
    ## iter  90 value 788.619633
    ## iter 100 value 788.611970
    ## final  value 788.611970 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1278.087379 
    ## iter  10 value 726.825875
    ## iter  20 value 525.047999
    ## iter  30 value 407.003955
    ## iter  40 value 359.236221
    ## iter  50 value 333.248630
    ## iter  60 value 314.622593
    ## iter  70 value 303.750809
    ## iter  80 value 299.412499
    ## iter  90 value 298.013808
    ## iter 100 value 297.429297
    ## final  value 297.429297 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1356.731203 
    ## iter  10 value 440.705235
    ## iter  20 value 277.140909
    ## iter  30 value 200.551803
    ## iter  40 value 170.194466
    ## iter  50 value 158.979197
    ## iter  60 value 152.470439
    ## iter  70 value 149.461963
    ## iter  80 value 147.434206
    ## iter  90 value 146.745754
    ## iter 100 value 146.667405
    ## final  value 146.667405 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1303.001776 
    ## iter  10 value 1141.604544
    ## iter  20 value 1122.488128
    ## iter  30 value 1104.293346
    ## iter  40 value 995.647249
    ## iter  50 value 928.042638
    ## iter  60 value 884.889246
    ## iter  70 value 863.154584
    ## iter  80 value 845.010047
    ## iter  90 value 829.305970
    ## iter 100 value 800.063816
    ## final  value 800.063816 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1255.375395 
    ## iter  10 value 558.017761
    ## iter  20 value 324.973597
    ## iter  30 value 246.363751
    ## iter  40 value 209.521964
    ## iter  50 value 189.906898
    ## iter  60 value 166.490269
    ## iter  70 value 156.184697
    ## iter  80 value 151.171607
    ## iter  90 value 149.595932
    ## iter 100 value 148.334596
    ## final  value 148.334596 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1389.061225 
    ## iter  10 value 575.544506
    ## iter  20 value 326.175880
    ## iter  30 value 242.572030
    ## iter  40 value 184.146808
    ## iter  50 value 146.685438
    ## iter  60 value 105.751633
    ## iter  70 value 76.867715
    ## iter  80 value 59.784244
    ## iter  90 value 47.179489
    ## iter 100 value 42.137853
    ## final  value 42.137853 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1242.131329 
    ## iter  10 value 845.392521
    ## iter  20 value 764.084611
    ## iter  30 value 681.344877
    ## iter  40 value 638.958275
    ## iter  50 value 570.987397
    ## iter  60 value 531.057207
    ## iter  70 value 511.846608
    ## iter  80 value 500.624240
    ## iter  90 value 477.650053
    ## iter 100 value 469.053955
    ## final  value 469.053955 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1328.444521 
    ## iter  10 value 732.094183
    ## iter  20 value 520.188297
    ## iter  30 value 453.430018
    ## iter  40 value 388.745493
    ## iter  50 value 350.462955
    ## iter  60 value 317.456305
    ## iter  70 value 290.777931
    ## iter  80 value 283.415305
    ## iter  90 value 279.059013
    ## iter 100 value 274.041401
    ## final  value 274.041401 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1305.402703 
    ## iter  10 value 347.216534
    ## iter  20 value 96.206118
    ## iter  30 value 40.286329
    ## iter  40 value 21.310631
    ## iter  50 value 16.800490
    ## iter  60 value 12.623102
    ## iter  70 value 10.119725
    ## iter  80 value 6.425021
    ## iter  90 value 3.365776
    ## iter 100 value 2.895440
    ## final  value 2.895440 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1211.656159 
    ## iter  10 value 1129.703133
    ## iter  20 value 983.146143
    ## iter  30 value 908.833241
    ## iter  40 value 873.365709
    ## iter  50 value 835.735007
    ## iter  60 value 808.656976
    ## iter  70 value 793.320670
    ## iter  80 value 784.729115
    ## iter  90 value 780.574397
    ## iter 100 value 778.542463
    ## final  value 778.542463 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1327.680286 
    ## iter  10 value 757.130036
    ## iter  20 value 494.805768
    ## iter  30 value 404.632393
    ## iter  40 value 344.357072
    ## iter  50 value 325.519353
    ## iter  60 value 317.303849
    ## iter  70 value 312.237931
    ## iter  80 value 309.528878
    ## iter  90 value 304.942528
    ## iter 100 value 303.657229
    ## final  value 303.657229 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1453.878708 
    ## iter  10 value 991.854983
    ## iter  20 value 596.523376
    ## iter  30 value 469.833283
    ## iter  40 value 395.287853
    ## iter  50 value 291.304952
    ## iter  60 value 227.429878
    ## iter  70 value 193.713877
    ## iter  80 value 178.791675
    ## iter  90 value 170.856981
    ## iter 100 value 161.737280
    ## final  value 161.737280 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1280.178270 
    ## iter  10 value 955.031434
    ## iter  20 value 825.379895
    ## iter  30 value 748.068315
    ## iter  40 value 696.131466
    ## iter  50 value 663.351675
    ## iter  60 value 639.430772
    ## iter  70 value 624.601675
    ## iter  80 value 619.472899
    ## iter  90 value 616.094133
    ## iter 100 value 609.709299
    ## final  value 609.709299 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1307.394160 
    ## iter  10 value 689.258266
    ## iter  20 value 360.001716
    ## iter  30 value 311.189238
    ## iter  40 value 287.241429
    ## iter  50 value 274.370357
    ## iter  60 value 255.576035
    ## iter  70 value 249.188648
    ## iter  80 value 247.672329
    ## iter  90 value 247.438828
    ## iter 100 value 245.656425
    ## final  value 245.656425 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1298.116991 
    ## iter  10 value 388.011541
    ## iter  20 value 158.731535
    ## iter  30 value 65.290593
    ## iter  40 value 42.333993
    ## iter  50 value 37.125575
    ## iter  60 value 29.752349
    ## iter  70 value 26.308767
    ## iter  80 value 25.823587
    ## iter  90 value 25.511876
    ## iter 100 value 23.254879
    ## final  value 23.254879 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1231.478125 
    ## iter  10 value 1044.288361
    ## iter  20 value 827.243722
    ## iter  30 value 767.940042
    ## iter  40 value 720.727620
    ## iter  50 value 686.330533
    ## iter  60 value 662.548874
    ## iter  70 value 646.944524
    ## iter  80 value 640.409093
    ## iter  90 value 638.697902
    ## iter 100 value 636.876291
    ## final  value 636.876291 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1257.528370 
    ## iter  10 value 700.823222
    ## iter  20 value 449.274629
    ## iter  30 value 398.927158
    ## iter  40 value 365.253813
    ## iter  50 value 328.842202
    ## iter  60 value 321.217691
    ## iter  70 value 315.581905
    ## iter  80 value 313.249528
    ## iter  90 value 310.763383
    ## iter 100 value 309.454378
    ## final  value 309.454378 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1304.617466 
    ## iter  10 value 475.518865
    ## iter  20 value 132.892005
    ## iter  30 value 88.499684
    ## iter  40 value 57.379689
    ## iter  50 value 41.312938
    ## iter  60 value 34.678488
    ## iter  70 value 32.363095
    ## iter  80 value 30.461095
    ## iter  90 value 30.224975
    ## iter 100 value 30.220777
    ## final  value 30.220777 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1266.471737 
    ## iter  10 value 934.618466
    ## iter  20 value 866.521511
    ## iter  30 value 824.980677
    ## iter  40 value 778.405217
    ## iter  50 value 757.239247
    ## iter  60 value 752.180363
    ## iter  70 value 749.818571
    ## iter  80 value 748.866630
    ## iter  90 value 748.735212
    ## iter 100 value 748.730893
    ## final  value 748.730893 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1395.663555 
    ## iter  10 value 777.171944
    ## iter  20 value 666.251794
    ## iter  30 value 592.782165
    ## iter  40 value 518.049470
    ## iter  50 value 450.179974
    ## iter  60 value 413.009994
    ## iter  70 value 390.173008
    ## iter  80 value 363.336181
    ## iter  90 value 320.596032
    ## iter 100 value 309.207023
    ## final  value 309.207023 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1243.730424 
    ## iter  10 value 522.994189
    ## iter  20 value 321.864770
    ## iter  30 value 207.577808
    ## iter  40 value 174.411192
    ## iter  50 value 161.331021
    ## iter  60 value 156.922814
    ## iter  70 value 155.459906
    ## iter  80 value 153.965457
    ## iter  90 value 148.675481
    ## iter 100 value 144.613423
    ## final  value 144.613423 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1257.786567 
    ## iter  10 value 930.544688
    ## iter  20 value 805.663813
    ## iter  30 value 769.649624
    ## iter  40 value 720.494789
    ## iter  50 value 690.998427
    ## iter  60 value 668.351499
    ## iter  70 value 666.557488
    ## iter  80 value 663.908671
    ## iter  90 value 663.100310
    ## iter 100 value 661.742136
    ## final  value 661.742136 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1299.723752 
    ## iter  10 value 1060.634255
    ## iter  20 value 965.607996
    ## iter  30 value 878.522367
    ## iter  40 value 734.461316
    ## iter  50 value 677.462154
    ## iter  60 value 661.472085
    ## iter  70 value 638.536322
    ## iter  80 value 624.726760
    ## iter  90 value 612.214057
    ## iter 100 value 607.806922
    ## final  value 607.806922 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1304.674730 
    ## iter  10 value 351.817048
    ## iter  20 value 137.852933
    ## iter  30 value 47.524354
    ## iter  40 value 17.548594
    ## iter  50 value 8.295119
    ## iter  60 value 7.379145
    ## iter  70 value 6.923015
    ## iter  80 value 6.564314
    ## iter  90 value 6.360232
    ## iter 100 value 6.150092
    ## final  value 6.150092 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1366.687977 
    ## iter  10 value 1104.154702
    ## iter  20 value 1095.367499
    ## iter  30 value 1095.294989
    ## final  value 1095.293335 
    ## converged
    ## # weights:  255
    ## initial  value 1381.189322 
    ## iter  10 value 833.941345
    ## iter  20 value 547.911747
    ## iter  30 value 471.812811
    ## iter  40 value 423.855922
    ## iter  50 value 398.555381
    ## iter  60 value 392.196041
    ## iter  70 value 384.388318
    ## iter  80 value 364.141928
    ## iter  90 value 339.941032
    ## iter 100 value 331.218479
    ## final  value 331.218479 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1253.485666 
    ## iter  10 value 361.826721
    ## iter  20 value 124.122307
    ## iter  30 value 47.211368
    ## iter  40 value 21.029031
    ## iter  50 value 13.682974
    ## iter  60 value 11.327873
    ## iter  70 value 8.206055
    ## iter  80 value 6.606091
    ## iter  90 value 6.460030
    ## iter 100 value 6.452648
    ## final  value 6.452648 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1188.186280 
    ## iter  10 value 1088.201659
    ## iter  20 value 966.277602
    ## iter  30 value 887.814352
    ## iter  40 value 851.313940
    ## iter  50 value 814.409608
    ## iter  60 value 787.782492
    ## iter  70 value 782.491422
    ## iter  80 value 781.751691
    ## iter  90 value 780.413684
    ## iter 100 value 779.116282
    ## final  value 779.116282 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1209.027870 
    ## iter  10 value 732.836755
    ## iter  20 value 511.221418
    ## iter  30 value 392.513016
    ## iter  40 value 342.309161
    ## iter  50 value 323.747957
    ## iter  60 value 310.176748
    ## iter  70 value 302.703960
    ## iter  80 value 299.976384
    ## iter  90 value 297.793354
    ## iter 100 value 296.527203
    ## final  value 296.527203 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1347.282342 
    ## iter  10 value 476.407772
    ## iter  20 value 285.992312
    ## iter  30 value 235.866026
    ## iter  40 value 196.042886
    ## iter  50 value 169.444026
    ## iter  60 value 158.859767
    ## iter  70 value 152.313093
    ## iter  80 value 150.487235
    ## iter  90 value 149.958724
    ## iter 100 value 149.641568
    ## final  value 149.641568 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1226.863282 
    ## iter  10 value 886.549732
    ## iter  20 value 790.438677
    ## iter  30 value 706.760212
    ## iter  40 value 646.005464
    ## iter  50 value 594.420760
    ## iter  60 value 576.922735
    ## iter  70 value 573.384437
    ## iter  80 value 571.463973
    ## iter  90 value 571.244639
    ## iter 100 value 569.888178
    ## final  value 569.888178 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1339.604965 
    ## iter  10 value 664.381063
    ## iter  20 value 419.467881
    ## iter  30 value 372.289588
    ## iter  40 value 307.593108
    ## iter  50 value 269.158146
    ## iter  60 value 237.336833
    ## iter  70 value 229.075955
    ## iter  80 value 223.128402
    ## iter  90 value 214.219890
    ## iter 100 value 210.253758
    ## final  value 210.253758 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1327.008962 
    ## iter  10 value 622.951807
    ## iter  20 value 317.475249
    ## iter  30 value 249.564478
    ## iter  40 value 166.685543
    ## iter  50 value 129.460711
    ## iter  60 value 104.154680
    ## iter  70 value 94.229711
    ## iter  80 value 91.134419
    ## iter  90 value 89.051267
    ## iter 100 value 84.605491
    ## final  value 84.605491 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1285.467477 
    ## iter  10 value 1109.125736
    ## iter  20 value 1097.315092
    ## final  value 1097.310504 
    ## converged
    ## # weights:  255
    ## initial  value 1279.410158 
    ## iter  10 value 836.506752
    ## iter  20 value 509.433524
    ## iter  30 value 446.765643
    ## iter  40 value 397.169357
    ## iter  50 value 354.360475
    ## iter  60 value 331.535146
    ## iter  70 value 303.485885
    ## iter  80 value 206.566020
    ## iter  90 value 154.918350
    ## iter 100 value 135.213114
    ## final  value 135.213114 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1380.065610 
    ## iter  10 value 534.602127
    ## iter  20 value 175.332407
    ## iter  30 value 66.779834
    ## iter  40 value 53.695560
    ## iter  50 value 43.781487
    ## iter  60 value 20.929649
    ## iter  70 value 12.876557
    ## iter  80 value 8.640017
    ## iter  90 value 8.085829
    ## iter 100 value 7.283368
    ## final  value 7.283368 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1250.524459 
    ## iter  10 value 1035.304655
    ## iter  20 value 872.365919
    ## iter  30 value 810.360538
    ## iter  40 value 784.672590
    ## iter  50 value 766.324178
    ## iter  60 value 753.505098
    ## iter  70 value 749.459256
    ## iter  80 value 748.792561
    ## iter  90 value 748.590408
    ## iter 100 value 748.556973
    ## final  value 748.556973 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1258.388449 
    ## iter  10 value 637.048874
    ## iter  20 value 449.614804
    ## iter  30 value 365.638868
    ## iter  40 value 341.156070
    ## iter  50 value 326.347474
    ## iter  60 value 314.244016
    ## iter  70 value 302.610955
    ## iter  80 value 292.030554
    ## iter  90 value 289.008906
    ## iter 100 value 287.775360
    ## final  value 287.775360 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1244.304475 
    ## iter  10 value 654.248426
    ## iter  20 value 423.545622
    ## iter  30 value 337.529141
    ## iter  40 value 293.617678
    ## iter  50 value 237.995464
    ## iter  60 value 212.607894
    ## iter  70 value 191.472673
    ## iter  80 value 170.152759
    ## iter  90 value 161.905757
    ## iter 100 value 153.136857
    ## final  value 153.136857 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1290.970948 
    ## iter  10 value 1005.271473
    ## iter  20 value 860.392626
    ## iter  30 value 785.527241
    ## iter  40 value 739.471770
    ## iter  50 value 693.697480
    ## iter  60 value 678.312822
    ## iter  70 value 671.362782
    ## iter  80 value 664.118808
    ## iter  90 value 640.969903
    ## iter 100 value 631.971609
    ## final  value 631.971609 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1224.367056 
    ## iter  10 value 583.711561
    ## iter  20 value 347.433465
    ## iter  30 value 194.789995
    ## iter  40 value 162.014640
    ## iter  50 value 142.662841
    ## iter  60 value 128.042784
    ## iter  70 value 118.327095
    ## iter  80 value 112.741895
    ## iter  90 value 109.297325
    ## iter 100 value 106.403625
    ## final  value 106.403625 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1341.931809 
    ## iter  10 value 947.653551
    ## iter  20 value 764.563403
    ## iter  30 value 618.951043
    ## iter  40 value 558.581617
    ## iter  50 value 509.583245
    ## iter  60 value 470.789404
    ## iter  70 value 418.996635
    ## iter  80 value 351.607844
    ## iter  90 value 258.605828
    ## iter 100 value 240.295207
    ## final  value 240.295207 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1242.000321 
    ## iter  10 value 1070.045602
    ## iter  20 value 944.696271
    ## iter  30 value 861.465556
    ## iter  40 value 829.019955
    ## iter  50 value 788.146232
    ## iter  60 value 743.319508
    ## iter  70 value 738.059709
    ## iter  80 value 731.671427
    ## iter  90 value 716.259252
    ## iter 100 value 708.912399
    ## final  value 708.912399 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1256.487855 
    ## iter  10 value 688.131982
    ## iter  20 value 477.658284
    ## iter  30 value 431.937844
    ## iter  40 value 386.100293
    ## iter  50 value 359.456729
    ## iter  60 value 339.389801
    ## iter  70 value 336.709083
    ## iter  80 value 333.516916
    ## iter  90 value 330.088024
    ## iter 100 value 321.628349
    ## final  value 321.628349 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1369.083676 
    ## iter  10 value 889.889286
    ## iter  20 value 505.062586
    ## iter  30 value 340.793780
    ## iter  40 value 180.238625
    ## iter  50 value 83.113070
    ## iter  60 value 70.095925
    ## iter  70 value 61.648983
    ## iter  80 value 58.338382
    ## iter  90 value 54.224575
    ## iter 100 value 45.007146
    ## final  value 45.007146 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1264.059384 
    ## iter  10 value 938.461604
    ## iter  20 value 869.399129
    ## iter  30 value 829.103085
    ## iter  40 value 806.846848
    ## iter  50 value 793.401305
    ## iter  60 value 780.087386
    ## iter  70 value 767.809083
    ## iter  80 value 762.049801
    ## iter  90 value 759.597402
    ## iter 100 value 759.233628
    ## final  value 759.233628 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1298.097216 
    ## iter  10 value 834.443535
    ## iter  20 value 584.053839
    ## iter  30 value 482.141613
    ## iter  40 value 437.720857
    ## iter  50 value 390.448826
    ## iter  60 value 348.663668
    ## iter  70 value 325.985548
    ## iter  80 value 316.010382
    ## iter  90 value 311.448179
    ## iter 100 value 309.296688
    ## final  value 309.296688 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1327.865250 
    ## iter  10 value 540.565195
    ## iter  20 value 306.674452
    ## iter  30 value 203.899513
    ## iter  40 value 169.167108
    ## iter  50 value 159.241339
    ## iter  60 value 156.519194
    ## iter  70 value 154.335688
    ## iter  80 value 149.068821
    ## iter  90 value 147.463997
    ## iter 100 value 147.028723
    ## final  value 147.028723 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1309.694030 
    ## iter  10 value 1026.275577
    ## iter  20 value 893.965268
    ## iter  30 value 850.675813
    ## iter  40 value 805.107562
    ## iter  50 value 750.130811
    ## iter  60 value 712.239400
    ## iter  70 value 704.185502
    ## iter  80 value 697.444656
    ## iter  90 value 695.040181
    ## iter 100 value 693.914823
    ## final  value 693.914823 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1346.090450 
    ## iter  10 value 1015.401217
    ## iter  20 value 668.825895
    ## iter  30 value 490.586634
    ## iter  40 value 367.440858
    ## iter  50 value 270.041408
    ## iter  60 value 197.841165
    ## iter  70 value 146.091572
    ## iter  80 value 130.116833
    ## iter  90 value 113.224183
    ## iter 100 value 101.056705
    ## final  value 101.056705 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1293.401147 
    ## iter  10 value 440.392706
    ## iter  20 value 146.808549
    ## iter  30 value 85.716219
    ## iter  40 value 63.473867
    ## iter  50 value 50.729049
    ## iter  60 value 45.242904
    ## iter  70 value 42.134360
    ## iter  80 value 36.533590
    ## iter  90 value 32.688252
    ## iter 100 value 29.170091
    ## final  value 29.170091 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1232.969393 
    ## iter  10 value 990.799458
    ## iter  20 value 788.799137
    ## iter  30 value 748.896410
    ## iter  40 value 726.623882
    ## iter  50 value 711.356972
    ## iter  60 value 701.922281
    ## iter  70 value 696.345735
    ## iter  80 value 694.078921
    ## iter  90 value 692.785777
    ## iter 100 value 691.097730
    ## final  value 691.097730 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1268.621484 
    ## iter  10 value 988.852565
    ## iter  20 value 750.985879
    ## iter  30 value 613.665534
    ## iter  40 value 508.105380
    ## iter  50 value 435.224522
    ## iter  60 value 329.067721
    ## iter  70 value 270.565028
    ## iter  80 value 228.514052
    ## iter  90 value 196.061915
    ## iter 100 value 170.201402
    ## final  value 170.201402 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1362.450579 
    ## iter  10 value 453.249257
    ## iter  20 value 184.539772
    ## iter  30 value 92.613730
    ## iter  40 value 58.413036
    ## iter  50 value 40.632186
    ## iter  60 value 27.511316
    ## iter  70 value 21.949400
    ## iter  80 value 19.914756
    ## iter  90 value 17.095275
    ## iter 100 value 14.977966
    ## final  value 14.977966 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1245.588461 
    ## iter  10 value 1024.088111
    ## iter  20 value 921.962546
    ## iter  30 value 861.862271
    ## iter  40 value 831.365304
    ## iter  50 value 813.291236
    ## iter  60 value 800.052158
    ## iter  70 value 791.708515
    ## iter  80 value 789.470409
    ## iter  90 value 788.721089
    ## iter 100 value 788.668261
    ## final  value 788.668261 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1316.735165 
    ## iter  10 value 625.841679
    ## iter  20 value 495.751389
    ## iter  30 value 402.189367
    ## iter  40 value 354.005961
    ## iter  50 value 339.530896
    ## iter  60 value 330.830194
    ## iter  70 value 322.031311
    ## iter  80 value 313.717020
    ## iter  90 value 302.635448
    ## iter 100 value 298.376144
    ## final  value 298.376144 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1355.382894 
    ## iter  10 value 532.548527
    ## iter  20 value 300.650606
    ## iter  30 value 229.186664
    ## iter  40 value 174.516310
    ## iter  50 value 161.353635
    ## iter  60 value 151.035451
    ## iter  70 value 143.167628
    ## iter  80 value 141.020364
    ## iter  90 value 139.323190
    ## iter 100 value 138.763218
    ## final  value 138.763218 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1247.081867 
    ## iter  10 value 886.562068
    ## iter  20 value 797.858215
    ## iter  30 value 765.883115
    ## iter  40 value 740.740825
    ## iter  50 value 682.569512
    ## iter  60 value 645.484406
    ## iter  70 value 600.484474
    ## iter  80 value 571.664793
    ## iter  90 value 560.594144
    ## iter 100 value 555.568999
    ## final  value 555.568999 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1319.246750 
    ## iter  10 value 695.557943
    ## iter  20 value 472.909699
    ## iter  30 value 342.755083
    ## iter  40 value 309.475160
    ## iter  50 value 282.364748
    ## iter  60 value 262.692110
    ## iter  70 value 258.916635
    ## iter  80 value 257.199602
    ## iter  90 value 251.380302
    ## iter 100 value 248.269834
    ## final  value 248.269834 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1452.527263 
    ## iter  10 value 1079.479433
    ## iter  20 value 716.118147
    ## iter  30 value 405.470750
    ## iter  40 value 147.907430
    ## iter  50 value 100.498855
    ## iter  60 value 71.859710
    ## iter  70 value 55.051135
    ## iter  80 value 45.430771
    ## iter  90 value 40.166730
    ## iter 100 value 34.787760
    ## final  value 34.787760 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1294.022264 
    ## iter  10 value 840.884067
    ## iter  20 value 797.850883
    ## iter  30 value 775.063577
    ## iter  40 value 749.905536
    ## iter  50 value 738.053635
    ## iter  60 value 734.807793
    ## iter  70 value 732.737440
    ## iter  80 value 731.987876
    ## iter  90 value 731.784132
    ## iter 100 value 731.705336
    ## final  value 731.705336 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1274.065605 
    ## iter  10 value 711.239154
    ## iter  20 value 359.611145
    ## iter  30 value 255.611974
    ## iter  40 value 195.563799
    ## iter  50 value 158.124150
    ## iter  60 value 132.575525
    ## iter  70 value 112.166588
    ## iter  80 value 102.046016
    ## iter  90 value 97.004081
    ## iter 100 value 93.732670
    ## final  value 93.732670 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1266.623762 
    ## iter  10 value 503.288438
    ## iter  20 value 147.956498
    ## iter  30 value 78.110852
    ## iter  40 value 34.422719
    ## iter  50 value 20.025670
    ## iter  60 value 10.409810
    ## iter  70 value 7.867690
    ## iter  80 value 6.559357
    ## iter  90 value 6.172187
    ## iter 100 value 5.872606
    ## final  value 5.872606 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1267.158357 
    ## iter  10 value 1107.630729
    ## iter  20 value 1078.938823
    ## iter  30 value 904.068938
    ## iter  40 value 845.694444
    ## iter  50 value 818.387616
    ## iter  60 value 793.219339
    ## iter  70 value 780.341392
    ## iter  80 value 770.370955
    ## iter  90 value 769.043874
    ## iter 100 value 768.900090
    ## final  value 768.900090 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1313.538552 
    ## iter  10 value 957.350692
    ## iter  20 value 839.429655
    ## iter  30 value 740.909957
    ## iter  40 value 620.498207
    ## iter  50 value 568.001783
    ## iter  60 value 524.621448
    ## iter  70 value 484.218971
    ## iter  80 value 409.538151
    ## iter  90 value 328.897243
    ## iter 100 value 309.862083
    ## final  value 309.862083 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1268.264335 
    ## iter  10 value 571.304355
    ## iter  20 value 375.796231
    ## iter  30 value 268.814882
    ## iter  40 value 238.255682
    ## iter  50 value 214.602827
    ## iter  60 value 187.646898
    ## iter  70 value 175.613966
    ## iter  80 value 166.903138
    ## iter  90 value 160.729556
    ## iter 100 value 156.781861
    ## final  value 156.781861 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1209.865036 
    ## iter  10 value 942.705961
    ## iter  20 value 826.412839
    ## iter  30 value 778.709424
    ## iter  40 value 753.733056
    ## iter  50 value 748.720414
    ## iter  60 value 745.684532
    ## iter  70 value 744.127441
    ## iter  80 value 742.647169
    ## iter  90 value 737.855797
    ## iter 100 value 736.320618
    ## final  value 736.320618 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1308.603197 
    ## iter  10 value 758.711639
    ## iter  20 value 399.939368
    ## iter  30 value 301.271569
    ## iter  40 value 246.588489
    ## iter  50 value 220.812502
    ## iter  60 value 205.701665
    ## iter  70 value 199.544085
    ## iter  80 value 196.785146
    ## iter  90 value 195.206257
    ## iter 100 value 194.362601
    ## final  value 194.362601 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1366.515692 
    ## iter  10 value 511.372531
    ## iter  20 value 136.982690
    ## iter  30 value 66.195600
    ## iter  40 value 44.994860
    ## iter  50 value 20.860107
    ## iter  60 value 10.826392
    ## iter  70 value 8.277685
    ## iter  80 value 7.389829
    ## iter  90 value 6.687854
    ## iter 100 value 5.144730
    ## final  value 5.144730 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1257.502651 
    ## iter  10 value 955.166725
    ## iter  20 value 830.825164
    ## iter  30 value 738.972008
    ## iter  40 value 667.950055
    ## iter  50 value 620.291574
    ## iter  60 value 576.587476
    ## iter  70 value 533.677030
    ## iter  80 value 512.287921
    ## iter  90 value 506.259006
    ## iter 100 value 499.853770
    ## final  value 499.853770 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1325.321307 
    ## iter  10 value 827.005123
    ## iter  20 value 598.921884
    ## iter  30 value 405.460553
    ## iter  40 value 297.532059
    ## iter  50 value 252.921044
    ## iter  60 value 220.453594
    ## iter  70 value 171.560998
    ## iter  80 value 144.395232
    ## iter  90 value 123.943640
    ## iter 100 value 104.716105
    ## final  value 104.716105 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1403.725673 
    ## iter  10 value 861.705208
    ## iter  20 value 425.135880
    ## iter  30 value 348.662149
    ## iter  40 value 286.144499
    ## iter  50 value 211.014462
    ## iter  60 value 176.175245
    ## iter  70 value 158.019190
    ## iter  80 value 134.147771
    ## iter  90 value 119.790546
    ## iter 100 value 103.128364
    ## final  value 103.128364 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1301.311215 
    ## iter  10 value 1040.824153
    ## iter  20 value 883.325128
    ## iter  30 value 823.710058
    ## iter  40 value 801.909833
    ## iter  50 value 789.205877
    ## iter  60 value 781.428749
    ## iter  70 value 774.215426
    ## iter  80 value 772.231756
    ## iter  90 value 769.865755
    ## iter 100 value 769.642366
    ## final  value 769.642366 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1195.064562 
    ## iter  10 value 646.581207
    ## iter  20 value 477.970177
    ## iter  30 value 408.002724
    ## iter  40 value 366.565347
    ## iter  50 value 346.487253
    ## iter  60 value 334.668224
    ## iter  70 value 320.386447
    ## iter  80 value 315.526449
    ## iter  90 value 313.220937
    ## iter 100 value 311.964527
    ## final  value 311.964527 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1356.078777 
    ## iter  10 value 543.029804
    ## iter  20 value 316.951760
    ## iter  30 value 243.840297
    ## iter  40 value 197.431530
    ## iter  50 value 174.411963
    ## iter  60 value 162.265294
    ## iter  70 value 157.978272
    ## iter  80 value 155.722853
    ## iter  90 value 154.411398
    ## iter 100 value 152.597454
    ## final  value 152.597454 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1281.329107 
    ## iter  10 value 932.802378
    ## iter  20 value 810.765721
    ## iter  30 value 762.724846
    ## iter  40 value 727.082623
    ## iter  50 value 716.062468
    ## iter  60 value 712.316622
    ## iter  70 value 711.290568
    ## iter  80 value 710.512043
    ## iter  90 value 709.045161
    ## iter 100 value 706.593355
    ## final  value 706.593355 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1307.517043 
    ## iter  10 value 603.473520
    ## iter  20 value 358.946166
    ## iter  30 value 303.801900
    ## iter  40 value 246.837297
    ## iter  50 value 225.544386
    ## iter  60 value 205.541414
    ## iter  70 value 191.323021
    ## iter  80 value 170.014772
    ## iter  90 value 156.914563
    ## iter 100 value 141.820051
    ## final  value 141.820051 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1398.393105 
    ## iter  10 value 920.279546
    ## iter  20 value 475.435636
    ## iter  30 value 183.829305
    ## iter  40 value 80.780435
    ## iter  50 value 44.334646
    ## iter  60 value 33.897605
    ## iter  70 value 30.337407
    ## iter  80 value 23.846113
    ## iter  90 value 18.276185
    ## iter 100 value 14.960309
    ## final  value 14.960309 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1175.856779 
    ## iter  10 value 886.627716
    ## iter  20 value 818.666245
    ## iter  30 value 801.289913
    ## iter  40 value 767.806407
    ## iter  50 value 711.738570
    ## iter  60 value 675.278822
    ## iter  70 value 650.134203
    ## iter  80 value 641.555076
    ## iter  90 value 640.366585
    ## iter 100 value 639.824513
    ## final  value 639.824513 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1233.216022 
    ## iter  10 value 554.391429
    ## iter  20 value 337.895965
    ## iter  30 value 265.787548
    ## iter  40 value 220.285932
    ## iter  50 value 185.580420
    ## iter  60 value 163.186350
    ## iter  70 value 150.638499
    ## iter  80 value 143.879042
    ## iter  90 value 137.557140
    ## iter 100 value 125.298100
    ## final  value 125.298100 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1252.480323 
    ## iter  10 value 509.266852
    ## iter  20 value 201.628176
    ## iter  30 value 102.459561
    ## iter  40 value 58.972194
    ## iter  50 value 50.484023
    ## iter  60 value 47.602929
    ## iter  70 value 45.198358
    ## iter  80 value 44.628860
    ## iter  90 value 44.387701
    ## iter 100 value 43.686320
    ## final  value 43.686320 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1247.319698 
    ## iter  10 value 973.606909
    ## iter  20 value 864.764267
    ## iter  30 value 834.152836
    ## iter  40 value 803.448834
    ## iter  50 value 791.108892
    ## iter  60 value 782.505942
    ## iter  70 value 777.930521
    ## iter  80 value 776.106738
    ## iter  90 value 775.811437
    ## iter 100 value 775.587271
    ## final  value 775.587271 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1343.711268 
    ## iter  10 value 997.953916
    ## iter  20 value 787.363439
    ## iter  30 value 639.170767
    ## iter  40 value 533.311930
    ## iter  50 value 442.739789
    ## iter  60 value 404.584653
    ## iter  70 value 368.271732
    ## iter  80 value 334.328897
    ## iter  90 value 316.854051
    ## iter 100 value 312.484594
    ## final  value 312.484594 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1365.011670 
    ## iter  10 value 397.795135
    ## iter  20 value 258.014603
    ## iter  30 value 197.239426
    ## iter  40 value 168.337390
    ## iter  50 value 161.011477
    ## iter  60 value 156.719918
    ## iter  70 value 154.865636
    ## iter  80 value 154.053434
    ## iter  90 value 153.630856
    ## iter 100 value 152.607672
    ## final  value 152.607672 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1200.537310 
    ## iter  10 value 863.013025
    ## iter  20 value 745.704000
    ## iter  30 value 653.310120
    ## iter  40 value 586.474734
    ## iter  50 value 537.383187
    ## iter  60 value 490.832329
    ## iter  70 value 452.460951
    ## iter  80 value 439.574025
    ## iter  90 value 434.906306
    ## iter 100 value 429.178179
    ## final  value 429.178179 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1370.158961 
    ## iter  10 value 622.767814
    ## iter  20 value 354.643364
    ## iter  30 value 271.906304
    ## iter  40 value 195.189796
    ## iter  50 value 135.763447
    ## iter  60 value 112.157412
    ## iter  70 value 90.966956
    ## iter  80 value 80.663521
    ## iter  90 value 74.881653
    ## iter 100 value 73.453382
    ## final  value 73.453382 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1257.940534 
    ## iter  10 value 552.854256
    ## iter  20 value 247.725896
    ## iter  30 value 192.032991
    ## iter  40 value 151.859161
    ## iter  50 value 114.875324
    ## iter  60 value 98.169762
    ## iter  70 value 82.063327
    ## iter  80 value 77.141967
    ## iter  90 value 75.964101
    ## iter 100 value 75.349516
    ## final  value 75.349516 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1211.076750 
    ## iter  10 value 914.179131
    ## iter  20 value 834.206926
    ## iter  30 value 812.372635
    ## iter  40 value 786.007610
    ## iter  50 value 776.599249
    ## iter  60 value 773.710158
    ## iter  70 value 772.982008
    ## iter  80 value 772.864566
    ## iter  90 value 772.790897
    ## iter 100 value 772.747565
    ## final  value 772.747565 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1245.287352 
    ## iter  10 value 523.395636
    ## iter  20 value 374.421231
    ## iter  30 value 293.943466
    ## iter  40 value 256.744594
    ## iter  50 value 228.998404
    ## iter  60 value 202.991912
    ## iter  70 value 189.389061
    ## iter  80 value 179.908298
    ## iter  90 value 174.851221
    ## iter 100 value 171.081407
    ## final  value 171.081407 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1375.499714 
    ## iter  10 value 467.136702
    ## iter  20 value 149.966968
    ## iter  30 value 115.445011
    ## iter  40 value 78.324745
    ## iter  50 value 60.970608
    ## iter  60 value 50.644463
    ## iter  70 value 46.002634
    ## iter  80 value 43.759854
    ## iter  90 value 40.175296
    ## iter 100 value 39.750029
    ## final  value 39.750029 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1290.739320 
    ## iter  10 value 1101.726032
    ## iter  20 value 974.824615
    ## iter  30 value 914.752162
    ## iter  40 value 898.122930
    ## iter  50 value 862.330671
    ## iter  60 value 834.561062
    ## iter  70 value 828.552498
    ## iter  80 value 826.425753
    ## iter  90 value 815.290305
    ## iter 100 value 813.770555
    ## final  value 813.770555 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1324.940174 
    ## iter  10 value 696.343576
    ## iter  20 value 543.903484
    ## iter  30 value 444.323853
    ## iter  40 value 408.306088
    ## iter  50 value 387.989641
    ## iter  60 value 371.456406
    ## iter  70 value 344.661339
    ## iter  80 value 331.819401
    ## iter  90 value 325.545997
    ## iter 100 value 320.561184
    ## final  value 320.561184 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1365.014746 
    ## iter  10 value 653.959963
    ## iter  20 value 328.572869
    ## iter  30 value 226.940427
    ## iter  40 value 179.270308
    ## iter  50 value 166.614896
    ## iter  60 value 160.299180
    ## iter  70 value 157.911284
    ## iter  80 value 156.903776
    ## iter  90 value 154.517520
    ## iter 100 value 153.951713
    ## final  value 153.951713 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1235.630343 
    ## iter  10 value 856.941704
    ## iter  20 value 766.685938
    ## iter  30 value 699.412472
    ## iter  40 value 655.477586
    ## iter  50 value 626.171708
    ## iter  60 value 605.471405
    ## iter  70 value 587.933351
    ## iter  80 value 583.469811
    ## iter  90 value 580.014138
    ## iter 100 value 577.013807
    ## final  value 577.013807 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1191.244819 
    ## iter  10 value 814.611877
    ## iter  20 value 588.618176
    ## iter  30 value 399.634682
    ## iter  40 value 296.408255
    ## iter  50 value 253.885367
    ## iter  60 value 219.789874
    ## iter  70 value 201.061242
    ## iter  80 value 188.607199
    ## iter  90 value 176.980925
    ## iter 100 value 173.579043
    ## final  value 173.579043 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1306.596369 
    ## iter  10 value 327.270414
    ## iter  20 value 107.536745
    ## iter  30 value 55.229320
    ## iter  40 value 31.510806
    ## iter  50 value 14.753883
    ## iter  60 value 11.906555
    ## iter  70 value 10.222157
    ## iter  80 value 8.205057
    ## iter  90 value 7.814816
    ## iter 100 value 7.275833
    ## final  value 7.275833 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1230.266590 
    ## iter  10 value 835.863702
    ## iter  20 value 757.369425
    ## iter  30 value 701.190112
    ## iter  40 value 639.847048
    ## iter  50 value 600.875734
    ## iter  60 value 580.669953
    ## iter  70 value 573.156410
    ## iter  80 value 567.717731
    ## iter  90 value 566.916663
    ## iter 100 value 566.806697
    ## final  value 566.806697 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1278.099060 
    ## iter  10 value 698.624001
    ## iter  20 value 402.745930
    ## iter  30 value 296.255746
    ## iter  40 value 253.012410
    ## iter  50 value 236.004252
    ## iter  60 value 222.671625
    ## iter  70 value 217.296354
    ## iter  80 value 215.693564
    ## iter  90 value 215.255134
    ## iter 100 value 215.057201
    ## final  value 215.057201 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1265.211482 
    ## iter  10 value 366.814299
    ## iter  20 value 124.865431
    ## iter  30 value 56.014908
    ## iter  40 value 30.380761
    ## iter  50 value 19.155784
    ## iter  60 value 14.882596
    ## iter  70 value 14.059561
    ## iter  80 value 13.138546
    ## iter  90 value 10.603321
    ## iter 100 value 8.880269
    ## final  value 8.880269 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1200.757796 
    ## iter  10 value 944.047946
    ## iter  20 value 873.865371
    ## iter  30 value 846.983951
    ## iter  40 value 828.405413
    ## iter  50 value 819.379344
    ## iter  60 value 813.842594
    ## iter  70 value 810.322452
    ## iter  80 value 799.268770
    ## iter  90 value 794.437725
    ## iter 100 value 794.046417
    ## final  value 794.046417 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1252.698674 
    ## iter  10 value 737.081155
    ## iter  20 value 617.133527
    ## iter  30 value 467.899251
    ## iter  40 value 392.151118
    ## iter  50 value 367.318928
    ## iter  60 value 348.282271
    ## iter  70 value 330.696653
    ## iter  80 value 321.179274
    ## iter  90 value 313.195812
    ## iter 100 value 310.427383
    ## final  value 310.427383 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1216.149499 
    ## iter  10 value 572.472796
    ## iter  20 value 304.751443
    ## iter  30 value 216.917112
    ## iter  40 value 176.988263
    ## iter  50 value 163.853555
    ## iter  60 value 159.690968
    ## iter  70 value 158.334655
    ## iter  80 value 154.913973
    ## iter  90 value 151.951115
    ## iter 100 value 150.763369
    ## final  value 150.763369 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1282.756868 
    ## iter  10 value 931.469595
    ## iter  20 value 818.929982
    ## iter  30 value 792.219175
    ## iter  40 value 771.095485
    ## iter  50 value 766.354336
    ## iter  60 value 765.342564
    ## iter  70 value 762.985012
    ## iter  80 value 761.712193
    ## iter  90 value 761.518844
    ## iter 100 value 761.122557
    ## final  value 761.122557 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1377.055894 
    ## iter  10 value 1097.330416
    ## iter  20 value 766.763339
    ## iter  30 value 622.426017
    ## iter  40 value 482.237913
    ## iter  50 value 417.459982
    ## iter  60 value 367.239531
    ## iter  70 value 339.181052
    ## iter  80 value 328.248404
    ## iter  90 value 321.448069
    ## iter 100 value 320.313332
    ## final  value 320.313332 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1296.967745 
    ## iter  10 value 398.343189
    ## iter  20 value 95.451301
    ## iter  30 value 38.485685
    ## iter  40 value 17.688975
    ## iter  50 value 9.225531
    ## iter  60 value 6.850155
    ## iter  70 value 5.549473
    ## iter  80 value 4.967567
    ## iter  90 value 4.603626
    ## iter 100 value 4.398690
    ## final  value 4.398690 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1290.241702 
    ## iter  10 value 940.657322
    ## iter  20 value 757.361765
    ## iter  30 value 669.336274
    ## iter  40 value 620.699659
    ## iter  50 value 577.528793
    ## iter  60 value 557.448045
    ## iter  70 value 554.424701
    ## iter  80 value 553.654249
    ## iter  90 value 552.298586
    ## iter 100 value 549.180514
    ## final  value 549.180514 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1380.889136 
    ## iter  10 value 627.143961
    ## iter  20 value 407.904876
    ## iter  30 value 323.406370
    ## iter  40 value 268.633469
    ## iter  50 value 222.520766
    ## iter  60 value 168.200595
    ## iter  70 value 144.184770
    ## iter  80 value 129.577148
    ## iter  90 value 111.612932
    ## iter 100 value 100.543301
    ## final  value 100.543301 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1222.209429 
    ## iter  10 value 533.112019
    ## iter  20 value 99.290716
    ## iter  30 value 33.622748
    ## iter  40 value 8.962953
    ## iter  50 value 2.077588
    ## iter  60 value 1.918920
    ## iter  70 value 1.912054
    ## iter  80 value 1.910689
    ## iter  90 value 1.910049
    ## iter 100 value 1.909972
    ## final  value 1.909972 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1333.223779 
    ## iter  10 value 880.297554
    ## iter  20 value 768.769762
    ## iter  30 value 745.058785
    ## iter  40 value 734.296360
    ## iter  50 value 729.826788
    ## iter  60 value 726.202050
    ## iter  70 value 723.184201
    ## iter  80 value 721.544699
    ## iter  90 value 721.095892
    ## iter 100 value 721.084459
    ## final  value 721.084459 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1247.019753 
    ## iter  10 value 779.047315
    ## iter  20 value 589.101064
    ## iter  30 value 496.978348
    ## iter  40 value 423.763947
    ## iter  50 value 352.667546
    ## iter  60 value 322.758818
    ## iter  70 value 304.954945
    ## iter  80 value 293.657884
    ## iter  90 value 289.677658
    ## iter 100 value 287.379191
    ## final  value 287.379191 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1261.869523 
    ## iter  10 value 586.276391
    ## iter  20 value 308.431573
    ## iter  30 value 205.986875
    ## iter  40 value 168.802776
    ## iter  50 value 159.123966
    ## iter  60 value 151.472910
    ## iter  70 value 148.689477
    ## iter  80 value 144.852230
    ## iter  90 value 142.772969
    ## iter 100 value 142.371199
    ## final  value 142.371199 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1351.680714 
    ## iter  10 value 1102.578910
    ## iter  20 value 1089.636822
    ## final  value 1089.612087 
    ## converged
    ## # weights:  255
    ## initial  value 1278.362591 
    ## iter  10 value 575.542711
    ## iter  20 value 309.772038
    ## iter  30 value 219.481501
    ## iter  40 value 177.862934
    ## iter  50 value 139.279621
    ## iter  60 value 101.805657
    ## iter  70 value 71.431997
    ## iter  80 value 59.671492
    ## iter  90 value 51.782306
    ## iter 100 value 43.282610
    ## final  value 43.282610 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1270.764603 
    ## iter  10 value 521.815635
    ## iter  20 value 249.749581
    ## iter  30 value 125.930832
    ## iter  40 value 64.985301
    ## iter  50 value 42.811252
    ## iter  60 value 31.621513
    ## iter  70 value 30.562547
    ## iter  80 value 29.634934
    ## iter  90 value 29.236725
    ## iter 100 value 29.056984
    ## final  value 29.056984 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1303.852721 
    ## iter  10 value 983.467183
    ## iter  20 value 814.356119
    ## iter  30 value 758.999696
    ## iter  40 value 708.313680
    ## iter  50 value 679.853117
    ## iter  60 value 660.638325
    ## iter  70 value 658.155288
    ## iter  80 value 656.341386
    ## iter  90 value 655.491374
    ## iter 100 value 655.079526
    ## final  value 655.079526 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1372.293297 
    ## iter  10 value 597.117494
    ## iter  20 value 321.398601
    ## iter  30 value 214.795730
    ## iter  40 value 165.112874
    ## iter  50 value 140.503477
    ## iter  60 value 125.151169
    ## iter  70 value 113.017842
    ## iter  80 value 108.486938
    ## iter  90 value 106.197758
    ## iter 100 value 105.210029
    ## final  value 105.210029 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1238.354305 
    ## iter  10 value 487.488863
    ## iter  20 value 186.929916
    ## iter  30 value 116.971079
    ## iter  40 value 50.651349
    ## iter  50 value 34.052476
    ## iter  60 value 23.150659
    ## iter  70 value 16.883429
    ## iter  80 value 10.364322
    ## iter  90 value 7.213762
    ## iter 100 value 5.671066
    ## final  value 5.671066 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1260.098803 
    ## iter  10 value 1002.013224
    ## iter  20 value 826.286547
    ## iter  30 value 791.725527
    ## iter  40 value 780.973635
    ## iter  50 value 765.786286
    ## iter  60 value 762.858614
    ## iter  70 value 761.697014
    ## iter  80 value 761.546649
    ## iter  90 value 761.497577
    ## final  value 761.494697 
    ## converged
    ## # weights:  255
    ## initial  value 1213.768071 
    ## iter  10 value 689.413307
    ## iter  20 value 543.319534
    ## iter  30 value 485.377109
    ## iter  40 value 423.370190
    ## iter  50 value 369.802823
    ## iter  60 value 345.102424
    ## iter  70 value 320.222548
    ## iter  80 value 305.803972
    ## iter  90 value 296.728281
    ## iter 100 value 292.598150
    ## final  value 292.598150 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1243.230057 
    ## iter  10 value 533.124790
    ## iter  20 value 349.601136
    ## iter  30 value 255.484284
    ## iter  40 value 199.206073
    ## iter  50 value 179.158007
    ## iter  60 value 158.172829
    ## iter  70 value 153.452229
    ## iter  80 value 151.880020
    ## iter  90 value 150.985888
    ## iter 100 value 150.401346
    ## final  value 150.401346 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1266.888883 
    ## iter  10 value 1122.466678
    ## iter  20 value 957.857727
    ## iter  30 value 831.757412
    ## iter  40 value 713.913886
    ## iter  50 value 678.340210
    ## iter  60 value 651.417390
    ## iter  70 value 643.804581
    ## iter  80 value 639.623364
    ## iter  90 value 637.775976
    ## iter 100 value 624.444172
    ## final  value 624.444172 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1304.583110 
    ## iter  10 value 609.110772
    ## iter  20 value 387.479480
    ## iter  30 value 231.136785
    ## iter  40 value 145.058093
    ## iter  50 value 109.505334
    ## iter  60 value 94.033618
    ## iter  70 value 87.047589
    ## iter  80 value 84.216327
    ## iter  90 value 80.141652
    ## iter 100 value 76.692217
    ## final  value 76.692217 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1346.125682 
    ## iter  10 value 430.199643
    ## iter  20 value 93.973522
    ## iter  30 value 31.148004
    ## iter  40 value 9.431532
    ## iter  50 value 2.303737
    ## iter  60 value 1.918043
    ## iter  70 value 1.731747
    ## iter  80 value 1.567495
    ## iter  90 value 1.441179
    ## iter 100 value 1.352798
    ## final  value 1.352798 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1203.977530 
    ## iter  10 value 949.335195
    ## iter  20 value 810.046716
    ## iter  30 value 737.626524
    ## iter  40 value 696.465490
    ## iter  50 value 654.290419
    ## iter  60 value 623.635374
    ## iter  70 value 605.574749
    ## iter  80 value 593.482759
    ## iter  90 value 589.225392
    ## iter 100 value 588.523673
    ## final  value 588.523673 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1290.940208 
    ## iter  10 value 637.856898
    ## iter  20 value 391.087749
    ## iter  30 value 325.881095
    ## iter  40 value 270.362769
    ## iter  50 value 241.281638
    ## iter  60 value 225.776483
    ## iter  70 value 198.832314
    ## iter  80 value 177.933684
    ## iter  90 value 171.925957
    ## iter 100 value 164.285617
    ## final  value 164.285617 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1367.059478 
    ## iter  10 value 429.061113
    ## iter  20 value 138.420788
    ## iter  30 value 83.471448
    ## iter  40 value 52.332812
    ## iter  50 value 38.424444
    ## iter  60 value 34.767979
    ## iter  70 value 32.517799
    ## iter  80 value 31.321854
    ## iter  90 value 23.063523
    ## iter 100 value 13.967453
    ## final  value 13.967453 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1329.142144 
    ## iter  10 value 1027.712996
    ## iter  20 value 838.839003
    ## iter  30 value 780.668260
    ## iter  40 value 760.579752
    ## iter  50 value 749.416037
    ## iter  60 value 744.822385
    ## iter  70 value 741.469010
    ## iter  80 value 741.028444
    ## iter  90 value 740.203631
    ## iter 100 value 740.132539
    ## final  value 740.132539 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1250.888807 
    ## iter  10 value 685.704437
    ## iter  20 value 469.480813
    ## iter  30 value 376.429683
    ## iter  40 value 349.895044
    ## iter  50 value 323.120804
    ## iter  60 value 309.121822
    ## iter  70 value 299.767675
    ## iter  80 value 296.419633
    ## iter  90 value 294.878744
    ## iter 100 value 294.247971
    ## final  value 294.247971 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1261.300625 
    ## iter  10 value 510.189896
    ## iter  20 value 268.915968
    ## iter  30 value 186.453852
    ## iter  40 value 158.463134
    ## iter  50 value 151.654769
    ## iter  60 value 146.686712
    ## iter  70 value 145.371669
    ## iter  80 value 144.674607
    ## iter  90 value 144.514657
    ## iter 100 value 144.447064
    ## final  value 144.447064 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1280.901040 
    ## iter  10 value 934.871222
    ## iter  20 value 761.625738
    ## iter  30 value 710.610381
    ## iter  40 value 677.041351
    ## iter  50 value 663.267631
    ## iter  60 value 653.428120
    ## iter  70 value 622.071530
    ## iter  80 value 594.120048
    ## iter  90 value 550.248455
    ## iter 100 value 537.552999
    ## final  value 537.552999 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1309.716858 
    ## iter  10 value 861.952003
    ## iter  20 value 496.422650
    ## iter  30 value 364.343712
    ## iter  40 value 294.920651
    ## iter  50 value 259.227494
    ## iter  60 value 220.392016
    ## iter  70 value 204.466932
    ## iter  80 value 198.655918
    ## iter  90 value 194.496890
    ## iter 100 value 188.518788
    ## final  value 188.518788 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1294.167304 
    ## iter  10 value 488.612837
    ## iter  20 value 144.760313
    ## iter  30 value 95.927948
    ## iter  40 value 70.743317
    ## iter  50 value 30.352687
    ## iter  60 value 17.958900
    ## iter  70 value 11.972141
    ## iter  80 value 9.064736
    ## iter  90 value 8.417290
    ## iter 100 value 6.347050
    ## final  value 6.347050 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1243.634748 
    ## iter  10 value 971.358526
    ## iter  20 value 897.108203
    ## iter  30 value 849.039368
    ## iter  40 value 800.516742
    ## iter  50 value 773.276512
    ## iter  60 value 746.136266
    ## iter  70 value 670.651782
    ## iter  80 value 642.839357
    ## iter  90 value 618.072167
    ## iter 100 value 611.612572
    ## final  value 611.612572 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1338.839876 
    ## iter  10 value 837.424304
    ## iter  20 value 645.830907
    ## iter  30 value 564.584214
    ## iter  40 value 483.318765
    ## iter  50 value 432.825752
    ## iter  60 value 393.331542
    ## iter  70 value 366.820670
    ## iter  80 value 324.799301
    ## iter  90 value 301.916413
    ## iter 100 value 264.380664
    ## final  value 264.380664 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1234.949348 
    ## iter  10 value 436.967580
    ## iter  20 value 129.615658
    ## iter  30 value 58.360450
    ## iter  40 value 36.057397
    ## iter  50 value 28.968574
    ## iter  60 value 25.999299
    ## iter  70 value 24.675679
    ## iter  80 value 23.731178
    ## iter  90 value 22.979386
    ## iter 100 value 22.811118
    ## final  value 22.811118 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1271.368953 
    ## iter  10 value 953.713734
    ## iter  20 value 894.524050
    ## iter  30 value 851.344974
    ## iter  40 value 831.787137
    ## iter  50 value 822.803999
    ## iter  60 value 817.141212
    ## iter  70 value 810.882204
    ## iter  80 value 809.621148
    ## iter  90 value 809.362167
    ## iter 100 value 809.337673
    ## final  value 809.337673 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1292.202027 
    ## iter  10 value 611.896935
    ## iter  20 value 456.778943
    ## iter  30 value 402.684230
    ## iter  40 value 366.738372
    ## iter  50 value 342.896130
    ## iter  60 value 333.628127
    ## iter  70 value 328.414899
    ## iter  80 value 321.176679
    ## iter  90 value 315.597317
    ## iter 100 value 313.542596
    ## final  value 313.542596 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1470.787870 
    ## iter  10 value 584.690313
    ## iter  20 value 274.301126
    ## iter  30 value 195.851465
    ## iter  40 value 163.334650
    ## iter  50 value 157.259526
    ## iter  60 value 153.770524
    ## iter  70 value 149.030073
    ## iter  80 value 144.935075
    ## iter  90 value 143.669041
    ## iter 100 value 143.538117
    ## final  value 143.538117 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1277.496111 
    ## iter  10 value 910.629308
    ## iter  20 value 833.139615
    ## iter  30 value 780.098717
    ## iter  40 value 733.662834
    ## iter  50 value 704.434569
    ## iter  60 value 665.257564
    ## iter  70 value 653.046649
    ## iter  80 value 643.241811
    ## iter  90 value 625.687266
    ## iter 100 value 620.104135
    ## final  value 620.104135 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1317.632550 
    ## iter  10 value 705.115485
    ## iter  20 value 376.294799
    ## iter  30 value 287.943413
    ## iter  40 value 241.360074
    ## iter  50 value 214.882505
    ## iter  60 value 187.152227
    ## iter  70 value 173.317660
    ## iter  80 value 163.694265
    ## iter  90 value 158.892451
    ## iter 100 value 157.425518
    ## final  value 157.425518 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1316.125676 
    ## iter  10 value 437.790885
    ## iter  20 value 153.825366
    ## iter  30 value 79.278897
    ## iter  40 value 38.514120
    ## iter  50 value 22.790474
    ## iter  60 value 13.140840
    ## iter  70 value 9.205450
    ## iter  80 value 4.045917
    ## iter  90 value 2.464712
    ## iter 100 value 2.201036
    ## final  value 2.201036 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1320.073020 
    ## iter  10 value 1112.213709
    ## iter  20 value 1102.288729
    ## iter  30 value 1102.158244
    ## final  value 1102.156091 
    ## converged
    ## # weights:  255
    ## initial  value 1197.905765 
    ## iter  10 value 731.366590
    ## iter  20 value 475.740300
    ## iter  30 value 346.903413
    ## iter  40 value 280.611915
    ## iter  50 value 252.124214
    ## iter  60 value 242.105489
    ## iter  70 value 235.942775
    ## iter  80 value 234.200908
    ## iter  90 value 233.637740
    ## iter 100 value 233.346972
    ## final  value 233.346972 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1368.608895 
    ## iter  10 value 463.979748
    ## iter  20 value 220.323785
    ## iter  30 value 143.870883
    ## iter  40 value 124.961003
    ## iter  50 value 91.341773
    ## iter  60 value 73.494406
    ## iter  70 value 57.875005
    ## iter  80 value 51.362449
    ## iter  90 value 45.921128
    ## iter 100 value 37.549271
    ## final  value 37.549271 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1278.969696 
    ## iter  10 value 931.358410
    ## iter  20 value 836.131170
    ## iter  30 value 785.466926
    ## iter  40 value 767.926462
    ## iter  50 value 759.765134
    ## iter  60 value 756.292772
    ## iter  70 value 756.075373
    ## iter  80 value 756.071119
    ## final  value 756.071049 
    ## converged
    ## # weights:  255
    ## initial  value 1267.080545 
    ## iter  10 value 736.793133
    ## iter  20 value 594.808859
    ## iter  30 value 502.802396
    ## iter  40 value 411.800559
    ## iter  50 value 362.911108
    ## iter  60 value 337.716738
    ## iter  70 value 317.013973
    ## iter  80 value 312.427480
    ## iter  90 value 311.120671
    ## iter 100 value 310.806486
    ## final  value 310.806486 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1530.571528 
    ## iter  10 value 951.697469
    ## iter  20 value 802.144816
    ## iter  30 value 737.959595
    ## iter  40 value 507.607599
    ## iter  50 value 394.903413
    ## iter  60 value 249.314386
    ## iter  70 value 198.129721
    ## iter  80 value 181.331719
    ## iter  90 value 168.247239
    ## iter 100 value 160.526778
    ## final  value 160.526778 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1238.663270 
    ## iter  10 value 980.076638
    ## iter  20 value 871.059625
    ## iter  30 value 848.004438
    ## iter  40 value 823.607847
    ## iter  50 value 773.091239
    ## iter  60 value 745.819887
    ## iter  70 value 738.933598
    ## iter  80 value 736.639218
    ## iter  90 value 735.034512
    ## iter 100 value 733.975507
    ## final  value 733.975507 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1296.521760 
    ## iter  10 value 640.631235
    ## iter  20 value 411.844579
    ## iter  30 value 326.442071
    ## iter  40 value 282.732122
    ## iter  50 value 252.382629
    ## iter  60 value 226.588839
    ## iter  70 value 206.027755
    ## iter  80 value 200.564595
    ## iter  90 value 199.000717
    ## iter 100 value 197.811571
    ## final  value 197.811571 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1290.603524 
    ## iter  10 value 516.051198
    ## iter  20 value 241.893165
    ## iter  30 value 125.094132
    ## iter  40 value 94.985773
    ## iter  50 value 57.482105
    ## iter  60 value 52.862213
    ## iter  70 value 51.389513
    ## iter  80 value 45.460831
    ## iter  90 value 42.346648
    ## iter 100 value 39.262047
    ## final  value 39.262047 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1264.184058 
    ## iter  10 value 1046.023767
    ## iter  20 value 830.756034
    ## iter  30 value 815.855846
    ## iter  40 value 805.172210
    ## iter  50 value 788.209387
    ## iter  60 value 780.364656
    ## iter  70 value 775.551860
    ## iter  80 value 774.437541
    ## iter  90 value 774.149404
    ## iter 100 value 774.045281
    ## final  value 774.045281 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1239.359175 
    ## iter  10 value 743.845590
    ## iter  20 value 441.009854
    ## iter  30 value 341.413818
    ## iter  40 value 294.463697
    ## iter  50 value 232.852228
    ## iter  60 value 200.949481
    ## iter  70 value 183.005262
    ## iter  80 value 169.612643
    ## iter  90 value 161.836514
    ## iter 100 value 154.591112
    ## final  value 154.591112 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1361.110979 
    ## iter  10 value 518.842125
    ## iter  20 value 198.276306
    ## iter  30 value 89.959251
    ## iter  40 value 55.131586
    ## iter  50 value 38.688914
    ## iter  60 value 25.001677
    ## iter  70 value 21.753953
    ## iter  80 value 19.612991
    ## iter  90 value 17.824612
    ## iter 100 value 14.053498
    ## final  value 14.053498 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1198.492891 
    ## iter  10 value 872.144718
    ## iter  20 value 821.618061
    ## iter  30 value 802.068211
    ## iter  40 value 787.152927
    ## iter  50 value 766.729075
    ## iter  60 value 757.554927
    ## iter  70 value 754.801589
    ## iter  80 value 753.281689
    ## iter  90 value 752.606147
    ## iter 100 value 752.528656
    ## final  value 752.528656 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1411.732111 
    ## iter  10 value 644.091795
    ## iter  20 value 436.381580
    ## iter  30 value 393.892437
    ## iter  40 value 365.845625
    ## iter  50 value 338.552231
    ## iter  60 value 328.144476
    ## iter  70 value 322.091191
    ## iter  80 value 315.539676
    ## iter  90 value 308.154887
    ## iter 100 value 304.415466
    ## final  value 304.415466 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1310.545630 
    ## iter  10 value 541.820244
    ## iter  20 value 370.342556
    ## iter  30 value 257.119781
    ## iter  40 value 203.888098
    ## iter  50 value 180.204557
    ## iter  60 value 167.433680
    ## iter  70 value 158.450689
    ## iter  80 value 154.732940
    ## iter  90 value 153.725762
    ## iter 100 value 152.956643
    ## final  value 152.956643 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1232.747572 
    ## iter  10 value 942.261380
    ## iter  20 value 842.113384
    ## iter  30 value 760.026566
    ## iter  40 value 713.342537
    ## iter  50 value 689.154475
    ## iter  60 value 673.474142
    ## iter  70 value 671.472123
    ## iter  80 value 671.033418
    ## iter  90 value 668.259844
    ## iter 100 value 667.503987
    ## final  value 667.503987 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1241.796247 
    ## iter  10 value 611.111181
    ## iter  20 value 461.773171
    ## iter  30 value 383.984467
    ## iter  40 value 334.536931
    ## iter  50 value 276.748388
    ## iter  60 value 248.399928
    ## iter  70 value 235.774598
    ## iter  80 value 233.720506
    ## iter  90 value 231.129376
    ## iter 100 value 227.491826
    ## final  value 227.491826 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1337.255436 
    ## iter  10 value 409.829902
    ## iter  20 value 105.565979
    ## iter  30 value 52.353197
    ## iter  40 value 18.190639
    ## iter  50 value 6.130289
    ## iter  60 value 3.898413
    ## iter  70 value 3.754674
    ## iter  80 value 3.570058
    ## iter  90 value 3.429334
    ## iter 100 value 3.310118
    ## final  value 3.310118 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1244.499761 
    ## iter  10 value 954.219752
    ## iter  20 value 805.029828
    ## iter  30 value 785.695100
    ## iter  40 value 758.240764
    ## iter  50 value 738.712539
    ## iter  60 value 736.346786
    ## iter  70 value 734.852597
    ## iter  80 value 733.973909
    ## iter  90 value 733.713623
    ## iter 100 value 733.600326
    ## final  value 733.600326 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1218.069230 
    ## iter  10 value 565.089301
    ## iter  20 value 381.500374
    ## iter  30 value 315.968566
    ## iter  40 value 278.164984
    ## iter  50 value 260.476260
    ## iter  60 value 249.277818
    ## iter  70 value 238.665030
    ## iter  80 value 229.523187
    ## iter  90 value 225.806504
    ## iter 100 value 223.587339
    ## final  value 223.587339 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1292.941218 
    ## iter  10 value 295.732414
    ## iter  20 value 82.702076
    ## iter  30 value 25.738544
    ## iter  40 value 16.669217
    ## iter  50 value 6.620744
    ## iter  60 value 0.932065
    ## iter  70 value 0.071973
    ## iter  80 value 0.001832
    ## iter  90 value 0.000445
    ## final  value 0.000088 
    ## converged
    ## # weights:  95
    ## initial  value 1202.813078 
    ## iter  10 value 974.079747
    ## iter  20 value 881.928026
    ## iter  30 value 848.532418
    ## iter  40 value 813.707210
    ## iter  50 value 802.942077
    ## iter  60 value 797.254007
    ## iter  70 value 774.300186
    ## iter  80 value 768.164743
    ## iter  90 value 765.891200
    ## iter 100 value 765.602983
    ## final  value 765.602983 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1433.619374 
    ## iter  10 value 816.679129
    ## iter  20 value 631.595601
    ## iter  30 value 507.527652
    ## iter  40 value 405.748894
    ## iter  50 value 356.564538
    ## iter  60 value 327.285478
    ## iter  70 value 314.410514
    ## iter  80 value 308.139122
    ## iter  90 value 305.625863
    ## iter 100 value 302.397753
    ## final  value 302.397753 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1268.897055 
    ## iter  10 value 644.341440
    ## iter  20 value 327.876458
    ## iter  30 value 202.672116
    ## iter  40 value 172.685649
    ## iter  50 value 162.163794
    ## iter  60 value 160.501115
    ## iter  70 value 159.076248
    ## iter  80 value 158.437271
    ## iter  90 value 156.526673
    ## iter 100 value 154.804797
    ## final  value 154.804797 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1220.061579 
    ## iter  10 value 837.278861
    ## iter  20 value 772.810146
    ## iter  30 value 712.876576
    ## iter  40 value 678.314786
    ## iter  50 value 654.472210
    ## iter  60 value 641.698792
    ## iter  70 value 626.878990
    ## iter  80 value 620.784103
    ## iter  90 value 616.549614
    ## iter 100 value 612.107632
    ## final  value 612.107632 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1295.853056 
    ## iter  10 value 716.376695
    ## iter  20 value 449.095040
    ## iter  30 value 387.781416
    ## iter  40 value 343.187108
    ## iter  50 value 264.806727
    ## iter  60 value 228.938989
    ## iter  70 value 204.310962
    ## iter  80 value 193.227961
    ## iter  90 value 190.885484
    ## iter 100 value 190.041630
    ## final  value 190.041630 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1354.416681 
    ## iter  10 value 733.951011
    ## iter  20 value 439.115196
    ## iter  30 value 274.614128
    ## iter  40 value 98.891357
    ## iter  50 value 45.371292
    ## iter  60 value 27.235265
    ## iter  70 value 15.632751
    ## iter  80 value 9.501601
    ## iter  90 value 6.155100
    ## iter 100 value 4.118183
    ## final  value 4.118183 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1360.526009 
    ## iter  10 value 1058.005366
    ## iter  20 value 1018.851530
    ## iter  30 value 946.236106
    ## iter  40 value 843.549556
    ## iter  50 value 765.997358
    ## iter  60 value 738.097493
    ## iter  70 value 724.049580
    ## iter  80 value 717.484248
    ## iter  90 value 702.911308
    ## iter 100 value 685.536553
    ## final  value 685.536553 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1256.663600 
    ## iter  10 value 767.430591
    ## iter  20 value 518.704569
    ## iter  30 value 423.097027
    ## iter  40 value 352.093539
    ## iter  50 value 305.148012
    ## iter  60 value 269.308476
    ## iter  70 value 242.702815
    ## iter  80 value 222.243523
    ## iter  90 value 218.750327
    ## iter 100 value 217.144033
    ## final  value 217.144033 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1364.450521 
    ## iter  10 value 504.745829
    ## iter  20 value 145.471173
    ## iter  30 value 63.195286
    ## iter  40 value 43.913953
    ## iter  50 value 30.076844
    ## iter  60 value 17.346986
    ## iter  70 value 9.310569
    ## iter  80 value 4.504378
    ## iter  90 value 4.164917
    ## iter 100 value 3.939122
    ## final  value 3.939122 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1134.437168 
    ## iter  10 value 896.586379
    ## iter  20 value 836.314995
    ## iter  30 value 809.098051
    ## iter  40 value 779.062533
    ## iter  50 value 769.447981
    ## iter  60 value 765.290864
    ## iter  70 value 760.032375
    ## iter  80 value 759.034777
    ## iter  90 value 758.192395
    ## iter 100 value 755.666145
    ## final  value 755.666145 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1241.765102 
    ## iter  10 value 699.957669
    ## iter  20 value 501.635657
    ## iter  30 value 414.017206
    ## iter  40 value 370.377618
    ## iter  50 value 345.926435
    ## iter  60 value 331.634431
    ## iter  70 value 321.640364
    ## iter  80 value 318.311093
    ## iter  90 value 317.072879
    ## iter 100 value 316.042065
    ## final  value 316.042065 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1272.765832 
    ## iter  10 value 548.722453
    ## iter  20 value 331.859340
    ## iter  30 value 239.187636
    ## iter  40 value 202.864481
    ## iter  50 value 189.238184
    ## iter  60 value 178.117012
    ## iter  70 value 170.541420
    ## iter  80 value 163.098145
    ## iter  90 value 161.572146
    ## iter 100 value 160.751887
    ## final  value 160.751887 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1251.082933 
    ## iter  10 value 1110.727919
    ## iter  20 value 1109.376880
    ## iter  30 value 1109.356080
    ## iter  40 value 1106.778698
    ## iter  50 value 1036.335738
    ## iter  60 value 970.188787
    ## iter  70 value 930.768530
    ## iter  80 value 891.035093
    ## iter  90 value 879.841953
    ## iter 100 value 871.776684
    ## final  value 871.776684 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1260.934601 
    ## iter  10 value 710.052232
    ## iter  20 value 386.292319
    ## iter  30 value 327.873490
    ## iter  40 value 291.086349
    ## iter  50 value 251.277063
    ## iter  60 value 223.755789
    ## iter  70 value 209.576871
    ## iter  80 value 182.839235
    ## iter  90 value 176.164239
    ## iter 100 value 169.863411
    ## final  value 169.863411 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1344.433222 
    ## iter  10 value 522.260179
    ## iter  20 value 231.949522
    ## iter  30 value 104.919935
    ## iter  40 value 61.587131
    ## iter  50 value 50.758815
    ## iter  60 value 47.927383
    ## iter  70 value 46.220742
    ## iter  80 value 42.185714
    ## iter  90 value 37.512169
    ## iter 100 value 33.919222
    ## final  value 33.919222 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1279.612700 
    ## iter  10 value 934.526395
    ## iter  20 value 802.927324
    ## iter  30 value 750.242087
    ## iter  40 value 710.273827
    ## iter  50 value 683.024223
    ## iter  60 value 661.924160
    ## iter  70 value 642.829083
    ## iter  80 value 615.782597
    ## iter  90 value 603.820212
    ## iter 100 value 595.842908
    ## final  value 595.842908 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1348.179726 
    ## iter  10 value 592.212187
    ## iter  20 value 310.671042
    ## iter  30 value 244.151282
    ## iter  40 value 207.706214
    ## iter  50 value 189.114260
    ## iter  60 value 174.735877
    ## iter  70 value 159.617172
    ## iter  80 value 151.137096
    ## iter  90 value 146.467201
    ## iter 100 value 144.417110
    ## final  value 144.417110 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1372.419469 
    ## iter  10 value 463.508871
    ## iter  20 value 100.632351
    ## iter  30 value 32.754998
    ## iter  40 value 16.630436
    ## iter  50 value 11.196001
    ## iter  60 value 7.880011
    ## iter  70 value 2.970170
    ## iter  80 value 2.656828
    ## iter  90 value 2.576700
    ## iter 100 value 2.542985
    ## final  value 2.542985 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1218.782299 
    ## iter  10 value 904.264382
    ## iter  20 value 824.651827
    ## iter  30 value 802.120937
    ## iter  40 value 777.051640
    ## iter  50 value 764.004760
    ## iter  60 value 757.925101
    ## iter  70 value 742.067501
    ## iter  80 value 735.252579
    ## iter  90 value 734.275099
    ## iter 100 value 734.187190
    ## final  value 734.187190 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1270.000222 
    ## iter  10 value 633.152000
    ## iter  20 value 511.315589
    ## iter  30 value 469.245126
    ## iter  40 value 443.220609
    ## iter  50 value 417.773724
    ## iter  60 value 397.362138
    ## iter  70 value 375.439411
    ## iter  80 value 352.120031
    ## iter  90 value 331.737913
    ## iter 100 value 312.612571
    ## final  value 312.612571 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1404.286740 
    ## iter  10 value 556.928877
    ## iter  20 value 312.559565
    ## iter  30 value 218.242961
    ## iter  40 value 171.650390
    ## iter  50 value 162.899569
    ## iter  60 value 155.078073
    ## iter  70 value 151.810265
    ## iter  80 value 150.772159
    ## iter  90 value 149.491908
    ## iter 100 value 147.894265
    ## final  value 147.894265 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1161.363015 
    ## iter  10 value 892.527620
    ## iter  20 value 836.510479
    ## iter  30 value 763.460715
    ## iter  40 value 708.901776
    ## iter  50 value 637.699594
    ## iter  60 value 603.724097
    ## iter  70 value 579.603333
    ## iter  80 value 561.560390
    ## iter  90 value 488.199800
    ## iter 100 value 464.537306
    ## final  value 464.537306 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1342.383073 
    ## iter  10 value 817.620009
    ## iter  20 value 516.505764
    ## iter  30 value 438.657092
    ## iter  40 value 394.811341
    ## iter  50 value 367.583578
    ## iter  60 value 353.751598
    ## iter  70 value 338.566509
    ## iter  80 value 329.022735
    ## iter  90 value 325.794594
    ## iter 100 value 323.948326
    ## final  value 323.948326 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1364.191892 
    ## iter  10 value 445.269548
    ## iter  20 value 125.987934
    ## iter  30 value 65.611342
    ## iter  40 value 49.580669
    ## iter  50 value 39.072615
    ## iter  60 value 33.565972
    ## iter  70 value 31.649318
    ## iter  80 value 30.198897
    ## iter  90 value 27.060437
    ## iter 100 value 23.951894
    ## final  value 23.951894 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1305.797573 
    ## iter  10 value 891.662277
    ## iter  20 value 763.140878
    ## iter  30 value 713.798138
    ## iter  40 value 660.039625
    ## iter  50 value 623.810938
    ## iter  60 value 615.888742
    ## iter  70 value 605.538940
    ## iter  80 value 601.860437
    ## iter  90 value 600.150533
    ## iter 100 value 599.555953
    ## final  value 599.555953 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1264.726186 
    ## iter  10 value 559.007104
    ## iter  20 value 263.574192
    ## iter  30 value 208.384787
    ## iter  40 value 177.821090
    ## iter  50 value 148.140399
    ## iter  60 value 129.492363
    ## iter  70 value 120.551909
    ## iter  80 value 114.639506
    ## iter  90 value 111.104010
    ## iter 100 value 106.790234
    ## final  value 106.790234 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1359.393360 
    ## iter  10 value 495.812477
    ## iter  20 value 156.258816
    ## iter  30 value 66.096444
    ## iter  40 value 40.151027
    ## iter  50 value 19.308069
    ## iter  60 value 14.671463
    ## iter  70 value 13.908232
    ## iter  80 value 13.506338
    ## iter  90 value 13.175382
    ## iter 100 value 13.097292
    ## final  value 13.097292 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1250.642391 
    ## iter  10 value 1120.290199
    ## iter  20 value 949.205639
    ## iter  30 value 843.679029
    ## iter  40 value 805.464545
    ## iter  50 value 778.844160
    ## iter  60 value 770.989766
    ## iter  70 value 764.573957
    ## iter  80 value 755.821551
    ## iter  90 value 754.304068
    ## iter 100 value 752.776109
    ## final  value 752.776109 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1250.454055 
    ## iter  10 value 767.701683
    ## iter  20 value 530.548957
    ## iter  30 value 431.018418
    ## iter  40 value 377.426928
    ## iter  50 value 349.275371
    ## iter  60 value 323.846912
    ## iter  70 value 315.072912
    ## iter  80 value 313.068357
    ## iter  90 value 312.386138
    ## iter 100 value 311.986493
    ## final  value 311.986493 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1288.773833 
    ## iter  10 value 490.260453
    ## iter  20 value 296.646497
    ## iter  30 value 222.527606
    ## iter  40 value 181.845285
    ## iter  50 value 163.657598
    ## iter  60 value 159.853573
    ## iter  70 value 159.090100
    ## iter  80 value 158.927428
    ## iter  90 value 158.857381
    ## iter 100 value 158.532737
    ## final  value 158.532737 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1345.336475 
    ## iter  10 value 1105.403829
    ## iter  20 value 1099.077614
    ## iter  30 value 1039.866588
    ## iter  40 value 852.112546
    ## iter  50 value 762.678540
    ## iter  60 value 718.748614
    ## iter  70 value 701.802578
    ## iter  80 value 697.268959
    ## iter  90 value 693.598029
    ## iter 100 value 692.211338
    ## final  value 692.211338 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1237.301811 
    ## iter  10 value 652.265686
    ## iter  20 value 357.393391
    ## iter  30 value 267.660200
    ## iter  40 value 211.427215
    ## iter  50 value 186.842476
    ## iter  60 value 171.661236
    ## iter  70 value 159.212254
    ## iter  80 value 155.221120
    ## iter  90 value 153.191421
    ## iter 100 value 149.367156
    ## final  value 149.367156 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1223.109019 
    ## iter  10 value 456.756044
    ## iter  20 value 156.453514
    ## iter  30 value 72.817165
    ## iter  40 value 35.237990
    ## iter  50 value 24.668625
    ## iter  60 value 15.439892
    ## iter  70 value 14.095563
    ## iter  80 value 13.750989
    ## iter  90 value 13.453831
    ## iter 100 value 13.262850
    ## final  value 13.262850 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1259.116469 
    ## iter  10 value 969.327058
    ## iter  20 value 828.035912
    ## iter  30 value 801.553968
    ## iter  40 value 781.252013
    ## iter  50 value 766.445147
    ## iter  60 value 762.678442
    ## iter  70 value 760.628051
    ## iter  80 value 759.757129
    ## iter  90 value 759.264263
    ## iter 100 value 759.116024
    ## final  value 759.116024 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1252.144992 
    ## iter  10 value 603.062894
    ## iter  20 value 339.088232
    ## iter  30 value 297.991597
    ## iter  40 value 265.197368
    ## iter  50 value 242.409371
    ## iter  60 value 232.321486
    ## iter  70 value 222.435262
    ## iter  80 value 219.079315
    ## iter  90 value 216.390240
    ## iter 100 value 210.966146
    ## final  value 210.966146 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1291.532958 
    ## iter  10 value 368.529311
    ## iter  20 value 98.036536
    ## iter  30 value 26.429980
    ## iter  40 value 7.971894
    ## iter  50 value 0.301162
    ## iter  60 value 0.017728
    ## iter  70 value 0.003507
    ## iter  80 value 0.001534
    ## iter  90 value 0.000717
    ## iter 100 value 0.000428
    ## final  value 0.000428 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1223.437637 
    ## iter  10 value 953.125753
    ## iter  20 value 883.777599
    ## iter  30 value 843.295788
    ## iter  40 value 832.306093
    ## iter  50 value 825.333847
    ## iter  60 value 816.261337
    ## iter  70 value 801.124843
    ## iter  80 value 788.601474
    ## iter  90 value 784.744719
    ## iter 100 value 783.473534
    ## final  value 783.473534 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1238.494856 
    ## iter  10 value 710.951950
    ## iter  20 value 522.876888
    ## iter  30 value 466.176082
    ## iter  40 value 406.514118
    ## iter  50 value 386.235414
    ## iter  60 value 374.011099
    ## iter  70 value 361.333945
    ## iter  80 value 345.579647
    ## iter  90 value 339.264071
    ## iter 100 value 325.126026
    ## final  value 325.126026 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1292.222254 
    ## iter  10 value 635.072414
    ## iter  20 value 432.676785
    ## iter  30 value 324.207684
    ## iter  40 value 255.795841
    ## iter  50 value 186.039949
    ## iter  60 value 166.933762
    ## iter  70 value 162.208190
    ## iter  80 value 159.272386
    ## iter  90 value 157.255781
    ## iter 100 value 156.114447
    ## final  value 156.114447 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1266.542422 
    ## iter  10 value 902.591397
    ## iter  20 value 828.111538
    ## iter  30 value 781.617594
    ## iter  40 value 740.050547
    ## iter  50 value 676.752274
    ## iter  60 value 635.482721
    ## iter  70 value 614.016498
    ## iter  80 value 608.082019
    ## iter  90 value 597.315134
    ## iter 100 value 595.441995
    ## final  value 595.441995 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1239.378609 
    ## iter  10 value 582.013715
    ## iter  20 value 297.336115
    ## iter  30 value 242.849419
    ## iter  40 value 205.458390
    ## iter  50 value 175.761277
    ## iter  60 value 160.510884
    ## iter  70 value 156.717147
    ## iter  80 value 155.014714
    ## iter  90 value 152.837424
    ## iter 100 value 151.701768
    ## final  value 151.701768 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1338.784942 
    ## iter  10 value 474.730768
    ## iter  20 value 163.928915
    ## iter  30 value 85.972475
    ## iter  40 value 55.338869
    ## iter  50 value 45.174538
    ## iter  60 value 37.958025
    ## iter  70 value 29.640638
    ## iter  80 value 23.286378
    ## iter  90 value 10.217516
    ## iter 100 value 8.186426
    ## final  value 8.186426 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1256.026327 
    ## iter  10 value 1115.328364
    ## iter  20 value 1112.910054
    ## iter  30 value 1079.069713
    ## iter  40 value 967.176132
    ## iter  50 value 939.734660
    ## iter  60 value 926.953840
    ## iter  70 value 912.102144
    ## iter  80 value 907.020591
    ## iter  90 value 904.057831
    ## iter 100 value 902.644623
    ## final  value 902.644623 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1207.573848 
    ## iter  10 value 592.296191
    ## iter  20 value 441.888753
    ## iter  30 value 360.388907
    ## iter  40 value 304.634695
    ## iter  50 value 267.271735
    ## iter  60 value 230.820523
    ## iter  70 value 214.770930
    ## iter  80 value 191.519074
    ## iter  90 value 176.426407
    ## iter 100 value 164.876660
    ## final  value 164.876660 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1247.679616 
    ## iter  10 value 543.386414
    ## iter  20 value 224.590844
    ## iter  30 value 143.085263
    ## iter  40 value 105.812560
    ## iter  50 value 94.878584
    ## iter  60 value 91.582388
    ## iter  70 value 88.294860
    ## iter  80 value 80.563353
    ## iter  90 value 70.817879
    ## iter 100 value 68.302267
    ## final  value 68.302267 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1228.696344 
    ## iter  10 value 882.541009
    ## iter  20 value 828.894468
    ## iter  30 value 802.427251
    ## iter  40 value 786.496548
    ## iter  50 value 776.550550
    ## iter  60 value 762.570544
    ## iter  70 value 757.587813
    ## iter  80 value 757.174271
    ## iter  90 value 757.151578
    ## final  value 757.150697 
    ## converged
    ## # weights:  255
    ## initial  value 1297.062292 
    ## iter  10 value 1030.973931
    ## iter  20 value 874.597005
    ## iter  30 value 742.770136
    ## iter  40 value 605.960885
    ## iter  50 value 538.431027
    ## iter  60 value 484.106513
    ## iter  70 value 407.378280
    ## iter  80 value 357.407652
    ## iter  90 value 332.925005
    ## iter 100 value 326.471241
    ## final  value 326.471241 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1353.538465 
    ## iter  10 value 530.907375
    ## iter  20 value 279.343813
    ## iter  30 value 199.049326
    ## iter  40 value 171.789588
    ## iter  50 value 161.626569
    ## iter  60 value 155.074466
    ## iter  70 value 152.646448
    ## iter  80 value 151.362899
    ## iter  90 value 150.807762
    ## iter 100 value 150.522979
    ## final  value 150.522979 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1230.807888 
    ## iter  10 value 944.139016
    ## iter  20 value 807.769716
    ## iter  30 value 750.120513
    ## iter  40 value 713.826720
    ## iter  50 value 651.443992
    ## iter  60 value 618.485345
    ## iter  70 value 609.836790
    ## iter  80 value 606.347986
    ## iter  90 value 605.459943
    ## iter 100 value 602.620698
    ## final  value 602.620698 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1272.146709 
    ## iter  10 value 700.900297
    ## iter  20 value 424.438943
    ## iter  30 value 361.565407
    ## iter  40 value 314.915625
    ## iter  50 value 269.589214
    ## iter  60 value 239.428685
    ## iter  70 value 215.919766
    ## iter  80 value 201.649741
    ## iter  90 value 194.200788
    ## iter 100 value 190.903807
    ## final  value 190.903807 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1277.979264 
    ## iter  10 value 309.746009
    ## iter  20 value 146.965700
    ## iter  30 value 79.181877
    ## iter  40 value 42.986485
    ## iter  50 value 29.979976
    ## iter  60 value 17.372271
    ## iter  70 value 11.270770
    ## iter  80 value 10.108344
    ## iter  90 value 9.577934
    ## iter 100 value 8.904796
    ## final  value 8.904796 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1216.742050 
    ## iter  10 value 1013.102652
    ## iter  20 value 796.779864
    ## iter  30 value 744.485661
    ## iter  40 value 713.700189
    ## iter  50 value 673.301524
    ## iter  60 value 650.129499
    ## iter  70 value 637.160031
    ## iter  80 value 626.166451
    ## iter  90 value 622.073726
    ## iter 100 value 619.792088
    ## final  value 619.792088 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1287.320211 
    ## iter  10 value 616.644952
    ## iter  20 value 394.701033
    ## iter  30 value 321.297567
    ## iter  40 value 272.436181
    ## iter  50 value 238.316003
    ## iter  60 value 204.434973
    ## iter  70 value 182.509756
    ## iter  80 value 168.240790
    ## iter  90 value 166.178471
    ## iter 100 value 164.870226
    ## final  value 164.870226 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1356.804205 
    ## iter  10 value 531.948171
    ## iter  20 value 149.000366
    ## iter  30 value 52.521878
    ## iter  40 value 20.944168
    ## iter  50 value 11.196062
    ## iter  60 value 1.961149
    ## iter  70 value 0.407623
    ## iter  80 value 0.098484
    ## iter  90 value 0.039517
    ## iter 100 value 0.015032
    ## final  value 0.015032 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1259.188673 
    ## iter  10 value 896.161184
    ## iter  20 value 824.833055
    ## iter  30 value 790.530269
    ## iter  40 value 770.374696
    ## iter  50 value 752.429251
    ## iter  60 value 742.031162
    ## iter  70 value 737.610534
    ## iter  80 value 736.506732
    ## iter  90 value 736.241807
    ## iter 100 value 736.221852
    ## final  value 736.221852 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1294.120789 
    ## iter  10 value 762.047362
    ## iter  20 value 623.497643
    ## iter  30 value 527.950203
    ## iter  40 value 468.027727
    ## iter  50 value 378.795552
    ## iter  60 value 340.936498
    ## iter  70 value 319.458313
    ## iter  80 value 306.474535
    ## iter  90 value 300.112259
    ## iter 100 value 296.674055
    ## final  value 296.674055 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1324.622960 
    ## iter  10 value 515.618502
    ## iter  20 value 283.486734
    ## iter  30 value 200.799011
    ## iter  40 value 164.280973
    ## iter  50 value 155.269874
    ## iter  60 value 152.739377
    ## iter  70 value 151.531976
    ## iter  80 value 151.136646
    ## iter  90 value 150.849849
    ## iter 100 value 150.631970
    ## final  value 150.631970 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1285.781616 
    ## iter  10 value 1103.978820
    ## iter  20 value 1098.670773
    ## iter  30 value 1098.649051
    ## iter  40 value 1034.711992
    ## iter  50 value 960.576739
    ## iter  60 value 925.624771
    ## iter  70 value 870.105208
    ## iter  80 value 852.783758
    ## iter  90 value 844.656628
    ## iter 100 value 841.675777
    ## final  value 841.675777 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1304.548127 
    ## iter  10 value 655.650738
    ## iter  20 value 428.205219
    ## iter  30 value 342.202286
    ## iter  40 value 309.153103
    ## iter  50 value 283.055394
    ## iter  60 value 268.640076
    ## iter  70 value 260.890499
    ## iter  80 value 257.601231
    ## iter  90 value 256.731905
    ## iter 100 value 255.579033
    ## final  value 255.579033 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1470.370990 
    ## iter  10 value 1113.997235
    ## iter  20 value 962.783377
    ## iter  30 value 802.718405
    ## iter  40 value 473.025407
    ## iter  50 value 258.969743
    ## iter  60 value 161.073418
    ## iter  70 value 130.636374
    ## iter  80 value 119.347941
    ## iter  90 value 109.888837
    ## iter 100 value 99.082949
    ## final  value 99.082949 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1197.505251 
    ## iter  10 value 831.061428
    ## iter  20 value 784.614171
    ## iter  30 value 737.171642
    ## iter  40 value 710.226063
    ## iter  50 value 688.782764
    ## iter  60 value 685.560897
    ## iter  70 value 685.221694
    ## iter  80 value 685.172787
    ## iter  90 value 685.126077
    ## iter 100 value 685.060121
    ## final  value 685.060121 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1301.091472 
    ## iter  10 value 585.864986
    ## iter  20 value 363.586505
    ## iter  30 value 297.903018
    ## iter  40 value 247.719377
    ## iter  50 value 211.161237
    ## iter  60 value 189.584427
    ## iter  70 value 174.510401
    ## iter  80 value 162.630224
    ## iter  90 value 156.272517
    ## iter 100 value 153.746368
    ## final  value 153.746368 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1248.476133 
    ## iter  10 value 379.825740
    ## iter  20 value 86.640194
    ## iter  30 value 41.404776
    ## iter  40 value 15.476255
    ## iter  50 value 3.373088
    ## iter  60 value 1.983808
    ## iter  70 value 1.920089
    ## iter  80 value 1.911769
    ## iter  90 value 1.910299
    ## iter 100 value 1.909671
    ## final  value 1.909671 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1238.969090 
    ## iter  10 value 956.137474
    ## iter  20 value 850.090132
    ## iter  30 value 809.552332
    ## iter  40 value 778.224352
    ## iter  50 value 765.068040
    ## iter  60 value 760.350988
    ## iter  70 value 757.231236
    ## iter  80 value 753.627375
    ## iter  90 value 750.064571
    ## iter 100 value 749.485379
    ## final  value 749.485379 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1353.052704 
    ## iter  10 value 674.798191
    ## iter  20 value 502.463366
    ## iter  30 value 442.240593
    ## iter  40 value 388.846091
    ## iter  50 value 367.540590
    ## iter  60 value 350.223179
    ## iter  70 value 337.452892
    ## iter  80 value 324.182405
    ## iter  90 value 319.697364
    ## iter 100 value 316.121483
    ## final  value 316.121483 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1312.566816 
    ## iter  10 value 512.225033
    ## iter  20 value 297.307961
    ## iter  30 value 201.174454
    ## iter  40 value 167.451071
    ## iter  50 value 157.000284
    ## iter  60 value 154.961748
    ## iter  70 value 153.661053
    ## iter  80 value 152.061640
    ## iter  90 value 151.236849
    ## iter 100 value 150.912494
    ## final  value 150.912494 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1232.597504 
    ## iter  10 value 870.879767
    ## iter  20 value 800.634310
    ## iter  30 value 760.221892
    ## iter  40 value 724.779781
    ## iter  50 value 703.372329
    ## iter  60 value 699.750804
    ## iter  70 value 699.332184
    ## iter  80 value 699.179340
    ## iter  90 value 698.662465
    ## iter 100 value 698.192353
    ## final  value 698.192353 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1300.604683 
    ## iter  10 value 664.195897
    ## iter  20 value 486.858389
    ## iter  30 value 378.734159
    ## iter  40 value 334.660869
    ## iter  50 value 286.132127
    ## iter  60 value 255.020280
    ## iter  70 value 237.188673
    ## iter  80 value 226.972346
    ## iter  90 value 210.635841
    ## iter 100 value 199.842801
    ## final  value 199.842801 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1255.352604 
    ## iter  10 value 413.206370
    ## iter  20 value 196.818133
    ## iter  30 value 99.874159
    ## iter  40 value 54.537038
    ## iter  50 value 40.533461
    ## iter  60 value 35.088354
    ## iter  70 value 33.264218
    ## iter  80 value 27.529453
    ## iter  90 value 26.513589
    ## iter 100 value 22.642461
    ## final  value 22.642461 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1237.971991 
    ## iter  10 value 975.073817
    ## iter  20 value 857.534545
    ## iter  30 value 843.935853
    ## iter  40 value 832.960079
    ## iter  50 value 827.379539
    ## iter  60 value 826.439590
    ## iter  70 value 825.647100
    ## iter  80 value 821.118209
    ## iter  90 value 820.552875
    ## iter 100 value 820.525570
    ## final  value 820.525570 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1293.900879 
    ## iter  10 value 732.787090
    ## iter  20 value 481.390615
    ## iter  30 value 378.612090
    ## iter  40 value 347.847876
    ## iter  50 value 292.945047
    ## iter  60 value 246.963145
    ## iter  70 value 197.265533
    ## iter  80 value 174.770065
    ## iter  90 value 142.380765
    ## iter 100 value 122.045272
    ## final  value 122.045272 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1437.471025 
    ## iter  10 value 515.875381
    ## iter  20 value 120.164857
    ## iter  30 value 50.264832
    ## iter  40 value 21.117425
    ## iter  50 value 6.626232
    ## iter  60 value 3.154668
    ## iter  70 value 2.822039
    ## iter  80 value 2.779954
    ## iter  90 value 2.775153
    ## iter 100 value 2.773922
    ## final  value 2.773922 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1181.132102 
    ## iter  10 value 908.858577
    ## iter  20 value 841.365549
    ## iter  30 value 801.948544
    ## iter  40 value 782.104508
    ## iter  50 value 774.661728
    ## iter  60 value 768.020705
    ## iter  70 value 766.197332
    ## iter  80 value 765.006903
    ## iter  90 value 764.932474
    ## final  value 764.932108 
    ## converged
    ## # weights:  255
    ## initial  value 1250.255858 
    ## iter  10 value 793.129360
    ## iter  20 value 668.450823
    ## iter  30 value 587.286756
    ## iter  40 value 532.233718
    ## iter  50 value 464.858439
    ## iter  60 value 418.640592
    ## iter  70 value 395.750194
    ## iter  80 value 366.176123
    ## iter  90 value 348.099268
    ## iter 100 value 331.668334
    ## final  value 331.668334 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1270.811599 
    ## iter  10 value 470.798153
    ## iter  20 value 284.566187
    ## iter  30 value 204.363506
    ## iter  40 value 171.370815
    ## iter  50 value 164.445336
    ## iter  60 value 161.547023
    ## iter  70 value 159.842341
    ## iter  80 value 158.792628
    ## iter  90 value 157.493918
    ## iter 100 value 156.018910
    ## final  value 156.018910 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1259.472055 
    ## iter  10 value 919.027379
    ## iter  20 value 811.706970
    ## iter  30 value 760.385111
    ## iter  40 value 714.488979
    ## iter  50 value 674.495558
    ## iter  60 value 639.666387
    ## iter  70 value 633.134521
    ## iter  80 value 626.095667
    ## iter  90 value 622.519309
    ## iter 100 value 621.962326
    ## final  value 621.962326 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1242.934783 
    ## iter  10 value 661.244275
    ## iter  20 value 349.294114
    ## iter  30 value 302.184892
    ## iter  40 value 282.277471
    ## iter  50 value 265.924590
    ## iter  60 value 257.829828
    ## iter  70 value 251.962111
    ## iter  80 value 243.554369
    ## iter  90 value 236.731306
    ## iter 100 value 234.398048
    ## final  value 234.398048 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1332.031649 
    ## iter  10 value 549.989510
    ## iter  20 value 233.059389
    ## iter  30 value 120.904553
    ## iter  40 value 81.216859
    ## iter  50 value 61.309905
    ## iter  60 value 55.400724
    ## iter  70 value 51.879475
    ## iter  80 value 45.769833
    ## iter  90 value 41.625113
    ## iter 100 value 40.502769
    ## final  value 40.502769 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1313.087766 
    ## iter  10 value 1050.743250
    ## iter  20 value 847.229565
    ## iter  30 value 810.682293
    ## iter  40 value 786.176183
    ## iter  50 value 768.301866
    ## iter  60 value 765.378808
    ## iter  70 value 761.880626
    ## iter  80 value 749.668035
    ## iter  90 value 744.356768
    ## iter 100 value 741.174943
    ## final  value 741.174943 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1265.033096 
    ## iter  10 value 637.307412
    ## iter  20 value 374.730043
    ## iter  30 value 267.480646
    ## iter  40 value 204.106655
    ## iter  50 value 168.288071
    ## iter  60 value 144.640285
    ## iter  70 value 132.579859
    ## iter  80 value 121.848913
    ## iter  90 value 108.195620
    ## iter 100 value 95.921437
    ## final  value 95.921437 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1379.521924 
    ## iter  10 value 558.404419
    ## iter  20 value 198.467889
    ## iter  30 value 120.760806
    ## iter  40 value 93.722650
    ## iter  50 value 66.338014
    ## iter  60 value 54.965165
    ## iter  70 value 52.673615
    ## iter  80 value 48.498770
    ## iter  90 value 43.102881
    ## iter 100 value 41.030037
    ## final  value 41.030037 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1311.288669 
    ## iter  10 value 995.449563
    ## iter  20 value 867.757872
    ## iter  30 value 834.092985
    ## iter  40 value 821.806870
    ## iter  50 value 806.969932
    ## iter  60 value 798.425764
    ## iter  70 value 796.375453
    ## iter  80 value 794.938987
    ## iter  90 value 789.656618
    ## iter 100 value 784.465179
    ## final  value 784.465179 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1284.281507 
    ## iter  10 value 683.724742
    ## iter  20 value 481.468864
    ## iter  30 value 404.684593
    ## iter  40 value 363.800878
    ## iter  50 value 333.510249
    ## iter  60 value 324.933640
    ## iter  70 value 319.558501
    ## iter  80 value 317.529897
    ## iter  90 value 313.329047
    ## iter 100 value 309.669610
    ## final  value 309.669610 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1299.523852 
    ## iter  10 value 677.457751
    ## iter  20 value 363.341579
    ## iter  30 value 241.513490
    ## iter  40 value 200.668694
    ## iter  50 value 174.950117
    ## iter  60 value 166.217581
    ## iter  70 value 157.380648
    ## iter  80 value 154.475610
    ## iter  90 value 154.033363
    ## iter 100 value 153.815900
    ## final  value 153.815900 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1296.675934 
    ## iter  10 value 977.226254
    ## iter  20 value 875.878463
    ## iter  30 value 824.228309
    ## iter  40 value 777.123069
    ## iter  50 value 757.666249
    ## iter  60 value 737.500382
    ## iter  70 value 727.256480
    ## iter  80 value 720.600018
    ## iter  90 value 719.760390
    ## iter 100 value 718.974998
    ## final  value 718.974998 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1285.669932 
    ## iter  10 value 843.086052
    ## iter  20 value 608.015058
    ## iter  30 value 512.975116
    ## iter  40 value 469.092129
    ## iter  50 value 417.733370
    ## iter  60 value 377.510372
    ## iter  70 value 360.350849
    ## iter  80 value 349.728446
    ## iter  90 value 345.845864
    ## iter 100 value 342.712745
    ## final  value 342.712745 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1309.202654 
    ## iter  10 value 474.849018
    ## iter  20 value 161.036838
    ## iter  30 value 73.996803
    ## iter  40 value 45.872042
    ## iter  50 value 32.459903
    ## iter  60 value 24.211569
    ## iter  70 value 21.858017
    ## iter  80 value 21.221854
    ## iter  90 value 13.212489
    ## iter 100 value 9.801844
    ## final  value 9.801844 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1347.376412 
    ## iter  10 value 953.034838
    ## iter  20 value 824.067002
    ## iter  30 value 708.064415
    ## iter  40 value 646.390562
    ## iter  50 value 592.987975
    ## iter  60 value 568.048420
    ## iter  70 value 564.246929
    ## iter  80 value 560.366210
    ## iter  90 value 559.599977
    ## iter 100 value 558.940472
    ## final  value 558.940472 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1313.032034 
    ## iter  10 value 669.693482
    ## iter  20 value 477.342422
    ## iter  30 value 391.290981
    ## iter  40 value 338.700118
    ## iter  50 value 309.211267
    ## iter  60 value 284.885570
    ## iter  70 value 274.384419
    ## iter  80 value 267.751460
    ## iter  90 value 263.839489
    ## iter 100 value 260.526199
    ## final  value 260.526199 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1290.762415 
    ## iter  10 value 688.605330
    ## iter  20 value 338.947384
    ## iter  30 value 254.648050
    ## iter  40 value 175.018585
    ## iter  50 value 129.833631
    ## iter  60 value 107.259733
    ## iter  70 value 85.073321
    ## iter  80 value 75.207463
    ## iter  90 value 70.264373
    ## iter 100 value 65.992570
    ## final  value 65.992570 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1254.842324 
    ## iter  10 value 918.573283
    ## iter  20 value 834.121926
    ## iter  30 value 808.748791
    ## iter  40 value 784.613092
    ## iter  50 value 775.285881
    ## iter  60 value 772.332689
    ## iter  70 value 768.176159
    ## iter  80 value 758.012243
    ## iter  90 value 748.361728
    ## iter 100 value 747.222571
    ## final  value 747.222571 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1357.067625 
    ## iter  10 value 1062.743644
    ## iter  20 value 869.644480
    ## iter  30 value 787.514033
    ## iter  40 value 611.506598
    ## iter  50 value 497.437348
    ## iter  60 value 399.078468
    ## iter  70 value 355.609466
    ## iter  80 value 328.504295
    ## iter  90 value 318.076173
    ## iter 100 value 309.662255
    ## final  value 309.662255 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1215.277879 
    ## iter  10 value 457.035124
    ## iter  20 value 306.225473
    ## iter  30 value 221.058776
    ## iter  40 value 168.287166
    ## iter  50 value 160.513361
    ## iter  60 value 156.965383
    ## iter  70 value 154.460582
    ## iter  80 value 153.334398
    ## iter  90 value 152.985708
    ## iter 100 value 152.719782
    ## final  value 152.719782 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1224.385841 
    ## iter  10 value 953.856049
    ## iter  20 value 808.988540
    ## iter  30 value 790.394879
    ## iter  40 value 769.996622
    ## iter  50 value 761.357889
    ## iter  60 value 757.221839
    ## iter  70 value 748.529278
    ## iter  80 value 737.282544
    ## iter  90 value 736.877917
    ## iter 100 value 736.738442
    ## final  value 736.738442 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1322.660792 
    ## iter  10 value 645.677925
    ## iter  20 value 437.206710
    ## iter  30 value 398.797922
    ## iter  40 value 340.343457
    ## iter  50 value 277.355631
    ## iter  60 value 199.318298
    ## iter  70 value 165.658794
    ## iter  80 value 141.285115
    ## iter  90 value 124.682051
    ## iter 100 value 118.239766
    ## final  value 118.239766 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1226.891934 
    ## iter  10 value 541.176651
    ## iter  20 value 328.061089
    ## iter  30 value 239.823362
    ## iter  40 value 164.637393
    ## iter  50 value 78.902974
    ## iter  60 value 50.076551
    ## iter  70 value 38.547733
    ## iter  80 value 30.817388
    ## iter  90 value 27.707895
    ## iter 100 value 26.238381
    ## final  value 26.238381 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1235.949450 
    ## iter  10 value 903.397754
    ## iter  20 value 778.733064
    ## iter  30 value 745.881793
    ## iter  40 value 713.261737
    ## iter  50 value 685.652866
    ## iter  60 value 665.660723
    ## iter  70 value 646.520240
    ## iter  80 value 636.230734
    ## iter  90 value 633.629547
    ## iter 100 value 632.100809
    ## final  value 632.100809 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1350.598157 
    ## iter  10 value 874.105685
    ## iter  20 value 754.662160
    ## iter  30 value 677.658882
    ## iter  40 value 620.047038
    ## iter  50 value 587.519769
    ## iter  60 value 563.369845
    ## iter  70 value 543.966388
    ## iter  80 value 487.101423
    ## iter  90 value 455.513808
    ## iter 100 value 416.837665
    ## final  value 416.837665 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1278.020147 
    ## iter  10 value 411.092321
    ## iter  20 value 189.973438
    ## iter  30 value 124.520499
    ## iter  40 value 86.201181
    ## iter  50 value 62.411106
    ## iter  60 value 46.533368
    ## iter  70 value 34.566157
    ## iter  80 value 26.252043
    ## iter  90 value 22.681590
    ## iter 100 value 19.591162
    ## final  value 19.591162 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1271.326246 
    ## iter  10 value 931.461318
    ## iter  20 value 845.086542
    ## iter  30 value 817.747322
    ## iter  40 value 803.557712
    ## iter  50 value 792.612011
    ## iter  60 value 774.839776
    ## iter  70 value 757.215389
    ## iter  80 value 755.142170
    ## iter  90 value 752.151265
    ## iter 100 value 745.426776
    ## final  value 745.426776 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1227.069478 
    ## iter  10 value 744.904317
    ## iter  20 value 511.347624
    ## iter  30 value 404.035305
    ## iter  40 value 358.941979
    ## iter  50 value 341.111964
    ## iter  60 value 331.162630
    ## iter  70 value 326.533558
    ## iter  80 value 318.650609
    ## iter  90 value 303.392930
    ## iter 100 value 292.523126
    ## final  value 292.523126 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1306.051612 
    ## iter  10 value 563.499161
    ## iter  20 value 279.688606
    ## iter  30 value 197.108232
    ## iter  40 value 157.643869
    ## iter  50 value 149.104343
    ## iter  60 value 145.400360
    ## iter  70 value 141.653156
    ## iter  80 value 140.814506
    ## iter  90 value 140.570158
    ## iter 100 value 140.478851
    ## final  value 140.478851 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1308.413456 
    ## iter  10 value 1016.715918
    ## iter  20 value 866.556835
    ## iter  30 value 820.193299
    ## iter  40 value 780.877270
    ## iter  50 value 755.005556
    ## iter  60 value 742.006268
    ## iter  70 value 738.390016
    ## iter  80 value 737.006387
    ## iter  90 value 736.329276
    ## iter 100 value 735.805693
    ## final  value 735.805693 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1281.726167 
    ## iter  10 value 636.603715
    ## iter  20 value 345.102317
    ## iter  30 value 282.256259
    ## iter  40 value 237.020276
    ## iter  50 value 202.129713
    ## iter  60 value 166.873755
    ## iter  70 value 143.130210
    ## iter  80 value 128.492250
    ## iter  90 value 125.743733
    ## iter 100 value 123.864503
    ## final  value 123.864503 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1328.233240 
    ## iter  10 value 480.412904
    ## iter  20 value 112.483774
    ## iter  30 value 58.597188
    ## iter  40 value 36.643595
    ## iter  50 value 16.352712
    ## iter  60 value 11.421017
    ## iter  70 value 9.038157
    ## iter  80 value 5.844734
    ## iter  90 value 4.661764
    ## iter 100 value 4.203290
    ## final  value 4.203290 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1264.864912 
    ## iter  10 value 864.149801
    ## iter  20 value 796.721158
    ## iter  30 value 766.299029
    ## iter  40 value 730.367234
    ## iter  50 value 711.732718
    ## iter  60 value 708.082270
    ## iter  70 value 706.365200
    ## iter  80 value 704.494312
    ## iter  90 value 702.027267
    ## iter 100 value 701.314905
    ## final  value 701.314905 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1354.336859 
    ## iter  10 value 650.806229
    ## iter  20 value 425.958889
    ## iter  30 value 345.361517
    ## iter  40 value 290.045791
    ## iter  50 value 249.263417
    ## iter  60 value 213.687871
    ## iter  70 value 176.490182
    ## iter  80 value 136.039797
    ## iter  90 value 109.597904
    ## iter 100 value 92.949629
    ## final  value 92.949629 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1250.400445 
    ## iter  10 value 487.959570
    ## iter  20 value 194.025088
    ## iter  30 value 157.742282
    ## iter  40 value 133.287123
    ## iter  50 value 109.612313
    ## iter  60 value 94.333641
    ## iter  70 value 75.602213
    ## iter  80 value 60.413167
    ## iter  90 value 46.210376
    ## iter 100 value 33.390966
    ## final  value 33.390966 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1321.758745 
    ## iter  10 value 1144.007325
    ## iter  20 value 1125.914511
    ## iter  30 value 980.912488
    ## iter  40 value 885.553589
    ## iter  50 value 843.935823
    ## iter  60 value 828.368278
    ## iter  70 value 820.858040
    ## iter  80 value 800.749338
    ## iter  90 value 787.113836
    ## iter 100 value 782.201634
    ## final  value 782.201634 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1370.336100 
    ## iter  10 value 693.932637
    ## iter  20 value 523.504275
    ## iter  30 value 453.024203
    ## iter  40 value 415.530641
    ## iter  50 value 404.101605
    ## iter  60 value 388.486319
    ## iter  70 value 365.527570
    ## iter  80 value 338.737507
    ## iter  90 value 317.567216
    ## iter 100 value 311.622407
    ## final  value 311.622407 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1338.538874 
    ## iter  10 value 614.041457
    ## iter  20 value 324.318901
    ## iter  30 value 231.114899
    ## iter  40 value 186.707675
    ## iter  50 value 166.378248
    ## iter  60 value 161.185744
    ## iter  70 value 155.773209
    ## iter  80 value 151.930479
    ## iter  90 value 148.420010
    ## iter 100 value 146.343998
    ## final  value 146.343998 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1274.351868 
    ## iter  10 value 843.265779
    ## iter  20 value 766.905933
    ## iter  30 value 706.017873
    ## iter  40 value 658.284375
    ## iter  50 value 615.927701
    ## iter  60 value 581.770322
    ## iter  70 value 548.085598
    ## iter  80 value 542.905614
    ## iter  90 value 533.993840
    ## iter 100 value 521.717105
    ## final  value 521.717105 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1233.086351 
    ## iter  10 value 536.840532
    ## iter  20 value 297.738416
    ## iter  30 value 225.117161
    ## iter  40 value 191.382667
    ## iter  50 value 172.295146
    ## iter  60 value 147.722633
    ## iter  70 value 130.258174
    ## iter  80 value 125.364129
    ## iter  90 value 123.537635
    ## iter 100 value 121.732686
    ## final  value 121.732686 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1271.110521 
    ## iter  10 value 528.017366
    ## iter  20 value 282.159870
    ## iter  30 value 192.687458
    ## iter  40 value 150.121412
    ## iter  50 value 129.268966
    ## iter  60 value 117.512806
    ## iter  70 value 97.987756
    ## iter  80 value 94.934576
    ## iter  90 value 92.334575
    ## iter 100 value 90.250503
    ## final  value 90.250503 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1272.586928 
    ## iter  10 value 1064.587986
    ## iter  20 value 935.801508
    ## iter  30 value 917.297734
    ## iter  40 value 900.487786
    ## iter  50 value 889.311600
    ## iter  60 value 878.837036
    ## iter  70 value 876.338527
    ## iter  80 value 874.356983
    ## iter  90 value 873.515377
    ## iter 100 value 873.248727
    ## final  value 873.248727 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1214.567615 
    ## iter  10 value 665.911638
    ## iter  20 value 368.705658
    ## iter  30 value 293.999827
    ## iter  40 value 233.010890
    ## iter  50 value 175.849921
    ## iter  60 value 147.406470
    ## iter  70 value 131.262913
    ## iter  80 value 120.158820
    ## iter  90 value 115.025622
    ## iter 100 value 104.904570
    ## final  value 104.904570 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1338.062516 
    ## iter  10 value 421.119671
    ## iter  20 value 94.035383
    ## iter  30 value 17.947418
    ## iter  40 value 4.503476
    ## iter  50 value 2.737116
    ## iter  60 value 2.522076
    ## iter  70 value 2.506051
    ## iter  80 value 2.503383
    ## iter  90 value 2.502441
    ## iter 100 value 2.502275
    ## final  value 2.502275 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1249.138773 
    ## iter  10 value 974.134892
    ## iter  20 value 856.704595
    ## iter  30 value 815.083596
    ## iter  40 value 795.985370
    ## iter  50 value 788.475186
    ## iter  60 value 771.710115
    ## iter  70 value 766.553453
    ## iter  80 value 765.913036
    ## iter  90 value 765.851120
    ## final  value 765.850511 
    ## converged
    ## # weights:  255
    ## initial  value 1249.135861 
    ## iter  10 value 583.529613
    ## iter  20 value 469.080268
    ## iter  30 value 375.238114
    ## iter  40 value 339.925410
    ## iter  50 value 316.529343
    ## iter  60 value 305.521154
    ## iter  70 value 299.542552
    ## iter  80 value 296.970027
    ## iter  90 value 294.318558
    ## iter 100 value 293.576319
    ## final  value 293.576319 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1307.764142 
    ## iter  10 value 489.162660
    ## iter  20 value 317.262058
    ## iter  30 value 231.133478
    ## iter  40 value 171.836566
    ## iter  50 value 161.204029
    ## iter  60 value 158.207304
    ## iter  70 value 157.035236
    ## iter  80 value 156.105835
    ## iter  90 value 153.249029
    ## iter 100 value 152.684796
    ## final  value 152.684796 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1310.590354 
    ## iter  10 value 1079.188563
    ## iter  20 value 936.292804
    ## iter  30 value 786.746459
    ## iter  40 value 732.888540
    ## iter  50 value 675.848607
    ## iter  60 value 649.958707
    ## iter  70 value 611.333152
    ## iter  80 value 598.165122
    ## iter  90 value 595.636789
    ## iter 100 value 592.867937
    ## final  value 592.867937 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1299.535493 
    ## iter  10 value 620.047009
    ## iter  20 value 384.146575
    ## iter  30 value 302.565645
    ## iter  40 value 265.499404
    ## iter  50 value 229.355302
    ## iter  60 value 189.991232
    ## iter  70 value 171.457568
    ## iter  80 value 157.247390
    ## iter  90 value 152.263258
    ## iter 100 value 147.416940
    ## final  value 147.416940 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1353.211166 
    ## iter  10 value 460.297277
    ## iter  20 value 162.029437
    ## iter  30 value 83.382526
    ## iter  40 value 48.404984
    ## iter  50 value 31.647905
    ## iter  60 value 22.133214
    ## iter  70 value 19.781795
    ## iter  80 value 16.102820
    ## iter  90 value 11.640307
    ## iter 100 value 7.318062
    ## final  value 7.318062 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1308.911624 
    ## iter  10 value 966.879078
    ## iter  20 value 867.923002
    ## iter  30 value 831.796159
    ## iter  40 value 787.120355
    ## iter  50 value 754.573452
    ## iter  60 value 740.933191
    ## iter  70 value 728.824438
    ## iter  80 value 711.744913
    ## iter  90 value 697.203586
    ## iter 100 value 685.440799
    ## final  value 685.440799 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1216.282754 
    ## iter  10 value 779.192708
    ## iter  20 value 350.419430
    ## iter  30 value 237.142638
    ## iter  40 value 210.284614
    ## iter  50 value 187.153659
    ## iter  60 value 174.336085
    ## iter  70 value 161.615629
    ## iter  80 value 156.656588
    ## iter  90 value 155.617425
    ## iter 100 value 154.130980
    ## final  value 154.130980 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1442.537157 
    ## iter  10 value 408.597904
    ## iter  20 value 142.872039
    ## iter  30 value 65.102317
    ## iter  40 value 44.203731
    ## iter  50 value 33.151368
    ## iter  60 value 22.005759
    ## iter  70 value 18.521381
    ## iter  80 value 15.034954
    ## iter  90 value 14.435138
    ## iter 100 value 13.940625
    ## final  value 13.940625 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1251.227649 
    ## iter  10 value 940.812491
    ## iter  20 value 885.945520
    ## iter  30 value 852.720036
    ## iter  40 value 826.359118
    ## iter  50 value 803.178341
    ## iter  60 value 783.196909
    ## iter  70 value 775.147892
    ## iter  80 value 773.794027
    ## iter  90 value 773.438318
    ## iter 100 value 773.408539
    ## final  value 773.408539 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1225.929134 
    ## iter  10 value 642.892811
    ## iter  20 value 498.784261
    ## iter  30 value 426.357090
    ## iter  40 value 390.407844
    ## iter  50 value 372.759738
    ## iter  60 value 354.978436
    ## iter  70 value 344.150053
    ## iter  80 value 340.265589
    ## iter  90 value 331.798954
    ## iter 100 value 320.721211
    ## final  value 320.721211 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1332.218693 
    ## iter  10 value 649.723498
    ## iter  20 value 372.253419
    ## iter  30 value 257.525188
    ## iter  40 value 189.994468
    ## iter  50 value 163.217430
    ## iter  60 value 151.371586
    ## iter  70 value 147.458231
    ## iter  80 value 145.501668
    ## iter  90 value 144.546548
    ## iter 100 value 143.625569
    ## final  value 143.625569 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1273.432950 
    ## iter  10 value 1122.846677
    ## iter  20 value 1116.190347
    ## iter  30 value 1047.875384
    ## iter  40 value 1016.859357
    ## iter  50 value 1001.148504
    ## iter  60 value 967.630566
    ## iter  70 value 921.992624
    ## iter  80 value 848.092769
    ## iter  90 value 809.037297
    ## iter 100 value 778.087786
    ## final  value 778.087786 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1256.046762 
    ## iter  10 value 550.942084
    ## iter  20 value 375.919833
    ## iter  30 value 290.095208
    ## iter  40 value 251.812591
    ## iter  50 value 217.246048
    ## iter  60 value 188.487798
    ## iter  70 value 155.424541
    ## iter  80 value 141.393365
    ## iter  90 value 136.265693
    ## iter 100 value 132.811159
    ## final  value 132.811159 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1314.670641 
    ## iter  10 value 360.032504
    ## iter  20 value 147.964623
    ## iter  30 value 53.499566
    ## iter  40 value 27.403280
    ## iter  50 value 8.009378
    ## iter  60 value 4.040172
    ## iter  70 value 3.472968
    ## iter  80 value 2.989815
    ## iter  90 value 2.634951
    ## iter 100 value 2.378428
    ## final  value 2.378428 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1246.860299 
    ## iter  10 value 1060.372810
    ## iter  20 value 916.363091
    ## iter  30 value 858.516481
    ## iter  40 value 823.634195
    ## iter  50 value 800.587186
    ## iter  60 value 780.030110
    ## iter  70 value 769.321601
    ## iter  80 value 764.152293
    ## iter  90 value 759.112873
    ## iter 100 value 754.575565
    ## final  value 754.575565 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1359.455739 
    ## iter  10 value 639.683402
    ## iter  20 value 313.221308
    ## iter  30 value 244.234149
    ## iter  40 value 190.740974
    ## iter  50 value 152.157089
    ## iter  60 value 133.527821
    ## iter  70 value 126.712200
    ## iter  80 value 124.934685
    ## iter  90 value 123.515577
    ## iter 100 value 121.140267
    ## final  value 121.140267 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1332.915566 
    ## iter  10 value 517.846010
    ## iter  20 value 271.497841
    ## iter  30 value 166.253969
    ## iter  40 value 115.334074
    ## iter  50 value 93.591990
    ## iter  60 value 87.187706
    ## iter  70 value 85.381527
    ## iter  80 value 82.525501
    ## iter  90 value 80.584123
    ## iter 100 value 79.673979
    ## final  value 79.673979 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1236.936972 
    ## iter  10 value 1081.266353
    ## iter  20 value 943.002287
    ## iter  30 value 870.591359
    ## iter  40 value 838.031141
    ## iter  50 value 827.461551
    ## iter  60 value 822.162056
    ## iter  70 value 820.035374
    ## iter  80 value 814.032254
    ## iter  90 value 800.829738
    ## iter 100 value 798.724169
    ## final  value 798.724169 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1251.438763 
    ## iter  10 value 803.631245
    ## iter  20 value 604.466796
    ## iter  30 value 472.719829
    ## iter  40 value 395.383860
    ## iter  50 value 362.545041
    ## iter  60 value 346.232690
    ## iter  70 value 337.927531
    ## iter  80 value 333.265909
    ## iter  90 value 329.873117
    ## iter 100 value 326.160906
    ## final  value 326.160906 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1319.117518 
    ## iter  10 value 526.487503
    ## iter  20 value 298.113100
    ## iter  30 value 212.215145
    ## iter  40 value 187.728362
    ## iter  50 value 172.991038
    ## iter  60 value 165.083779
    ## iter  70 value 161.693177
    ## iter  80 value 157.238826
    ## iter  90 value 155.549460
    ## iter 100 value 155.125233
    ## final  value 155.125233 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1246.109608 
    ## iter  10 value 1013.320898
    ## iter  20 value 926.075164
    ## iter  30 value 871.069501
    ## iter  40 value 835.874972
    ## iter  50 value 810.281739
    ## iter  60 value 805.851835
    ## iter  70 value 805.317924
    ## iter  80 value 804.542941
    ## iter  90 value 804.044337
    ## iter 100 value 803.427680
    ## final  value 803.427680 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1420.479738 
    ## iter  10 value 1054.022285
    ## iter  20 value 744.253887
    ## iter  30 value 508.948477
    ## iter  40 value 388.315429
    ## iter  50 value 300.768030
    ## iter  60 value 236.742982
    ## iter  70 value 210.700750
    ## iter  80 value 198.573930
    ## iter  90 value 179.244529
    ## iter 100 value 169.408562
    ## final  value 169.408562 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1457.204372 
    ## iter  10 value 371.435951
    ## iter  20 value 111.195998
    ## iter  30 value 55.520078
    ## iter  40 value 31.015290
    ## iter  50 value 12.500449
    ## iter  60 value 7.724749
    ## iter  70 value 2.952317
    ## iter  80 value 2.332961
    ## iter  90 value 2.146246
    ## iter 100 value 1.976882
    ## final  value 1.976882 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1250.938204 
    ## iter  10 value 957.867454
    ## iter  20 value 802.802945
    ## iter  30 value 755.253865
    ## iter  40 value 717.399998
    ## iter  50 value 686.899362
    ## iter  60 value 682.326684
    ## iter  70 value 682.292906
    ## iter  80 value 682.289093
    ## iter  90 value 682.272114
    ## iter 100 value 682.237438
    ## final  value 682.237438 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1257.141856 
    ## iter  10 value 723.573445
    ## iter  20 value 364.185762
    ## iter  30 value 300.338896
    ## iter  40 value 259.685746
    ## iter  50 value 200.962334
    ## iter  60 value 171.805623
    ## iter  70 value 147.791346
    ## iter  80 value 129.799449
    ## iter  90 value 115.888759
    ## iter 100 value 107.662897
    ## final  value 107.662897 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1340.584151 
    ## iter  10 value 435.912376
    ## iter  20 value 87.568112
    ## iter  30 value 25.745937
    ## iter  40 value 15.421904
    ## iter  50 value 11.592401
    ## iter  60 value 4.296303
    ## iter  70 value 3.852894
    ## iter  80 value 3.819061
    ## iter  90 value 3.739371
    ## iter 100 value 3.049997
    ## final  value 3.049997 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1266.258652 
    ## iter  10 value 1050.804541
    ## iter  20 value 874.490797
    ## iter  30 value 844.103496
    ## iter  40 value 821.442048
    ## iter  50 value 810.683596
    ## iter  60 value 792.662500
    ## iter  70 value 785.803012
    ## iter  80 value 782.722625
    ## iter  90 value 781.633522
    ## iter 100 value 781.571968
    ## final  value 781.571968 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1276.615122 
    ## iter  10 value 812.024956
    ## iter  20 value 618.196521
    ## iter  30 value 460.226904
    ## iter  40 value 412.987662
    ## iter  50 value 389.515694
    ## iter  60 value 358.572638
    ## iter  70 value 343.557062
    ## iter  80 value 334.555673
    ## iter  90 value 329.772794
    ## iter 100 value 326.083893
    ## final  value 326.083893 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1230.707014 
    ## iter  10 value 552.431092
    ## iter  20 value 315.971067
    ## iter  30 value 222.711401
    ## iter  40 value 179.162445
    ## iter  50 value 166.915242
    ## iter  60 value 160.960754
    ## iter  70 value 158.744143
    ## iter  80 value 158.202089
    ## iter  90 value 157.981082
    ## iter 100 value 157.897015
    ## final  value 157.897015 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1256.118573 
    ## iter  10 value 1002.170353
    ## iter  20 value 879.269041
    ## iter  30 value 823.495849
    ## iter  40 value 780.885503
    ## iter  50 value 748.370082
    ## iter  60 value 713.331700
    ## iter  70 value 692.361701
    ## iter  80 value 687.823478
    ## iter  90 value 677.171417
    ## iter 100 value 670.463518
    ## final  value 670.463518 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1274.323520 
    ## iter  10 value 753.988074
    ## iter  20 value 456.267427
    ## iter  30 value 368.870877
    ## iter  40 value 323.555245
    ## iter  50 value 266.554620
    ## iter  60 value 228.831943
    ## iter  70 value 203.193878
    ## iter  80 value 194.809856
    ## iter  90 value 191.934868
    ## iter 100 value 189.692900
    ## final  value 189.692900 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1312.012013 
    ## iter  10 value 327.202353
    ## iter  20 value 67.294649
    ## iter  30 value 22.759206
    ## iter  40 value 7.355195
    ## iter  50 value 3.839895
    ## iter  60 value 3.651557
    ## iter  70 value 3.499744
    ## iter  80 value 3.409251
    ## iter  90 value 3.337758
    ## iter 100 value 3.283133
    ## final  value 3.283133 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1276.140527 
    ## iter  10 value 831.607023
    ## iter  20 value 789.571816
    ## iter  30 value 749.065747
    ## iter  40 value 671.414789
    ## iter  50 value 618.320530
    ## iter  60 value 566.381406
    ## iter  70 value 505.191745
    ## iter  80 value 477.309123
    ## iter  90 value 467.434093
    ## iter 100 value 465.548555
    ## final  value 465.548555 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1224.266856 
    ## iter  10 value 571.051618
    ## iter  20 value 341.390437
    ## iter  30 value 227.377950
    ## iter  40 value 174.946243
    ## iter  50 value 148.110942
    ## iter  60 value 129.523612
    ## iter  70 value 116.742924
    ## iter  80 value 102.972814
    ## iter  90 value 86.640307
    ## iter 100 value 77.726991
    ## final  value 77.726991 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1259.781565 
    ## iter  10 value 321.227322
    ## iter  20 value 150.875599
    ## iter  30 value 87.752303
    ## iter  40 value 68.745306
    ## iter  50 value 51.766930
    ## iter  60 value 41.640128
    ## iter  70 value 29.522621
    ## iter  80 value 21.063984
    ## iter  90 value 17.684691
    ## iter 100 value 16.188647
    ## final  value 16.188647 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1196.891273 
    ## iter  10 value 897.619941
    ## iter  20 value 847.525436
    ## iter  30 value 818.965684
    ## iter  40 value 790.123073
    ## iter  50 value 768.716315
    ## iter  60 value 748.804418
    ## iter  70 value 741.623515
    ## iter  80 value 739.121950
    ## iter  90 value 738.743530
    ## iter 100 value 738.688891
    ## final  value 738.688891 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1241.944022 
    ## iter  10 value 743.674821
    ## iter  20 value 460.096577
    ## iter  30 value 388.927079
    ## iter  40 value 354.484094
    ## iter  50 value 332.627216
    ## iter  60 value 323.211081
    ## iter  70 value 315.340080
    ## iter  80 value 309.949581
    ## iter  90 value 307.041940
    ## iter 100 value 305.393467
    ## final  value 305.393467 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1249.351404 
    ## iter  10 value 548.356735
    ## iter  20 value 305.968635
    ## iter  30 value 225.297638
    ## iter  40 value 179.159566
    ## iter  50 value 171.561686
    ## iter  60 value 163.410737
    ## iter  70 value 160.869388
    ## iter  80 value 159.882470
    ## iter  90 value 158.922338
    ## iter 100 value 156.693963
    ## final  value 156.693963 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1260.819716 
    ## iter  10 value 885.015583
    ## iter  20 value 776.228164
    ## iter  30 value 728.898491
    ## iter  40 value 700.862502
    ## iter  50 value 681.936725
    ## iter  60 value 657.086880
    ## iter  70 value 655.464019
    ## iter  80 value 655.084270
    ## iter  90 value 654.745868
    ## iter 100 value 654.021904
    ## final  value 654.021904 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1239.133977 
    ## iter  10 value 615.068616
    ## iter  20 value 331.177695
    ## iter  30 value 251.088489
    ## iter  40 value 194.247669
    ## iter  50 value 155.272359
    ## iter  60 value 124.726231
    ## iter  70 value 104.478420
    ## iter  80 value 93.066477
    ## iter  90 value 90.856324
    ## iter 100 value 89.370047
    ## final  value 89.370047 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1310.361503 
    ## iter  10 value 438.031540
    ## iter  20 value 152.020900
    ## iter  30 value 65.969745
    ## iter  40 value 35.681182
    ## iter  50 value 20.537675
    ## iter  60 value 13.175057
    ## iter  70 value 12.286707
    ## iter  80 value 11.035052
    ## iter  90 value 8.505689
    ## iter 100 value 6.910809
    ## final  value 6.910809 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1284.034496 
    ## iter  10 value 948.109091
    ## iter  20 value 835.797854
    ## iter  30 value 780.572095
    ## iter  40 value 726.822593
    ## iter  50 value 702.409531
    ## iter  60 value 683.838444
    ## iter  70 value 675.309858
    ## iter  80 value 667.088928
    ## iter  90 value 658.790944
    ## iter 100 value 655.797067
    ## final  value 655.797067 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1284.268545 
    ## iter  10 value 635.279462
    ## iter  20 value 368.646777
    ## iter  30 value 272.410132
    ## iter  40 value 205.855520
    ## iter  50 value 166.297991
    ## iter  60 value 140.285190
    ## iter  70 value 120.005848
    ## iter  80 value 98.433134
    ## iter  90 value 87.052110
    ## iter 100 value 81.419722
    ## final  value 81.419722 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1470.205967 
    ## iter  10 value 733.204476
    ## iter  20 value 290.033998
    ## iter  30 value 175.677824
    ## iter  40 value 127.524539
    ## iter  50 value 71.625646
    ## iter  60 value 50.472556
    ## iter  70 value 35.515212
    ## iter  80 value 24.531981
    ## iter  90 value 17.991116
    ## iter 100 value 16.055989
    ## final  value 16.055989 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1253.873200 
    ## iter  10 value 935.829718
    ## iter  20 value 854.944001
    ## iter  30 value 824.171198
    ## iter  40 value 808.992437
    ## iter  50 value 802.537629
    ## iter  60 value 792.037868
    ## iter  70 value 779.740685
    ## iter  80 value 772.780477
    ## iter  90 value 772.273346
    ## iter 100 value 771.834829
    ## final  value 771.834829 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1264.748073 
    ## iter  10 value 779.867023
    ## iter  20 value 585.822941
    ## iter  30 value 513.597333
    ## iter  40 value 458.443165
    ## iter  50 value 424.554502
    ## iter  60 value 397.393452
    ## iter  70 value 373.782425
    ## iter  80 value 360.383027
    ## iter  90 value 349.766162
    ## iter 100 value 343.817773
    ## final  value 343.817773 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1293.990089 
    ## iter  10 value 544.375869
    ## iter  20 value 302.006823
    ## iter  30 value 204.560060
    ## iter  40 value 178.113953
    ## iter  50 value 164.735303
    ## iter  60 value 157.342352
    ## iter  70 value 154.216170
    ## iter  80 value 151.897714
    ## iter  90 value 150.933615
    ## iter 100 value 150.555059
    ## final  value 150.555059 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1263.304155 
    ## iter  10 value 962.244272
    ## iter  20 value 772.854669
    ## iter  30 value 720.282356
    ## iter  40 value 669.064953
    ## iter  50 value 634.498985
    ## iter  60 value 612.906696
    ## iter  70 value 596.204033
    ## iter  80 value 589.456302
    ## iter  90 value 583.223104
    ## iter 100 value 579.939562
    ## final  value 579.939562 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1319.189007 
    ## iter  10 value 690.050263
    ## iter  20 value 391.631582
    ## iter  30 value 290.433726
    ## iter  40 value 222.646196
    ## iter  50 value 187.815907
    ## iter  60 value 147.949652
    ## iter  70 value 118.318210
    ## iter  80 value 89.633802
    ## iter  90 value 70.374334
    ## iter 100 value 49.948684
    ## final  value 49.948684 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1392.645250 
    ## iter  10 value 403.476660
    ## iter  20 value 149.142896
    ## iter  30 value 86.322646
    ## iter  40 value 74.111396
    ## iter  50 value 58.249886
    ## iter  60 value 50.136863
    ## iter  70 value 43.379247
    ## iter  80 value 38.746475
    ## iter  90 value 38.242837
    ## iter 100 value 37.348900
    ## final  value 37.348900 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1249.630821 
    ## iter  10 value 967.077399
    ## iter  20 value 830.241373
    ## iter  30 value 770.764264
    ## iter  40 value 729.168688
    ## iter  50 value 706.770867
    ## iter  60 value 700.791228
    ## iter  70 value 698.219171
    ## iter  80 value 687.731621
    ## iter  90 value 671.231475
    ## iter 100 value 615.903279
    ## final  value 615.903279 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1254.594555 
    ## iter  10 value 776.925253
    ## iter  20 value 392.072922
    ## iter  30 value 292.267342
    ## iter  40 value 239.266406
    ## iter  50 value 204.500102
    ## iter  60 value 189.671947
    ## iter  70 value 176.838251
    ## iter  80 value 171.788514
    ## iter  90 value 168.591350
    ## iter 100 value 167.091708
    ## final  value 167.091708 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1274.101834 
    ## iter  10 value 589.526088
    ## iter  20 value 192.154052
    ## iter  30 value 84.963706
    ## iter  40 value 46.761987
    ## iter  50 value 39.817324
    ## iter  60 value 35.457041
    ## iter  70 value 32.613327
    ## iter  80 value 29.434098
    ## iter  90 value 25.905061
    ## iter 100 value 22.632656
    ## final  value 22.632656 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1280.695429 
    ## iter  10 value 1121.529973
    ## iter  20 value 963.615884
    ## iter  30 value 855.972956
    ## iter  40 value 835.321299
    ## iter  50 value 818.504920
    ## iter  60 value 800.065789
    ## iter  70 value 777.506210
    ## iter  80 value 765.602918
    ## iter  90 value 761.327442
    ## iter 100 value 759.963088
    ## final  value 759.963088 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1276.264167 
    ## iter  10 value 808.839913
    ## iter  20 value 673.874938
    ## iter  30 value 540.936044
    ## iter  40 value 456.017829
    ## iter  50 value 429.238597
    ## iter  60 value 397.411725
    ## iter  70 value 360.922849
    ## iter  80 value 341.260428
    ## iter  90 value 329.573046
    ## iter 100 value 324.800056
    ## final  value 324.800056 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1339.833499 
    ## iter  10 value 655.424246
    ## iter  20 value 349.053949
    ## iter  30 value 224.580450
    ## iter  40 value 186.653188
    ## iter  50 value 168.858631
    ## iter  60 value 160.718469
    ## iter  70 value 158.287836
    ## iter  80 value 156.343237
    ## iter  90 value 152.555172
    ## iter 100 value 149.883073
    ## final  value 149.883073 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1278.128940 
    ## iter  10 value 1046.376352
    ## iter  20 value 854.890709
    ## iter  30 value 795.058498
    ## iter  40 value 747.631769
    ## iter  50 value 712.230226
    ## iter  60 value 698.918375
    ## iter  70 value 688.831936
    ## iter  80 value 685.971277
    ## iter  90 value 685.641861
    ## iter 100 value 684.964101
    ## final  value 684.964101 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1199.852763 
    ## iter  10 value 562.548521
    ## iter  20 value 360.712283
    ## iter  30 value 249.695906
    ## iter  40 value 181.469649
    ## iter  50 value 135.587229
    ## iter  60 value 100.275908
    ## iter  70 value 83.182210
    ## iter  80 value 72.696053
    ## iter  90 value 69.611785
    ## iter 100 value 68.755512
    ## final  value 68.755512 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1246.138084 
    ## iter  10 value 424.511359
    ## iter  20 value 126.254436
    ## iter  30 value 77.426677
    ## iter  40 value 52.692254
    ## iter  50 value 34.144084
    ## iter  60 value 26.187300
    ## iter  70 value 24.580309
    ## iter  80 value 22.875143
    ## iter  90 value 19.956720
    ## iter 100 value 19.699270
    ## final  value 19.699270 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1282.822815 
    ## iter  10 value 863.780792
    ## iter  20 value 780.044247
    ## iter  30 value 748.795585
    ## iter  40 value 715.842140
    ## iter  50 value 693.753079
    ## iter  60 value 683.164625
    ## iter  70 value 681.538122
    ## iter  80 value 680.829856
    ## iter  90 value 680.381941
    ## iter 100 value 680.110090
    ## final  value 680.110090 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1287.188950 
    ## iter  10 value 681.562100
    ## iter  20 value 405.383009
    ## iter  30 value 337.283202
    ## iter  40 value 271.811144
    ## iter  50 value 187.719649
    ## iter  60 value 155.856994
    ## iter  70 value 140.296757
    ## iter  80 value 136.406152
    ## iter  90 value 130.966095
    ## iter 100 value 127.024343
    ## final  value 127.024343 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1233.794715 
    ## iter  10 value 447.702195
    ## iter  20 value 197.350475
    ## iter  30 value 142.760597
    ## iter  40 value 123.131208
    ## iter  50 value 104.662545
    ## iter  60 value 83.797535
    ## iter  70 value 55.601935
    ## iter  80 value 32.161581
    ## iter  90 value 25.950503
    ## iter 100 value 18.173383
    ## final  value 18.173383 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1308.252035 
    ## iter  10 value 913.948941
    ## iter  20 value 841.437747
    ## iter  30 value 801.524156
    ## iter  40 value 778.648023
    ## iter  50 value 771.596091
    ## iter  60 value 766.366925
    ## iter  70 value 761.784578
    ## iter  80 value 760.881295
    ## iter  90 value 760.816190
    ## iter 100 value 760.814110
    ## final  value 760.814110 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1362.148536 
    ## iter  10 value 695.858972
    ## iter  20 value 482.041768
    ## iter  30 value 412.816854
    ## iter  40 value 371.911354
    ## iter  50 value 339.450511
    ## iter  60 value 321.802492
    ## iter  70 value 315.517019
    ## iter  80 value 312.213083
    ## iter  90 value 311.370834
    ## iter 100 value 308.715767
    ## final  value 308.715767 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1306.764304 
    ## iter  10 value 482.347756
    ## iter  20 value 266.594007
    ## iter  30 value 193.493309
    ## iter  40 value 174.267824
    ## iter  50 value 162.420448
    ## iter  60 value 159.201501
    ## iter  70 value 157.041767
    ## iter  80 value 156.343479
    ## iter  90 value 155.918644
    ## iter 100 value 154.264975
    ## final  value 154.264975 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1340.160627 
    ## iter  10 value 984.868094
    ## iter  20 value 901.390135
    ## iter  30 value 847.439715
    ## iter  40 value 807.620250
    ## iter  50 value 789.250381
    ## iter  60 value 782.304583
    ## iter  70 value 773.437862
    ## iter  80 value 771.456352
    ## iter  90 value 767.895214
    ## iter 100 value 766.536394
    ## final  value 766.536394 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1263.572174 
    ## iter  10 value 551.247147
    ## iter  20 value 379.263151
    ## iter  30 value 286.407363
    ## iter  40 value 245.929826
    ## iter  50 value 207.911013
    ## iter  60 value 174.117975
    ## iter  70 value 154.102000
    ## iter  80 value 133.992362
    ## iter  90 value 119.247217
    ## iter 100 value 111.566165
    ## final  value 111.566165 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1236.199334 
    ## iter  10 value 393.351860
    ## iter  20 value 120.891624
    ## iter  30 value 71.337297
    ## iter  40 value 48.110778
    ## iter  50 value 34.989879
    ## iter  60 value 30.058356
    ## iter  70 value 28.929359
    ## iter  80 value 27.828605
    ## iter  90 value 26.557832
    ## iter 100 value 25.961222
    ## final  value 25.961222 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1204.535444 
    ## iter  10 value 915.289226
    ## iter  20 value 787.394292
    ## iter  30 value 752.355447
    ## iter  40 value 726.480941
    ## iter  50 value 690.507040
    ## iter  60 value 674.601601
    ## iter  70 value 659.437110
    ## iter  80 value 652.830682
    ## iter  90 value 649.215680
    ## iter 100 value 646.305030
    ## final  value 646.305030 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1360.623270 
    ## iter  10 value 1150.280911
    ## iter  20 value 1097.723158
    ## iter  30 value 1001.362722
    ## iter  40 value 598.510518
    ## iter  50 value 523.167207
    ## iter  60 value 447.639454
    ## iter  70 value 386.171001
    ## iter  80 value 236.896973
    ## iter  90 value 155.162301
    ## iter 100 value 116.209751
    ## final  value 116.209751 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1422.487699 
    ## iter  10 value 560.303572
    ## iter  20 value 214.223465
    ## iter  30 value 116.581801
    ## iter  40 value 77.693036
    ## iter  50 value 39.316972
    ## iter  60 value 20.955471
    ## iter  70 value 9.716903
    ## iter  80 value 7.289987
    ## iter  90 value 6.847978
    ## iter 100 value 6.785142
    ## final  value 6.785142 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1290.887318 
    ## iter  10 value 939.171008
    ## iter  20 value 857.371082
    ## iter  30 value 816.040718
    ## iter  40 value 781.849994
    ## iter  50 value 763.413413
    ## iter  60 value 752.639506
    ## iter  70 value 746.319840
    ## iter  80 value 743.761092
    ## iter  90 value 742.247466
    ## iter 100 value 742.044530
    ## final  value 742.044530 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1496.864478 
    ## iter  10 value 1122.958107
    ## iter  20 value 926.343014
    ## iter  30 value 708.985470
    ## iter  40 value 585.856635
    ## iter  50 value 524.616193
    ## iter  60 value 418.146121
    ## iter  70 value 357.594207
    ## iter  80 value 332.022628
    ## iter  90 value 318.303993
    ## iter 100 value 312.520423
    ## final  value 312.520423 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1278.707328 
    ## iter  10 value 564.699342
    ## iter  20 value 311.631013
    ## iter  30 value 231.094258
    ## iter  40 value 183.139890
    ## iter  50 value 166.446539
    ## iter  60 value 158.988163
    ## iter  70 value 155.969173
    ## iter  80 value 154.138984
    ## iter  90 value 153.335052
    ## iter 100 value 152.030485
    ## final  value 152.030485 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1281.981872 
    ## iter  10 value 919.795870
    ## iter  20 value 779.417379
    ## iter  30 value 741.021111
    ## iter  40 value 688.429412
    ## iter  50 value 649.088490
    ## iter  60 value 620.730959
    ## iter  70 value 607.093782
    ## iter  80 value 602.415614
    ## iter  90 value 593.833898
    ## iter 100 value 589.893059
    ## final  value 589.893059 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1195.215661 
    ## iter  10 value 693.248939
    ## iter  20 value 460.670462
    ## iter  30 value 392.337601
    ## iter  40 value 353.455480
    ## iter  50 value 322.455210
    ## iter  60 value 303.058231
    ## iter  70 value 294.749978
    ## iter  80 value 289.961257
    ## iter  90 value 287.447111
    ## iter 100 value 285.402836
    ## final  value 285.402836 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1312.453370 
    ## iter  10 value 374.502707
    ## iter  20 value 129.156754
    ## iter  30 value 58.157549
    ## iter  40 value 35.014232
    ## iter  50 value 27.355591
    ## iter  60 value 23.842432
    ## iter  70 value 22.493244
    ## iter  80 value 20.374266
    ## iter  90 value 13.040913
    ## iter 100 value 7.467712
    ## final  value 7.467712 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1301.439427 
    ## iter  10 value 908.205780
    ## iter  20 value 827.516439
    ## iter  30 value 774.699157
    ## iter  40 value 704.750434
    ## iter  50 value 654.971294
    ## iter  60 value 620.214960
    ## iter  70 value 594.470375
    ## iter  80 value 589.468889
    ## iter  90 value 585.452034
    ## iter 100 value 582.420794
    ## final  value 582.420794 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1275.901657 
    ## iter  10 value 719.674715
    ## iter  20 value 601.277139
    ## iter  30 value 501.704921
    ## iter  40 value 435.697059
    ## iter  50 value 397.170179
    ## iter  60 value 371.932593
    ## iter  70 value 343.898732
    ## iter  80 value 336.705877
    ## iter  90 value 335.554354
    ## iter 100 value 334.987979
    ## final  value 334.987979 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1207.255414 
    ## iter  10 value 386.701267
    ## iter  20 value 125.638357
    ## iter  30 value 63.151688
    ## iter  40 value 28.113267
    ## iter  50 value 16.961777
    ## iter  60 value 11.524576
    ## iter  70 value 9.042052
    ## iter  80 value 6.768589
    ## iter  90 value 5.490506
    ## iter 100 value 5.006701
    ## final  value 5.006701 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1237.974274 
    ## iter  10 value 1044.018373
    ## iter  20 value 906.437736
    ## iter  30 value 874.208386
    ## iter  40 value 845.089192
    ## iter  50 value 825.212941
    ## iter  60 value 812.148823
    ## iter  70 value 804.070526
    ## iter  80 value 801.110406
    ## iter  90 value 799.835655
    ## iter 100 value 799.701190
    ## final  value 799.701190 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1353.110601 
    ## iter  10 value 1005.237864
    ## iter  20 value 891.840746
    ## iter  30 value 850.460796
    ## iter  40 value 818.615914
    ## iter  50 value 768.386153
    ## iter  60 value 584.819970
    ## iter  70 value 471.740526
    ## iter  80 value 407.163200
    ## iter  90 value 377.920412
    ## iter 100 value 361.599190
    ## final  value 361.599190 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1279.183338 
    ## iter  10 value 583.795644
    ## iter  20 value 312.877187
    ## iter  30 value 224.769384
    ## iter  40 value 188.347365
    ## iter  50 value 168.752647
    ## iter  60 value 161.662082
    ## iter  70 value 158.769695
    ## iter  80 value 157.549873
    ## iter  90 value 154.347248
    ## iter 100 value 150.016514
    ## final  value 150.016514 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1317.877933 
    ## iter  10 value 1134.693011
    ## iter  20 value 974.411188
    ## iter  30 value 810.278617
    ## iter  40 value 753.919805
    ## iter  50 value 717.577878
    ## iter  60 value 682.122554
    ## iter  70 value 672.914291
    ## iter  80 value 671.384794
    ## iter  90 value 670.818427
    ## iter 100 value 670.373447
    ## final  value 670.373447 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1236.037117 
    ## iter  10 value 720.375602
    ## iter  20 value 510.321532
    ## iter  30 value 372.457959
    ## iter  40 value 322.845046
    ## iter  50 value 287.854460
    ## iter  60 value 255.159496
    ## iter  70 value 241.938871
    ## iter  80 value 231.080941
    ## iter  90 value 226.124352
    ## iter 100 value 223.657884
    ## final  value 223.657884 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1341.080920 
    ## iter  10 value 440.894834
    ## iter  20 value 115.267583
    ## iter  30 value 31.693744
    ## iter  40 value 21.614769
    ## iter  50 value 13.969612
    ## iter  60 value 11.947014
    ## iter  70 value 9.007411
    ## iter  80 value 7.193868
    ## iter  90 value 6.633915
    ## iter 100 value 5.355950
    ## final  value 5.355950 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1291.914811 
    ## iter  10 value 1060.937850
    ## iter  20 value 864.492457
    ## iter  30 value 765.800927
    ## iter  40 value 724.109655
    ## iter  50 value 702.118180
    ## iter  60 value 681.927601
    ## iter  70 value 672.178994
    ## iter  80 value 668.673886
    ## iter  90 value 665.916567
    ## iter 100 value 663.235143
    ## final  value 663.235143 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1277.655051 
    ## iter  10 value 645.170345
    ## iter  20 value 334.007449
    ## iter  30 value 226.391895
    ## iter  40 value 147.168590
    ## iter  50 value 103.862100
    ## iter  60 value 69.133578
    ## iter  70 value 49.610037
    ## iter  80 value 39.326152
    ## iter  90 value 31.745222
    ## iter 100 value 29.929849
    ## final  value 29.929849 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1401.586024 
    ## iter  10 value 1079.751438
    ## iter  20 value 796.018621
    ## iter  30 value 424.548991
    ## iter  40 value 188.271782
    ## iter  50 value 103.413519
    ## iter  60 value 61.625574
    ## iter  70 value 41.923292
    ## iter  80 value 30.303165
    ## iter  90 value 23.325645
    ## iter 100 value 18.488674
    ## final  value 18.488674 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1278.605172 
    ## iter  10 value 919.735436
    ## iter  20 value 849.797607
    ## iter  30 value 800.421583
    ## iter  40 value 775.639550
    ## iter  50 value 768.410847
    ## iter  60 value 766.075706
    ## iter  70 value 764.528433
    ## iter  80 value 763.913118
    ## iter  90 value 763.897570
    ## final  value 763.897519 
    ## converged
    ## # weights:  255
    ## initial  value 1419.277103 
    ## iter  10 value 1136.285590
    ## iter  20 value 931.781637
    ## iter  30 value 659.910019
    ## iter  40 value 458.641554
    ## iter  50 value 382.236080
    ## iter  60 value 336.935144
    ## iter  70 value 315.564687
    ## iter  80 value 303.457673
    ## iter  90 value 297.891854
    ## iter 100 value 291.416345
    ## final  value 291.416345 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1263.869312 
    ## iter  10 value 546.544020
    ## iter  20 value 321.771502
    ## iter  30 value 238.576872
    ## iter  40 value 199.428754
    ## iter  50 value 182.971856
    ## iter  60 value 166.305781
    ## iter  70 value 157.963720
    ## iter  80 value 154.679860
    ## iter  90 value 153.426325
    ## iter 100 value 151.236541
    ## final  value 151.236541 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1207.074375 
    ## iter  10 value 902.774210
    ## iter  20 value 808.127168
    ## iter  30 value 798.668864
    ## iter  40 value 776.171669
    ## iter  50 value 741.030185
    ## iter  60 value 723.087949
    ## iter  70 value 714.212002
    ## iter  80 value 713.429914
    ## iter  90 value 713.106257
    ## iter 100 value 712.457569
    ## final  value 712.457569 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1225.227389 
    ## iter  10 value 722.070219
    ## iter  20 value 439.592329
    ## iter  30 value 341.298112
    ## iter  40 value 277.626111
    ## iter  50 value 250.703258
    ## iter  60 value 217.828302
    ## iter  70 value 207.539155
    ## iter  80 value 201.567075
    ## iter  90 value 199.578773
    ## iter 100 value 196.813757
    ## final  value 196.813757 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1363.282810 
    ## iter  10 value 463.100537
    ## iter  20 value 172.963130
    ## iter  30 value 68.028672
    ## iter  40 value 40.538281
    ## iter  50 value 20.086216
    ## iter  60 value 10.876882
    ## iter  70 value 10.182227
    ## iter  80 value 9.400269
    ## iter  90 value 8.058333
    ## iter 100 value 7.626327
    ## final  value 7.626327 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1222.895530 
    ## iter  10 value 979.247211
    ## iter  20 value 820.433244
    ## iter  30 value 786.746126
    ## iter  40 value 756.712121
    ## iter  50 value 728.866298
    ## iter  60 value 716.433867
    ## iter  70 value 712.526970
    ## iter  80 value 710.878606
    ## iter  90 value 710.454887
    ## iter 100 value 710.333282
    ## final  value 710.333282 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1273.411190 
    ## iter  10 value 961.669226
    ## iter  20 value 649.820742
    ## iter  30 value 545.368270
    ## iter  40 value 406.956200
    ## iter  50 value 334.120368
    ## iter  60 value 277.412542
    ## iter  70 value 247.509130
    ## iter  80 value 231.231811
    ## iter  90 value 215.994437
    ## iter 100 value 207.910951
    ## final  value 207.910951 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1291.320102 
    ## iter  10 value 511.273005
    ## iter  20 value 126.732132
    ## iter  30 value 56.573258
    ## iter  40 value 37.332095
    ## iter  50 value 15.018491
    ## iter  60 value 8.448734
    ## iter  70 value 6.227122
    ## iter  80 value 5.779728
    ## iter  90 value 5.501865
    ## iter 100 value 3.887766
    ## final  value 3.887766 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1187.767048 
    ## iter  10 value 916.706800
    ## iter  20 value 848.953772
    ## iter  30 value 815.978104
    ## iter  40 value 802.008006
    ## iter  50 value 795.144527
    ## iter  60 value 783.864333
    ## iter  70 value 780.034587
    ## iter  80 value 779.365716
    ## iter  90 value 777.911365
    ## iter 100 value 777.865750
    ## final  value 777.865750 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1305.661150 
    ## iter  10 value 715.047478
    ## iter  20 value 589.558823
    ## iter  30 value 540.523306
    ## iter  40 value 486.676224
    ## iter  50 value 421.129695
    ## iter  60 value 373.772944
    ## iter  70 value 351.799175
    ## iter  80 value 326.458455
    ## iter  90 value 317.266078
    ## iter 100 value 314.993710
    ## final  value 314.993710 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1291.704019 
    ## iter  10 value 482.258606
    ## iter  20 value 249.338480
    ## iter  30 value 189.757851
    ## iter  40 value 165.563952
    ## iter  50 value 157.711442
    ## iter  60 value 152.955017
    ## iter  70 value 150.184397
    ## iter  80 value 149.726217
    ## iter  90 value 149.585250
    ## iter 100 value 149.235637
    ## final  value 149.235637 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1264.066890 
    ## iter  10 value 898.419692
    ## iter  20 value 820.061664
    ## iter  30 value 779.155159
    ## iter  40 value 740.336220
    ## iter  50 value 705.508613
    ## iter  60 value 694.671876
    ## iter  70 value 691.486475
    ## iter  80 value 688.450643
    ## iter  90 value 684.416170
    ## iter 100 value 677.550647
    ## final  value 677.550647 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1308.394747 
    ## iter  10 value 741.634487
    ## iter  20 value 596.260518
    ## iter  30 value 471.072553
    ## iter  40 value 405.364272
    ## iter  50 value 353.715360
    ## iter  60 value 317.308729
    ## iter  70 value 292.872933
    ## iter  80 value 221.767274
    ## iter  90 value 193.833497
    ## iter 100 value 177.310473
    ## final  value 177.310473 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1408.209140 
    ## iter  10 value 429.636505
    ## iter  20 value 97.614646
    ## iter  30 value 38.814592
    ## iter  40 value 18.687467
    ## iter  50 value 12.905458
    ## iter  60 value 7.470975
    ## iter  70 value 6.095912
    ## iter  80 value 2.735689
    ## iter  90 value 1.963940
    ## iter 100 value 1.570482
    ## final  value 1.570482 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1260.703123 
    ## iter  10 value 992.398170
    ## iter  20 value 900.944420
    ## iter  30 value 873.188145
    ## iter  40 value 860.502695
    ## iter  50 value 853.777129
    ## iter  60 value 850.817226
    ## iter  70 value 850.139323
    ## iter  80 value 849.992573
    ## iter  90 value 849.955280
    ## iter 100 value 849.861766
    ## final  value 849.861766 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1253.821602 
    ## iter  10 value 701.383358
    ## iter  20 value 395.709858
    ## iter  30 value 267.874119
    ## iter  40 value 185.627544
    ## iter  50 value 142.716645
    ## iter  60 value 123.768127
    ## iter  70 value 108.417622
    ## iter  80 value 96.298593
    ## iter  90 value 87.163852
    ## iter 100 value 84.140374
    ## final  value 84.140374 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1360.670303 
    ## iter  10 value 480.516656
    ## iter  20 value 187.067891
    ## iter  30 value 100.204158
    ## iter  40 value 69.921546
    ## iter  50 value 37.132131
    ## iter  60 value 23.436521
    ## iter  70 value 15.058032
    ## iter  80 value 7.166230
    ## iter  90 value 5.263761
    ## iter 100 value 5.157342
    ## final  value 5.157342 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1305.611600 
    ## iter  10 value 957.704805
    ## iter  20 value 873.491043
    ## iter  30 value 834.381145
    ## iter  40 value 814.028853
    ## iter  50 value 793.700118
    ## iter  60 value 773.512747
    ## iter  70 value 767.094226
    ## iter  80 value 760.979131
    ## iter  90 value 758.543762
    ## iter 100 value 757.821681
    ## final  value 757.821681 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1170.309607 
    ## iter  10 value 646.172242
    ## iter  20 value 519.904555
    ## iter  30 value 432.055663
    ## iter  40 value 375.011201
    ## iter  50 value 351.447880
    ## iter  60 value 338.088113
    ## iter  70 value 329.377387
    ## iter  80 value 324.463144
    ## iter  90 value 320.127460
    ## iter 100 value 315.413676
    ## final  value 315.413676 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1330.092419 
    ## iter  10 value 550.232475
    ## iter  20 value 303.948015
    ## iter  30 value 196.852147
    ## iter  40 value 172.322382
    ## iter  50 value 162.250931
    ## iter  60 value 158.825509
    ## iter  70 value 157.149121
    ## iter  80 value 153.405916
    ## iter  90 value 148.060252
    ## iter 100 value 144.411858
    ## final  value 144.411858 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1284.530138 
    ## iter  10 value 1114.494528
    ## iter  20 value 1070.595734
    ## iter  30 value 935.119275
    ## iter  40 value 887.208774
    ## iter  50 value 868.092567
    ## iter  60 value 854.331584
    ## iter  70 value 851.689749
    ## iter  80 value 851.212535
    ## iter  90 value 850.760454
    ## iter 100 value 849.497254
    ## final  value 849.497254 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1337.184735 
    ## iter  10 value 1027.896498
    ## iter  20 value 847.935687
    ## iter  30 value 684.445245
    ## iter  40 value 562.752881
    ## iter  50 value 540.014958
    ## iter  60 value 509.146974
    ## iter  70 value 492.732324
    ## iter  80 value 484.292945
    ## iter  90 value 474.694505
    ## iter 100 value 471.553494
    ## final  value 471.553494 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1451.710187 
    ## iter  10 value 445.924102
    ## iter  20 value 88.275307
    ## iter  30 value 19.843520
    ## iter  40 value 5.525976
    ## iter  50 value 1.701608
    ## iter  60 value 1.326866
    ## iter  70 value 1.273884
    ## iter  80 value 1.212685
    ## iter  90 value 1.135593
    ## iter 100 value 1.037143
    ## final  value 1.037143 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1205.150645 
    ## iter  10 value 890.936237
    ## iter  20 value 825.434321
    ## iter  30 value 804.629226
    ## iter  40 value 781.614996
    ## iter  50 value 773.918804
    ## iter  60 value 771.638232
    ## iter  70 value 771.426747
    ## iter  80 value 771.324825
    ## iter  90 value 771.308589
    ## iter 100 value 771.302622
    ## final  value 771.302622 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1251.640283 
    ## iter  10 value 707.188953
    ## iter  20 value 358.875062
    ## iter  30 value 297.859650
    ## iter  40 value 244.320956
    ## iter  50 value 201.652607
    ## iter  60 value 179.714577
    ## iter  70 value 168.807338
    ## iter  80 value 165.379827
    ## iter  90 value 163.984364
    ## iter 100 value 163.242852
    ## final  value 163.242852 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1252.193649 
    ## iter  10 value 419.836833
    ## iter  20 value 128.958469
    ## iter  30 value 52.990258
    ## iter  40 value 26.840298
    ## iter  50 value 7.386080
    ## iter  60 value 1.360239
    ## iter  70 value 0.094622
    ## iter  80 value 0.010893
    ## iter  90 value 0.004345
    ## iter 100 value 0.001621
    ## final  value 0.001621 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1219.438686 
    ## iter  10 value 961.830973
    ## iter  20 value 865.111063
    ## iter  30 value 840.795043
    ## iter  40 value 822.675212
    ## iter  50 value 800.922206
    ## iter  60 value 790.797758
    ## iter  70 value 787.753595
    ## iter  80 value 786.695645
    ## iter  90 value 786.658800
    ## iter 100 value 786.627313
    ## final  value 786.627313 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1337.271942 
    ## iter  10 value 806.129339
    ## iter  20 value 661.911983
    ## iter  30 value 533.246839
    ## iter  40 value 448.297590
    ## iter  50 value 403.822756
    ## iter  60 value 373.687410
    ## iter  70 value 347.969449
    ## iter  80 value 338.944087
    ## iter  90 value 332.100667
    ## iter 100 value 327.553314
    ## final  value 327.553314 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1271.118182 
    ## iter  10 value 578.461143
    ## iter  20 value 310.956565
    ## iter  30 value 197.833982
    ## iter  40 value 177.513302
    ## iter  50 value 168.408509
    ## iter  60 value 164.104096
    ## iter  70 value 163.177336
    ## iter  80 value 163.079797
    ## iter  90 value 163.035254
    ## iter 100 value 163.020072
    ## final  value 163.020072 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1288.488499 
    ## iter  10 value 1014.540763
    ## iter  20 value 871.360952
    ## iter  30 value 798.500409
    ## iter  40 value 759.528950
    ## iter  50 value 735.809833
    ## iter  60 value 724.679210
    ## iter  70 value 721.167257
    ## iter  80 value 718.485427
    ## iter  90 value 717.069616
    ## iter 100 value 714.416941
    ## final  value 714.416941 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1384.927700 
    ## iter  10 value 1065.650202
    ## iter  20 value 764.173769
    ## iter  30 value 484.766242
    ## iter  40 value 329.269453
    ## iter  50 value 255.515491
    ## iter  60 value 218.970455
    ## iter  70 value 140.273834
    ## iter  80 value 96.411954
    ## iter  90 value 64.092550
    ## iter 100 value 56.215762
    ## final  value 56.215762 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1280.919012 
    ## iter  10 value 348.743805
    ## iter  20 value 134.856484
    ## iter  30 value 80.584517
    ## iter  40 value 49.348633
    ## iter  50 value 34.060703
    ## iter  60 value 20.046629
    ## iter  70 value 17.314883
    ## iter  80 value 12.677090
    ## iter  90 value 8.487867
    ## iter 100 value 7.425707
    ## final  value 7.425707 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1233.388009 
    ## iter  10 value 922.068809
    ## iter  20 value 817.923080
    ## iter  30 value 763.550021
    ## iter  40 value 710.608509
    ## iter  50 value 678.962638
    ## iter  60 value 670.692932
    ## iter  70 value 668.695506
    ## iter  80 value 666.710068
    ## iter  90 value 662.850469
    ## iter 100 value 659.469924
    ## final  value 659.469924 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1201.606675 
    ## iter  10 value 583.378952
    ## iter  20 value 397.861900
    ## iter  30 value 330.249208
    ## iter  40 value 281.363539
    ## iter  50 value 248.230613
    ## iter  60 value 231.745359
    ## iter  70 value 216.305270
    ## iter  80 value 211.608803
    ## iter  90 value 204.560621
    ## iter 100 value 199.389466
    ## final  value 199.389466 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1254.815881 
    ## iter  10 value 445.008433
    ## iter  20 value 126.365617
    ## iter  30 value 45.490210
    ## iter  40 value 27.152078
    ## iter  50 value 24.544668
    ## iter  60 value 20.652526
    ## iter  70 value 19.149252
    ## iter  80 value 17.529084
    ## iter  90 value 15.098686
    ## iter 100 value 13.806935
    ## final  value 13.806935 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1224.074483 
    ## iter  10 value 1102.245737
    ## iter  20 value 967.962302
    ## iter  30 value 856.507004
    ## iter  40 value 827.611713
    ## iter  50 value 815.304721
    ## iter  60 value 800.882120
    ## iter  70 value 768.604783
    ## iter  80 value 761.455499
    ## iter  90 value 749.876466
    ## iter 100 value 746.483518
    ## final  value 746.483518 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1405.478191 
    ## iter  10 value 905.720462
    ## iter  20 value 791.653189
    ## iter  30 value 649.643618
    ## iter  40 value 586.415176
    ## iter  50 value 544.718272
    ## iter  60 value 517.789118
    ## iter  70 value 471.212957
    ## iter  80 value 408.382091
    ## iter  90 value 368.055981
    ## iter 100 value 342.625468
    ## final  value 342.625468 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1293.388401 
    ## iter  10 value 503.323181
    ## iter  20 value 298.122042
    ## iter  30 value 227.901282
    ## iter  40 value 181.490599
    ## iter  50 value 172.780819
    ## iter  60 value 164.974371
    ## iter  70 value 162.763945
    ## iter  80 value 160.161919
    ## iter  90 value 158.529920
    ## iter 100 value 157.917468
    ## final  value 157.917468 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1306.295932 
    ## iter  10 value 974.937518
    ## iter  20 value 866.531685
    ## iter  30 value 768.763728
    ## iter  40 value 706.394478
    ## iter  50 value 661.139465
    ## iter  60 value 633.996902
    ## iter  70 value 613.838351
    ## iter  80 value 603.975675
    ## iter  90 value 596.763759
    ## iter 100 value 592.614952
    ## final  value 592.614952 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1134.946733 
    ## iter  10 value 551.973623
    ## iter  20 value 376.689933
    ## iter  30 value 310.765940
    ## iter  40 value 275.097559
    ## iter  50 value 255.305833
    ## iter  60 value 245.454531
    ## iter  70 value 239.876164
    ## iter  80 value 234.370919
    ## iter  90 value 231.661521
    ## iter 100 value 229.956122
    ## final  value 229.956122 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1210.148727 
    ## iter  10 value 289.766358
    ## iter  20 value 102.707439
    ## iter  30 value 26.354709
    ## iter  40 value 7.406693
    ## iter  50 value 4.395481
    ## iter  60 value 4.153944
    ## iter  70 value 3.956705
    ## iter  80 value 3.778130
    ## iter  90 value 3.646691
    ## iter 100 value 3.561775
    ## final  value 3.561775 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1301.584246 
    ## iter  10 value 994.774049
    ## iter  20 value 884.105773
    ## iter  30 value 833.899754
    ## iter  40 value 808.649434
    ## iter  50 value 774.929140
    ## iter  60 value 753.472284
    ## iter  70 value 736.173289
    ## iter  80 value 728.050716
    ## iter  90 value 725.216428
    ## iter 100 value 721.550971
    ## final  value 721.550971 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1396.338755 
    ## iter  10 value 869.728513
    ## iter  20 value 656.744471
    ## iter  30 value 474.705672
    ## iter  40 value 407.516686
    ## iter  50 value 364.263829
    ## iter  60 value 307.687228
    ## iter  70 value 241.206248
    ## iter  80 value 208.452147
    ## iter  90 value 197.389669
    ## iter 100 value 187.830678
    ## final  value 187.830678 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1357.983426 
    ## iter  10 value 509.814510
    ## iter  20 value 248.765687
    ## iter  30 value 117.075393
    ## iter  40 value 54.263251
    ## iter  50 value 32.767459
    ## iter  60 value 26.613579
    ## iter  70 value 25.629279
    ## iter  80 value 25.560124
    ## iter  90 value 25.550595
    ## iter 100 value 25.548606
    ## final  value 25.548606 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1295.045027 
    ## iter  10 value 846.718429
    ## iter  20 value 818.219592
    ## iter  30 value 805.637464
    ## iter  40 value 787.455363
    ## iter  50 value 767.524411
    ## iter  60 value 757.779191
    ## iter  70 value 755.750796
    ## iter  80 value 754.988298
    ## iter  90 value 754.911772
    ## iter 100 value 754.907652
    ## final  value 754.907652 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1212.758624 
    ## iter  10 value 851.515057
    ## iter  20 value 547.015620
    ## iter  30 value 447.809580
    ## iter  40 value 365.163037
    ## iter  50 value 339.106873
    ## iter  60 value 326.485584
    ## iter  70 value 320.886032
    ## iter  80 value 317.034596
    ## iter  90 value 314.742486
    ## iter 100 value 312.795286
    ## final  value 312.795286 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1312.902529 
    ## iter  10 value 635.840745
    ## iter  20 value 370.901013
    ## iter  30 value 243.381197
    ## iter  40 value 195.203450
    ## iter  50 value 179.847563
    ## iter  60 value 170.307469
    ## iter  70 value 163.319856
    ## iter  80 value 161.406635
    ## iter  90 value 160.367703
    ## iter 100 value 159.323424
    ## final  value 159.323424 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1259.355235 
    ## iter  10 value 1079.248329
    ## iter  20 value 909.064878
    ## iter  30 value 791.413485
    ## iter  40 value 761.705025
    ## iter  50 value 747.137374
    ## iter  60 value 743.423072
    ## iter  70 value 740.443355
    ## iter  80 value 735.316746
    ## iter  90 value 732.202409
    ## iter 100 value 713.006843
    ## final  value 713.006843 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1200.466408 
    ## iter  10 value 623.990379
    ## iter  20 value 360.915066
    ## iter  30 value 284.656219
    ## iter  40 value 228.472980
    ## iter  50 value 196.106423
    ## iter  60 value 175.336782
    ## iter  70 value 158.330648
    ## iter  80 value 152.117601
    ## iter  90 value 145.479866
    ## iter 100 value 142.195658
    ## final  value 142.195658 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1350.470175 
    ## iter  10 value 509.739929
    ## iter  20 value 165.537500
    ## iter  30 value 73.648726
    ## iter  40 value 26.471876
    ## iter  50 value 9.772193
    ## iter  60 value 6.423580
    ## iter  70 value 6.072453
    ## iter  80 value 5.785270
    ## iter  90 value 5.578319
    ## iter 100 value 5.365215
    ## final  value 5.365215 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1356.166396 
    ## iter  10 value 1097.380473
    ## iter  20 value 1093.271855
    ## iter  30 value 1093.256115
    ## final  value 1093.256041 
    ## converged
    ## # weights:  255
    ## initial  value 1313.742801 
    ## iter  10 value 768.604621
    ## iter  20 value 443.585921
    ## iter  30 value 359.633423
    ## iter  40 value 314.299842
    ## iter  50 value 280.718220
    ## iter  60 value 247.103259
    ## iter  70 value 222.422692
    ## iter  80 value 211.219936
    ## iter  90 value 204.197469
    ## iter 100 value 194.727075
    ## final  value 194.727075 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1136.046763 
    ## iter  10 value 509.815835
    ## iter  20 value 170.539117
    ## iter  30 value 102.956909
    ## iter  40 value 60.847063
    ## iter  50 value 28.147284
    ## iter  60 value 13.728595
    ## iter  70 value 10.844296
    ## iter  80 value 9.242556
    ## iter  90 value 8.023460
    ## iter 100 value 5.807356
    ## final  value 5.807356 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1226.607760 
    ## iter  10 value 1007.034467
    ## iter  20 value 889.434078
    ## iter  30 value 834.859874
    ## iter  40 value 809.314599
    ## iter  50 value 789.913429
    ## iter  60 value 783.059384
    ## iter  70 value 779.688777
    ## iter  80 value 764.397615
    ## iter  90 value 753.768273
    ## iter 100 value 750.108933
    ## final  value 750.108933 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1275.397094 
    ## iter  10 value 723.620316
    ## iter  20 value 495.842218
    ## iter  30 value 396.564998
    ## iter  40 value 347.345805
    ## iter  50 value 325.718337
    ## iter  60 value 308.837499
    ## iter  70 value 304.357275
    ## iter  80 value 296.242600
    ## iter  90 value 292.546428
    ## iter 100 value 291.294277
    ## final  value 291.294277 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1352.991130 
    ## iter  10 value 583.749419
    ## iter  20 value 341.278775
    ## iter  30 value 193.656901
    ## iter  40 value 162.533891
    ## iter  50 value 155.811248
    ## iter  60 value 153.232414
    ## iter  70 value 151.509816
    ## iter  80 value 151.024871
    ## iter  90 value 150.698884
    ## iter 100 value 150.316217
    ## final  value 150.316217 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1196.452302 
    ## iter  10 value 873.075081
    ## iter  20 value 762.194758
    ## iter  30 value 709.778743
    ## iter  40 value 663.760336
    ## iter  50 value 631.992839
    ## iter  60 value 618.080520
    ## iter  70 value 614.338998
    ## iter  80 value 611.810947
    ## iter  90 value 611.298481
    ## iter 100 value 611.218674
    ## final  value 611.218674 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1255.462453 
    ## iter  10 value 663.134706
    ## iter  20 value 412.250517
    ## iter  30 value 357.503171
    ## iter  40 value 312.832747
    ## iter  50 value 277.520233
    ## iter  60 value 262.289448
    ## iter  70 value 255.632005
    ## iter  80 value 252.241965
    ## iter  90 value 248.731261
    ## iter 100 value 243.006653
    ## final  value 243.006653 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1258.940963 
    ## iter  10 value 347.101943
    ## iter  20 value 121.580515
    ## iter  30 value 74.874636
    ## iter  40 value 61.121326
    ## iter  50 value 45.825439
    ## iter  60 value 25.218950
    ## iter  70 value 6.289314
    ## iter  80 value 4.593265
    ## iter  90 value 4.235543
    ## iter 100 value 4.044314
    ## final  value 4.044314 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1348.208876 
    ## iter  10 value 1011.037498
    ## iter  20 value 826.527599
    ## iter  30 value 762.063521
    ## iter  40 value 703.960382
    ## iter  50 value 679.313560
    ## iter  60 value 654.880776
    ## iter  70 value 644.130189
    ## iter  80 value 638.337469
    ## iter  90 value 633.955449
    ## iter 100 value 632.778991
    ## final  value 632.778991 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1292.495293 
    ## iter  10 value 693.926261
    ## iter  20 value 391.426411
    ## iter  30 value 314.079131
    ## iter  40 value 291.959723
    ## iter  50 value 277.588103
    ## iter  60 value 264.353348
    ## iter  70 value 253.211964
    ## iter  80 value 250.355415
    ## iter  90 value 249.847425
    ## iter 100 value 248.964204
    ## final  value 248.964204 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1292.191371 
    ## iter  10 value 393.262173
    ## iter  20 value 104.894277
    ## iter  30 value 43.208471
    ## iter  40 value 12.484580
    ## iter  50 value 1.778963
    ## iter  60 value 0.103890
    ## iter  70 value 0.018328
    ## iter  80 value 0.002272
    ## iter  90 value 0.000527
    ## iter 100 value 0.000316
    ## final  value 0.000316 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1275.650393 
    ## iter  10 value 942.884802
    ## iter  20 value 848.506699
    ## iter  30 value 798.998732
    ## iter  40 value 765.943748
    ## iter  50 value 754.937554
    ## iter  60 value 752.839234
    ## iter  70 value 752.446440
    ## iter  80 value 752.354109
    ## iter  90 value 752.352303
    ## iter  90 value 752.352299
    ## iter  90 value 752.352299
    ## final  value 752.352299 
    ## converged
    ## # weights:  255
    ## initial  value 1267.474031 
    ## iter  10 value 648.906562
    ## iter  20 value 531.089909
    ## iter  30 value 474.473379
    ## iter  40 value 400.628223
    ## iter  50 value 349.798401
    ## iter  60 value 332.781364
    ## iter  70 value 318.575931
    ## iter  80 value 310.488268
    ## iter  90 value 308.030429
    ## iter 100 value 307.260222
    ## final  value 307.260222 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1253.408774 
    ## iter  10 value 592.699420
    ## iter  20 value 301.871803
    ## iter  30 value 217.507710
    ## iter  40 value 184.902576
    ## iter  50 value 169.884010
    ## iter  60 value 155.609474
    ## iter  70 value 150.928973
    ## iter  80 value 144.983756
    ## iter  90 value 142.524890
    ## iter 100 value 140.181384
    ## final  value 140.181384 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1210.052319 
    ## iter  10 value 841.492095
    ## iter  20 value 775.179722
    ## iter  30 value 727.021276
    ## iter  40 value 667.029850
    ## iter  50 value 626.524247
    ## iter  60 value 598.593344
    ## iter  70 value 577.697578
    ## iter  80 value 558.535646
    ## iter  90 value 554.228404
    ## iter 100 value 551.257936
    ## final  value 551.257936 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1265.515862 
    ## iter  10 value 649.134048
    ## iter  20 value 410.709957
    ## iter  30 value 307.651469
    ## iter  40 value 230.307943
    ## iter  50 value 191.014702
    ## iter  60 value 165.002252
    ## iter  70 value 148.053168
    ## iter  80 value 139.409496
    ## iter  90 value 128.444723
    ## iter 100 value 124.440015
    ## final  value 124.440015 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1373.983344 
    ## iter  10 value 549.147600
    ## iter  20 value 237.728131
    ## iter  30 value 134.118191
    ## iter  40 value 109.911537
    ## iter  50 value 92.205071
    ## iter  60 value 77.021880
    ## iter  70 value 63.932201
    ## iter  80 value 54.575378
    ## iter  90 value 36.278774
    ## iter 100 value 33.450298
    ## final  value 33.450298 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1240.050995 
    ## iter  10 value 1010.056568
    ## iter  20 value 849.056447
    ## iter  30 value 810.781516
    ## iter  40 value 749.230584
    ## iter  50 value 720.704897
    ## iter  60 value 704.090516
    ## iter  70 value 700.150078
    ## iter  80 value 693.987908
    ## iter  90 value 687.825030
    ## iter 100 value 686.817519
    ## final  value 686.817519 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1240.341527 
    ## iter  10 value 680.849540
    ## iter  20 value 468.758176
    ## iter  30 value 417.659523
    ## iter  40 value 370.216020
    ## iter  50 value 327.573860
    ## iter  60 value 261.033898
    ## iter  70 value 242.085183
    ## iter  80 value 231.069766
    ## iter  90 value 218.647388
    ## iter 100 value 212.320520
    ## final  value 212.320520 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1272.157659 
    ## iter  10 value 530.639156
    ## iter  20 value 212.629709
    ## iter  30 value 104.611761
    ## iter  40 value 56.779976
    ## iter  50 value 34.194968
    ## iter  60 value 25.438578
    ## iter  70 value 18.586963
    ## iter  80 value 16.551493
    ## iter  90 value 16.230131
    ## iter 100 value 12.566613
    ## final  value 12.566613 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1286.222761 
    ## iter  10 value 947.564331
    ## iter  20 value 867.711716
    ## iter  30 value 835.961409
    ## iter  40 value 800.548805
    ## iter  50 value 791.828064
    ## iter  60 value 787.942646
    ## iter  70 value 781.882402
    ## iter  80 value 769.858216
    ## iter  90 value 765.120452
    ## iter 100 value 763.747193
    ## final  value 763.747193 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1331.541238 
    ## iter  10 value 938.770093
    ## iter  20 value 653.164087
    ## iter  30 value 496.904167
    ## iter  40 value 443.937995
    ## iter  50 value 398.580735
    ## iter  60 value 373.514274
    ## iter  70 value 353.373649
    ## iter  80 value 335.698096
    ## iter  90 value 322.005497
    ## iter 100 value 314.868184
    ## final  value 314.868184 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1313.347947 
    ## iter  10 value 617.756674
    ## iter  20 value 366.871670
    ## iter  30 value 241.427284
    ## iter  40 value 180.426277
    ## iter  50 value 161.397226
    ## iter  60 value 154.917465
    ## iter  70 value 152.240609
    ## iter  80 value 150.959038
    ## iter  90 value 148.841639
    ## iter 100 value 147.845576
    ## final  value 147.845576 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1217.479693 
    ## iter  10 value 844.890687
    ## iter  20 value 744.404475
    ## iter  30 value 688.471834
    ## iter  40 value 644.576009
    ## iter  50 value 613.562178
    ## iter  60 value 602.031638
    ## iter  70 value 598.000779
    ## iter  80 value 593.235065
    ## iter  90 value 592.419279
    ## iter 100 value 590.125722
    ## final  value 590.125722 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1396.880124 
    ## iter  10 value 1128.974836
    ## iter  20 value 860.769681
    ## iter  30 value 611.508762
    ## iter  40 value 442.665237
    ## iter  50 value 391.278801
    ## iter  60 value 367.883083
    ## iter  70 value 335.079344
    ## iter  80 value 297.654835
    ## iter  90 value 251.656595
    ## iter 100 value 224.168105
    ## final  value 224.168105 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1288.503742 
    ## iter  10 value 506.845073
    ## iter  20 value 214.950286
    ## iter  30 value 148.158456
    ## iter  40 value 55.600964
    ## iter  50 value 28.443877
    ## iter  60 value 17.163573
    ## iter  70 value 11.925469
    ## iter  80 value 10.400305
    ## iter  90 value 9.476589
    ## iter 100 value 8.874271
    ## final  value 8.874271 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1219.748391 
    ## iter  10 value 934.227978
    ## iter  20 value 823.758697
    ## iter  30 value 809.022009
    ## iter  40 value 797.990351
    ## iter  50 value 797.787804
    ## final  value 797.778290 
    ## converged
    ## # weights:  255
    ## initial  value 1273.108647 
    ## iter  10 value 784.726619
    ## iter  20 value 536.666381
    ## iter  30 value 451.464798
    ## iter  40 value 356.671294
    ## iter  50 value 322.006412
    ## iter  60 value 297.657305
    ## iter  70 value 255.790777
    ## iter  80 value 247.212505
    ## iter  90 value 238.518203
    ## iter 100 value 230.797922
    ## final  value 230.797922 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1332.904956 
    ## iter  10 value 476.341525
    ## iter  20 value 90.855461
    ## iter  30 value 29.351283
    ## iter  40 value 11.909796
    ## iter  50 value 5.712133
    ## iter  60 value 2.874603
    ## iter  70 value 1.956090
    ## iter  80 value 1.924794
    ## iter  90 value 1.914873
    ## iter 100 value 1.912086
    ## final  value 1.912086 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1254.386261 
    ## iter  10 value 924.340610
    ## iter  20 value 841.566588
    ## iter  30 value 808.975394
    ## iter  40 value 792.558222
    ## iter  50 value 786.135225
    ## iter  60 value 782.258765
    ## iter  70 value 777.421284
    ## iter  80 value 775.709136
    ## iter  90 value 775.419706
    ## iter 100 value 774.287406
    ## final  value 774.287406 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1318.785449 
    ## iter  10 value 780.461112
    ## iter  20 value 517.645303
    ## iter  30 value 420.363176
    ## iter  40 value 368.808468
    ## iter  50 value 343.087488
    ## iter  60 value 329.867339
    ## iter  70 value 322.376199
    ## iter  80 value 318.245002
    ## iter  90 value 315.113858
    ## iter 100 value 313.925889
    ## final  value 313.925889 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1231.501071 
    ## iter  10 value 518.588212
    ## iter  20 value 263.230522
    ## iter  30 value 196.205820
    ## iter  40 value 170.205433
    ## iter  50 value 162.417238
    ## iter  60 value 156.804940
    ## iter  70 value 152.835709
    ## iter  80 value 152.105935
    ## iter  90 value 151.776312
    ## iter 100 value 151.357255
    ## final  value 151.357255 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1319.344051 
    ## iter  10 value 965.153409
    ## iter  20 value 820.594215
    ## iter  30 value 790.310542
    ## iter  40 value 774.415821
    ## iter  50 value 767.555550
    ## iter  60 value 755.777796
    ## iter  70 value 745.437291
    ## iter  80 value 739.393985
    ## iter  90 value 732.106359
    ## iter 100 value 712.671699
    ## final  value 712.671699 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1348.463121 
    ## iter  10 value 792.426592
    ## iter  20 value 403.950776
    ## iter  30 value 314.316410
    ## iter  40 value 278.908086
    ## iter  50 value 238.270500
    ## iter  60 value 218.104738
    ## iter  70 value 204.898412
    ## iter  80 value 198.489926
    ## iter  90 value 196.135252
    ## iter 100 value 195.050900
    ## final  value 195.050900 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1309.225501 
    ## iter  10 value 411.917397
    ## iter  20 value 109.880600
    ## iter  30 value 46.199914
    ## iter  40 value 30.122705
    ## iter  50 value 12.067615
    ## iter  60 value 6.008032
    ## iter  70 value 4.621435
    ## iter  80 value 4.314149
    ## iter  90 value 4.078529
    ## iter 100 value 3.946000
    ## final  value 3.946000 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1303.483914 
    ## iter  10 value 966.571177
    ## iter  20 value 810.573888
    ## iter  30 value 787.994510
    ## iter  40 value 769.851705
    ## iter  50 value 747.322967
    ## iter  60 value 738.902433
    ## iter  70 value 737.489425
    ## iter  80 value 736.629906
    ## iter  90 value 735.759582
    ## iter 100 value 735.272902
    ## final  value 735.272902 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1253.768944 
    ## iter  10 value 665.155956
    ## iter  20 value 437.232384
    ## iter  30 value 337.521221
    ## iter  40 value 273.564721
    ## iter  50 value 236.848703
    ## iter  60 value 204.645418
    ## iter  70 value 176.690124
    ## iter  80 value 158.566720
    ## iter  90 value 142.719514
    ## iter 100 value 135.956176
    ## final  value 135.956176 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1198.711936 
    ## iter  10 value 449.104383
    ## iter  20 value 185.980547
    ## iter  30 value 129.349556
    ## iter  40 value 91.613161
    ## iter  50 value 60.014620
    ## iter  60 value 53.432339
    ## iter  70 value 48.756824
    ## iter  80 value 41.916132
    ## iter  90 value 35.520070
    ## iter 100 value 29.475885
    ## final  value 29.475885 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1248.441482 
    ## iter  10 value 877.130725
    ## iter  20 value 819.386064
    ## iter  30 value 791.859140
    ## iter  40 value 778.371821
    ## iter  50 value 766.813778
    ## iter  60 value 762.592984
    ## iter  70 value 759.393710
    ## iter  80 value 758.904015
    ## iter  90 value 758.696729
    ## final  value 758.696236 
    ## converged
    ## # weights:  255
    ## initial  value 1246.951442 
    ## iter  10 value 616.140989
    ## iter  20 value 473.429670
    ## iter  30 value 394.898069
    ## iter  40 value 350.218510
    ## iter  50 value 335.701274
    ## iter  60 value 329.271408
    ## iter  70 value 319.555630
    ## iter  80 value 315.494480
    ## iter  90 value 314.154158
    ## iter 100 value 313.376847
    ## final  value 313.376847 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1298.326552 
    ## iter  10 value 495.551427
    ## iter  20 value 297.288457
    ## iter  30 value 236.081452
    ## iter  40 value 187.619659
    ## iter  50 value 174.347045
    ## iter  60 value 169.255569
    ## iter  70 value 166.596663
    ## iter  80 value 165.273898
    ## iter  90 value 164.661545
    ## iter 100 value 164.179025
    ## final  value 164.179025 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1206.759012 
    ## iter  10 value 899.927805
    ## iter  20 value 800.535575
    ## iter  30 value 777.556357
    ## iter  40 value 740.829264
    ## iter  50 value 715.669095
    ## iter  60 value 701.183209
    ## iter  70 value 696.850933
    ## iter  80 value 694.797003
    ## iter  90 value 690.175053
    ## iter 100 value 666.392841
    ## final  value 666.392841 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1275.641396 
    ## iter  10 value 566.568369
    ## iter  20 value 343.304786
    ## iter  30 value 269.953627
    ## iter  40 value 216.957274
    ## iter  50 value 174.906713
    ## iter  60 value 149.985047
    ## iter  70 value 126.382745
    ## iter  80 value 113.713533
    ## iter  90 value 109.877283
    ## iter 100 value 106.032153
    ## final  value 106.032153 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1279.759633 
    ## iter  10 value 506.218104
    ## iter  20 value 140.017951
    ## iter  30 value 56.305038
    ## iter  40 value 36.480755
    ## iter  50 value 19.235894
    ## iter  60 value 12.412830
    ## iter  70 value 11.443929
    ## iter  80 value 10.329738
    ## iter  90 value 8.258823
    ## iter 100 value 7.820812
    ## final  value 7.820812 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1269.326602 
    ## iter  10 value 1090.587798
    ## iter  20 value 1082.607984
    ## iter  30 value 851.100845
    ## iter  40 value 767.486067
    ## iter  50 value 733.062690
    ## iter  60 value 706.368738
    ## iter  70 value 695.225703
    ## iter  80 value 674.440300
    ## iter  90 value 635.411129
    ## iter 100 value 583.954188
    ## final  value 583.954188 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1359.667927 
    ## iter  10 value 716.093437
    ## iter  20 value 566.982635
    ## iter  30 value 503.953295
    ## iter  40 value 459.961144
    ## iter  50 value 430.046369
    ## iter  60 value 411.065285
    ## iter  70 value 385.285736
    ## iter  80 value 374.150449
    ## iter  90 value 368.964782
    ## iter 100 value 366.902571
    ## final  value 366.902571 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1281.149424 
    ## iter  10 value 355.214173
    ## iter  20 value 142.780309
    ## iter  30 value 100.428907
    ## iter  40 value 64.093467
    ## iter  50 value 48.218444
    ## iter  60 value 44.387447
    ## iter  70 value 39.505593
    ## iter  80 value 32.481292
    ## iter  90 value 24.830968
    ## iter 100 value 15.837494
    ## final  value 15.837494 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1181.058573 
    ## iter  10 value 864.652003
    ## iter  20 value 812.163844
    ## iter  30 value 792.195348
    ## iter  40 value 769.175022
    ## iter  50 value 742.254828
    ## iter  60 value 734.096586
    ## iter  70 value 732.045588
    ## iter  80 value 730.225995
    ## iter  90 value 729.743660
    ## iter 100 value 728.492047
    ## final  value 728.492047 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1261.270114 
    ## iter  10 value 691.064366
    ## iter  20 value 494.635376
    ## iter  30 value 423.338846
    ## iter  40 value 385.859876
    ## iter  50 value 358.441158
    ## iter  60 value 341.580966
    ## iter  70 value 333.460688
    ## iter  80 value 327.019331
    ## iter  90 value 321.142691
    ## iter 100 value 313.719869
    ## final  value 313.719869 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1278.759085 
    ## iter  10 value 568.751294
    ## iter  20 value 332.134805
    ## iter  30 value 257.750552
    ## iter  40 value 209.684928
    ## iter  50 value 195.572948
    ## iter  60 value 178.299951
    ## iter  70 value 167.169846
    ## iter  80 value 164.306078
    ## iter  90 value 163.445544
    ## iter 100 value 163.140969
    ## final  value 163.140969 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1244.417859 
    ## iter  10 value 828.873198
    ## iter  20 value 739.227430
    ## iter  30 value 677.248139
    ## iter  40 value 608.038291
    ## iter  50 value 573.259180
    ## iter  60 value 538.373194
    ## iter  70 value 497.202237
    ## iter  80 value 487.209937
    ## iter  90 value 486.095726
    ## iter 100 value 485.647057
    ## final  value 485.647057 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1236.213383 
    ## iter  10 value 680.783705
    ## iter  20 value 444.136745
    ## iter  30 value 370.086069
    ## iter  40 value 355.342845
    ## iter  50 value 337.423740
    ## iter  60 value 329.625336
    ## iter  70 value 317.978607
    ## iter  80 value 316.829158
    ## iter  90 value 314.205356
    ## iter 100 value 313.504473
    ## final  value 313.504473 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1268.034080 
    ## iter  10 value 328.466566
    ## iter  20 value 108.661328
    ## iter  30 value 61.292722
    ## iter  40 value 35.774221
    ## iter  50 value 24.983144
    ## iter  60 value 21.914606
    ## iter  70 value 20.078806
    ## iter  80 value 18.085071
    ## iter  90 value 17.022597
    ## iter 100 value 16.721797
    ## final  value 16.721797 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1320.582541 
    ## iter  10 value 947.847625
    ## iter  20 value 836.294752
    ## iter  30 value 753.201094
    ## iter  40 value 693.942099
    ## iter  50 value 642.702077
    ## iter  60 value 627.277520
    ## iter  70 value 621.397892
    ## iter  80 value 615.872010
    ## iter  90 value 608.197389
    ## iter 100 value 603.629449
    ## final  value 603.629449 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1381.800060 
    ## iter  10 value 872.600138
    ## iter  20 value 762.682161
    ## iter  30 value 675.452709
    ## iter  40 value 615.852862
    ## iter  50 value 473.468877
    ## iter  60 value 383.795504
    ## iter  70 value 322.579780
    ## iter  80 value 294.442035
    ## iter  90 value 267.324939
    ## iter 100 value 250.142853
    ## final  value 250.142853 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1211.911021 
    ## iter  10 value 396.921140
    ## iter  20 value 161.923665
    ## iter  30 value 77.193103
    ## iter  40 value 46.389381
    ## iter  50 value 26.183524
    ## iter  60 value 14.625503
    ## iter  70 value 8.488861
    ## iter  80 value 6.573033
    ## iter  90 value 5.915839
    ## iter 100 value 5.089190
    ## final  value 5.089190 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1275.597378 
    ## iter  10 value 1130.020565
    ## iter  20 value 1123.010633
    ## iter  30 value 928.384478
    ## iter  40 value 857.287743
    ## iter  50 value 828.594438
    ## iter  60 value 800.888166
    ## iter  70 value 787.222165
    ## iter  80 value 783.983129
    ## iter  90 value 778.871327
    ## iter 100 value 773.062341
    ## final  value 773.062341 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1272.301093 
    ## iter  10 value 729.425215
    ## iter  20 value 513.495640
    ## iter  30 value 384.225215
    ## iter  40 value 346.909156
    ## iter  50 value 332.762239
    ## iter  60 value 324.808062
    ## iter  70 value 320.766024
    ## iter  80 value 319.181596
    ## iter  90 value 316.718252
    ## iter 100 value 309.155707
    ## final  value 309.155707 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1315.216545 
    ## iter  10 value 648.973231
    ## iter  20 value 372.712674
    ## iter  30 value 268.314618
    ## iter  40 value 203.177773
    ## iter  50 value 174.701024
    ## iter  60 value 166.060501
    ## iter  70 value 163.177050
    ## iter  80 value 161.097639
    ## iter  90 value 158.546846
    ## iter 100 value 157.527094
    ## final  value 157.527094 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1264.711981 
    ## iter  10 value 901.045185
    ## iter  20 value 779.532521
    ## iter  30 value 707.137461
    ## iter  40 value 646.966345
    ## iter  50 value 597.681019
    ## iter  60 value 547.008386
    ## iter  70 value 526.980536
    ## iter  80 value 525.644455
    ## iter  90 value 524.595608
    ## iter 100 value 523.307287
    ## final  value 523.307287 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1354.878426 
    ## iter  10 value 881.639388
    ## iter  20 value 727.152660
    ## iter  30 value 651.796186
    ## iter  40 value 585.821221
    ## iter  50 value 532.765205
    ## iter  60 value 502.432015
    ## iter  70 value 476.490348
    ## iter  80 value 455.864758
    ## iter  90 value 451.925026
    ## iter 100 value 443.235208
    ## final  value 443.235208 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1368.032112 
    ## iter  10 value 585.348747
    ## iter  20 value 223.947602
    ## iter  30 value 115.644881
    ## iter  40 value 75.320691
    ## iter  50 value 54.116605
    ## iter  60 value 42.579401
    ## iter  70 value 37.804050
    ## iter  80 value 34.375418
    ## iter  90 value 31.353647
    ## iter 100 value 27.706399
    ## final  value 27.706399 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1223.487655 
    ## iter  10 value 1088.894167
    ## iter  20 value 844.160822
    ## iter  30 value 741.106885
    ## iter  40 value 686.725682
    ## iter  50 value 656.517811
    ## iter  60 value 635.323304
    ## iter  70 value 631.891389
    ## iter  80 value 625.006882
    ## iter  90 value 623.279924
    ## iter 100 value 617.433754
    ## final  value 617.433754 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1222.721568 
    ## iter  10 value 698.483819
    ## iter  20 value 498.645419
    ## iter  30 value 412.813837
    ## iter  40 value 364.226187
    ## iter  50 value 347.188085
    ## iter  60 value 328.564014
    ## iter  70 value 318.699648
    ## iter  80 value 307.332555
    ## iter  90 value 294.866030
    ## iter 100 value 289.355162
    ## final  value 289.355162 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1266.388236 
    ## iter  10 value 489.092270
    ## iter  20 value 172.389576
    ## iter  30 value 116.433579
    ## iter  40 value 75.612266
    ## iter  50 value 52.068284
    ## iter  60 value 43.310724
    ## iter  70 value 40.578184
    ## iter  80 value 40.357713
    ## iter  90 value 40.336564
    ## iter 100 value 40.332105
    ## final  value 40.332105 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1258.728830 
    ## iter  10 value 1119.511175
    ## iter  20 value 928.304228
    ## iter  30 value 881.280675
    ## iter  40 value 857.523582
    ## iter  50 value 837.347731
    ## iter  60 value 814.642699
    ## iter  70 value 790.517329
    ## iter  80 value 783.939927
    ## iter  90 value 782.475560
    ## iter 100 value 781.723104
    ## final  value 781.723104 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1312.603705 
    ## iter  10 value 1127.083118
    ## iter  20 value 854.816542
    ## iter  30 value 640.211825
    ## iter  40 value 538.833136
    ## iter  50 value 433.802964
    ## iter  60 value 379.145776
    ## iter  70 value 351.142251
    ## iter  80 value 325.132059
    ## iter  90 value 314.507377
    ## iter 100 value 307.860215
    ## final  value 307.860215 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1362.814381 
    ## iter  10 value 503.427226
    ## iter  20 value 316.503416
    ## iter  30 value 223.008049
    ## iter  40 value 169.532465
    ## iter  50 value 160.068416
    ## iter  60 value 155.874061
    ## iter  70 value 153.389100
    ## iter  80 value 151.912038
    ## iter  90 value 150.755356
    ## iter 100 value 148.272213
    ## final  value 148.272213 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1270.006754 
    ## iter  10 value 937.717393
    ## iter  20 value 769.787201
    ## iter  30 value 695.166443
    ## iter  40 value 636.845153
    ## iter  50 value 570.413242
    ## iter  60 value 529.758953
    ## iter  70 value 507.296477
    ## iter  80 value 498.096025
    ## iter  90 value 492.635194
    ## iter 100 value 486.626195
    ## final  value 486.626195 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1253.720043 
    ## iter  10 value 638.931917
    ## iter  20 value 288.022386
    ## iter  30 value 182.248568
    ## iter  40 value 119.352412
    ## iter  50 value 85.169979
    ## iter  60 value 68.299878
    ## iter  70 value 56.347228
    ## iter  80 value 45.976989
    ## iter  90 value 43.922809
    ## iter 100 value 42.571799
    ## final  value 42.571799 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1251.145213 
    ## iter  10 value 405.573940
    ## iter  20 value 260.225257
    ## iter  30 value 187.990076
    ## iter  40 value 130.707338
    ## iter  50 value 113.613547
    ## iter  60 value 106.181320
    ## iter  70 value 99.304446
    ## iter  80 value 96.650558
    ## iter  90 value 88.857075
    ## iter 100 value 84.167549
    ## final  value 84.167549 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1281.797080 
    ## iter  10 value 1107.150042
    ## iter  20 value 1006.822324
    ## iter  30 value 991.965614
    ## iter  40 value 942.347794
    ## iter  50 value 864.345315
    ## iter  60 value 791.635539
    ## iter  70 value 724.949181
    ## iter  80 value 673.371500
    ## iter  90 value 600.688650
    ## iter 100 value 544.461288
    ## final  value 544.461288 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1215.844948 
    ## iter  10 value 604.815140
    ## iter  20 value 388.929517
    ## iter  30 value 305.024817
    ## iter  40 value 241.655120
    ## iter  50 value 202.580121
    ## iter  60 value 176.504724
    ## iter  70 value 152.404266
    ## iter  80 value 119.778910
    ## iter  90 value 103.018221
    ## iter 100 value 97.250332
    ## final  value 97.250332 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1389.184931 
    ## iter  10 value 548.640684
    ## iter  20 value 135.564738
    ## iter  30 value 47.361071
    ## iter  40 value 30.238867
    ## iter  50 value 16.875414
    ## iter  60 value 13.002651
    ## iter  70 value 11.685121
    ## iter  80 value 8.750991
    ## iter  90 value 7.856949
    ## iter 100 value 7.248775
    ## final  value 7.248775 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1318.196648 
    ## iter  10 value 1064.262692
    ## iter  20 value 988.185906
    ## iter  30 value 925.070836
    ## iter  40 value 856.973600
    ## iter  50 value 824.407808
    ## iter  60 value 800.417488
    ## iter  70 value 789.308757
    ## iter  80 value 781.910693
    ## iter  90 value 774.402155
    ## iter 100 value 771.402814
    ## final  value 771.402814 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1235.897782 
    ## iter  10 value 668.458769
    ## iter  20 value 530.753300
    ## iter  30 value 483.087628
    ## iter  40 value 441.149001
    ## iter  50 value 387.456204
    ## iter  60 value 346.789023
    ## iter  70 value 320.048014
    ## iter  80 value 310.066727
    ## iter  90 value 304.781845
    ## iter 100 value 302.971299
    ## final  value 302.971299 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1278.070596 
    ## iter  10 value 632.096456
    ## iter  20 value 311.691852
    ## iter  30 value 217.642934
    ## iter  40 value 170.743753
    ## iter  50 value 162.921901
    ## iter  60 value 159.361999
    ## iter  70 value 157.851891
    ## iter  80 value 157.180228
    ## iter  90 value 156.678741
    ## iter 100 value 154.040687
    ## final  value 154.040687 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1307.761731 
    ## iter  10 value 804.096701
    ## iter  20 value 727.807354
    ## iter  30 value 685.412884
    ## iter  40 value 656.445943
    ## iter  50 value 632.696626
    ## iter  60 value 603.570570
    ## iter  70 value 597.619367
    ## iter  80 value 594.699276
    ## iter  90 value 593.834197
    ## iter 100 value 593.686318
    ## final  value 593.686318 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1345.598108 
    ## iter  10 value 656.418184
    ## iter  20 value 407.563724
    ## iter  30 value 343.565652
    ## iter  40 value 289.351137
    ## iter  50 value 256.016203
    ## iter  60 value 222.080400
    ## iter  70 value 202.265474
    ## iter  80 value 199.264517
    ## iter  90 value 197.998503
    ## iter 100 value 196.616060
    ## final  value 196.616060 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1301.205542 
    ## iter  10 value 394.174801
    ## iter  20 value 171.982596
    ## iter  30 value 89.051103
    ## iter  40 value 49.190768
    ## iter  50 value 34.295618
    ## iter  60 value 26.874710
    ## iter  70 value 21.174133
    ## iter  80 value 19.555317
    ## iter  90 value 18.595197
    ## iter 100 value 17.778535
    ## final  value 17.778535 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1215.298520 
    ## iter  10 value 881.890114
    ## iter  20 value 818.260097
    ## iter  30 value 795.642314
    ## iter  40 value 774.248515
    ## iter  50 value 767.155957
    ## iter  60 value 766.836153
    ## iter  70 value 766.803986
    ## iter  80 value 766.798553
    ## iter  90 value 766.797664
    ## iter 100 value 766.797290
    ## final  value 766.797290 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1229.774865 
    ## iter  10 value 553.829305
    ## iter  20 value 280.999591
    ## iter  30 value 199.551081
    ## iter  40 value 159.503451
    ## iter  50 value 126.335685
    ## iter  60 value 108.015117
    ## iter  70 value 97.598737
    ## iter  80 value 84.446332
    ## iter  90 value 76.079413
    ## iter 100 value 70.994540
    ## final  value 70.994540 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1299.399914 
    ## iter  10 value 537.154899
    ## iter  20 value 284.717752
    ## iter  30 value 188.311436
    ## iter  40 value 150.149359
    ## iter  50 value 126.363810
    ## iter  60 value 109.277153
    ## iter  70 value 99.228342
    ## iter  80 value 90.437716
    ## iter  90 value 71.985665
    ## iter 100 value 60.439976
    ## final  value 60.439976 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1291.818587 
    ## iter  10 value 1156.088797
    ## iter  20 value 1134.663423
    ## iter  30 value 983.230141
    ## iter  40 value 912.755526
    ## iter  50 value 878.939896
    ## iter  60 value 859.749145
    ## iter  70 value 843.503217
    ## iter  80 value 825.786687
    ## iter  90 value 805.719250
    ## iter 100 value 800.645839
    ## final  value 800.645839 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1240.326111 
    ## iter  10 value 656.136044
    ## iter  20 value 509.245139
    ## iter  30 value 423.154379
    ## iter  40 value 388.524987
    ## iter  50 value 358.524182
    ## iter  60 value 337.782511
    ## iter  70 value 324.923681
    ## iter  80 value 319.397051
    ## iter  90 value 318.450724
    ## iter 100 value 317.461046
    ## final  value 317.461046 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1333.519294 
    ## iter  10 value 564.714124
    ## iter  20 value 325.933895
    ## iter  30 value 214.527916
    ## iter  40 value 177.803379
    ## iter  50 value 164.400699
    ## iter  60 value 154.117872
    ## iter  70 value 150.203699
    ## iter  80 value 148.281827
    ## iter  90 value 147.726552
    ## iter 100 value 147.233966
    ## final  value 147.233966 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1251.232050 
    ## iter  10 value 952.116919
    ## iter  20 value 867.299211
    ## iter  30 value 836.332227
    ## iter  40 value 819.729674
    ## iter  50 value 802.722040
    ## iter  60 value 796.048826
    ## iter  70 value 784.359536
    ## iter  80 value 779.078775
    ## iter  90 value 772.452880
    ## iter 100 value 771.582332
    ## final  value 771.582332 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1365.398124 
    ## iter  10 value 526.876032
    ## iter  20 value 337.017794
    ## iter  30 value 275.877467
    ## iter  40 value 224.974634
    ## iter  50 value 178.118965
    ## iter  60 value 147.083325
    ## iter  70 value 128.306061
    ## iter  80 value 107.862119
    ## iter  90 value 95.895698
    ## iter 100 value 94.607830
    ## final  value 94.607830 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1358.878508 
    ## iter  10 value 330.402608
    ## iter  20 value 118.064387
    ## iter  30 value 58.639937
    ## iter  40 value 27.721689
    ## iter  50 value 17.918808
    ## iter  60 value 11.968494
    ## iter  70 value 10.399488
    ## iter  80 value 7.311370
    ## iter  90 value 3.951149
    ## iter 100 value 3.259953
    ## final  value 3.259953 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1357.175378 
    ## iter  10 value 1110.708233
    ## iter  20 value 1103.120745
    ## final  value 1103.079427 
    ## converged
    ## # weights:  255
    ## initial  value 1226.069529 
    ## iter  10 value 637.031901
    ## iter  20 value 341.946013
    ## iter  30 value 275.310229
    ## iter  40 value 232.701052
    ## iter  50 value 198.421124
    ## iter  60 value 171.612680
    ## iter  70 value 144.551362
    ## iter  80 value 123.796735
    ## iter  90 value 115.368291
    ## iter 100 value 108.238796
    ## final  value 108.238796 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1341.474441 
    ## iter  10 value 672.322893
    ## iter  20 value 317.169543
    ## iter  30 value 245.996060
    ## iter  40 value 184.307335
    ## iter  50 value 130.682535
    ## iter  60 value 111.563676
    ## iter  70 value 99.097329
    ## iter  80 value 85.391014
    ## iter  90 value 80.762768
    ## iter 100 value 76.985199
    ## final  value 76.985199 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1321.967340 
    ## iter  10 value 934.001999
    ## iter  20 value 826.034379
    ## iter  30 value 794.299152
    ## iter  40 value 784.627256
    ## iter  50 value 779.144180
    ## iter  60 value 773.512806
    ## iter  70 value 767.826239
    ## iter  80 value 765.230689
    ## iter  90 value 754.376412
    ## iter 100 value 751.867153
    ## final  value 751.867153 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1276.557219 
    ## iter  10 value 748.490151
    ## iter  20 value 515.991832
    ## iter  30 value 456.022659
    ## iter  40 value 427.323022
    ## iter  50 value 396.837486
    ## iter  60 value 378.597372
    ## iter  70 value 358.991355
    ## iter  80 value 338.999744
    ## iter  90 value 324.346693
    ## iter 100 value 319.239415
    ## final  value 319.239415 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1326.134280 
    ## iter  10 value 547.877089
    ## iter  20 value 291.353484
    ## iter  30 value 214.780220
    ## iter  40 value 181.403665
    ## iter  50 value 168.947535
    ## iter  60 value 165.270247
    ## iter  70 value 162.639297
    ## iter  80 value 160.264249
    ## iter  90 value 159.914917
    ## iter 100 value 159.740244
    ## final  value 159.740244 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1307.447878 
    ## iter  10 value 947.206476
    ## iter  20 value 754.913052
    ## iter  30 value 679.545952
    ## iter  40 value 611.900741
    ## iter  50 value 540.263979
    ## iter  60 value 481.571031
    ## iter  70 value 443.143272
    ## iter  80 value 422.101635
    ## iter  90 value 413.840028
    ## iter 100 value 404.682318
    ## final  value 404.682318 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1259.774684 
    ## iter  10 value 641.917473
    ## iter  20 value 365.687504
    ## iter  30 value 309.538614
    ## iter  40 value 259.179747
    ## iter  50 value 217.026849
    ## iter  60 value 183.863359
    ## iter  70 value 166.782868
    ## iter  80 value 152.105813
    ## iter  90 value 141.617269
    ## iter 100 value 127.985677
    ## final  value 127.985677 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1277.745980 
    ## iter  10 value 449.324833
    ## iter  20 value 196.457065
    ## iter  30 value 144.481717
    ## iter  40 value 112.406765
    ## iter  50 value 87.593176
    ## iter  60 value 65.918220
    ## iter  70 value 57.095620
    ## iter  80 value 54.593244
    ## iter  90 value 53.757275
    ## iter 100 value 50.930256
    ## final  value 50.930256 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1279.185449 
    ## iter  10 value 1023.938752
    ## iter  20 value 887.972540
    ## iter  30 value 798.949512
    ## iter  40 value 746.008738
    ## iter  50 value 698.774961
    ## iter  60 value 671.409437
    ## iter  70 value 661.127678
    ## iter  80 value 625.864727
    ## iter  90 value 591.584655
    ## iter 100 value 571.329994
    ## final  value 571.329994 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1339.931488 
    ## iter  10 value 688.882466
    ## iter  20 value 581.633536
    ## iter  30 value 461.729894
    ## iter  40 value 412.259477
    ## iter  50 value 369.194407
    ## iter  60 value 343.170248
    ## iter  70 value 316.509548
    ## iter  80 value 285.922405
    ## iter  90 value 255.272747
    ## iter 100 value 228.040020
    ## final  value 228.040020 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1406.762302 
    ## iter  10 value 470.398759
    ## iter  20 value 216.587184
    ## iter  30 value 132.273107
    ## iter  40 value 88.729095
    ## iter  50 value 61.619497
    ## iter  60 value 49.125134
    ## iter  70 value 46.608979
    ## iter  80 value 42.055395
    ## iter  90 value 38.780723
    ## iter 100 value 36.655667
    ## final  value 36.655667 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1270.557704 
    ## iter  10 value 903.171524
    ## iter  20 value 854.351379
    ## iter  30 value 820.928498
    ## iter  40 value 800.398729
    ## iter  50 value 785.102890
    ## iter  60 value 776.835427
    ## iter  70 value 771.295259
    ## iter  80 value 767.059196
    ## iter  90 value 766.549943
    ## iter 100 value 766.524180
    ## final  value 766.524180 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1313.241939 
    ## iter  10 value 787.031326
    ## iter  20 value 496.303073
    ## iter  30 value 402.336475
    ## iter  40 value 369.530500
    ## iter  50 value 352.031379
    ## iter  60 value 340.121519
    ## iter  70 value 331.408375
    ## iter  80 value 324.984911
    ## iter  90 value 319.849761
    ## iter 100 value 316.178917
    ## final  value 316.178917 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1281.103580 
    ## iter  10 value 532.290959
    ## iter  20 value 290.209973
    ## iter  30 value 227.605850
    ## iter  40 value 194.299678
    ## iter  50 value 166.103390
    ## iter  60 value 158.364921
    ## iter  70 value 154.604268
    ## iter  80 value 150.954270
    ## iter  90 value 148.164058
    ## iter 100 value 147.135212
    ## final  value 147.135212 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1309.762152 
    ## iter  10 value 899.387391
    ## iter  20 value 792.551426
    ## iter  30 value 716.278060
    ## iter  40 value 671.332927
    ## iter  50 value 631.956058
    ## iter  60 value 615.969988
    ## iter  70 value 608.340405
    ## iter  80 value 606.825562
    ## iter  90 value 606.394045
    ## iter 100 value 605.191919
    ## final  value 605.191919 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1248.791463 
    ## iter  10 value 585.891240
    ## iter  20 value 308.768005
    ## iter  30 value 223.316389
    ## iter  40 value 165.756957
    ## iter  50 value 134.723101
    ## iter  60 value 117.972265
    ## iter  70 value 106.666227
    ## iter  80 value 98.960577
    ## iter  90 value 93.044113
    ## iter 100 value 90.046973
    ## final  value 90.046973 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1220.927301 
    ## iter  10 value 297.630211
    ## iter  20 value 92.167776
    ## iter  30 value 60.444734
    ## iter  40 value 46.469139
    ## iter  50 value 36.884363
    ## iter  60 value 33.760950
    ## iter  70 value 31.288498
    ## iter  80 value 29.048859
    ## iter  90 value 28.362028
    ## iter 100 value 27.792547
    ## final  value 27.792547 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1272.537193 
    ## iter  10 value 1100.711776
    ## iter  20 value 964.259295
    ## iter  30 value 826.318074
    ## iter  40 value 782.925842
    ## iter  50 value 762.259269
    ## iter  60 value 754.923240
    ## iter  70 value 752.153614
    ## iter  80 value 747.073154
    ## iter  90 value 746.281905
    ## iter 100 value 746.034674
    ## final  value 746.034674 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1321.740161 
    ## iter  10 value 816.027653
    ## iter  20 value 479.433729
    ## iter  30 value 415.179657
    ## iter  40 value 350.326406
    ## iter  50 value 275.408483
    ## iter  60 value 220.954505
    ## iter  70 value 185.620041
    ## iter  80 value 159.364009
    ## iter  90 value 145.651040
    ## iter 100 value 135.835092
    ## final  value 135.835092 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1326.143575 
    ## iter  10 value 634.162913
    ## iter  20 value 373.016805
    ## iter  30 value 201.083834
    ## iter  40 value 124.530632
    ## iter  50 value 104.558892
    ## iter  60 value 98.166598
    ## iter  70 value 87.415070
    ## iter  80 value 84.419761
    ## iter  90 value 80.820606
    ## iter 100 value 76.912366
    ## final  value 76.912366 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1221.712917 
    ## iter  10 value 893.517325
    ## iter  20 value 840.147754
    ## iter  30 value 811.041177
    ## iter  40 value 785.898725
    ## iter  50 value 773.937285
    ## iter  60 value 765.104307
    ## iter  70 value 752.223396
    ## iter  80 value 749.472672
    ## iter  90 value 747.525413
    ## iter 100 value 747.302169
    ## final  value 747.302169 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1376.320424 
    ## iter  10 value 755.974034
    ## iter  20 value 641.924519
    ## iter  30 value 582.340904
    ## iter  40 value 531.697336
    ## iter  50 value 435.664498
    ## iter  60 value 383.578212
    ## iter  70 value 347.927659
    ## iter  80 value 322.064452
    ## iter  90 value 307.954246
    ## iter 100 value 300.823314
    ## final  value 300.823314 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1267.677108 
    ## iter  10 value 479.766922
    ## iter  20 value 265.593395
    ## iter  30 value 196.337296
    ## iter  40 value 177.492405
    ## iter  50 value 162.160445
    ## iter  60 value 152.352052
    ## iter  70 value 148.816476
    ## iter  80 value 147.267018
    ## iter  90 value 146.710522
    ## iter 100 value 146.481229
    ## final  value 146.481229 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1199.430343 
    ## iter  10 value 940.528500
    ## iter  20 value 790.692476
    ## iter  30 value 734.363999
    ## iter  40 value 696.647324
    ## iter  50 value 662.437885
    ## iter  60 value 652.063351
    ## iter  70 value 650.832427
    ## iter  80 value 650.373564
    ## iter  90 value 650.041703
    ## iter 100 value 649.734314
    ## final  value 649.734314 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1167.380391 
    ## iter  10 value 682.226093
    ## iter  20 value 352.104380
    ## iter  30 value 257.816648
    ## iter  40 value 184.216355
    ## iter  50 value 152.113452
    ## iter  60 value 136.897974
    ## iter  70 value 132.317106
    ## iter  80 value 131.522160
    ## iter  90 value 131.052273
    ## iter 100 value 130.383509
    ## final  value 130.383509 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1263.956353 
    ## iter  10 value 401.628055
    ## iter  20 value 113.355191
    ## iter  30 value 72.760792
    ## iter  40 value 56.247480
    ## iter  50 value 49.279690
    ## iter  60 value 38.916682
    ## iter  70 value 28.888686
    ## iter  80 value 16.558287
    ## iter  90 value 6.601822
    ## iter 100 value 4.144747
    ## final  value 4.144747 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1266.929748 
    ## iter  10 value 918.819422
    ## iter  20 value 790.329442
    ## iter  30 value 769.604108
    ## iter  40 value 741.996875
    ## iter  50 value 725.545975
    ## iter  60 value 722.055355
    ## iter  70 value 719.760669
    ## iter  80 value 715.502428
    ## iter  90 value 713.809412
    ## iter 100 value 713.075739
    ## final  value 713.075739 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1286.242000 
    ## iter  10 value 658.022372
    ## iter  20 value 385.721578
    ## iter  30 value 319.321362
    ## iter  40 value 285.718196
    ## iter  50 value 269.434226
    ## iter  60 value 257.232557
    ## iter  70 value 245.449107
    ## iter  80 value 235.349674
    ## iter  90 value 226.453138
    ## iter 100 value 220.925453
    ## final  value 220.925453 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1261.497607 
    ## iter  10 value 484.788441
    ## iter  20 value 149.996343
    ## iter  30 value 93.202944
    ## iter  40 value 71.849146
    ## iter  50 value 50.702066
    ## iter  60 value 36.948322
    ## iter  70 value 29.634788
    ## iter  80 value 21.949591
    ## iter  90 value 15.027283
    ## iter 100 value 9.993216
    ## final  value 9.993216 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1273.246812 
    ## iter  10 value 922.449623
    ## iter  20 value 879.441448
    ## iter  30 value 813.200102
    ## iter  40 value 792.195958
    ## iter  50 value 765.877673
    ## iter  60 value 749.116327
    ## iter  70 value 745.066645
    ## iter  80 value 744.796732
    ## iter  90 value 744.767691
    ## final  value 744.767346 
    ## converged
    ## # weights:  255
    ## initial  value 1391.582584 
    ## iter  10 value 757.868897
    ## iter  20 value 550.143531
    ## iter  30 value 451.764709
    ## iter  40 value 400.484084
    ## iter  50 value 384.756631
    ## iter  60 value 371.079408
    ## iter  70 value 359.747094
    ## iter  80 value 342.239195
    ## iter  90 value 334.886592
    ## iter 100 value 332.562494
    ## final  value 332.562494 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1283.882545 
    ## iter  10 value 486.331571
    ## iter  20 value 246.181005
    ## iter  30 value 204.759881
    ## iter  40 value 181.218250
    ## iter  50 value 174.659530
    ## iter  60 value 168.102336
    ## iter  70 value 165.220975
    ## iter  80 value 162.169841
    ## iter  90 value 159.460787
    ## iter 100 value 157.685543
    ## final  value 157.685543 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1249.961046 
    ## iter  10 value 1091.928637
    ## iter  20 value 873.601054
    ## iter  30 value 748.928549
    ## iter  40 value 707.117702
    ## iter  50 value 693.908019
    ## iter  60 value 683.150561
    ## iter  70 value 676.804016
    ## iter  80 value 674.960647
    ## iter  90 value 673.387383
    ## iter 100 value 670.830034
    ## final  value 670.830034 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1345.352253 
    ## iter  10 value 582.753491
    ## iter  20 value 366.547255
    ## iter  30 value 303.180827
    ## iter  40 value 213.286419
    ## iter  50 value 155.561018
    ## iter  60 value 125.937433
    ## iter  70 value 102.334038
    ## iter  80 value 81.243331
    ## iter  90 value 65.055478
    ## iter 100 value 56.879292
    ## final  value 56.879292 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1329.262621 
    ## iter  10 value 501.858197
    ## iter  20 value 116.877812
    ## iter  30 value 41.422679
    ## iter  40 value 17.199710
    ## iter  50 value 7.605395
    ## iter  60 value 3.133873
    ## iter  70 value 2.529855
    ## iter  80 value 2.253875
    ## iter  90 value 2.161360
    ## iter 100 value 2.015356
    ## final  value 2.015356 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1252.700798 
    ## iter  10 value 834.750632
    ## iter  20 value 770.337092
    ## iter  30 value 733.490342
    ## iter  40 value 693.403314
    ## iter  50 value 683.954754
    ## iter  60 value 668.498601
    ## iter  70 value 657.768551
    ## iter  80 value 651.996313
    ## iter  90 value 650.378861
    ## iter 100 value 648.819546
    ## final  value 648.819546 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1327.010226 
    ## iter  10 value 1050.203631
    ## iter  20 value 843.057910
    ## iter  30 value 761.641470
    ## iter  40 value 674.104182
    ## iter  50 value 552.483228
    ## iter  60 value 493.450036
    ## iter  70 value 477.717589
    ## iter  80 value 437.785478
    ## iter  90 value 338.597772
    ## iter 100 value 308.055939
    ## final  value 308.055939 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1309.786318 
    ## iter  10 value 456.742811
    ## iter  20 value 183.128186
    ## iter  30 value 92.937658
    ## iter  40 value 60.567269
    ## iter  50 value 47.930641
    ## iter  60 value 38.624803
    ## iter  70 value 34.098391
    ## iter  80 value 33.673551
    ## iter  90 value 33.646746
    ## iter 100 value 33.642216
    ## final  value 33.642216 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1279.587662 
    ## iter  10 value 1001.602752
    ## iter  20 value 886.453622
    ## iter  30 value 780.868656
    ## iter  40 value 747.713273
    ## iter  50 value 730.160087
    ## iter  60 value 724.547194
    ## iter  70 value 723.003726
    ## iter  80 value 722.655431
    ## iter  90 value 722.645669
    ## final  value 722.645620 
    ## converged
    ## # weights:  255
    ## initial  value 1237.855313 
    ## iter  10 value 717.786784
    ## iter  20 value 505.912625
    ## iter  30 value 392.982450
    ## iter  40 value 344.108304
    ## iter  50 value 328.639121
    ## iter  60 value 318.226850
    ## iter  70 value 309.543964
    ## iter  80 value 304.578185
    ## iter  90 value 299.522655
    ## iter 100 value 294.251519
    ## final  value 294.251519 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1284.012945 
    ## iter  10 value 479.133670
    ## iter  20 value 277.283854
    ## iter  30 value 203.213850
    ## iter  40 value 177.840125
    ## iter  50 value 164.933238
    ## iter  60 value 158.511137
    ## iter  70 value 154.247466
    ## iter  80 value 152.198580
    ## iter  90 value 151.258131
    ## iter 100 value 150.207698
    ## final  value 150.207698 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1265.696271 
    ## iter  10 value 898.544831
    ## iter  20 value 784.635458
    ## iter  30 value 753.697327
    ## iter  40 value 729.630786
    ## iter  50 value 718.188788
    ## iter  60 value 709.863492
    ## iter  70 value 708.288306
    ## iter  80 value 706.826658
    ## iter  90 value 706.333103
    ## iter 100 value 706.017478
    ## final  value 706.017478 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1315.237268 
    ## iter  10 value 887.568722
    ## iter  20 value 600.028017
    ## iter  30 value 521.816084
    ## iter  40 value 459.428144
    ## iter  50 value 420.060632
    ## iter  60 value 387.539733
    ## iter  70 value 352.408945
    ## iter  80 value 324.668004
    ## iter  90 value 302.504241
    ## iter 100 value 285.559458
    ## final  value 285.559458 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1220.839559 
    ## iter  10 value 466.237576
    ## iter  20 value 147.070520
    ## iter  30 value 52.632066
    ## iter  40 value 23.376808
    ## iter  50 value 10.390539
    ## iter  60 value 6.412176
    ## iter  70 value 5.982784
    ## iter  80 value 5.711948
    ## iter  90 value 5.478539
    ## iter 100 value 5.357803
    ## final  value 5.357803 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1256.126907 
    ## iter  10 value 1063.086142
    ## iter  20 value 973.321528
    ## iter  30 value 889.746935
    ## iter  40 value 875.153782
    ## iter  50 value 851.230055
    ## iter  60 value 838.405736
    ## iter  70 value 829.340553
    ## iter  80 value 822.592462
    ## iter  90 value 816.379928
    ## iter 100 value 813.793440
    ## final  value 813.793440 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1497.384798 
    ## iter  10 value 1075.601435
    ## iter  20 value 1065.925481
    ## iter  30 value 1065.822405
    ## iter  40 value 959.275168
    ## iter  50 value 789.367227
    ## iter  60 value 757.294661
    ## iter  70 value 730.839848
    ## iter  80 value 716.245259
    ## iter  90 value 706.042261
    ## iter 100 value 648.662386
    ## final  value 648.662386 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1346.893259 
    ## iter  10 value 510.604757
    ## iter  20 value 170.418055
    ## iter  30 value 82.012035
    ## iter  40 value 41.859423
    ## iter  50 value 21.454935
    ## iter  60 value 13.887191
    ## iter  70 value 12.518702
    ## iter  80 value 11.528243
    ## iter  90 value 11.053772
    ## iter 100 value 11.034398
    ## final  value 11.034398 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1283.688728 
    ## iter  10 value 948.086296
    ## iter  20 value 815.483215
    ## iter  30 value 774.020460
    ## iter  40 value 747.713864
    ## iter  50 value 738.041851
    ## iter  60 value 734.007730
    ## iter  70 value 729.777709
    ## iter  80 value 728.754546
    ## iter  90 value 720.553116
    ## iter 100 value 711.618665
    ## final  value 711.618665 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1385.111225 
    ## iter  10 value 726.211511
    ## iter  20 value 523.795056
    ## iter  30 value 438.253096
    ## iter  40 value 394.621699
    ## iter  50 value 376.894423
    ## iter  60 value 362.692106
    ## iter  70 value 350.112960
    ## iter  80 value 339.108952
    ## iter  90 value 329.917390
    ## iter 100 value 323.311664
    ## final  value 323.311664 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1242.594738 
    ## iter  10 value 612.170242
    ## iter  20 value 327.425366
    ## iter  30 value 221.356575
    ## iter  40 value 177.748949
    ## iter  50 value 171.524275
    ## iter  60 value 170.191023
    ## iter  70 value 169.080079
    ## iter  80 value 166.101117
    ## iter  90 value 164.650422
    ## iter 100 value 162.274483
    ## final  value 162.274483 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1341.756093 
    ## iter  10 value 1071.356813
    ## iter  20 value 1065.860524
    ## iter  30 value 1065.826198
    ## final  value 1065.824672 
    ## converged
    ## # weights:  255
    ## initial  value 1347.616745 
    ## iter  10 value 561.110268
    ## iter  20 value 266.285791
    ## iter  30 value 201.760221
    ## iter  40 value 163.512765
    ## iter  50 value 142.454891
    ## iter  60 value 133.610435
    ## iter  70 value 130.240243
    ## iter  80 value 128.470927
    ## iter  90 value 127.587368
    ## iter 100 value 125.116018
    ## final  value 125.116018 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1376.766692 
    ## iter  10 value 353.540058
    ## iter  20 value 88.188455
    ## iter  30 value 37.343598
    ## iter  40 value 17.433078
    ## iter  50 value 11.496893
    ## iter  60 value 8.148716
    ## iter  70 value 7.361212
    ## iter  80 value 6.892476
    ## iter  90 value 4.925369
    ## iter 100 value 3.747150
    ## final  value 3.747150 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1358.522109 
    ## iter  10 value 1097.916642
    ## iter  20 value 1085.806795
    ## iter  30 value 1085.791134
    ## final  value 1085.791031 
    ## converged
    ## # weights:  255
    ## initial  value 1280.394804 
    ## iter  10 value 441.835968
    ## iter  20 value 285.196761
    ## iter  30 value 213.825131
    ## iter  40 value 177.580757
    ## iter  50 value 157.226798
    ## iter  60 value 143.759727
    ## iter  70 value 133.688539
    ## iter  80 value 129.481321
    ## iter  90 value 127.912391
    ## iter 100 value 126.154573
    ## final  value 126.154573 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1373.690724 
    ## iter  10 value 584.585520
    ## iter  20 value 315.799590
    ## iter  30 value 203.909482
    ## iter  40 value 157.162474
    ## iter  50 value 122.019264
    ## iter  60 value 96.009621
    ## iter  70 value 75.954889
    ## iter  80 value 64.314246
    ## iter  90 value 51.777970
    ## iter 100 value 43.565265
    ## final  value 43.565265 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1194.382428 
    ## iter  10 value 904.855925
    ## iter  20 value 832.768982
    ## iter  30 value 806.251247
    ## iter  40 value 778.628430
    ## iter  50 value 752.807232
    ## iter  60 value 741.715609
    ## iter  70 value 737.194237
    ## iter  80 value 735.946514
    ## iter  90 value 735.513229
    ## iter 100 value 735.450876
    ## final  value 735.450876 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1336.863852 
    ## iter  10 value 677.255380
    ## iter  20 value 497.549282
    ## iter  30 value 448.815639
    ## iter  40 value 387.753126
    ## iter  50 value 346.258683
    ## iter  60 value 327.893539
    ## iter  70 value 317.531280
    ## iter  80 value 312.704643
    ## iter  90 value 309.929344
    ## iter 100 value 307.958764
    ## final  value 307.958764 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1378.962351 
    ## iter  10 value 546.814467
    ## iter  20 value 292.510658
    ## iter  30 value 211.495413
    ## iter  40 value 179.610799
    ## iter  50 value 164.558116
    ## iter  60 value 160.337977
    ## iter  70 value 158.505741
    ## iter  80 value 156.550701
    ## iter  90 value 154.647765
    ## iter 100 value 153.758215
    ## final  value 153.758215 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1321.844679 
    ## iter  10 value 1016.839342
    ## iter  20 value 809.174151
    ## iter  30 value 783.131715
    ## iter  40 value 749.625861
    ## iter  50 value 707.656531
    ## iter  60 value 687.642148
    ## iter  70 value 663.448487
    ## iter  80 value 659.770106
    ## iter  90 value 658.318051
    ## iter 100 value 657.955155
    ## final  value 657.955155 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1332.946434 
    ## iter  10 value 932.299741
    ## iter  20 value 772.340826
    ## iter  30 value 606.870455
    ## iter  40 value 509.036650
    ## iter  50 value 430.868725
    ## iter  60 value 343.960555
    ## iter  70 value 321.103655
    ## iter  80 value 308.250067
    ## iter  90 value 301.265595
    ## iter 100 value 295.993733
    ## final  value 295.993733 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1319.163561 
    ## iter  10 value 402.922241
    ## iter  20 value 130.771145
    ## iter  30 value 72.497548
    ## iter  40 value 34.324145
    ## iter  50 value 20.380630
    ## iter  60 value 10.299569
    ## iter  70 value 9.341608
    ## iter  80 value 8.738702
    ## iter  90 value 8.407663
    ## iter 100 value 8.149503
    ## final  value 8.149503 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1254.091848 
    ## iter  10 value 885.594419
    ## iter  20 value 781.127429
    ## iter  30 value 722.021419
    ## iter  40 value 674.573025
    ## iter  50 value 631.588656
    ## iter  60 value 600.331543
    ## iter  70 value 590.272569
    ## iter  80 value 582.452132
    ## iter  90 value 579.749302
    ## iter 100 value 578.569951
    ## final  value 578.569951 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1300.697920 
    ## iter  10 value 796.306958
    ## iter  20 value 511.336276
    ## iter  30 value 414.215068
    ## iter  40 value 316.234890
    ## iter  50 value 269.408137
    ## iter  60 value 230.659719
    ## iter  70 value 207.246656
    ## iter  80 value 197.269188
    ## iter  90 value 193.213104
    ## iter 100 value 184.398937
    ## final  value 184.398937 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1341.427752 
    ## iter  10 value 1150.598206
    ## iter  20 value 601.980677
    ## iter  30 value 275.029480
    ## iter  40 value 123.130352
    ## iter  50 value 58.818708
    ## iter  60 value 37.949260
    ## iter  70 value 27.180847
    ## iter  80 value 19.038751
    ## iter  90 value 15.461289
    ## iter 100 value 13.371333
    ## final  value 13.371333 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1212.650666 
    ## iter  10 value 1069.521126
    ## iter  20 value 979.490049
    ## iter  30 value 956.031561
    ## iter  40 value 875.010077
    ## iter  50 value 832.335985
    ## iter  60 value 814.990904
    ## iter  70 value 807.051090
    ## iter  80 value 785.999162
    ## iter  90 value 770.294712
    ## iter 100 value 765.780585
    ## final  value 765.780585 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1309.965437 
    ## iter  10 value 907.842480
    ## iter  20 value 584.723737
    ## iter  30 value 467.191805
    ## iter  40 value 406.009593
    ## iter  50 value 367.099587
    ## iter  60 value 329.805258
    ## iter  70 value 316.723927
    ## iter  80 value 308.488908
    ## iter  90 value 304.020165
    ## iter 100 value 299.229094
    ## final  value 299.229094 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1260.865769 
    ## iter  10 value 487.007762
    ## iter  20 value 315.102514
    ## iter  30 value 228.422488
    ## iter  40 value 207.146948
    ## iter  50 value 179.132309
    ## iter  60 value 163.510642
    ## iter  70 value 160.509924
    ## iter  80 value 152.385243
    ## iter  90 value 147.772467
    ## iter 100 value 147.226111
    ## final  value 147.226111 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1292.142234 
    ## iter  10 value 952.449207
    ## iter  20 value 861.570667
    ## iter  30 value 825.386138
    ## iter  40 value 791.159362
    ## iter  50 value 771.455136
    ## iter  60 value 765.603050
    ## iter  70 value 765.359565
    ## iter  80 value 765.145997
    ## iter  90 value 765.015763
    ## iter 100 value 764.396497
    ## final  value 764.396497 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1364.838565 
    ## iter  10 value 840.290305
    ## iter  20 value 662.669262
    ## iter  30 value 558.597984
    ## iter  40 value 474.297474
    ## iter  50 value 414.907029
    ## iter  60 value 367.065884
    ## iter  70 value 336.522627
    ## iter  80 value 313.893623
    ## iter  90 value 291.391265
    ## iter 100 value 274.686445
    ## final  value 274.686445 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1310.656037 
    ## iter  10 value 435.080375
    ## iter  20 value 194.623182
    ## iter  30 value 98.625041
    ## iter  40 value 17.251443
    ## iter  50 value 5.554707
    ## iter  60 value 4.501041
    ## iter  70 value 4.316397
    ## iter  80 value 4.156391
    ## iter  90 value 4.040729
    ## iter 100 value 3.922804
    ## final  value 3.922804 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1219.219135 
    ## iter  10 value 874.664248
    ## iter  20 value 810.674070
    ## iter  30 value 792.127847
    ## iter  40 value 761.024446
    ## iter  50 value 729.316352
    ## iter  60 value 724.265696
    ## iter  70 value 721.157274
    ## iter  80 value 718.860250
    ## iter  90 value 713.272805
    ## iter 100 value 710.788560
    ## final  value 710.788560 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1326.040160 
    ## iter  10 value 629.783809
    ## iter  20 value 295.819790
    ## iter  30 value 218.527110
    ## iter  40 value 171.353985
    ## iter  50 value 144.205275
    ## iter  60 value 126.718022
    ## iter  70 value 114.473876
    ## iter  80 value 110.157992
    ## iter  90 value 107.069354
    ## iter 100 value 104.637976
    ## final  value 104.637976 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1274.575790 
    ## iter  10 value 405.973889
    ## iter  20 value 132.082176
    ## iter  30 value 51.324929
    ## iter  40 value 31.190243
    ## iter  50 value 23.690959
    ## iter  60 value 20.706620
    ## iter  70 value 20.459663
    ## iter  80 value 20.441850
    ## iter  90 value 20.439026
    ## iter 100 value 20.438376
    ## final  value 20.438376 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1312.153692 
    ## iter  10 value 957.735349
    ## iter  20 value 828.863850
    ## iter  30 value 795.055664
    ## iter  40 value 784.358946
    ## iter  50 value 778.150372
    ## iter  60 value 767.631007
    ## iter  70 value 759.343521
    ## iter  80 value 758.446191
    ## iter  90 value 757.780095
    ## iter 100 value 757.719218
    ## final  value 757.719218 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1422.186929 
    ## iter  10 value 739.839910
    ## iter  20 value 538.335075
    ## iter  30 value 428.505766
    ## iter  40 value 384.565419
    ## iter  50 value 368.175443
    ## iter  60 value 355.591535
    ## iter  70 value 345.110291
    ## iter  80 value 338.250392
    ## iter  90 value 329.166877
    ## iter 100 value 324.480521
    ## final  value 324.480521 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1269.510457 
    ## iter  10 value 609.391922
    ## iter  20 value 341.568760
    ## iter  30 value 226.679628
    ## iter  40 value 185.589940
    ## iter  50 value 170.047641
    ## iter  60 value 164.453134
    ## iter  70 value 161.713380
    ## iter  80 value 161.041531
    ## iter  90 value 158.814622
    ## iter 100 value 157.505637
    ## final  value 157.505637 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1340.263219 
    ## iter  10 value 1124.254324
    ## iter  20 value 1114.890242
    ## iter  30 value 1114.883229
    ## iter  40 value 1114.881845
    ## final  value 1114.881676 
    ## converged
    ## # weights:  255
    ## initial  value 1242.217498 
    ## iter  10 value 668.048569
    ## iter  20 value 348.574823
    ## iter  30 value 265.232222
    ## iter  40 value 216.939347
    ## iter  50 value 173.416795
    ## iter  60 value 146.167934
    ## iter  70 value 139.457048
    ## iter  80 value 133.478514
    ## iter  90 value 129.930275
    ## iter 100 value 127.075023
    ## final  value 127.075023 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1228.648816 
    ## iter  10 value 355.092626
    ## iter  20 value 104.798105
    ## iter  30 value 53.207172
    ## iter  40 value 28.289847
    ## iter  50 value 15.976823
    ## iter  60 value 11.394829
    ## iter  70 value 6.231227
    ## iter  80 value 3.953085
    ## iter  90 value 2.870525
    ## iter 100 value 2.533473
    ## final  value 2.533473 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1201.746295 
    ## iter  10 value 915.193841
    ## iter  20 value 835.885196
    ## iter  30 value 793.228956
    ## iter  40 value 742.250843
    ## iter  50 value 698.975819
    ## iter  60 value 664.284058
    ## iter  70 value 637.208039
    ## iter  80 value 627.160246
    ## iter  90 value 626.198880
    ## iter 100 value 625.943837
    ## final  value 625.943837 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1339.562683 
    ## iter  10 value 866.863971
    ## iter  20 value 792.026772
    ## iter  30 value 758.980216
    ## iter  40 value 723.578741
    ## iter  50 value 707.547717
    ## iter  60 value 691.322798
    ## iter  70 value 681.987746
    ## iter  80 value 676.908392
    ## iter  90 value 674.340944
    ## iter 100 value 673.641835
    ## final  value 673.641835 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1340.085813 
    ## iter  10 value 538.948143
    ## iter  20 value 211.877824
    ## iter  30 value 128.189229
    ## iter  40 value 98.760333
    ## iter  50 value 68.469166
    ## iter  60 value 51.200020
    ## iter  70 value 41.516536
    ## iter  80 value 32.917373
    ## iter  90 value 30.379306
    ## iter 100 value 29.602909
    ## final  value 29.602909 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1324.672729 
    ## iter  10 value 1046.302119
    ## iter  20 value 925.623548
    ## iter  30 value 830.956142
    ## iter  40 value 798.669851
    ## iter  50 value 783.836058
    ## iter  60 value 777.650552
    ## iter  70 value 772.454916
    ## iter  80 value 770.565927
    ## iter  90 value 770.026312
    ## iter 100 value 769.992262
    ## final  value 769.992262 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1297.620533 
    ## iter  10 value 654.380840
    ## iter  20 value 497.750561
    ## iter  30 value 445.960609
    ## iter  40 value 409.572486
    ## iter  50 value 372.808264
    ## iter  60 value 349.189847
    ## iter  70 value 339.155126
    ## iter  80 value 332.154142
    ## iter  90 value 328.758746
    ## iter 100 value 325.864778
    ## final  value 325.864778 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1406.127513 
    ## iter  10 value 1168.885570
    ## iter  20 value 895.562464
    ## iter  30 value 577.336365
    ## iter  40 value 354.132516
    ## iter  50 value 235.787271
    ## iter  60 value 174.900990
    ## iter  70 value 160.554713
    ## iter  80 value 157.738578
    ## iter  90 value 156.149832
    ## iter 100 value 155.766555
    ## final  value 155.766555 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1221.537697 
    ## iter  10 value 890.353091
    ## iter  20 value 781.574391
    ## iter  30 value 701.396524
    ## iter  40 value 636.952205
    ## iter  50 value 585.372667
    ## iter  60 value 559.543940
    ## iter  70 value 546.744853
    ## iter  80 value 543.328310
    ## iter  90 value 541.562582
    ## iter 100 value 540.080152
    ## final  value 540.080152 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1405.504846 
    ## iter  10 value 842.326625
    ## iter  20 value 509.176969
    ## iter  30 value 375.301476
    ## iter  40 value 246.220945
    ## iter  50 value 177.219508
    ## iter  60 value 140.632077
    ## iter  70 value 119.274352
    ## iter  80 value 104.109493
    ## iter  90 value 95.660969
    ## iter 100 value 93.710686
    ## final  value 93.710686 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1285.312227 
    ## iter  10 value 379.836422
    ## iter  20 value 140.640899
    ## iter  30 value 59.592132
    ## iter  40 value 20.654551
    ## iter  50 value 6.464412
    ## iter  60 value 3.940758
    ## iter  70 value 3.634113
    ## iter  80 value 3.429314
    ## iter  90 value 3.259352
    ## iter 100 value 3.145067
    ## final  value 3.145067 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1263.024001 
    ## iter  10 value 848.781115
    ## iter  20 value 777.903412
    ## iter  30 value 741.292601
    ## iter  40 value 706.595146
    ## iter  50 value 683.348765
    ## iter  60 value 665.556558
    ## iter  70 value 658.026295
    ## iter  80 value 655.434036
    ## iter  90 value 654.421237
    ## iter 100 value 653.590584
    ## final  value 653.590584 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1206.644978 
    ## iter  10 value 556.495568
    ## iter  20 value 339.902644
    ## iter  30 value 232.619487
    ## iter  40 value 181.810835
    ## iter  50 value 139.058095
    ## iter  60 value 108.986708
    ## iter  70 value 90.611221
    ## iter  80 value 73.518692
    ## iter  90 value 65.645242
    ## iter 100 value 64.681651
    ## final  value 64.681651 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1265.083489 
    ## iter  10 value 476.432340
    ## iter  20 value 166.792639
    ## iter  30 value 68.728610
    ## iter  40 value 22.965304
    ## iter  50 value 10.368425
    ## iter  60 value 4.023436
    ## iter  70 value 0.366800
    ## iter  80 value 0.023099
    ## iter  90 value 0.005988
    ## iter 100 value 0.001569
    ## final  value 0.001569 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1296.917065 
    ## iter  10 value 1127.684411
    ## iter  20 value 959.313034
    ## iter  30 value 833.854111
    ## iter  40 value 791.647107
    ## iter  50 value 765.412519
    ## iter  60 value 759.947819
    ## iter  70 value 756.248323
    ## iter  80 value 755.208861
    ## iter  90 value 750.536805
    ## iter 100 value 749.008486
    ## final  value 749.008486 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1351.309999 
    ## iter  10 value 911.414853
    ## iter  20 value 701.417846
    ## iter  30 value 607.204555
    ## iter  40 value 560.007444
    ## iter  50 value 478.999486
    ## iter  60 value 409.155918
    ## iter  70 value 352.871986
    ## iter  80 value 335.086329
    ## iter  90 value 325.181333
    ## iter 100 value 307.069384
    ## final  value 307.069384 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1220.153180 
    ## iter  10 value 429.950384
    ## iter  20 value 280.256721
    ## iter  30 value 225.226226
    ## iter  40 value 189.559892
    ## iter  50 value 164.329817
    ## iter  60 value 153.978642
    ## iter  70 value 150.812722
    ## iter  80 value 149.035681
    ## iter  90 value 148.249062
    ## iter 100 value 147.600362
    ## final  value 147.600362 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1267.002200 
    ## iter  10 value 894.599380
    ## iter  20 value 764.449793
    ## iter  30 value 703.291375
    ## iter  40 value 658.999115
    ## iter  50 value 645.285005
    ## iter  60 value 642.612190
    ## iter  70 value 638.915410
    ## iter  80 value 637.671024
    ## iter  90 value 637.101865
    ## iter 100 value 636.719986
    ## final  value 636.719986 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1320.207987 
    ## iter  10 value 823.139726
    ## iter  20 value 770.507175
    ## iter  30 value 700.567392
    ## iter  40 value 661.864911
    ## iter  50 value 613.460796
    ## iter  60 value 569.551058
    ## iter  70 value 518.037981
    ## iter  80 value 482.494739
    ## iter  90 value 472.733812
    ## iter 100 value 469.167447
    ## final  value 469.167447 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1260.195291 
    ## iter  10 value 421.021636
    ## iter  20 value 79.808324
    ## iter  30 value 19.795853
    ## iter  40 value 8.758916
    ## iter  50 value 8.030156
    ## iter  60 value 7.652602
    ## iter  70 value 7.317311
    ## iter  80 value 6.686313
    ## iter  90 value 6.230610
    ## iter 100 value 4.970171
    ## final  value 4.970171 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1316.599775 
    ## iter  10 value 819.461845
    ## iter  20 value 735.592933
    ## iter  30 value 695.749728
    ## iter  40 value 649.878825
    ## iter  50 value 613.986998
    ## iter  60 value 585.320185
    ## iter  70 value 568.202020
    ## iter  80 value 563.502135
    ## iter  90 value 559.108852
    ## iter 100 value 555.384858
    ## final  value 555.384858 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1247.940578 
    ## iter  10 value 606.008081
    ## iter  20 value 396.359441
    ## iter  30 value 319.543979
    ## iter  40 value 272.220855
    ## iter  50 value 242.917271
    ## iter  60 value 222.576716
    ## iter  70 value 203.875196
    ## iter  80 value 194.866245
    ## iter  90 value 186.773616
    ## iter 100 value 184.915734
    ## final  value 184.915734 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1186.808728 
    ## iter  10 value 504.540090
    ## iter  20 value 153.513817
    ## iter  30 value 101.439504
    ## iter  40 value 74.511469
    ## iter  50 value 61.049458
    ## iter  60 value 55.298303
    ## iter  70 value 49.417120
    ## iter  80 value 41.813753
    ## iter  90 value 40.979068
    ## iter 100 value 40.246024
    ## final  value 40.246024 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1279.957458 
    ## iter  10 value 981.310214
    ## iter  20 value 892.301019
    ## iter  30 value 846.739178
    ## iter  40 value 815.736191
    ## iter  50 value 800.311785
    ## iter  60 value 791.885361
    ## iter  70 value 786.293860
    ## iter  80 value 785.482629
    ## iter  90 value 785.455445
    ## final  value 785.454985 
    ## converged
    ## # weights:  255
    ## initial  value 1467.151115 
    ## iter  10 value 780.320521
    ## iter  20 value 589.185615
    ## iter  30 value 467.657377
    ## iter  40 value 422.116894
    ## iter  50 value 402.118934
    ## iter  60 value 381.770966
    ## iter  70 value 361.766849
    ## iter  80 value 353.790240
    ## iter  90 value 347.657966
    ## iter 100 value 342.910584
    ## final  value 342.910584 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1193.723809 
    ## iter  10 value 483.562716
    ## iter  20 value 273.018662
    ## iter  30 value 193.280398
    ## iter  40 value 171.019692
    ## iter  50 value 164.791510
    ## iter  60 value 160.862870
    ## iter  70 value 157.861962
    ## iter  80 value 156.601721
    ## iter  90 value 156.189426
    ## iter 100 value 156.021864
    ## final  value 156.021864 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1317.193925 
    ## iter  10 value 852.781371
    ## iter  20 value 775.814766
    ## iter  30 value 730.274016
    ## iter  40 value 713.086685
    ## iter  50 value 700.820687
    ## iter  60 value 688.778502
    ## iter  70 value 680.746186
    ## iter  80 value 673.691004
    ## iter  90 value 663.653761
    ## iter 100 value 648.958908
    ## final  value 648.958908 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1378.794528 
    ## iter  10 value 750.320427
    ## iter  20 value 507.411152
    ## iter  30 value 368.451711
    ## iter  40 value 296.856155
    ## iter  50 value 247.158247
    ## iter  60 value 202.206311
    ## iter  70 value 155.639498
    ## iter  80 value 129.780028
    ## iter  90 value 106.609648
    ## iter 100 value 91.237863
    ## final  value 91.237863 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1182.997566 
    ## iter  10 value 385.732474
    ## iter  20 value 142.548756
    ## iter  30 value 83.253257
    ## iter  40 value 48.359929
    ## iter  50 value 34.642443
    ## iter  60 value 29.082779
    ## iter  70 value 27.907523
    ## iter  80 value 24.731736
    ## iter  90 value 13.796239
    ## iter 100 value 7.076152
    ## final  value 7.076152 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1327.727247 
    ## iter  10 value 1014.356290
    ## iter  20 value 870.896728
    ## iter  30 value 771.865671
    ## iter  40 value 717.240656
    ## iter  50 value 695.620703
    ## iter  60 value 684.402560
    ## iter  70 value 675.644228
    ## iter  80 value 665.650564
    ## iter  90 value 663.579666
    ## iter 100 value 662.250392
    ## final  value 662.250392 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1261.732542 
    ## iter  10 value 685.934633
    ## iter  20 value 343.845530
    ## iter  30 value 231.974811
    ## iter  40 value 204.633706
    ## iter  50 value 181.940280
    ## iter  60 value 172.644829
    ## iter  70 value 160.638758
    ## iter  80 value 153.737644
    ## iter  90 value 145.991955
    ## iter 100 value 141.171952
    ## final  value 141.171952 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1402.791623 
    ## iter  10 value 400.354148
    ## iter  20 value 130.499413
    ## iter  30 value 49.566930
    ## iter  40 value 28.731971
    ## iter  50 value 22.609040
    ## iter  60 value 20.978995
    ## iter  70 value 19.646527
    ## iter  80 value 17.825819
    ## iter  90 value 15.285706
    ## iter 100 value 14.660794
    ## final  value 14.660794 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1231.655936 
    ## iter  10 value 992.332947
    ## iter  20 value 829.862867
    ## iter  30 value 770.574401
    ## iter  40 value 748.025692
    ## iter  50 value 735.172288
    ## iter  60 value 728.682431
    ## iter  70 value 726.303591
    ## iter  80 value 725.155063
    ## iter  90 value 725.102885
    ## final  value 725.102673 
    ## converged
    ## # weights:  255
    ## initial  value 1392.600213 
    ## iter  10 value 694.175182
    ## iter  20 value 494.373242
    ## iter  30 value 431.544559
    ## iter  40 value 376.898160
    ## iter  50 value 347.366444
    ## iter  60 value 324.837964
    ## iter  70 value 313.891894
    ## iter  80 value 305.841293
    ## iter  90 value 303.386619
    ## iter 100 value 301.202303
    ## final  value 301.202303 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1360.108844 
    ## iter  10 value 637.413662
    ## iter  20 value 323.051708
    ## iter  30 value 239.407868
    ## iter  40 value 206.581153
    ## iter  50 value 176.014202
    ## iter  60 value 160.588697
    ## iter  70 value 158.471539
    ## iter  80 value 158.175014
    ## iter  90 value 158.025443
    ## iter 100 value 157.918880
    ## final  value 157.918880 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1343.663560 
    ## iter  10 value 1016.081079
    ## iter  20 value 792.159487
    ## iter  30 value 727.227147
    ## iter  40 value 689.620175
    ## iter  50 value 658.526738
    ## iter  60 value 638.772852
    ## iter  70 value 613.017198
    ## iter  80 value 608.564071
    ## iter  90 value 602.830941
    ## iter 100 value 594.099047
    ## final  value 594.099047 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1282.171517 
    ## iter  10 value 543.710410
    ## iter  20 value 384.095730
    ## iter  30 value 328.388169
    ## iter  40 value 291.690279
    ## iter  50 value 263.895139
    ## iter  60 value 235.885370
    ## iter  70 value 215.609793
    ## iter  80 value 199.166160
    ## iter  90 value 189.065858
    ## iter 100 value 182.411136
    ## final  value 182.411136 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1414.251869 
    ## iter  10 value 425.021172
    ## iter  20 value 139.177302
    ## iter  30 value 78.494247
    ## iter  40 value 54.176926
    ## iter  50 value 42.042800
    ## iter  60 value 40.735922
    ## iter  70 value 38.421287
    ## iter  80 value 36.143393
    ## iter  90 value 35.063177
    ## iter 100 value 33.674432
    ## final  value 33.674432 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1276.840347 
    ## iter  10 value 931.774623
    ## iter  20 value 781.972998
    ## iter  30 value 759.183050
    ## iter  40 value 700.393569
    ## iter  50 value 656.179694
    ## iter  60 value 630.940464
    ## iter  70 value 622.815239
    ## iter  80 value 620.334580
    ## iter  90 value 619.493574
    ## iter 100 value 619.280932
    ## final  value 619.280932 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1202.549853 
    ## iter  10 value 620.268903
    ## iter  20 value 422.742803
    ## iter  30 value 362.389746
    ## iter  40 value 312.788820
    ## iter  50 value 252.482148
    ## iter  60 value 200.739132
    ## iter  70 value 158.547183
    ## iter  80 value 138.103672
    ## iter  90 value 131.162979
    ## iter 100 value 125.961937
    ## final  value 125.961937 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1430.737605 
    ## iter  10 value 467.176936
    ## iter  20 value 168.090901
    ## iter  30 value 96.811285
    ## iter  40 value 61.832480
    ## iter  50 value 44.337457
    ## iter  60 value 33.296295
    ## iter  70 value 26.746254
    ## iter  80 value 24.630960
    ## iter  90 value 24.338340
    ## iter 100 value 24.085095
    ## final  value 24.085095 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1383.187738 
    ## iter  10 value 991.338570
    ## iter  20 value 928.396404
    ## iter  30 value 844.905914
    ## iter  40 value 784.653986
    ## iter  50 value 757.261157
    ## iter  60 value 748.504713
    ## iter  70 value 744.342908
    ## iter  80 value 742.996161
    ## iter  90 value 742.601631
    ## iter 100 value 742.552091
    ## final  value 742.552091 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1202.517902 
    ## iter  10 value 661.425399
    ## iter  20 value 453.298833
    ## iter  30 value 396.525296
    ## iter  40 value 355.447748
    ## iter  50 value 328.777446
    ## iter  60 value 316.673368
    ## iter  70 value 307.169044
    ## iter  80 value 302.026017
    ## iter  90 value 299.374699
    ## iter 100 value 296.489993
    ## final  value 296.489993 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1362.495686 
    ## iter  10 value 511.229113
    ## iter  20 value 303.865674
    ## iter  30 value 217.840981
    ## iter  40 value 182.426297
    ## iter  50 value 171.306470
    ## iter  60 value 163.648459
    ## iter  70 value 160.842736
    ## iter  80 value 159.668784
    ## iter  90 value 159.384811
    ## iter 100 value 159.040863
    ## final  value 159.040863 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1257.887225 
    ## iter  10 value 1097.096874
    ## iter  20 value 1094.194972
    ## iter  30 value 1094.156621
    ## final  value 1094.156569 
    ## converged
    ## # weights:  255
    ## initial  value 1241.296907 
    ## iter  10 value 650.071176
    ## iter  20 value 456.881294
    ## iter  30 value 399.512948
    ## iter  40 value 358.094666
    ## iter  50 value 333.815287
    ## iter  60 value 316.318720
    ## iter  70 value 309.093548
    ## iter  80 value 300.859832
    ## iter  90 value 297.227478
    ## iter 100 value 296.000483
    ## final  value 296.000483 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1333.765931 
    ## iter  10 value 1028.400297
    ## iter  20 value 671.092780
    ## iter  30 value 513.393260
    ## iter  40 value 429.025374
    ## iter  50 value 288.902886
    ## iter  60 value 241.077166
    ## iter  70 value 217.959957
    ## iter  80 value 201.994543
    ## iter  90 value 182.917939
    ## iter 100 value 145.287308
    ## final  value 145.287308 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1254.413128 
    ## iter  10 value 953.509629
    ## iter  20 value 831.251584
    ## iter  30 value 726.286559
    ## iter  40 value 626.039534
    ## iter  50 value 563.809501
    ## iter  60 value 505.500346
    ## iter  70 value 476.415893
    ## iter  80 value 445.566497
    ## iter  90 value 424.239150
    ## iter 100 value 409.883600
    ## final  value 409.883600 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1357.636338 
    ## iter  10 value 1108.343878
    ## iter  20 value 814.315704
    ## iter  30 value 717.462333
    ## iter  40 value 664.093147
    ## iter  50 value 636.469379
    ## iter  60 value 611.400371
    ## iter  70 value 600.075574
    ## iter  80 value 588.438144
    ## iter  90 value 585.057387
    ## iter 100 value 583.413709
    ## final  value 583.413709 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1290.989620 
    ## iter  10 value 354.944412
    ## iter  20 value 169.007445
    ## iter  30 value 78.822968
    ## iter  40 value 52.138035
    ## iter  50 value 35.569282
    ## iter  60 value 29.822197
    ## iter  70 value 21.069997
    ## iter  80 value 14.536164
    ## iter  90 value 14.126229
    ## iter 100 value 14.044878
    ## final  value 14.044878 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1270.230853 
    ## iter  10 value 1114.210537
    ## iter  20 value 871.213580
    ## iter  30 value 817.430357
    ## iter  40 value 792.189436
    ## iter  50 value 762.994578
    ## iter  60 value 747.462499
    ## iter  70 value 742.091861
    ## iter  80 value 741.763646
    ## iter  90 value 741.749060
    ## iter 100 value 741.744269
    ## final  value 741.744269 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1395.743369 
    ## iter  10 value 863.588783
    ## iter  20 value 681.021532
    ## iter  30 value 593.391709
    ## iter  40 value 545.139521
    ## iter  50 value 478.190752
    ## iter  60 value 432.446114
    ## iter  70 value 378.442562
    ## iter  80 value 334.348805
    ## iter  90 value 319.048681
    ## iter 100 value 307.837283
    ## final  value 307.837283 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1278.887268 
    ## iter  10 value 647.424325
    ## iter  20 value 401.816597
    ## iter  30 value 296.487151
    ## iter  40 value 220.553457
    ## iter  50 value 174.613951
    ## iter  60 value 160.613159
    ## iter  70 value 152.793712
    ## iter  80 value 149.110016
    ## iter  90 value 147.509508
    ## iter 100 value 147.018869
    ## final  value 147.018869 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1293.598629 
    ## iter  10 value 911.441116
    ## iter  20 value 746.135622
    ## iter  30 value 678.981195
    ## iter  40 value 635.537554
    ## iter  50 value 607.153111
    ## iter  60 value 571.270919
    ## iter  70 value 555.581293
    ## iter  80 value 554.104833
    ## iter  90 value 548.398284
    ## iter 100 value 545.543221
    ## final  value 545.543221 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1302.301612 
    ## iter  10 value 1088.848214
    ## iter  20 value 669.204696
    ## iter  30 value 490.384080
    ## iter  40 value 380.275505
    ## iter  50 value 322.611263
    ## iter  60 value 295.345243
    ## iter  70 value 270.384674
    ## iter  80 value 259.286533
    ## iter  90 value 252.103231
    ## iter 100 value 249.813363
    ## final  value 249.813363 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1295.576759 
    ## iter  10 value 457.654191
    ## iter  20 value 146.148105
    ## iter  30 value 59.812373
    ## iter  40 value 27.774180
    ## iter  50 value 16.060666
    ## iter  60 value 11.005598
    ## iter  70 value 10.593579
    ## iter  80 value 8.650625
    ## iter  90 value 6.351768
    ## iter 100 value 5.311221
    ## final  value 5.311221 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1236.897168 
    ## iter  10 value 935.552774
    ## iter  20 value 818.749032
    ## iter  30 value 770.518180
    ## iter  40 value 736.093155
    ## iter  50 value 706.204247
    ## iter  60 value 689.347532
    ## iter  70 value 681.654393
    ## iter  80 value 679.397034
    ## iter  90 value 677.954910
    ## iter 100 value 677.472088
    ## final  value 677.472088 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1385.512275 
    ## iter  10 value 1076.405508
    ## iter  20 value 731.447890
    ## iter  30 value 515.283932
    ## iter  40 value 428.233127
    ## iter  50 value 390.999756
    ## iter  60 value 353.968619
    ## iter  70 value 331.089408
    ## iter  80 value 321.732343
    ## iter  90 value 316.141209
    ## iter 100 value 313.328052
    ## final  value 313.328052 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1284.177312 
    ## iter  10 value 470.451639
    ## iter  20 value 114.549000
    ## iter  30 value 48.190132
    ## iter  40 value 29.466961
    ## iter  50 value 21.754278
    ## iter  60 value 15.260006
    ## iter  70 value 12.823302
    ## iter  80 value 11.466274
    ## iter  90 value 10.752337
    ## iter 100 value 10.327547
    ## final  value 10.327547 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1257.184558 
    ## iter  10 value 987.745529
    ## iter  20 value 926.471678
    ## iter  30 value 872.304914
    ## iter  40 value 834.421952
    ## iter  50 value 820.071703
    ## iter  60 value 810.505334
    ## iter  70 value 809.298079
    ## iter  80 value 809.108007
    ## iter  90 value 808.428509
    ## iter 100 value 807.952478
    ## final  value 807.952478 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1289.549961 
    ## iter  10 value 788.588756
    ## iter  20 value 528.710904
    ## iter  30 value 447.866729
    ## iter  40 value 400.074313
    ## iter  50 value 373.482417
    ## iter  60 value 355.401144
    ## iter  70 value 344.056253
    ## iter  80 value 332.301679
    ## iter  90 value 328.085290
    ## iter 100 value 324.153792
    ## final  value 324.153792 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1372.169438 
    ## iter  10 value 475.315858
    ## iter  20 value 255.214150
    ## iter  30 value 190.871258
    ## iter  40 value 161.110216
    ## iter  50 value 149.493848
    ## iter  60 value 147.300327
    ## iter  70 value 145.761224
    ## iter  80 value 143.560295
    ## iter  90 value 142.385167
    ## iter 100 value 140.831817
    ## final  value 140.831817 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1285.781937 
    ## iter  10 value 970.651584
    ## iter  20 value 933.506925
    ## iter  30 value 897.446777
    ## iter  40 value 846.093753
    ## iter  50 value 821.344758
    ## iter  60 value 805.257602
    ## iter  70 value 789.647759
    ## iter  80 value 785.182768
    ## iter  90 value 783.569551
    ## iter 100 value 782.671486
    ## final  value 782.671486 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1308.485905 
    ## iter  10 value 810.288537
    ## iter  20 value 651.047124
    ## iter  30 value 512.813750
    ## iter  40 value 430.807411
    ## iter  50 value 372.995617
    ## iter  60 value 337.218745
    ## iter  70 value 304.463436
    ## iter  80 value 277.956264
    ## iter  90 value 253.360971
    ## iter 100 value 245.013711
    ## final  value 245.013711 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1296.719252 
    ## iter  10 value 438.696103
    ## iter  20 value 103.818532
    ## iter  30 value 44.657730
    ## iter  40 value 20.524794
    ## iter  50 value 10.385458
    ## iter  60 value 5.815675
    ## iter  70 value 4.153344
    ## iter  80 value 2.554824
    ## iter  90 value 2.028066
    ## iter 100 value 1.902003
    ## final  value 1.902003 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1285.617879 
    ## iter  10 value 884.425831
    ## iter  20 value 809.005854
    ## iter  30 value 761.555588
    ## iter  40 value 714.048936
    ## iter  50 value 673.933606
    ## iter  60 value 640.297061
    ## iter  70 value 629.726973
    ## iter  80 value 626.678277
    ## iter  90 value 624.809271
    ## iter 100 value 623.553290
    ## final  value 623.553290 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1296.187072 
    ## iter  10 value 745.823926
    ## iter  20 value 615.672664
    ## iter  30 value 534.841940
    ## iter  40 value 457.234248
    ## iter  50 value 412.404395
    ## iter  60 value 352.391725
    ## iter  70 value 284.680348
    ## iter  80 value 245.446714
    ## iter  90 value 217.505703
    ## iter 100 value 197.685211
    ## final  value 197.685211 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1256.019752 
    ## iter  10 value 292.591056
    ## iter  20 value 108.900887
    ## iter  30 value 60.197469
    ## iter  40 value 37.409577
    ## iter  50 value 27.812882
    ## iter  60 value 24.887887
    ## iter  70 value 21.487937
    ## iter  80 value 21.171199
    ## iter  90 value 21.128703
    ## iter 100 value 21.121730
    ## final  value 21.121730 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1258.597184 
    ## iter  10 value 1085.073120
    ## iter  20 value 912.483832
    ## iter  30 value 873.726781
    ## iter  40 value 848.770187
    ## iter  50 value 818.052228
    ## iter  60 value 803.536147
    ## iter  70 value 797.209182
    ## iter  80 value 794.053946
    ## iter  90 value 793.847144
    ## iter 100 value 793.754821
    ## final  value 793.754821 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1295.466079 
    ## iter  10 value 817.340986
    ## iter  20 value 533.292744
    ## iter  30 value 421.515948
    ## iter  40 value 381.136759
    ## iter  50 value 350.886516
    ## iter  60 value 331.555679
    ## iter  70 value 323.611798
    ## iter  80 value 320.945895
    ## iter  90 value 320.262036
    ## iter 100 value 320.029297
    ## final  value 320.029297 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1291.606938 
    ## iter  10 value 456.681476
    ## iter  20 value 301.880729
    ## iter  30 value 227.409988
    ## iter  40 value 188.685237
    ## iter  50 value 170.490720
    ## iter  60 value 162.128460
    ## iter  70 value 159.055530
    ## iter  80 value 156.393908
    ## iter  90 value 154.567757
    ## iter 100 value 153.379700
    ## final  value 153.379700 
    ## stopped after 100 iterations
    ## # weights:  95
    ## initial  value 1354.473248 
    ## iter  10 value 1140.998178
    ## iter  20 value 1138.162362
    ## iter  30 value 904.801556
    ## iter  40 value 836.751015
    ## iter  50 value 825.962092
    ## iter  60 value 818.235429
    ## iter  70 value 803.586570
    ## iter  80 value 798.169977
    ## iter  90 value 774.567732
    ## iter 100 value 764.267819
    ## final  value 764.267819 
    ## stopped after 100 iterations
    ## # weights:  255
    ## initial  value 1263.824124 
    ## iter  10 value 610.619240
    ## iter  20 value 303.001889
    ## iter  30 value 245.357559
    ## iter  40 value 198.530863
    ## iter  50 value 163.597472
    ## iter  60 value 138.307082
    ## iter  70 value 120.915596
    ## iter  80 value 115.535512
    ## iter  90 value 114.620113
    ## iter 100 value 111.082187
    ## final  value 111.082187 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1328.725397 
    ## iter  10 value 431.322516
    ## iter  20 value 120.374001
    ## iter  30 value 41.171273
    ## iter  40 value 18.233307
    ## iter  50 value 9.375753
    ## iter  60 value 4.812299
    ## iter  70 value 3.553036
    ## iter  80 value 3.034787
    ## iter  90 value 2.748945
    ## iter 100 value 2.472090
    ## final  value 2.472090 
    ## stopped after 100 iterations
    ## # weights:  415
    ## initial  value 1267.341143 
    ## iter  10 value 629.141734
    ## iter  20 value 339.154391
    ## iter  30 value 242.935124
    ## iter  40 value 192.772628
    ## iter  50 value 180.058032
    ## iter  60 value 175.761891
    ## iter  70 value 174.167145
    ## iter  80 value 172.776000
    ## iter  90 value 168.555613
    ## iter 100 value 163.534520
    ## final  value 163.534520 
    ## stopped after 100 iterations

``` r
# train the GBM model
set.seed(7)
modelGbm <- train(Class~., data=soyTrain, method="gbm", metric="accuracy", trControl=cv_control, verbose=FALSE)
```

    ## Warning in train.default(x, y, weights = w, ...): The metric "accuracy" was
    ## not in the result set. Accuracy will be used instead.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 64: roots2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 47: ext.decay2 has no variation.

    ## Warning in (function (x, y, offset = NULL, misc = NULL, distribution =
    ## "bernoulli", : variable 53: fruit.pods2 has no variation.

``` r
# train the SVM model
set.seed(7)
modelSvm <- train(Class~., data=soyTrain, method="svmRadial", metric="accuracy", trControl=cv_control, verbose=FALSE)
```

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in train.default(x, y, weights = w, ...): The metric "accuracy" was
    ## not in the result set. Accuracy will be used instead.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

    ## Warning in .local(x, ...): Variable(s) `' constant. Cannot scale data.

``` r
# collect resamples
results <- resamples(list(NN=modelNN, GBM=modelGbm, SVM=modelSvm))
# summarize the distributions
summary(results)
```

    ## 
    ## Call:
    ## summary.resamples(object = results)
    ## 
    ## Models: NN, GBM, SVM 
    ## Number of resamples: 100 
    ## 
    ## Accuracy 
    ##          Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
    ## NN  0.8152866 0.8816165 0.8962399 0.8949980 0.9099075 0.9526627    0
    ## GBM 0.8470588 0.8916814 0.9041882 0.9044028 0.9175856 0.9629630    0
    ## SVM 0.8713450 0.9041504 0.9147034 0.9134670 0.9261518 0.9540230    0
    ## 
    ## Kappa 
    ##          Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
    ## NN  0.7952327 0.8663320 0.8837667 0.8823115 0.8982047 0.9466456    0
    ## GBM 0.8259157 0.8786572 0.8925447 0.8927855 0.9079758 0.9583226    0
    ## SVM 0.8562184 0.8930030 0.9038380 0.9029037 0.9172614 0.9479626    0

``` r
# boxplots of results
bwplot(results)
```

![](sm10_supervisedLearning_files/figure-markdown_github/unnamed-chunk-19-1.png)

> You can finetune the parameters for different models in Caret using tuneGrid.

### Additional metrics for evaluation

We can also consider other metrics to assess the performance of a model. This becomes particularly important when, for example, your your test set is imbalanced. In that scenario, evaluating the accuracy of a model might not be the best indication that your classifier works well. In fact, it may be biased by the over-represented class in your test set.

#### Kappa-score

Overall accuracy is a misleading statistic in case of unbalanced datasets. The kappa statistic overcomes thsi by taking the expected error into account.

#### ROC, Precision & Recall

Receiver-Operator Characteristic (ROC) Curves can be used to characterize the performance of our model in a binary classification setting.
For binary classification, we might also be interested in precision and recall, i.e. a metric of our True Positive and True Negative rates.

``` r
binary_preds = as.integer(y_pred_cv)-1
binary_true = as.integer(soyTest$Class)-1
precision <- sum(binary_preds & binary_true) / sum(binary_preds)
recall <- sum(binary_preds & binary_true) / sum(binary_true)
```

Precision tells us how often, when we predict a class 'y', do we get it correct.
Recall tells us how often, out of all our 'y' instances, do we predict them correctly.

``` r
#Precision    
print(precision)
```

    ## [1] 0.1560284

``` r
#Recall  
print(recall)
```

    ## [1] 0.1405751

In case of multi-class classification, another metric that comes in handy is the F1-score. This is the harmonic mean of precision and recall. A high F1-score will tell you that you are quite precise and sensitive in your prediction of *all* classes.

``` r
Fmeasure <- 2 * precision * recall / (precision + recall)
```

> In the models we fit for our multiclass problem, what was the Kappa statistic on our test set? What was the Accuracy?
> Which model selection appraoch worked better, the bootstrapped SVM or the repeated 5-CV approach? How do we assess this (on the test set, or from the cross-validation results?)

Discussion points
-----------------

How good do you think a classifier will have to be to be clinically relevant? What level of specificity or sensitivity do you think is "good enough"? The author's data set has half lymph-node positive and half negative. The incidence of lymph-node-positive breast cancer is about 33% at the time of diagnosis (according to \[<http://seer.cancer.gov/statfacts/html/breast.html>\]). How does that affect your opinion about the classifier?

### Take-home exercises (Restated from above)

Recall our 6-fold cross validation on the binary dataset.

**Exercise 1**: perform 100 runs of this CV before selecting a model to test! Add at least on model to the list of models, e.g., use genes with a p-val threshold &lt; cutoff.

**Exercise 2**: Use AUC as a criteria to select a model based on the training data! Tip: extract the predicted probabilities from each method and use the roc function in ROCR.

References
==========

1: Varma S, Simon R: Bias in error estimation when using cross-validation for model selection. BMC Bioinformatics 2006, 7:91.

2: Statnikov A, Aliferis CF, Tsamardinos I, Hardin D, Levy S: A comprehensive evaluation of multicategory classification methods for microarray gene expression cancer diagnosis. Bioinformatics 2005, 21:631-643.
