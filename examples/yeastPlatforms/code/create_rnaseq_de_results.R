runEdgerAnalysis <- function(counts, groups){
    # Use the `edgeR` package to test for differentially expressed genes between two groups.
    #
    # Args:
    #   counts : (data.frame) Data frame of counts with genes as rows and samples a columns.
    #   groups : (factor) Grouping of samples.
    #
    # Results:
    #   results : (data.frame)
    
    library('edgeR')
    
    # Build data object
    dat <- DGEList(counts=counts, group=groups)
    
    # Build design matrix
    des <- model.matrix(~groups)
    
    # Normalise data
    dat <- estimateGLMCommonDisp(dat, des)
    
    dat <- estimateGLMTrendedDisp(dat, des)
    
    dat <- estimateGLMTagwiseDisp(dat, des)
    
    # GLM fitting
    fit <- glmFit(dat, des)
    
    lrt <- glmLRT(fit, coef=2)
    
    # Extract results
    top.tags <- topTags(lrt, adjust.method="BH", n=Inf)
    
    top.table <- top.tags$table
    
    # Store results in standard format data frame.
    results <- data.frame(gene.id = I(row.names(top.table)),
                          p.value = top.table$PValue,
                          q.value = top.table$FDR,
                          log.fc = top.table$logFC,
                          test.stat = top.table$LR)
    return(results)    
}

runLimmaAnalysis <- function(counts, groups){
    # Use the `limma` package to test for differentially expressed genes between two groups.
    #
    # Args:
    #   counts : (data.frame) Data frame of counts with genes as rows and samples a columns.
    #   groups : (factor) Grouping of samples.
    #
    # Results:
    #   results : (data.frame)
    
    library("edgeR")
    
    library("limma")
    
    # Build design matrix
    des <- model.matrix(~groups)
    
    # Compute normalisation
    norm.factors <- calcNormFactors(counts)
    
    lib.size <- norm.factors * colSums(counts)
    
    # Convert counts to log2-cpm
    dat <- voom(counts, des, lib.size=lib.size)
    
    # Fit model
    fit <- lmFit(dat, des)
    
    fit <- eBayes(fit)
    
    # Extract all results
    top.table <- topTable(fit, coef=2, n=Inf, p.value=Inf)
    
    # Put results in a new data frame
    results <- with(top.table, data.frame(gene.id = I(ID),                                     
                                          p.value = P.Value,
                                          q.value = adj.P.Val,
                                          log.fc = logFC,
                                          test.stat = t))
    
    return(results)
}

runDESeqAnalysis <- function(counts, groups){
  # Use the `DESeq` package to test for differentially expressed genes between two groups.
  #
  # Args:
  #   counts : (data.frame) Data frame of counts with genes as rows and samples a columns.
  #   groups : (factor) Grouping of samples.
  #
  # Results:
  #   results : (data.frame)
  
  library(DESeq)
  
  # Build data object
  dat <- newCountDataSet(counts, groups)
  
  # Normalise data
  dat <- estimateSizeFactors(dat)
  dat <- estimateDispersions(dat)
  
  #Negative Binomial fit
  results <- nbinomTest(dat, levels(groups)[1], levels(groups)[2])
  
  #Sort results by ascending p.value and descending log.fc
  results<-results[order(results$pval,-results$log2FoldChange),]
  
  # Store results in standard format data frame
  results <- data.frame(gene.id = I(results$id),
                        p.value = results$pval,
                        q.value = results$padj,
                        log.fc = results$log2FoldChange)#,
                        #test.stat = results$LR) #there is not such output value
  
  return(results)
}


#-------------------------------------------------------------
# DEA
## working directory assumed to be the parent of the code directory where this
## file lies; at the time of writing, that means the 'yeastPlatforms' directory

#setwd('/home/andrew/code/stat540/rmd/homework/homework_02')

counts <- list(stampy=read.table('data/stampy.counts.tsv', header=TRUE, row.names=1),
               tophat=read.table('data/tophat.counts.tsv', header=TRUE, row.names=1))

groups <- factor(c(rep("b", 3), c(rep("c", 3))))

edger.results <- lapply(counts, runEdgerAnalysis, groups=groups)

limma.results <- lapply(counts, runLimmaAnalysis, groups=groups)

deseq.results <- lapply(counts, runDESeqAnalysis, groups=groups)

for(aligner in names(edger.results)){
  write.table(edger.results[[aligner]], paste("results/", aligner, ".edger.results.tsv", sep=""), row.names=TRUE, col.names=NA, sep="\t")
}

for(aligner in names(limma.results)){
  write.table(limma.results[[aligner]], paste("results/", aligner, ".limma.results.tsv", sep=""), row.names=TRUE, col.names=NA, sep="\t")
}

for(aligner in names(deseq.results)){
  write.table(deseq.results[[aligner]], paste("results/", aligner, ".deseq.results.tsv", sep=""), row.names=TRUE, col.names=NA, sep="\t")
}