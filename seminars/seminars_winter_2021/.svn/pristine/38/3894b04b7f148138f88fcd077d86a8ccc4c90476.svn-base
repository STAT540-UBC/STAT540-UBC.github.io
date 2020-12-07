
# lets assess differential expression for all genes across all dev_stage in the wildtype condition only

# filter for wildtype data
wildTypeSamples <- samplesMetadata %>% filter(genotype == "wt")

wildTypeSamples # confirm only wildtype data exist

# expressionMatrix %>% transformGeneExpressionMatrix()

# reusable function for pulling expression data for given samples
getExpressionForSamples <- function(sampleIds, expressionMatrix) {
  # use gene column as row name
  dataFrame <- expressionMatrix %>% 
    as.data.frame() %>% 
    column_to_rownames("gene")
  return(dataFrame[sampleIds])
}

# use the wildTypeSamples to pull out the wildtype expression data from the expression matrix
wildTypeExpressionMatrix <- getExpressionForSamples(wildTypeSamples$sample_id, expressionMatrix)

head(wildTypeExpressionMatrix) # now we have a data frame containing all expression data for the wildtype samples

# fit the linear model

# construct the design matrix for fitting the linear model
designMatrix <- model.matrix(~dev_stage, wildTypeSamples)

# fit using limma

limmaFit <- lmFit(wildTypeExpressionMatrix, designMatrix)
limmaFitEb <- eBayes(limmaFit)

topTable(limmaFitEb)


# this gives us anything that is different from 
# topTable(linearFitEb,
#          number = Inf, 
#          adjust.method = "fdr", 
#          p.value = 0.05)





########################
### CHALLENGE #######
########################

## Find what's changed from stage p6 to p10

contrastMatrix <- makeContrasts(
  p10vsp6 = dev_stageP10 - dev_stageP6,
  fourweeksVsP10 = dev_stage4_weeks - dev_stageP10,
  levels = designMatrix
)

contrastFit <- contrasts.fit(limmaFit, contrastMatrix)
contrastFitEb <- eBayes(contrastFit)

topTable(contrastFitEb)


cutoff <- 1e-04
contrastCounts <- decideTests(contrastFitEb, p.value = cutoff, method = "global")
summary(contrastCounts)


topTable(contrastFitEb, number = Inf) %>% 
  as_tibble()



xxx <- summary(contrastCounts)





########################
### CHALLENGE #######
########################

## Find whatever changes up to P6 and then hold steady until 4_weeks

bonusContrast <- makeContrasts(
  p2 = dev_stageP2,
  p6vsp2 = dev_stageP6 - dev_stageP2,
  p10vsp6 = dev_stageP10 - dev_stageP6,
  fourweeksVsP10 = dev_stage4_weeks - dev_stageP10,
  levels = designMatrix
)


contrastFit <- contrasts.fit(limmaFit, bonusContrast)
contrastFitEb <- eBayes(contrastFit)
contrastFitEb %>% topTable(number = Inf) %>% as.data.frame() %>% rownames_to_column("gene") %>% View()

contrastCounts <- decideTests(contrastFitEb, p.value = 0.05, method = "global")
summary(contrastCounts)

contrastChanges <- contrastCounts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene")

hits <- contrastChanges %>% filter(p2 != 0, p6vsp2 != 0, p10vsp6 == 0, fourweeksVsP10 == 0)
hits$gene

contrastHits <- contrastFitEb %>% topTable(number = Inf) %>% 
  as_tibble() %>% 
  rownames_to_column("gene")

contrastHits %>% 
  filter(gene %in% hits$gene) %>% 
  arrange(adj.P.Val)




expressionData <- expressionMatrix %>% filter(gene == "1451555_at") %>% transformGeneExpressionMatrix()

expressionData %>% 
  right_join(samplesMetadata %>% filter(genotype == "wt")) %>% 
  ggplot(aes(x = dev_stage, y = expression)) +
  geom_point() +
  geom_smooth()


