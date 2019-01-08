cutoff <- 1e-06
changeDirections <- decideTests(interactionFit, p.value = cutoff, method = "global") %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  as_tibble()

# look for down regulation across developmental stages in wt but up regulation in genotypeNrlKO 
hits <- changeDirections %>% filter(dev_stage4_weeks < 0, `genotypeNrlKO:dev_stage4_weeks` > 0)

hits <- changeDirections %>% filter(dev_stage4_weeks > 0)

# -------

library(EGAD)

# compute multifunction scores for genes
annotations <- GO.mouse[c("name", "GO")]
annotations <- make_annotations(annotations, unique(annotations[["name"]]), unique(annotations[["GO"]]))
multiFuncScores <- calculate_multifunc(annotations)
# geneMultifuncScores %>%
#   remove_rownames() %>%
#   write_csv("data/gene_multifunction_scores.csv")

# library(tidyverse)
# library(mouse4302.db)

# explore correlations with dea scores and mulitfunction

multiFuncScores <- read_csv("data/gene_multifunction_scores.csv")

# get the probe mappings
# mappings <- mouse4302ALIAS2PROBE %>% as.list()
# probes <- character(0)
# genes <- character(0)
# for (gene in names(mappings)) {
#   currGeneProbes <- mappings[[gene]]
#   probes <- c(probes, currGeneProbes)
#   genes <- c(genes, rep(gene, length(currGeneProbes)))
# }
# mappingsDf <- tibble(probe = probes, gene = genes)
# mappingsDf %>% write_csv("data/probe_gene_mappings.csv")

probeMappings <- read_csv("data/probe_gene_mappings.csv")


deaResult <- interactionFit %>% topTable(n = Inf) %>% rownames_to_column("gene") %>% as_tibble()

deaGeneResult <- deaResult %>% 
  dplyr::select(probe = gene, everything()) %>% 
  left_join(probeMappings, by = "probe") %>% 
  dplyr::select(probe, gene, everything())

uniqueProbes <- deaGeneResult %>% 
  group_by(probe) %>% 
  summarise(num_genes = n()) %>% 
  filter(num_genes == 1)

result <- deaGeneResult %>% 
  filter(probe %in% uniqueProbes$probe) %>% 
  dplyr::select(-probe) %>% 
  na.omit() %>% 
  group_by(gene) %>% 
  summarize_all(.funs = mean)


result %>% 
  inner_join(dplyr::select(multiFuncScores, gene = Gene, everything()), by = "gene") %>% 
  ggplot(aes(x = -log10(P.Value), y = MF.score)) +
  geom_point()



# ======= try with expressions from gemma

experimentDesign <- read_csv("data/design.csv")
expressionData <- read_csv("data/expression_data.csv") %>% dplyr::select(gene = GeneSymbol, everything())
expressionData <- expressionData[c("gene", experimentDesign$Bioassay)] %>% 
  group_by(gene) %>% 
  summarize_all(.funs = mean) %>% 
  ungroup() %>% 
  na.omit()

multiFuncScores <- read_csv("data/gene_multifunction_scores.csv") %>% 
  dplyr::select(gene = Gene, everything())


# expressionData %>% inner_join(multiFuncScores, by = "gene")




# DEA

targetSamples <- experimentDesign
targetSamples$timepoint <- factor(targetSamples$timepoint, level = c("E16", "P2", "P6", "P10", "4_weeks"))
targetSamples$genotype <- factor(targetSamples$genotype, levels = c("wild_type_genotype", "Nrl_[mouse]_neural_retina_leucine_zipper_gene_|_Homozygous_negative_|"))
targetSamples <- targetSamples %>% arrange(genotype, timepoint)
targetExpressions <- expressionData[c("gene", targetSamples$Bioassay)] %>% as.data.frame() %>% column_to_rownames("gene")

targetDesignMatrix <- model.matrix(~genotype, targetSamples)


# fit linear model
targetFit <- lmFit(targetExpressions, targetDesignMatrix)

# run ebayes to calculate moderated t-statistics
targetFitEbayes <- eBayes(targetFit)

targetDeaResult <- topTable(targetFitEbayes, n = Inf) %>% rownames_to_column("gene") %>% as_tibble()


# assess correlation with multifunctionality

multiFunCor <- targetDeaResult %>% 
  inner_join(multiFuncScores, by = "gene") %>% 
  mutate(diff_rank = rank(P.Value))



multiFunCor %>% 
  ggplot(aes(x = abs(logFC), y = MF.score)) + 
  geom_point(alpha = 0.1) +
  ggtitle(paste0("r = ", cor(abs(multiFunCor$logFC), multiFunCor$MF.score, method = "spearman")))




# plot genes

expressionDataForGenes <- targetExpressions %>% 
  rownames_to_column("gene") %>% 
  filter(gene == "Fam135a") %>%
  transformGeneExpressionMatrix() %>% 
  left_join(targetSamples %>% dplyr::select(sample_id = Bioassay, everything()), id = "sample_id", by = "sample_id")

expressionDataForGenes %>% 
  ggplot(aes(x = genotype, y = expression)) +
  geom_point() +
  geom_jitter() +
  stat_summary(aes(y = expression, group=1), fun.y = mean, geom="line") +
  facet_wrap(~gene)



# ====================================================================================================================
# ====================================================================================================================
# ====================================================================================================================
# ====================================================================================================================



# TRY HBV DATASET

setwd("~/ws/STAT540-instructors-only/seminars/seminars_winter_2017/seminar9/")

hbvDesign <- read_csv("data/GSE38941_design.csv") 
hbvExpressions <- read_csv("data/GSE38941_expressions.csv") %>% dplyr::select(gene = GeneSymbol, everything())

hbvDesign$disease <- factor(hbvDesign$disease, levels = c("reference_subject_role", "Acute_hepatic_failure_|_hepatitis_B_|"))

targetSamples <- hbvDesign %>% arrange(disease)

targetMatrix <- model.matrix(~disease, targetSamples)

# targetExpressions <- hbvExpressions[c("gene", targetSamples$Bioassay)] %>% 
#   group_by(gene) %>% 
#   summarize_all(.funs = mean) %>% 
#   ungroup() %>% 
#   na.omit() %>% 
#   as.data.frame() %>% 
#   column_to_rownames("gene")

targetExpressions <- hbvExpressions[c("Probe", targetSamples$Bioassay)] %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  column_to_rownames("Probe")


# fit linear model
targetFit <- lmFit(targetExpressions, targetMatrix)

# run ebayes to calculate moderated t-statistics
targetFitEbayes <- eBayes(targetFit)

# get result
targetDeaResult <- topTable(targetFitEbayes, n = Inf) %>% rownames_to_column("gene") %>% as_tibble()

# get probe mapping to map probes to genes
probeMappings <- read_csv("data/GSE38941_expressions.csv") %>% dplyr::select(probe = Probe, gene = GeneSymbol)

# join probe and gene names
targetDeaResult %>% 
  select(probe = gene, logFC) %>% 
  left_join(probeMappings, by = "probe") %>% 
  arrange(desc(abs(logFC))) %>% 
  group_by(gene) %>%
  summarize(logFC = head(logFC, 1)) -> deaSummaryData

# targetDeaResult %>% 
#   select(probe = gene, logFC) %>% 
#   left_join(probeMappings, by = "probe") %>% 
#   mutate(fdsgsdfds = (logFC / abs(logFC)))

# deaSummaryData %>% write_csv("data/ranked_gene_list.csv")

# ==================== correlating with multifunction


# get multifunctionality score for human GO
annotations <- GO.human[c("name", "GO")]
annotations <- make_annotations(annotations, unique(annotations[["name"]]), unique(annotations[["GO"]]))
multiFuncScores <- calculate_multifunc(annotations)

multiFuncScores <- multiFuncScores %>% as_tibble() %>% dplyr::select(gene = Gene, everything())

# merge dea result & mf scores
xx <- targetDeaResult %>% 
  select(Probe = gene, everything()) %>% 
  left_join(hbvExpressions %>% select(gene, Probe), by = "Probe") %>% 
  left_join(multiFuncScores, by = "gene") %>% 
  na.omit()

xx %>% 
  ggplot(aes(x = abs(logFC), y = log(MF.score))) + 
  geom_point(alpha = 0.1) +
  ggtitle(paste0("r = ", round(cor(abs(xx$logFC), log10(xx$MF.score), method = "spearman"), 2)))

xx %>% 
  ggplot(aes(x = -log10(P.Value), y = log(MF.score))) + 
  geom_point(alpha = 0.1) +
  ggtitle(paste0("r = ", round(cor(-log10(xx$P.Value), log10(xx$MF.score), method = "spearman"), 2)))


xx %>% 
  ggplot(aes(x = abs(t), y = log(MF.score))) + 
  geom_point(alpha = 0.1) +
  ggtitle(paste0("r = ", round(cor(xx$t, log10(xx$MF.score), method = "spearman"), 2)))


# ====================

ermineInput <- targetDeaResult %>% 
  arrange(desc(abs(logFC))) %>% 
  select(gene) %>% 
  mutate(rank = seq(1, nrow(targetDeaResult))) %>% 
  mutate(rank = (1 - rank/max(rank))) %>% 
  as.data.frame() %>% 
  column_to_rownames("gene")

roc(scores = ermineInput, 
    scoreColumn = 1, 
    bigIsBetter = TRUE,
    annotation = "GPL570", 
    aspects = "B",
    geneSetDescription = "GO.xml.gz") -> rocResult

precRecall(scores = ermineInput, 
           scoreColumn = 1, 
           bigIsBetter = TRUE,
           annotation = "GPL570", 
           aspects = "B",
           geneSetDescription = "GO.xml.gz") -> prResult

prResult$results %>% ggplot(aes(x = Pval, y = MFPvalue)) + geom_point()

prResult$results$Pval %>% cor(prResult$results$MFPvalue)

prResult$results %>%
  select(ID, Pval, MFPvalue) %>% 
  melt() %>% 
  as_tibble() %>% 
  ggplot(aes(x = variable, y = value)) +
  geom_point() +
  geom_line(aes(group = ID))



