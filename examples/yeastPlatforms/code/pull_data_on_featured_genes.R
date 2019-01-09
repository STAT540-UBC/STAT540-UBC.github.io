## support for Homework 2 Q4b part ii
## pick some interesting hits and non-hits and revisit "raw" data

library(VennDiagram) # venn.diagram() for a sanity check

## working directory assumed to be the homework_02 directory
array.limma.results <- read.table('results/array-dea-results.tsv',
                                  header=TRUE, row.names=1, as.is=c('probe.id', 'gene.id'))
stampy.limma.results <- read.table('results/stampy.limma.results.tsv', header=TRUE,
                                   row.names=1, as.is=c('gene.id'))
stampy.edger.results <- read.table('results/stampy.edger.results.tsv', header=TRUE,
                                   row.names=1, as.is=c('gene.id'))

counts <- read.table('data/stampy.counts.tsv', header=TRUE, row.names=1)
str(counts)
(groups <- factor(substr(names(counts), start = 1, stop = 1)))
(groups <- factor(ifelse(groups == 'b', 'batch', 'chemostat')))

array.dat <- read.table('data/array-AFTER-SWAP-FIX.tsv', header=TRUE, row.names=1)
str(array.dat)

## useful in many places:
## genes that are represented in the array data by exactly one probe
str(arraySingles <- array.limma.results$gene.id)
table(foo <- table(arraySingles))
##    1    2    3 
## 5611   44    2 
str(arraySingles <- names(foo[foo == 1])) #5611 such genes/probes

## useful in many places:
## genes that are represented in the array and RNA-Seq data 
str(bothPlatforms <- intersect(array.limma.results$gene.id, stampy.limma.results$gene.id)) # 5585

## useful in many places:
## intersection of the above 
str(happyGenes <- intersect(arraySingles, bothPlatforms)) # 5543

## function for keeping hits and returning gene.id and rank
## based on function getSignificantGenes in homework 2 solution
jFun <- function(x, max.q.value) {
  genes <- subset(x, q.value < max.q.value, select = c(gene.id, q.value))
  genes$rank <- rank(genes$q.value)
  genes$q.value <- NULL
  return(genes)
}

max.q.value <- 1e-5
edger.de.ranks <- jFun(stampy.edger.results, max.q.value)
limma.de.ranks <- jFun(stampy.limma.results, max.q.value)
array.de.ranks <- jFun(array.limma.results, max.q.value)

de.ranks <- merge(edger.de.ranks, limma.de.ranks, by = "gene.id",
                  suffixes = c(".edR", ".liV"), all = TRUE)
de.ranks <- merge(de.ranks, array.de.ranks, by = "gene.id", all = TRUE)
names(de.ranks) <- c(names(de.ranks)[1:3], "rank.liA")
de.ranks$avgRank <- apply(de.ranks[ , c("rank.edR", "rank.liV", "rank.liA")],
                          1, mean, na.rm = TRUE)

## sanity check these counts against venn diagram
vennTmp <- list(Array=array.de.ranks$gene.id,
                edgeR=edger.de.ranks$gene.id, limma=limma.de.ranks$gene.id)
plot.new()
venn.plot <- venn.diagram(vennTmp, filename=NULL, height=600, width=600,
                          fill=c('red', 'blue', 'yellow'), force.unique = TRUE)
grid.draw(venn.plot)

with(de.ranks, table(is.na(rank.edR), is.na(rank.liV), is.na(rank.liA)))
## up to discrepancies caused by non-uniqueness ... match

## marshal the hits from various regions of the Venn diagram
hits <- with(de.ranks,
             by(de.ranks,
                list(ifelse(is.na(rank.liA), "xx", "liA"),
                     ifelse(is.na(rank.edR), "xx", "edR"),
                     ifelse(is.na(rank.liV), "xx", "liV")),
                function(x) return(x$gene.id[order(x$avgRank)])))
tmp <- levels(interaction(dimnames(hits)))
dim(hits) <- NULL
names(hits) <- tmp
str(hits)

## retain only intersection hits with happyGenes
hits <- lapply(hits, function(x) intersect(x, happyGenes))
str(hits)

## marshall non hits and filter:
## require presence in array and RNA-Seq data
## eliminate gene.ids that appear more than once in array data
str(hitGenes <- de.ranks$gene.id)
str(nonHitGenes <- setdiff(happyGenes, hitGenes))

## did I accomplish my goal?
all(nonHitGenes %in% stampy.limma.results$gene.id)
all(nonHitGenes %in% array.limma.results$gene.id)
table(table(array.limma.results$gene.id[array.limma.results$gene.id %in% nonHitGenes]))
## yes

## select actual genes to feature
str(hits)
str(nonHitGenes)
set.seed(234)
featGenes <- sapply(hits, function(x) if(length(x) > 0) x[1])
featGenes <- unlist(featGenes)
featGenes <- data.frame(hitType = factor(names(featGenes), levels = names(featGenes)),
                        gene.id = I(featGenes))
featGenes <- rbind(featGenes,
                   data.frame(hitType = I("xx.xx.xx"),
                              gene.id = sample(nonHitGenes, size = nrow(featGenes)))) #GCF: I can't get Andy's sample with this seed
tmp <- with(featGenes, paste(gene.id, hitType))
featGenes$id <- factor(tmp, levels = tmp)

plotDat <- data.frame(gene.id = factor(rep(featGenes$gene.id, each = ncol(counts)),
                                       levels = featGenes$gene.id),
                      hitType = factor(rep(featGenes$hitType, each = ncol(counts)),
                                       levels = featGenes$hitType),
                      id = factor(rep(featGenes$id, each = ncol(counts)),
                                  levels = featGenes$id),
                      cond = groups,
                      count = as.vector(t(counts[featGenes$gene.id, ])))
probes <- array.limma.results$probe.id[match(featGenes$gene.id,
                                             array.limma.results$gene.id)]
plotDat$arrExp <- as.vector(t(array.dat[probes, ]))
plotDat$status <- ifelse(plotDat$hitType == "xx.xx.xx", "not DE by any method",
                         "DE by >= 1 method")

str(plotDat)

write.table(plotDat, "data/featGenesExpressionData.tsv", row.names=TRUE, col.names=NA)  #GCF: I did not change Andy's output
save(plotDat, file = "data/featGenesExpressionData.robj")
dput(plotDat, file = "data/featGenesExpressionData_DPUT.txt")

#--------------------------------------------
# GCF: a simpler plotDat for HW2-Q3 in 2014

# Get array data for the genes are explored in Q2-e

featureMe <- c("YDR384C", "YDR345C")

(featureCounts <- counts[featureMe,])
featureDat <- data.frame(gene.id = factor(rep(rownames(featureCounts),
                                              ncol(featureCounts))),
                         cond = factor(rep(groups, each = nrow(featureCounts))),
                         rnaseqLogCounts = log2(unlist(featureCounts)))

probes<-array.results[match(featureMe,array.results$gene.id),c("probe.id","gene.id")]

(featureArray <- array.dat[probes$probe.id,])
featureDat.array <- data.frame(probe.id = factor(rep(rownames(featureArray),
                                                  ncol(featureArray))),
                            gene.id = factor(rep(probes$gene.id,
                                                 ncol(featureArray))),
                            cond = factor(rep(groups, each = nrow(featureArray))),
                            expression = unlist(featureArray))

# As I constructed both data.frames using the same structure, it is fine to add a vector of expression data to featureDat

featureDat$arrayExp <- featureDat.array$expression
featureDat$probe.id <- featureDat.array$probe.id

# Add genes identified by only one platform and boring genes 

# Array only: "YGL209W" from Andy's list

# EdgeR only:
edg.only <- hits[[6]][1:20]
goodGene <- edg.only[which(apply(counts[edg.only,],1,function(x)(min(x)>0))==1)]
goodProbe<-array.results$probe.id[match(goodGene,array.results$gene.id)]
array.dat[goodProbe,]  # from this plot, I've selected "1776330_at" matching to gene "YCL042W"

# Boring: "YBL025W" from Andy's list

# Create data.frame
interestMe <- c("YGL209W","YCL042W", "YBL025W")

(interestCounts <- counts[interestMe,])
interestDat <- data.frame(gene.id = factor(rep(rownames(interestCounts),
                                              ncol(interestCounts))),
                         cond = factor(rep(groups, each = nrow(interestCounts))),
                         rnaseqLogCounts = log2(unlist(interestCounts)))

interestProbes<-array.results[match(interestMe,array.results$gene.id),c("probe.id","gene.id")]

(interestArray <- array.dat[interestProbes$probe.id,])
interestDat.array <- data.frame(probe.id = factor(rep(rownames(interestArray),
                                                     ncol(interestArray))),
                               gene.id = factor(rep(interestProbes$gene.id,
                                                    ncol(interestArray))),
                               cond = factor(rep(groups, each = nrow(interestArray))),
                               expression = unlist(interestArray))

# As I constructed both data.frames using the same structure, it is fine to add a vector of expression data to featureDat

interestDat$arrayExp <- interestDat.array$expression
interestDat$probe.id <- interestDat.array$probe.id

plotDat <- rbind(featureDat,interestDat)
plotDat <- plotDat[order(plotDat$gene.id),]

write.table(plotDat, "data/featGenesData-q3.tsv", row.names=TRUE, col.names=NA)  
save(plotDat, file = "data/featGenesData-q3.robj")
dput(plotDat, file = "data/featGenesData-q3-DPUT.txt")

plotDat <- dget("data/featGenesData-q3-DPUT.txt") # creates plotDat
dePlot <- stripplot(gene.id ~ log.count + arrayExp, plotDat,
                    groups = cond,
                    auto.key = list(columns = 2, x = 0.05, y = 0.95, corner = c(0, 0.5)),
                    outer = TRUE, scales = list(x = list(relation = "free")),
                    xlab = "expression")

