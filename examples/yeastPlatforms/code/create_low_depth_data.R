## working directory assumed to be the parent of the code directory where this
## file lies; at the time of writing, that means the 'yeastPlatforms' directory

library(biomaRt)
library(Rsamtools)
library(GenomicAlignments)


ensembl <- useMart('ensembl')

yeast.ensembl <- useDataset("scerevisiae_gene_ensembl", ensembl)

annotation.fields <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")

gene.annotations <- getBM(annotation.fields, mart=yeast.ensembl)

bam.files <- Sys.glob("examples/yeastPlatforms/data/bam/*.bam")

sample.ids <- basename(bam.files)

sample.ids <- sub(".bam", "", sample.ids)

names(bam.files) <- sample.ids

range.data <- IRanges(start=gene.annotations$start_position, end=gene.annotations$end_position)

genomic.range.data <- GRanges(gene.annotations$chromosome_name, range.data, gene=gene.annotations$ensembl_gene_id)

bam.dat <- lapply(bam.files, readGAlignments) 

overlaps <- sapply(bam.dat, summarizeOverlaps, features=genomic.range.data, mode="IntersectionStrict", ignore.strand=TRUE)

counts <- sapply(overlaps, function(x) assays(x)$counts)

row.names(counts) <- gene.annotations$ensembl_gene_id

groups <- factor(c(rep("b", 3), c(rep("c", 3))))

write.table(counts, 'examples/yeastPlatforms/data/stampy.low.counts.tsv', row.name=TRUE, col.names=NA, sep="\t")
