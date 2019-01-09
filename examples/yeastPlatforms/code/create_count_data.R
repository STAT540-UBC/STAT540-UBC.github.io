## working directory assumed to be the parent of the code directory where this
## file lies; at the time of writing, that means the 'yeastPlatforms' directory

#setwd('/home/andrew/code/stat540/rmd/homework/homework_02')

counts <- read.table('data/SRS307298-rnaseq-RAW.tsv', header=TRUE, row.names=1)

stampy.counts <- counts[, grep('stampy', names(counts))]

tophat.counts <- counts[, grep('tophat', names(counts))]

names(stampy.counts) <- sub("stampy.", "", names(stampy.counts))

names(tophat.counts) <- sub("tophat.", "", names(tophat.counts))

stampy.counts.filtered <- stampy.counts[rowSums(stampy.counts) > 0, ]

tophat.counts.filtered <- tophat.counts[rowSums(tophat.counts) > 0, ]

write.table(stampy.counts.filtered, 'data/stampy.counts.tsv', row.names=TRUE, col.names=NA, sep="\t")

write.table(tophat.counts.filtered, 'data/tophat.counts.tsv', row.names=TRUE, col.names=NA, sep="\t")