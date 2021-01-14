# Feedback and Proposed Changes
* Run entirely within R, e.g;: 

```
library(Rsamtools)
library(Rsubread)
library(RNAseqData.HNRNPC.bam.chr14)

# get file location (within RNAseqData.HNRNPC.bam.chr14 package files) of one 
# control replicate of HELA RNA-seq, subsetted to Chr 14
bamfile <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[grepl("ERR127306",
                                                RNAseqData.HNRNPC.bam.chr14_BAMFILES)]
                                                
# copy bam file to wd
file.copy(bamfile, "hela")

# view basic info of bam
bamfile <- BamFile("hela.bam")
bamfile
seqinfo(bamfile)
quickBamFlagSummary(bamfile)
countBam(bamfile) 

# read in and look at 10 alignments
yieldSize(bamfile) <- 10
open(bamfile)
scanBam(bamfile)[[1]]$seq

# close those 10 alignments and go back to entire bam
close(bamfile)
yieldSize(bamfile) <- NA

# sort and index entire bam file
sortBam("hela.bam", destination="hela_sorted")
indexBam("hela.bam")

# count with subread
counts <- featureCounts("hela.bam",
                        annot.inbuilt = "hg19",
                        isPairedEnd = TRUE)
                        
# how many genes had nonzero counts? (remember this is just chr14)
sum(counts$counts>0)

# write to file
write.table(data.frame(counts$annotation[,c("GeneID","Length")],
                       counts$counts,
                       stringsAsFactors=FALSE),
            file="hela_counts.txt", quote=FALSE,sep="\t", row.names=FALSE)
```