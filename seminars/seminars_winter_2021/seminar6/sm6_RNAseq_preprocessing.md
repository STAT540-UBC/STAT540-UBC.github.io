Seminar 6: RNA-seq preprocessing: BAM files to count tables
================
Keegan Korthauer
(4 February 2021)

# Overview

This seminar focuses on a portion of the upstream preprocessing of
RNA-seq: going from aligned reads (BAM files) to count tables that will
be used in differential expression analysis. Note that we are skipping
the step of going from raw reads (fastq files) to aligned reads (BAM
files), as this is outside the scope of this course. Refer to the
lecture materials on High-dimensional genomics assays & data for details
and resources.

In order to use RNA-Seq data for the purposes of differential expression
analysis, we need to obtain an expression value for each gene (or
transcript). The digital nature of RNA-seq data allows us to use the
number of reads (an integer count) which align to a specific feature
(gene or transcript) as its expression value. In this seminar, we will
explore how to use BAM (or SAM) files to generate a count table of reads
aligning to the human transcriptome.

By the end of this seminar, you should be able to:

-   use the
    [`Rsamtools`](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html)
    package to **sort and index a BAM file**
-   use the
    [`Rsubread`](https://bioconductor.org/packages/release/bioc/html/Rsubread.html)
    package to **count the number of reads that align to each
    gene/transcript**

# BAM/SAM - Aligned Sequence data format

SAM (sequence alignment map) and BAM (the binary version of SAM) files
are the preferred output of most alignment tools. SAM and BAM files
carry the same information. The only difference is that SAM is human
readable, meaning that you can visually inspect each of the reads. You
can learn more about the SAM file format
[here](https://en.wikipedia.org/wiki/SAM_(file_format)).

The
[`Rsamtools`](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html)
package is one of many software tools that have been developed for
working with SAM/BAM alignment files - it is an R implementation of the
widely used command line tool [Samtools](http://www.htslib.org/). The
[`Rsubread`](https://bioconductor.org/packages/release/bioc/html/Rsubread.html)
package is a tool that summarizes SAM/BAM files into count tables - it
is an R implementation of the command line program
[Subread](http://subread.sourceforge.net/) which has tools for alignment
and quantification. The tool we’ll use within Subread is the
[featureCounts](https://academic.oup.com/bioinformatics/article/30/7/923/232889)
tool.

Although you could use the standalone command line versions of both of
these tools, to make things as convenient as possible, we’ll use the R
implementations and run today’s seminar entirely within R.

# HeLa cell line transcription

The alignment BAM file we will be working with was created by aligning
raw reads from an RNA-seq experiment of HNRNPC gene knockdown cell line
and control HeLa cells ([Zarnack et
al. 2012](http://europepmc.org/article/MED/23374342)). We will
specifically work with one of the control HeLa cell samples. [HeLa
cells](https://en.wikipedia.org/wiki/HeLa) are immortalized cervical
cancer cells derived taken from cancer patient Henrietta Lacks in 1951.
Henrietta Lacks was a 31-year-old African-American mother of five and a
patient at Johns Hopkins Hospital. She died from cervical cancer on
October 4, 1951. Her cells were taken from her without her knowledge or
consent, and the HeLa cell line is now the oldest and most used immortal
cell line in scientific research. To read more about Henrietta’s legacy
and the ethical issues involved in procuring and using her cells, read
[The Immortal Life of Henrietta
Lacks](http://rebeccaskloot.com/the-immortal-life/).

In this tutorial, we are working with an already aligned BAM file, but
the raw reads (fastq files) from this experiment are available in the
[European Nucleotide
Archive](https://www.ebi.ac.uk/ena/browser/view/PRJEB3048).

# RNA-seq upstream analysis pipeline

In this tutorial we will only be dealing with a portion of the RNA-seq
analysis pipeline. Here is an overview of the entire pipeline (image
source: [Yang & Kim
2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4742321/)). Note that
this figure is simplistic in its depiction of downstream analysis -
differential expression is **not** the only downstream analysis in
typical experiments. For example, we’ll learn about things like gene set
enrichment analysis later in the course.

![Typical workflow for RNA sequencing data analysis (Yang & Kim
2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4742321/bin/gni-13-119-g001.jpg)

As you can see, there are many steps to go from the raw reads to the
aligned bam files (red to light blue). These have been done for us in
this example, as we are making use of the aligned data in the
Bioconductor package
[RNAseqData.HNRNPC.bam.chr14](http://bioconductor.org/packages/release/data/experiment/html/RNAseqData.HNRNPC.bam.chr14.html).
Briefly, this package has preprocessed the data from [Zarnack et
al. 2012](http://europepmc.org/article/MED/23374342). The reads
(paired-end) were aligned to the full human genome (version hg19) with
TopHat2. In addition, the data have been subset to only include reads
mapped to Chromosome 14. This will reduce the computational time for our
tutorial compared to running on the entire genome.

**In this tutorial, we are only going to deal with the dark blue box:
“Expression Quantification”. Specifically, we are going to learn how to
convert a SAM/BAM alignment file to a count table that can be used for
differential expression analysis. Along the way, we will learn more
about the SAM and BAM files.**

# Loading libraries

Before we begin, we’ll need to load the necessary R packages for this
seminar. It is likely you don’t already have all of these installed on
your machine - if this is true, you’ll need to run this code chunk to
install them (currently set to `eval = FALSE`):

``` r
library(BiocManager)
install("Rsamtools")
install("Rsubread")
install("RNAseqData.HNRNPC.bam.chr14")
```

After they are successfully installed, we’ll load them for use in our
session.

``` r
library(Rsamtools)
library(Rsubread)
library(RNAseqData.HNRNPC.bam.chr14)
```

# Viewing BAM files

The first thing we are going to do is look at the alignment file, so we
can learn more about the structure of BAM/SAM files. Since we’re not
downloading the BAM file directly, but instead using a BAM filea
included in the `RNAseqData.HNRNPC.bam.chr14` package, we first need to
find out the file location of where the BAM file is stored. From the
package documentation, we learn that the file names are stored in the
object `RNAseqData.HNRNPC.bam.chr14_BAMFILES`.

``` r
RNAseqData.HNRNPC.bam.chr14_BAMFILES
```

    ##                                                                                                                ERR127306 
    ## "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/RNAseqData.HNRNPC.bam.chr14/extdata/ERR127306_chr14.bam" 
    ##                                                                                                                ERR127307 
    ## "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/RNAseqData.HNRNPC.bam.chr14/extdata/ERR127307_chr14.bam" 
    ##                                                                                                                ERR127308 
    ## "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/RNAseqData.HNRNPC.bam.chr14/extdata/ERR127308_chr14.bam" 
    ##                                                                                                                ERR127309 
    ## "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/RNAseqData.HNRNPC.bam.chr14/extdata/ERR127309_chr14.bam" 
    ##                                                                                                                ERR127302 
    ## "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/RNAseqData.HNRNPC.bam.chr14/extdata/ERR127302_chr14.bam" 
    ##                                                                                                                ERR127303 
    ## "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/RNAseqData.HNRNPC.bam.chr14/extdata/ERR127303_chr14.bam" 
    ##                                                                                                                ERR127304 
    ## "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/RNAseqData.HNRNPC.bam.chr14/extdata/ERR127304_chr14.bam" 
    ##                                                                                                                ERR127305 
    ## "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/RNAseqData.HNRNPC.bam.chr14/extdata/ERR127305_chr14.bam"

Here, we’ll of one \# control replicate of HELA RNA-seq, subsetted to
Chr 14

This file is not human readable, so we are going to convert it into a
SAM file, which we can read.

Now that we’re in the correct directory, we can convert the BAM file
into a SAM file, so we can view it.

    samtools view -o hela.sam ../stat540/HeLa.bam 

The *-o* flag signifies that hela.sam is the desired output file. The
BAM file I want to analyze is located in the common folder (stat540),
which is one level above my current folder. The characters “..” mean
“one level above my current folder”. So in “../stat540/hela.bam”, I’m
telling the computer to look for the file HeLa.bam in a folder called
stat540 that is one level above my current directory.

Now if you type:

    ls

you should see that you now have a new file: hela.sam.

If you type:

    less hela.sam

The contents of the file will print to console.

![HeLa head](figures/hela_head.jpg)

You can see that each line corresponds to a single read, and each field
refers to a different descriptor of the read. For a detailed explanation
of what each field means, see
[here](http://biobits.org/samtools_primer.html#UnderstandingtheSAMFormat).
Understanding what each field means is not of immediate relavence to the
following steps. However, it is always good to have a general
understanding of the structure of the files you are working with.

Type “q” to return to the command prompt.

\#Sorting and indexing BAM files

In order to be able to do things like extract reads that align to a
certain chromosome, or that map to a certain genomic region, we need to
generate a BAM file index. To do this, we first have to sort the file.
We can do this using the following command to sort our bam file by
chromosome and positions of each read:

    samtools sort ../stat540/HeLa.bam -o hela_sorted.bam

Once the file has been sorted, we can generate an index file.

    samtools index hela_sorted.bam

*Note: you must always sort before indexing a file*

# Generating a count table

Let’s take stock of where we are: we have a BAM file, a SAM file and a
BAI (BAM index) file. It is hard to do much analysis with what we have,
because even though we know which chromosome our reads map to, we don’t
know which genes they map to.

A variety of different tools can be used to generate a count table - i.e
a table of how many reads align to each gene. HTSeq and RSEM are two
examples of such programs. Each tool works optimally with the output of
different aligners, so it is best to read up on how the BAM files you
are working with were processed, and to use the feature counting tools.

Today, we are going to use HTSeq, as the BAM files we are using were
generated using the STAR tool, which generates BAM files that are
compatible with HTSeq. Because we are trying to correlate each read with
a transcript ID, we need an alignment file that contains information on
the genomic coordinates of each transcript. This alignment file already
exists in your repository, and is called gencode.gtf.

    htseq-count -s yes -a 10 hela.sam ../stat540/gencode.gtf > hela_count.txt

Let’s break this down.

The *-s* flag refers to whether the reads are stranded. In this case, we
know they are, so we write “yes”. Certain alignment
procedures/sequencing protocols may require that you specify “reverse”
instead for a stranded experiment. The difference in the number of
recognized features here should be drastic if this is specified
incorrectly The *-a* flag refers to the minimum quality of the read we
want to allow. Our quality score filter is set at 10, the maximum score
is 255. The SAM and GTF files are self-explanatory. The “&gt;” operator
means that we want to print the output of the commands on the left to
the file on the right.

The resulting file should look like this:

``` r
HTseq_output <- read.table("hela_count.txt", sep = "\t")
head(HTseq_output)
```

    ##                   V1 V2
    ## 1 ENSG00000000003.14  0
    ## 2  ENSG00000000005.5  0
    ## 3 ENSG00000000419.12  0
    ## 4 ENSG00000000457.13  0
    ## 5 ENSG00000000460.16 16
    ## 6 ENSG00000000938.12  0

This can be used as input to a differential expression tool such as
edgeR or Deseq2.

# Take home exercise

Now we have a table that summarizes the number of reads that align to
each transcript. The is the perfect input to use for algorithms like
edgeR, whose algorithm requires reads be in integer counts, not
normalized by library size. Other tools, however, may require that reads
are normalized first.

Your take home exercise is to convert the existing count data into
Transcripts per Million (CPM), as defined by Bo Li et al in [this
paper](https://academic.oup.com/bioinformatics/article/26/4/493/243395/RNA-Seq-gene-expression-estimation-with-read).

Compare your output with the file in the repository for this course.

This take home exercise will not be graded.

# Relevant Tutorial for RNA-Seq to Differential Expression Pipelines

[Griffifth Lab RNA-Seq
Tutorial](https://github.com/griffithlab/rnaseq_tutorial/wiki)
