Seminar 6: RNA-seq preprocessing - BAM files to count tables
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
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
```

# Viewing BAM files

The first thing we are going to do is look at the alignment file, so we
can learn more about the structure of BAM/SAM files. Since we’re not
downloading the BAM file directly, but instead using a BAM file included
in the `RNAseqData.HNRNPC.bam.chr14` package files, we first need to
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

We can see we have 8 BAM files. We’ll use sample ERR127306, which we can
determine from the [data
repository](https://www.ebi.ac.uk/ena/browser/view/PRJEB3048) is a
control HeLa sample.

``` r
bamfile <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[grepl("ERR127306",
                                                RNAseqData.HNRNPC.bam.chr14_BAMFILES)]
```

This file is not human readable, so we are going to convert it into a
SAM file, which we can read. The `asSam` function converts BAM files to
SAM files. We’ll save a SAM file to the current working directory.

``` r
asSam(bamfile, destination = "hela")
```

    ## [1] "hela.sam"

We have created the file “hela.sam”. We’ll take a peek at the top of the
header:

``` r
cat(system("head hela.sam", intern = TRUE), sep = '\n')
```

    ## @HD  VN:1.0  SO:coordinate
    ## @SQ  SN:chr1 LN:249250621
    ## @SQ  SN:chr10    LN:135534747
    ## @SQ  SN:chr11    LN:135006516
    ## @SQ  SN:chr11_gl000202_random    LN:40103
    ## @SQ  SN:chr12    LN:133851895
    ## @SQ  SN:chr13    LN:115169878
    ## @SQ  SN:chr14    LN:107349540
    ## @SQ  SN:chr15    LN:102531392
    ## @SQ  SN:chr16    LN:90354753

This contains information on the format and contents of the file. And at
the last few lines of the file, which illustrates some reads and their
mapping location, along with lots of extra mysterious-looking
information:

``` r
cat(system("tail hela.sam", intern = TRUE), sep = '\n')
```

    ## ERR127306.10965452   163 chr14   106981278   50  72M =   106981361   155 AAGAAACTTACAGTCATGATGGAAAGGGCAGCAAACACGTACCTCTTCACATGGCTGCAGGAGAGAGAAATG    +++++++++++*++++++++++++++++*+++++++++++++++++++)++++*+*+&+++++*++++*+)*    AS:i:-3 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:69G2   YT:Z:UU NH:i:1
    ## ERR127306.10965452   83  chr14   106981361   50  72M =   106981278   -155    GGGAGGCCCTTTATATAACCATCAGGTCTTGTGAGAACTCACTCACTAATAGGATAAAAGCATGGAGAGAAC    ##+)+++%+*&+++++++*++)+++++++++++++++)+++*+++++++)++*+%++++*++++++)+++++    AS:i:-3 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:68A3   YT:Z:UU NH:i:1
    ## ERR127306.11919172   163 chr14   106986423   50  72M =   106986528   177 CGGTGGTGTCCTCAGTGGCCCCTGGTGTCCTGAGCATCCCCTGGTGTCCTGAGTGCCCCCTGGTGGTTGTGA    ++**+++&+++)++)*++++++()*"''''++(+***)++*&**&)'')&'%&%'&'&''!$""%$"$$!!!    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:72 YT:Z:UU NH:i:1
    ## ERR127306.11919172   83  chr14   106986528   50  72M =   106986423   -177    CCTGAGAGCCCCTTGCTGTCCTGAGCACCTCCTGGTGTTCTGAGCGCCCTCTGGTGTTCTGATCACTCTCTG    &'&)%*&&)(*****+*+**+*+*++))++*++++++++++++++++**+(+++++++++++++++++++++    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:72 YT:Z:UU NH:i:1
    ## ERR127306.3567919    99  chr14   106989680   50  72M =   106989780   172 CAACTTTTATTTCTTAAACACAAGACATTCCAATGAGAAAGCTGTTCTCAGGTGAGCTGTCGAGCAGGGAGG    ++++++++*+++++++*++++++)++++++++++++++++++++++++++++)+&+)++*++*))*+++)+&    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:72 YT:Z:UU NH:i:1
    ## ERR127306.3567919    147 chr14   106989780   50  72M =   106989680   -172    GAAATTGCTGAAACTTGAAGACCAAGGCCACCTCTGAGGGGCAGAGATCCACCTATGAGTACATCACATCAG    (+++++++&+++)'*+++++*++++++*+*+++++++++++++++++)++++++*+++++++++++++++++    AS:i:-3 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:62G9   YT:Z:UU NH:i:1
    ## ERR127306.21510817   99  chr14   106994763   50  72M =   106994819   128 CAAAGCTGGATGTGTCTAGTGTTTTTATCAGAACCCACTTTCCGTAATAAGAGCATGTGTGGTTTTGCTGCC    ++++++++++++++++++*+)++++++++++++++++++++++++++++*+*++*++)+*+(*++***+**+    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:72 YT:Z:UU NH:i:1
    ## ERR127306.21510817   147 chr14   106994819   50  72M =   106994763   -128    GTGTGGTTTTGCTGCCCTCCAGCACTCTTCTGAAAATATGGAGAGAACTAGGATCCAGGCACATTAATTTTC    +%++*+*(**++++**+&++*++)+++++++++++++++++++++++*++++++++++++++++++++++++    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:72 YT:Z:UU NH:i:1
    ## ERR127306.661203 163 chr14   107003080   50  72M =   107003171   163 AAGGAACCCTTGAACTCCCTTGGCGACATGTACTCCTACTAGCACTGTGGCATTATGGTCCTCTGCCTCAAG    +++++++++++++++++++++++++++++++++++*++*+*')++++)*+++++(+*()+'*++*+((**')    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:72 YT:Z:UU NH:i:1
    ## ERR127306.661203 83  chr14   107003171   50  72M =   107003080   -163    CATGACTTGATGGCTGGAACAAATACATTTAGAGATTTTACCTCCAATACTAGCCTTTGCCATACAGTATTT    +(&+*)+)+++++*+++*+*++++*&+&++++++)+'+)(+)+++++++++*++++++++++*+++++++++    AS:i:-3 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:46G25  YT:Z:UU NH:i:1

You can see that each line corresponds to a single read, and each field
refers to a different descriptor of the read. For a detailed explanation
of what each field means, see [this
document](https://samtools.github.io/hts-specs/SAMv1.pdf). Understanding
what each field means is not of immediate relevence to the following
steps. However, it is always good to have a general understanding of the
structure of the files you are working with.

# Sorting and indexing BAM files

In order to be able to do things like extract reads that align to a
certain chromosome, or that map to a certain genomic region, we need to
have a BAM file index. In this example dataset, one has already been
generated for us. But for illustration we’ll review how to create one.
To do this, we first have to sort the file. We can do this using the
following command to sort our BAM file by chromosome and positions of
each read:

``` r
sortBam(bamfile, destination="hela_sorted")
```

    ## [1] "hela_sorted.bam"

This generated a file ‘hela\_sorted.bam’ in our current working
directory. Once this file has been sorted, we can generate an index
file.

``` r
indexBam("hela_sorted.bam")
```

    ##       hela_sorted.bam 
    ## "hela_sorted.bam.bai"

*Note: you must always sort before indexing a file*

# Viewing basic information about our BAM file

We don’t want to have to look through our entire file line by line.
Let’s instead use some of the handy parsing functions in `RSamtools` to
view some basic summary info about our file (note these functions can
use the non-human readable BAM file, which is more efficient than the
human-readable SAM file.

First, we create a special reference to the BAM file, which will be
passed to the parsing functions (this contains more information than
just the path to the BAM file - also points to the location of the BAM
index file):

``` r
# view basic info of bam
bamFile <- BamFile(bamfile)
bamFile
```

    ## class: BamFile 
    ## path: /Library/Frameworks/R.framework/Versions/4.0/Resour.../ERR127306_chr14.bam
    ## index: /Library/Frameworks/R.framework/Versions/4.0/R.../ERR127306_chr14.bam.bai
    ## isOpen: FALSE 
    ## yieldSize: NA 
    ## obeyQname: FALSE 
    ## asMates: FALSE 
    ## qnamePrefixEnd: NA 
    ## qnameSuffixStart: NA

We’ll pass this special reference to `countBam`, which will by default
tell us how many reads and bases are contained in our BAM file.

``` r
countBam(bamFile) 
```

    ##   space start end width                file records nucleotides
    ## 1    NA    NA  NA    NA ERR127306_chr14.bam  800484    57634848

We can see we have 800484 reads, which amounts to a total of about 5.76
million nucleotides. From this output, can you determine how long are
our reads?

Next, we’ll use the `quickBamFlagSummary` function to look at the
different ‘flags’ for our reads.

``` r
quickBamFlagSummary(bamFile)
```

    ##                                 group |    nb of |    nb of | mean / max
    ##                                    of |  records |   unique | records per
    ##                               records | in group |   QNAMEs | unique QNAME
    ## All records........................ A |   800484 |   393300 | 2.04 / 10
    ##   o template has single segment.... S |        0 |        0 |   NA / NA
    ##   o template has multiple segments. M |   800484 |   393300 | 2.04 / 10
    ##       - first segment.............. F |   400242 |   393300 | 1.02 / 5
    ##       - last segment............... L |   400242 |   393300 | 1.02 / 5
    ##       - other segment.............. O |        0 |        0 |   NA / NA
    ## 
    ## Note that (S, M) is a partitioning of A, and (F, L, O) is a partitioning of M.
    ## Indentation reflects this.
    ## 
    ## Details for group M:
    ##   o record is mapped.............. M1 |   800484 |   393300 | 2.04 / 10
    ##       - primary alignment......... M2 |   779228 |   389614 |    2 / 2
    ##       - secondary alignment....... M3 |    21256 |    10052 | 2.11 / 8
    ##   o record is unmapped............ M4 |        0 |        0 |   NA / NA
    ## 
    ## Details for group F:
    ##   o record is mapped.............. F1 |   400242 |   393300 | 1.02 / 5
    ##       - primary alignment......... F2 |   389614 |   389614 |    1 / 1
    ##       - secondary alignment....... F3 |    10628 |    10052 | 1.06 / 4
    ##   o record is unmapped............ F4 |        0 |        0 |   NA / NA
    ## 
    ## Details for group L:
    ##   o record is mapped.............. L1 |   400242 |   393300 | 1.02 / 5
    ##       - primary alignment......... L2 |   389614 |   389614 |    1 / 1
    ##       - secondary alignment....... L3 |    10628 |    10052 | 1.06 / 4
    ##   o record is unmapped............ L4 |        0 |        0 |   NA / NA

These ‘flags’ can be used to get detailed counts of certain types of
reads. As a quick example, we might want to grab the count of reads that
on the minus strand, we could do the following:

``` r
params <- ScanBamParam(flag = scanBamFlag(isMinusStrand = TRUE))
countBam(bamFile, param = params)
```

    ##   space start end width                file records nucleotides
    ## 1    NA    NA  NA    NA ERR127306_chr14.bam  400242    28817424

There’s a ton more possibilities - explore the documentation of the
`RSamtools` package (in particular the `scanBamFlag` and `scanBamParam`
functions), as well as the [meaning of SAM
flags](https://www.samformat.info/sam-format-flag) if you’d like to
learn more. These can be really helpful for quality control of our
alignment - to make sure we have what we generally expect (e.g. if we
performed stranded sequencing, we expect that roughly half of the reads
map to each strand). We don’t need to dive too deep into the flags if we
are satisfied with our alignment and just want to proceed to get our
count table, however.

# Generating a count table

Let’s take stock of where we are: we have a BAM file, a SAM file and a
BAI (BAM index) file. It is hard to do much analysis with what we have,
because even though we know which chromosome our reads map to, we don’t
know which *genes* (or *transcripts*) they map to.

A variety of different tools can be used to generate a count table - i.e
a table of how many reads align to each gene. `HTSeq`, `RSEM`, and
`featureCounts` are examples of such programs that are command line
tools. Each tool works optimally with the output of different aligners,
so it is best to read up on how the BAM files you are working with were
processed, and to use the feature counting tools.

Today, we are going to use `featureCounts` as implemented in the
`RSubread` package, as the BAM files we are using were generated using
the Tophat2 aligner, which generates BAM files that are compatible with
`featureCounts`. Because we are trying to associate each read with a
gene ID, we need an *alignment file* (such as a GTF file) that contains
information on the genomic coordinates of each transcript.

This type of alignment file can be downloaded from a database such as
[Gencode](https://www.gencodegenes.org/), but conveniently, `RSubread`
has some RefSeq annotations built-in, including for mm10, mm9, hg38, and
hg19. It is critical that you use the **same genome version as the
genome version used to align raw reads to your BAM file**. Our BAM file
was generated by aligning reads to hg19, so that’s what we’ll use. These
built-in annotations allow for summarizing counts over genes (the
default, and also the option we’ll choose) or exons. Note that external
annotation files can be used for other genomes and to summarize over
transcript/isoform annotations.

Now we’re ready to apply the `featureCounts` function that will do the
counting of reads over genes. Note that this step can take a while if we
have a large BAM file. We specify `isPairedEnd = TRUE` since as we
mentioned above our sequencing protocol generated paired end reads. We
also specify `strandSpecific = 1` to indicate that our sequencing
protocol was stranded.

``` r
counts <- featureCounts(bamfile,
                        annot.inbuilt = "hg19",
                        isPairedEnd = TRUE,
                        strandSpecific = 1)
```

    ## NCBI RefSeq annotation for hg19 (build 37.2) is used.
    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.4.2
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 1 BAM file                                       ||
    ## ||                                                                            ||
    ## ||                           ERR127306_chr14.bam                              ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : inbuilt (hg19)                                   ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 1                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : not counted                                      ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file hg19_RefSeq_exon.txt ...                              ||
    ## ||    Features : 225074                                                       ||
    ## ||    Meta-features : 25702                                                   ||
    ## ||    Chromosomes/contigs : 52                                                ||
    ## ||                                                                            ||
    ## || Process BAM file ERR127306_chr14.bam...                                    ||
    ## ||    Strand specific : stranded                                              ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 400242                                               ||
    ## ||    Successfully assigned alignments : 171267 (42.8%)                       ||
    ## ||    Running time : 0.02 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//

You can change other options, such as the `minMQS`, the minimum mapping
quality score to include in the counts. To learn more about different
parameters, see the help fule for the `featureCounts` function.

Let’s take a peek at the output.

``` r
str(counts)
```

    ## List of 4
    ##  $ counts    : int [1:25702, 1] 0 0 0 0 0 0 0 0 0 0 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:25702] "653635" "100422834" "645520" "79501" ...
    ##   .. ..$ : chr "ERR127306_chr14.bam"
    ##  $ annotation:'data.frame':  25702 obs. of  6 variables:
    ##   ..$ GeneID: int [1:25702] 653635 100422834 645520 79501 729737 100507658 100132287 100288646 729759 100131754 ...
    ##   ..$ Chr   : chr [1:25702] "chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1" "chr1" "chr1;chr1;chr1" "chr1" ...
    ##   ..$ Start : chr [1:25702] "14362;14970;15796;16607;16858;17233;17606;17915;18268;24738;29321" "30366" "34611;35277;35721" "69091" ...
    ##   ..$ End   : chr [1:25702] "14829;15038;15947;16765;17055;17368;17742;18061;18366;24891;29370" "30503" "35174;35481;36081" "70008" ...
    ##   ..$ Strand: chr [1:25702] "-;-;-;-;-;-;-;-;-;-;-" "+" "-;-;-" "+" ...
    ##   ..$ Length: int [1:25702] 1769 138 1130 918 3402 793 4370 472 939 1019 ...
    ##  $ targets   : chr "ERR127306_chr14.bam"
    ##  $ stat      :'data.frame':  14 obs. of  2 variables:
    ##   ..$ Status             : chr [1:14] "Assigned" "Unassigned_Unmapped" "Unassigned_Read_Type" "Unassigned_Singleton" ...
    ##   ..$ ERR127306_chr14.bam: int [1:14] 171267 0 0 0 0 0 0 0 0 0 ...

It looks like the counts are stored in the `counts` slot and the gene
annotation is stored in the `annotation` slot. There are also a couple
more slots for experiment metadata (`targets` has sample name, and
`stat` contains some stats about our BAM file).

Let’s take a peek at the count table, subsetted for genes on chromosome
14 (since our bam file was subsetted to reads mapping to chr 14).

``` r
chr14 <- which(grepl("chr14", counts$annotation$Chr))
str(chr14)
```

    ##  int [1:1044] 6746 6747 6748 6749 6750 6751 6752 6753 6754 6755 ...

``` r
head(counts$counts[chr14,])
```

    ## 100132257    440153    404785 100506092 100506303 100506350 
    ##         0         0         0         0         4         1

``` r
sum(counts$counts[chr14,] > 0)
```

    ## [1] 499

``` r
sum(counts$counts[-chr14,] > 0)
```

    ## [1] 0

Looks like there are 1044 genes on chromosome 14, 499 of which have
non-zero expression. And as we expected, no genes on any other
chromosomes have any counts.

Note that the gene IDs are RefSeq IDs. If we’d like to convert them to
another ID, such as Hugo symbols or Ensembl IDs, we can use things like
the
[`biomaRt`](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
package, or annotation packages specific to the organism we are dealing
with (e.g. for human:
[`org.Hs.eg.db`](https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html))

Let’s look at a density of expression of the chr14 genes
(log-transformed).

``` r
data.frame(count = counts$counts[chr14,]) %>%
  ggplot(aes(x = count + 1)) +
    geom_density() +
    scale_x_log10()
```

<img src="sm6_RNAseq_preprocessing_files/figure-gfm/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

We can save the count table to a text file for importing into future R
sessions for downstream analysis (such as with
[`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
or
[`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html).
We also include a column for gene length, since we might want to use
this to obtain quantities like RPKM.

``` r
write.table(data.frame(counts$annotation[,c("GeneID","Length")],
                       counts$counts),
            file = "hela_counts.txt", quote = FALSE, sep = "\t", row.names = FALSE)
```

Lastly, we’ll cleanup our directory by removing the SAM/BAM/index files
we generated earlier, and the text file of counts (comment the last line
if you’d like to keep the text file).

``` r
file.remove("hela.sam")
file.remove("hela_sorted.bam")
file.remove("hela_sorted.bam.bai")
file.remove("hela_counts.txt")
```

# Take home exercise

Now we have a table that summarizes the number of reads that align to
each transcript. The is the perfect input to use for methods mentioned
above like
[`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
and
[`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html),
whose algorithm requires reads be in integer counts, not normalized by
library size. For visualization purposes, however, it may be helpful to
use reads that are normalized first.

Your take home exercise is to convert the existing count data into
Transcripts per Million (CPM), as defined by Bo Li et al in [this
paper](https://academic.oup.com/bioinformatics/article/26/4/493/243395/RNA-Seq-gene-expression-estimation-with-read).
Hint: check out the `cpm` and `rpkm` functions in the `edgeR` package.

This take home exercise will not be graded.
