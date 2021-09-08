README
================
Adriana Salcedo

## Dependencies

``` r
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(argparse)
library(plyr)
library(httr)
library(jsonlite)
source("helper_functions.R")
```

## Introduction

The R script main.R is the main interface for variant annotation, and it will interface with helper functions. It will first call process\_inputs to process the command line inputs and generate the VCF-class object. It accepts three arguments:

``` bash
Rscript main.R --vcf_file --output_file --chunk_size
```

It will then pass the arguments onto run\_annotation which accepts a VCF object, and two optional parameters. The parameter chunk\_size determines the number of variants that will be annotated at once. It thus allows users to prioritize memory usage, or efficiency depending on their use case. The default annotates the entire VCF at once. The second parameter variant\_effects is an optional data frame ranking variant effect severity. It should have two columns, one titled "impact" which has the variant impact label and one labelled "effect" which gives an ascending numeric ranking of the variant severity (so the most benign variant impacts will have the smallest effect and the most severe variants will have the largest). Otherwise the default is used which can be seen as follows:

``` r
default_variant_effect()
```

    ##           impact effect
    ## 1         coding      1
    ## 2     synonymous      2
    ## 3   intergenenic      3
    ## 4         intron      4
    ## 5       threeUTR      5
    ## 6        fiveUTR      6
    ## 7       promoter      7
    ## 8     spliceSite      8
    ## 9  nonsynonymous      9
    ## 10    frameshift     10
    ## 11      nonsense     11

A user can also generate a VCF-class object from file within R and annotate it with run\_annotation directly.

``` r
vcf <- readVcf("Challenge_data.vcf")
annotation <- run_annotation(vcf)
```

## Output

``` r
head(annotation)
```

    ##   chr   start        qual ref alt num_alt_alleles total_depth reference_depth
    ## 1   1  931393 2.17938e-13   G   T               1        4124            4029
    ## 2   1  935222 1.68667e+04   C   A               1        1134             480
    ## 3   1 1277533 2.81686e+04   T   C               1         786               0
    ## 4   1 1284490 6.30056e+03   G   A               1         228               0
    ## 5   1 1571850 0.00000e+00   G   A               1        4055            3961
    ## 6   1 1572579 0.00000e+00   A   G               1        3456            3430
    ##   n_alt_reads alt_read_percent variant_type              ExAC_AF        impact
    ## 1          95           2.3036          snp                 <NA>          <NA>
    ## 2         652          57.4956          snp    0.661110826852231 nonsynonymous
    ## 3         786              100          snp    0.997571146075144      promoter
    ## 4         228              100          snp    0.897299369568921      promoter
    ## 5          94           2.3181          snp 1.11375938342281e-05        intron
    ## 6          26           0.7523          snp  0.00217689117420759        intron

The output annotation is a data frame with the following columns:

chr: Chromosome name

start: Start of the variant

qual: Variant quality Phred-score

ref: the reference allele

alt: the alternate allele(s)

num\_alt\_alleles: the number of alternate alleles

total\_depth: the total depth of coverage at that site

reference depth: the number of reads supporting the reference allele

n\_alt\_reads: the number of reads supporting each alternate allele. If there is more than one alternate allele the number of reads is listed for each in the same order as they appear in the alt column.

alt\_read\_percent: the percent of reads supporting each alternate allele. If there is more than one alternate allele the percent of reads for each is listed for each in the same order as they appear in the alt column.

variant type: the variant class.

ExAC\_AF: The frequency of the variant in the ExAc database. If there is more than one alternate allele the allele frequencies for each is listed for each in the same order as they appear in the alt column.

variant type: the variant class.

impact: The impact of the variant as predicted from is location relating to known genes. If there were multiple possible impacts the most severe is listed (as determined by the variant\_effects dataframe).

The annotation will be writen to the user defined output file or the default based on the original vcf name \*\_annotation.txt.
