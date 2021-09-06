library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(stringr)
library(argparse)
library(stringr)
library(plyr)

source("helper_functions.R")

#read input arguments and obtain VCF object and output file name
tmp <-  process_inputs()
vcf <- tmp[[1]]
output_file <- tmp[[2]]
samples <- tmp[[3]]






