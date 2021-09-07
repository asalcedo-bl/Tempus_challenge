library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(argparse)
library(plyr)
library(httr)
library(jsonlite)


source("helper_functions.R")

#read input arguments and obtain VCF object and output file name
tmp <-  process_inputs()
vcf <- tmp[[1]]
output_file <- tmp[[2]]
chunk_size <- tmp[[3]]

annotation <- run_annotation(vcf, chunk_size=chunk_size)

write.table(annotation, file=output_file, sep='\t', quotes=TRUE, row.names=FALSE)



