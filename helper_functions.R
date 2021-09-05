
process_inputs <- function(){
	# reads a vcf from a system file and outputs a list with a VCF object and a character vector of the outputfile
	# vcf_file is a character vertor length (1) with the vcf file to read
	# output_file is a character vertor length (1) with the name of the output annotation file
	# genome is the genome version
	parser <- ArgumentParser()
	parser$add_argument('--vcf_file', type = 'character', default = NULL)
	parser$add_argument('--output_file', type = 'character', default = NULL)
	parser$add_argument('--genome', type = 'character', default = "hg19")
	args <- parser$parse_args()

	vcf <- readVcf(file=args$vcf_file, genome=args$genome)

	if(is.null(args$output_file)){
		output_file <- gsub(".vcf", "_annotation.txt", args$vcf_file)
	}else{
		output_file <- args$output_file
	}
	return(list(vcf=vcf, output_file=output_file))
}


