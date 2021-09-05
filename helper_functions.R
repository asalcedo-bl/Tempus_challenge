
process_inputs <- function(){
	# reads a vcf from a system file and outputs a list with a VCF object and a character vector of the outputfile
	# vcf_file is a character vertor length (1) with the vcf file to read
	# output_file is a character vertor length (1) with the name of the output annotation file
	# genome is the genome version
	# samples is a character vector with the samples for which to output annotations
	parser <- ArgumentParser()
	parser$add_argument('--vcf_file', type = 'character', default = NULL)
	parser$add_argument('--output_file', type = 'character', default = NULL)
	parser$add_argument('--genome', type = 'character', default = "hg19")
	parser$add_argument('--samples', type = 'character', default = "all")	
	args <- parser$parse_args()

	vcf <- 	tryCatch(
				{
					readVcf(file=args$vcf_file, genome=args$genome)
				},
				error =  function(cond) {
						    		print("unable to read VCF file")
						    		msg <- conditionMessage(cond)
						            stop(msg)         
									return(NULL)
						},
				finally = 1)

	if(is.null(args$output_file)){
	# if there is no output file generate one based on the vcf name	 
		output_file <- gsub(".vcf", "_annotation.txt", args$vcf_file)
	}else{
		output_file <- args$output_file
	}
	
	# gather samples to output
	all_samples <- samples(header(vcf))
	if( 'all' == args$samples ){
		samples <- all_samples
	}else{
	# check that all of the requested samples are in the vcf
		a_ply(args$samples, 1, function(x) if(!(x %in% all_samples)) print( paste(x, "not found in vcf file")))
		samples <- all_samples[all_samples %in% args$samples]
	}

	return(list(vcf=vcf, output_file=output_file, samples=samples))
}


