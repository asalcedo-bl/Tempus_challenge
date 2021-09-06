
process_inputs <- function(){
	# reads a vcf from a system file and outputs a list with a VCF object, a character vector of the outputfile, and a list of samples to analyze
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
						    		print("Unable to read VCF file")
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
		a_ply(args$samples, 1, function(x) if(!(x %in% all_samples)) warning( paste(x, "not found in vcf file")))
		samples <- all_samples[all_samples %in% args$samples]
	}

	if(length(samples) == 0){
		stop("No samples specified. Specify valid sample names for which to annotate variants or 'all' for to output annotation for all samples")
	}

	return(list(vcf=vcf, output_file=output_file, samples=samples))
}

collapse_alleles <- function(x, separator=","){
	# collapses a vector of alleles and outputs a character vector of length(1) with all alleles separated by a user defined separator
	# x is a character vector of alleles
	# separator is a character vector of length 1 with the character to be used to separate alternative alleles
	collapsed_alleles <- ifelse(length(x) == 1, as.character(x), paste(x, collapse=separator))

	return(collapsed_alleles)
}


annotate_variants <- function(vcf, samples){
	vcf_info <- info(vcf)
	total_depth <- vcf_info$DP
	total_alt_dp <- vcf_info$AO
	total_alt_dp_collapsed <- laply(total_alt_dp, function(x) ifelse(length(x) > 1, collapse_alleles(x), as.character(x)))

	total_alt_read_percent <- round(total_alt_dp/total_depth*100,2)
	total_alt_read_percent_collapsed <- laply(total_alt_read_percent, function(x) ifelse(length(x) > 1, collapse_alleles(x), as.character(x)))

	return(data.frame(total_depth=total_depth, n_alt_reads=total_alt_dp_collapsed, total_alt_read_percent=total_alt_read_percent_collapsed))	
}


run_annotation <- function(vcf, samples, chunk_size=1000){
	alt_alleles <- laply(alt(vcf), collapse_alleles)
	basic_annotation <- data.frame(chr=as.character(seqnames(vcf)), start=start(vcf), id=rownames(info(vcf)), qual=qual(vcf), ref=as.character(ref(vcf)), alt=alt_alleles, reference_dp=info(vcf)$RO)
	

	range_start <- seq(1,length(vcf), by=chunk_size)
	range_end <- seq((chunk_size+1), length(vcf),by=chunk_size)
	range_end <- c(range_end, length(vcf))
	vcf_ranges <- cbind(range_start, range_end)

	# source("helper_functions.R")
	variant_annotation <- adply(vcf_ranges, 1, .id='id',function(x) annotate_variants(vcf[x[1]:x[2]]))
}

