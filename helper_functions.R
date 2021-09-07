
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
	parser$add_argument('--chunk_size', type = 'integer', default = 1000)	
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

annotate_ExAC_AF <- function(basic_annotation, vcf, alt_alleles){

	# generate the basic ExAC query format necessary
	ExAC_query <- paste0(basic_annotation$chr,"-", basic_annotation$start,"-",basic_annotation$ref,"-",basic_annotation$alt)
	# if a variant has more than one alternate allele separate each into a separate query
	multi_allele <- sapply(grep(',',ExAC_query), function(i){
															site <- basic_annotation[i,];
															paste0(site$chr,"-", site$start,"-",site$ref,"-", as.character(alt_alleles[[i]]))
															})
	multi_allele <- unlist(multi_allele)
	ExAC_query <- c(ExAC_query[grep(',',ExAC_query, invert=TRUE)], multi_allele)
	# query ExAC and format as json
	ExAC_results <- httr::POST(url="http://exac.hms.harvard.edu/rest/bulk/variant", body=jsonlite::toJSON(as.character(ExAC_query)), encode = "json")
	ExAC_json <- content(ExAC_results)
	allele_freq <- ldply(ExAC_json, .id='exac_variant_id', function(x) ifelse(!is.null(x$variant$allele_freq), x$variant$allele_freq, NA))
	colnames(allele_freq)[2] <- "ExAC_AF"
	if(all(is.na(allele_freq$ExAC_AF))){
		stop("Something is wrong, all ExAC allele frequencies are missing")
	}
	
	allele_freq$loc_id <- gsub("-","_",gsub("-[ACTG]*-.*","",allele_freq$exac_variant_id))
	# reorder the data frame to reflect the order of the original query so that allele order will also be preserved
	allele_freq <- allele_freq[match(ExAC_query, allele_freq$exac_variant_id),]
	# consolidate allele frequencies for variants with more than a single alternate allele
	if(any(allele_freq$exac_variant_id %in% multi_allele)){
		multi_allele_freq <- ddply(allele_freq[which(allele_freq$exac_variant_id %in% multi_allele),], .(loc_id), function(x){
																															return(paste(x$ExAC_AF, collapse=','))
																															})
		colnames(multi_allele_freq)[2] <- "ExAC_AF"
		allele_freq <- rbind.fill(allele_freq[!(allele_freq$exac_variant_id %in% multi_allele),], multi_allele_freq)
	}
	if(nrow(allele_freq) != length(vcf)){
		stop("ExAC annotation length does not match vcf length")
	}
	return(allele_freq)

}

annotate_variants <- function(vcf, samples){
	alt_alleles <- alt(vcf)
	# there can be multiple alleles so collapse them into a single character value 
	alt_alleles_collapsed <- laply(alt_alleles, collapse_alleles)
	basic_annotation <- data.frame(chr=as.character(seqnames(vcf)), start=start(vcf), loc_id=rownames(info(vcf)), qual=qual(vcf), ref=as.character(ref(vcf)), alt=alt_alleles_collapsed, reference_dp=info(vcf)$RO)
	
	vcf_info <- info(vcf)
	total_depth <- vcf_info$DP
	total_alt_dp <- vcf_info$AO
	total_alt_dp_collapsed <- laply(total_alt_dp, function(x) ifelse(length(x) > 1, collapse_alleles(x), as.character(x)))

	total_alt_read_percent <- round(total_alt_dp/total_depth*100,2)
	total_alt_read_percent_collapsed <- laply(total_alt_read_percent, function(x) ifelse(length(x) > 1, collapse_alleles(x), as.character(x)))

	allele_freq <- annotate_ExAC_AF(basic_annotation, vcf, alt_alleles)
	
	return(data.frame(total_depth=total_depth, n_alt_reads=total_alt_dp_collapsed, total_alt_read_percent=total_alt_read_percent_collapsed))	
}


run_annotation <- function(vcf, samples, chunk_size=1000){
	
	n_variants <- length(vcf)
	if(chunk_size < Inf){
		#users may have different memory constraints depending on the file size and system 
		range_start <- seq(1,n_variants, by=chunk_size)
		range_end <- seq((chunk_size+1), n_variants,by=chunk_size)
		range_end <- c(range_end, n_variants)
		vcf_ranges <- cbind(range_start, range_end)
	}else{
		vcf_ranges <- matrix(c(1, n_variants), nrow=1, ncol=2)
	}
	# source("helper_functions.R")
	variant_annotation <- adply(vcf_ranges, 1, .id='id',function(x) annotate_variants(vcf[x[1]:x[2]]))
}

