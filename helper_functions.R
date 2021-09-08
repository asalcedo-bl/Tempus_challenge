
process_inputs <- function(){
	# reads a vcf from a system file and processes it with the VariantAnnotation packahe
	# returns a list with a VCF object, a character vector of the outputfile, and a list of samples to analyze
	# vcf_file is a character vertor length (1) with the vcf file to read
	# output_file is a character vertor length (1) with the name of the output annotation file
	# genome is the genome version
	# chunk size is the number of records in the VCF to analyze at once. chunk_size of Inf (the default) will analyze the entire VCF at once. Users that are more constrained by memory may benefit from a smaller chunk_size. 
	parser <- ArgumentParser()
	parser$add_argument('--vcf_file', type = 'character', default = NULL)
	parser$add_argument('--output_file', type = 'character', default = NULL)
	parser$add_argument('--genome', type = 'character', default = "hg19")	
	parser$add_argument('--chunk_size', type = 'integer', default = NULL)	
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

	return(list(vcf=vcf, output_file=output_file, chunk_size=args$chunk_size))
}



run_annotation <- function(vcf, chunk_size=NULL, variant_effects=NULL){
	# annotates a user provided VCF object by breaking it up into user defined chunks
	# returns a data frame with the annotation for each variant. The dataframe contains the following annotations for each variant: the chromosome name, the start of the variant, its quality, the reference and alternate alleles, the number of alternate alleles, the number of reads supporting the reference and alternate alleles, the percent of reads supporting the alternate allele, the variant type, its predicted functional effect, and its frequency in ExAC
	# vcf is the VCF-class object (see VariantAnnotation) to annotate
	# chunk size is the number of records in the VCF to analyze at once. chunk_size of Inf (the default) will analyze the entire VCF at once. Users that are more constrained by memory may benefit from a smaller chunk_size. 	
	# variant effects is an optional data frame with the ranking of variant severity. It should have two columns, one titled "impact" which has the variant impact label and one labelled "effect" which gives an ascending numeric ranking of the variant severity (so the most benign variant impacts will have the smallest effect and the most severe variants will have the largest)

	n_variants <- length(vcf)
	if( !is.null(chunk_size)){
		#users may have different memory constraints depending on the file size and system 
		range_start <- seq(1,n_variants, by=chunk_size)
		range_end <- seq((chunk_size+1), n_variants,by=chunk_size)
		range_end <- c(range_end, n_variants)
		vcf_ranges <- cbind(range_start, range_end)
	}else{
		vcf_ranges <- matrix(c(1, n_variants), nrow=1, ncol=2)
	}
	
	variant_annotation <- adply(vcf_ranges, 1, .id='id',function(x) annotate_variants(vcf[x[1]:x[2]], variant_effects=variant_effects))
	variant_annotation$id <- NULL
	if(nrow(variant_annotation) != length(vcf)){
		warn("number of rows in annotation do not match the original VCF")
	}
	return(variant_annotation)
}


annotate_variants <- function(vcf, variant_effects=NULL){
	# annotates a user provided vcf with descriptive features
	# returns a dataframe with the following information for each variant: the chromosome name, the start of the variant, its quality, the reference and alternate alleles, the number of alternate alleles, the number of reads supporting the reference and alternate alleles, the percent of reads supporting the alternate allele, the variant type, its predicted functional effect, and its frequency in ExAC
	# vcf is the VCF object to annotate
	# variant effects is an optional data frame with the ranking of variant severity. It should have two columns, one titled "impact" which has the variant impact label and one labelled "effect" which gives an ascending numeric ranking of the variant severity (so the most benign variant impacts will have the smallest effect and the most severe variants will have the largest)
	
	alt_alleles <- alt(vcf)
	# there can be multiple alleles so collapse them into a single character value 
	alt_alleles_collapsed <- laply(alt_alleles, collapse_alleles)
	annotation <- data.frame(chr=as.character(seqnames(vcf)), start=start(vcf),  qual=qual(vcf), ref=as.character(ref(vcf)), alt=alt_alleles_collapsed)
	loc_id <- paste0(annotation$chr, "_", annotation$start)
	annotation$loc_id <- loc_id
	# extract depth information from the info fields
	vcf_info <- info(vcf)
	annotation$num_alt_alleles <- vcf_info$NUMALT
	annotation$total_depth <- vcf_info$DP
	annotation$reference_depth <- vcf_info$RO

	total_alt_dp <- vcf_info$AO
	# there can be multiple alleles so collapse depths into a single character value 
	annotation$n_alt_reads <- laply(total_alt_dp, function(x) ifelse(length(x) > 1, collapse_alleles(x), as.character(x)))

	alt_read_percent <- round(total_alt_dp/annotation$total_depth*100,4)
	annotation$alt_read_percent <- laply(alt_read_percent, function(x) ifelse(length(x) > 1, collapse_alleles(x), as.character(x)))
	annotation$variant_type <- laply(vcf_info$TYPE, function(x) ifelse(length(x) > 1, collapse_alleles(unlist(x)), as.character(unlist(x))))

	variant_effects <- extract_variant_effects(vcf, variant_effects=variant_effects)
	allele_freq <- annotate_ExAC_AF(annotation, vcf)

	annotation <- merge(annotation, allele_freq[,c("loc_id", "ExAC_AF")], all.x=TRUE, by="loc_id")
	annotation <- merge(annotation, variant_effects, all.x=TRUE, by="loc_id")
	# put the annotation in the same order as the original vcf
	annotation <- annotation[match(loc_id, annotation$loc_id),]
	annotation$loc_id <- NULL
	return(annotation)
}

collapse_alleles <- function(x, separator=","){
	# collapses a vector of alleles 
	# returns a character vector of length(1) with all alleles separated by a user defined separator
	# x is a character vector of alleles
	# separator is a character vector of length 1 with the character to be used to separate alternative alleles
	collapsed_alleles <- ifelse(length(x) == 1, as.character(x), paste(x, collapse=separator))

	return(collapsed_alleles)
}


default_variant_effect <- function(){
	# outputs the default variant effect dataframe
	return(data.frame(impact = c("coding","synonymous", "intergenenic", "intron", "threeUTR", "fiveUTR", "promoter",  "spliceSite", "nonsynonymous", "frameshift", "nonsense"),
										  effect = c(1:11)
										))
}

extract_variant_effects <- function(vcf, variant_effects=NULL){
	# annotates each variant in a VCF object with its functional effect based on its sequence position relative to genes
	# returns a two column data frame with the locus ID and the most deleterious predicted impact for that variant
	# vcf is the VCF object to annotate
	# variant effects is an optional data frame with the ranking of variant severity. It should have two columns, one titled "impact" which has the variant impact label and one labelled "effect" which gives an ascending numeric ranking of the variant severity (so the most benign variant impacts will have the smallest effect and the most severe variants will have the largest)
	
	# change chromosome names to be compatible with UCSC
	vcf_query <- keepStandardChromosomes(vcf)
	vcf_query <- dropSeqlevels(vcf_query,"MT")
	seqlevelsStyle(vcf_query) <- "UCSC"
	seqlevels(vcf_query)

	all_var <- locateVariants(vcf_query, TxDb.Hsapiens.UCSC.hg19.knownGene, AllVariants())
	coding_var <- predictCoding(vcf_query,TxDb.Hsapiens.UCSC.hg19.knownGene,BSgenome.Hsapiens.UCSC.hg19 )

	# extract the relevant columns from each annotation 
	all_var_df <- unique(data.frame(id=rownames(mcols(all_var)), impact=as.character(mcols(all_var)[,1]),stringsAsFactors=FALSE))
	rm(all_var)
	coding_var_df <-  unique(data.frame(id=rownames(mcols(coding_var)), impact=as.character(mcols(coding_var)[,"CONSEQUENCE"]),stringsAsFactors=FALSE))
	rm(coding_var)

	# if any variants are coding they will be more impactful
	all_var_df <- all_var_df[!(all_var_df$id %in% coding_var_df[coding_var_df$impact != "synonymous",]$id),]
	# generic "coding" annotations are less descriptive than those from predictCoding so remove them
	all_var_df <- all_var_df[!(all_var_df$impact == "coding" & all_var_df$id %in% coding_var_df$id),]

	all_var_df <- rbind(all_var_df, coding_var_df)
	# extract a single id for each locus
	all_var_df$loc_id <- gsub(":","_", gsub("_[ACTG].*","",all_var_df$id))

	if(is.null(variant_effects)){
		variant_effects <- default_variant_effect()
	}

	# for each locus output the more deleterious consequence as defined by the variant_effects dataframe
	all_var_df <- merge(all_var_df, variant_effects, by="impact", all.x=TRUE)
	variant_summary <- ddply(all_var_df, .(loc_id), function(x){
																	return(x[which.max(x$effect),])
																	})

	return(variant_summary[,c("loc_id", "impact")])
}


annotate_ExAC_AF <- function(annotation, vcf){
	# annotates each variant in a VCF object with its allele frequency in the ExAC database by querying the ExAC API. Requires a functioning internet connection.
	# returns a three column data frame with the locus ID, the variant ID in the ExAC data base, and the allele frequency for each variant in ExAC
	# annotation is a data frame with at least four columns: chr which is a character vector with the chromosome id, start which is a numeric vector of variant start positions, ref which is a character vector with the reference allele, and alt which is a character vectpr with the alternate allele
	# vcf is the VCF object to annotate

	# generate the basic ExAC query format necessary
	ExAC_query <- paste0(annotation$chr,"-", annotation$start,"-",annotation$ref,"-",annotation$alt)	
	# if a variant has more than one alternate allele separate each into a separate query
	alt_alleles <- alt(vcf)
	multi_allele <- sapply(grep(',',ExAC_query), function(i){
															site <- annotation[i,];
															paste0(site$chr,"-", site$start,"-",site$ref,"-", as.character(alt_alleles[[i]]))
															})
	multi_allele <- unlist(multi_allele)
	ExAC_query <- c(ExAC_query[grep(',',ExAC_query, invert=TRUE)], multi_allele)
	# query ExAC, format as json, and extract the allele frequencies
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
