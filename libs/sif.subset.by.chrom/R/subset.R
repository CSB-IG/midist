# Reads the filename assuming it correspond to a chromosome index for your SIF.
#
# Returns a data frame containing the chromosome names for the genes in your SIF.
load.index <- function(filename){
	idx = data.table::fread(input = filename,
		data.table = FALSE,
		colClasses = c("character", "character"))
	idx$V1 = as.character(idx$V1)
	idx$V2 = as.character(idx$V2)
	return(as.data.frame(idx))
}

#
#Subgraphing functions
#
known.chrom <- function(chrom_index) {
	return(!is.na(chrom_index$V1) & !is.na(chrom_index$V2))
}

is.intra.chromosomic <- function(chrom_index){
	return(
		known.chrom(chrom_index)
		& (chrom_index$V1 == chrom_index$V2)
	)
}

is.inter.chromosomic <- function(chrom_index){
	return(
		known.chrom(chrom_index)
		& chrom_index$V1 != chrom_index$V2
	)
}

is.in.chrom <- function(chrom_index, chrom_name) {
	return(
		known.chrom(chrom_index)
		& (chrom_index$V1 == chrom_name | chrom_index$V2 == chrom_name)
	)
}

is.in.chroms <- function(chrom_index, ...) {
	return(
		known.chrom(chrom_index)
		& (chrom_index$V1 %in% c(...) | chrom_index$V2 %in% c(...))
	)
}

is.in.chrom.pair <- function(chrom_index, c1, c2) {
	return(
		known.chrom(chrom_index)
		& ((chrom_index$V1 == c1 & chrom_index$V2 == c2)
		| (chrom_index$V1 == c2 & chrom_index$V2 == c1))
	)
}
