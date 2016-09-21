# this function generates a dictionary of gene, chromosome
# for lookup using gene.chrom
#
# genes_in_chrom: data frame containing fields
#       hgnc_symbol, chromosome_name
#
gene.chrom.dict <- function(genes_in_chrom) {

	# filter only valid gene/chromosome pairs
	genes_in_valid_chroms = valid.chrom.info(genes_in_chrom)

	# convert to dictionary
	gene_dict = hash::hash(genes_in_valid_chroms$hgnc_symbol, genes_in_valid_chroms$chromosome_name)
	return(gene_dict)
}

# this function returns only the info where chromosome is simple
valid.chrom.info <- function(chrom_info) {
	valid_chroms = c(1:23, "X", "Y")
	valid_chrom_info <- chrom_info[chrom_info$chromosome_name %in% valid_chroms,]
	return(valid_chrom_info)
}

# Get chromosome info for given genes from bioMart.
#
# genes: vector of gene names (string).
#
chrom.info.biomart <- function(genes) {
	# find the genes on biomart
	mart = biomaRt::useMart(
		biomart="ENSEMBL_MART_ENSEMBL",
		dataset="hsapiens_gene_ensembl",
		host="www.ensembl.org"
	)
	genes_in_chrom <- biomaRt::getBM(
		attributes = c("hgnc_symbol", "chromosome_name"),
		filters = "hgnc_symbol",
		values = genes,
		mart = mart
	)
	return(genes_in_chrom)
}

# this function gets chromosome name for given gene
gene.chrom <- function(gene_name, gene_dict) {
	if (hash::has.key(gene_name, gene_dict)) { return(gene_dict[[gene_name]]) }
	else { return(NA) }
}

# Takes a list of gene names and a dict containing gene, chromosome.
#
# Returns a vector containing a chromosome name for each gene name on genes.
index.chromosome <- function(genes,
chrom_dict = gene.chrom.dict(chrom.info.biomart(genes)),
mc.cores = 1) {
	chrom_index <- parallel::mcmapply(
		FUN=gene.chrom,
		genes,
		MoreArgs = list(chrom_dict),
		USE.NAMES=FALSE,
		mc.cores=mc.cores
	)
	return(chrom_index)
}
