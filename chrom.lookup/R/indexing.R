# this function generates a dictionary of gene, chromosome
# for lookup using gene.chrom
gene.chrom.dict <- function(genes) {
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

	# filter only valid gene/chromosome pairs
	valid_chroms = c(1:23, "X", "Y")
	genes_in_valid_chroms <- genes_in_chrom[genes_in_chrom$chromosome_name%in%valid_chroms,]

	# convert to dictionary
	gene_dict = hash::hash(genes_in_valid_chroms$hgnc_symbol, genes_in_valid_chroms$chromosome_name)
	return(gene_dict)
}

# this function gets chromosome name for given gene
gene.chrom <- function(gene_name, gene_dict) {
	if (hash::has.key(gene_name, gene_dict)) { return(gene_dict[[gene_name]]) }
	else { return(NA) }
}

# Takes a list of gene names and a dict containing gene, chromosome.
#
# Returns a vector containing a chromosome name for each gene name on genes.
index.chromosome <- function(genes, chrom_dict = gene.chrom.dict(genes)) {
	chrom_index <- parallel::mcmapply(
		FUN=gene.chrom,
		genes,
		MoreArgs = list(chrom_dict),
		USE.NAMES=FALSE,
		mc.cores=16
	)
	return(chrom_index)
}
