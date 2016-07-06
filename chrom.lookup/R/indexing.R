# this function generates a dictionary of gene, chromosome
# for lookup using gene.chrom
gene.chrom.dict <- function(sif) {
	# get all the genes in the sif file
	genes<-unique(c(unique(sif$V1), unique(sif$V3)))

	# find the genes on biomart
	mart = biomaRt::useMart(
		biomart="ENSEMBL_MART_ENSEMBL",
		dataset="hsapiens_gene_ensembl",
		host="www.ensembl.org"
	)
	genes_in_chrom <- biomart::getBM(
		attributes = c("hgnc_symbol", "chromosome_name"),
		filters = "hgnc_symbol",
		values = genes,
		mart = mart
	)

	# filter only valid gene/chromosome pairs
	valid_chroms = c(1:23, "X", "Y")
	genes_in_valid_chroms <- genes_in_chrom[genes_in_chrom$chromosome_name%in%valid_chroms,]
	gene_dict = hash::hash(genes_in_valid_chroms$hgnc_symbol, genes_in_valid_chroms$chromosome_name)

	return(gene_dict)
}

# this function gets chromosome name for given gene
gene.chrom<-function(gene_name, gene_dict){
	if(has.key(gene_name, gene_dict)) { return(gene_dict[[gene_name]]) }
	else { return(NA) }
}

# Takes a sif and a dict containing gene, chromosome.
#
# Returns an index of crom_a, chrom_b in SIF.
index.chromosome <- function(sif, chrom_dict=gene.chrom.dict(sif)){
	return(
		as.matrix(
			cbind(
				parallel::mcmapply(FUN=gene.chrom, sif$V1, MoreArgs = list(chrom_dict), USE.NAMES=FALSE, mc.cores=16),
				parallel::mcmapply(FUN=gene.chrom, sif$V3, MoreArgs = list(chrom_dict), USE.NAMES=FALSE, mc.cores=16)
			),
			dimnames = NULL
		)
	)
}
