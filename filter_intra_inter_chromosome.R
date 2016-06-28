#!/usr/bin/Rscript

# This script should either:
#
# - split a network gene interaction in intra chromosomic interactions and inter chromosomic interations, or
#
# - plot the distribution of MI in intra chromosomic interactions and inter chromosomic interactions.

library("biomaRt")

# TODO:
# - should we use a network as input for an interactive script?
# - should we receive a name for the graph as an input?

# FIXME: we did not loaded "genes" in this script,
#        this script should run by it's own.

#get genes in networks
genes<-unique(
	c(unique(sanos_sif$V1),
	unique(sanos_sif$V3),
	unique(enfermos_sif$V1),
	unique(enfermos_sif$V3))
)

mart=useMart(
  biomart="ENSEMBL_MART_ENSEMBL",
  dataset="hsapiens_gene_ensembl",
  host="www.ensembl.org"
)


# Make a list of genes in chromosome For each chromosome.
#
# This should be done outside of the in_same_chrom checker.
#
#  - If we built and saved the entire list from biomaRt,
#    we could use it as input and save bandwidth and processing.
#
#  - It is important however to be able to update the list at will.


genes_in_chrom <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

valid_chroms = c(1:23, "X", "Y")

genes_in_valid_chroms <- genes_in_chrom[genes_in_chrom$chromosome_name%in%valid_chroms,]

#offline mode dictionary
#write.table(genes_in_valid_chroms,
#           file = "gene_chromosome_dictionary_july2016.txt",
#           quote = FALSE, sep = "\t", row.names = FALSE)
#


# not_in_chrom<-genes[!(genes%in%genes_in_valid_chroms)]

# TODO: Add a function to make an index of intra chromosome relationships.

# We could then map the network interactions through the "in_same_chrom" function
# to get the vector of intra chromosomic conections.
# Inter chromosomic connections are the negation on that vector.

# this function gets chromosome name for given gene
gene.chrom<-function(gene_name){
	stopifnot(is.character(gene_name))
	# get the symbol name in the dictionary
	chrom = genes_in_valid_chroms$chromosome_name[genes_in_valid_chroms$hgnc_symbol == gene_name]
  if(length(chrom) < 1){
    return(NA)
  }
  else
	return(chrom)
}


# This function takes an index and finds whether both genes are in the same chromosome.
# check whether any of the genes are in some chrom
in.same.chrom <- function (A, B) {
	chrom_a = gene.chrom(A)
	chrom_b = gene.chrom(B)

	# TODO: Can we have an N/A or other error condition?
	if (is.na(chrom_a) || is.na(chrom_b)) {
		return(NA)
	}

	if (chrom_a == chrom_b) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

# takes sif. Returns TRUE for each row that has intrachromosomic interaction
# FALSE if interchromosomic interaction
is.intra.chromosomic <- function(sif){
  a = sif[,1]
  b = sif[,3]
  mapply(FUN = in.same.chrom, A = a, B = b, USE.NAMES = FALSE)
}

# Then we can pass this network through the mi_distribution_plots to get intra and extra.
#
# ¿Save the intra and extra networks?
# ¿Save the intra_chromosomic_index?
# ¿Repeat Ourselves (not DRY)?
