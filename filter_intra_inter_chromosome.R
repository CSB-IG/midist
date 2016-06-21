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

# not_in_chrom<-genes[!(genes%in%genes_in_valid_chroms)]

# TODO: Add a function to make an index of intra chromosome relationships.

# We could then map the network interactions through the "in_same_chrom" function
# to get the vector of intra chromosomic conections.
# Inter chromosomic connections are the negation on that vector.

# This function takes an index and finds whether both genes are in the same chromosome.
function in_same_chrom(A, B) {
	# check whether any of the genes are in some chrom
	for c in valid_chroms:
		# test A in chrom c
		a_in_chrom = ??
		# test B in chrom c
		b_in_chrom = ??
		# if a_in_chrom and b_in_chrom, return true (both genes in the same chrom)
		# if a_in_chrom or b_in_chrom, return false (we found one and the other is not in the same chrom)
		# else, keep looking

	# We did not found any of the genes in the chroms
	# ¿should we return NA or raise an error?
}

# Map genes_in_valid_chroms trough the in_same_chrom function
intra_chromosomic_index <- ?map(genes_in_valid_chroms, in_same_chrom(x$V1, x$V3)?

# Then we can pass this network through the mi_distribution_plots to get intra and extra.
#
# ¿Save the intra and extra networks?
# ¿Save the intra_chromosomic_index?
# ¿Repeat Ourselves (not DRY)?
