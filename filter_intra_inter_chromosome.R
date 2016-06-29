#!/usr/bin/Rscript

# This script should either:
#
# - split a network gene interaction in intra chromosomic interactions and inter chromosomic interations, or
#
# - plot the distribution of MI in intra chromosomic interactions and inter chromosomic interactions.

library("argparse")
library("biomaRt")
library("data.table")
library("hash")
library("parallel")

# this function gets chromosome name for given gene
gene.chrom<-function(gene_name, gene_dict){
	if(has.key(gene_name, gene_dict)) { return(gene_dict[[gene_name]]) }
	else { return(NA) }
}

gene.chrom.dict <- function(sif) {
	# get all the genes in the sif file
	genes<-unique(c(unique(sif$V1), unique(sif$V3)))

	# find the genes on biomart
	mart=useMart(
	  biomart="ENSEMBL_MART_ENSEMBL",
	  dataset="hsapiens_gene_ensembl",
	  host="www.ensembl.org"
	)
	genes_in_chrom <- getBM(
	  attributes = c("hgnc_symbol", "chromosome_name"),
	  filters = "hgnc_symbol",
	  values = genes,
	  mart = mart
	)

	# filter only valid gene/chromosome pairs
	valid_chroms = c(1:23, "X", "Y")
	genes_in_valid_chroms <- genes_in_chrom[genes_in_chrom$chromosome_name%in%valid_chroms,]
	gene_dict = hash(genes_in_valid_chroms$hgnc_symbol, genes_in_valid_chroms$chromosome_name)

	return(gene_dict)
}

# This function takes an index and finds whether both genes are in the same chromosome.
# check whether any of the genes are in some chrom
in.same.chrom <- function (A, B, dict) {
	chrom_a = gene.chrom(A, dict)
	chrom_b = gene.chrom(B, dict)

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
index.intra.chromosomic <- function(sif, dict){

	a = sif[,1]
	b = sif[,3]
	mcmapply(FUN = in.same.chrom, A = a, B = b, USE.NAMES = FALSE, MoreArgs = list(dict), mc.cores = 16)
}

# MAIN

parser <- ArgumentParser(description = "Split a sif file in intra/inter chromosomic interactions.")
parser$add_argument(
	"file",
	metavar = 'file',
	type = 'character',
	help = "The SIF file to process."
)
args <- parser$parse_args()

sif <- fread(input = args$file, data.table = FALSE)

#offline mode dictionary
#write.table(genes_in_valid_chroms,
#           file = "gene_chromosome_dictionary_july2016.txt",
#           quote = FALSE, sep = "\t", row.names = FALSE)
#


# not_in_chrom<-genes[!(genes%in%genes_in_valid_chroms)]

chrom_dict = gene.chrom.dict(sif)
index_intra <- index.intra.chromosomic(sif, chrom_dict)
intra_sif <- sif[index_intra,]
inter_sif <- sif[!index_intra,]

intra_fname = gsub(".(txt|sif)", "_intra.sif", basename(args$file))
inter_fname = gsub(".(txt|sif)", "_inter.sif", basename(args$file))

write.table(
	intra_sif,
	file = intra_fname,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
)

write.table(
	inter_sif,
	file = inter_fname,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
)
