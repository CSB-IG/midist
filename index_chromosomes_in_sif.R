#!/usr/bin/Rscript

# This script takes a SIF file containting an interaction network
# and makes an index of chromosomes where the involved genes belong.

# Input is of the form:

# gene_1 <tab> M.I. value <tab> gene_2

# The output index is of the form:

# chrom_for_gene_1 <tab> chrom_for_gene_2

library("argparse")
library("chrom.lookup")
library("data.table")
library("hash")
# MAIN

parser <- ArgumentParser(description = "Split a sif file in intra/inter chromosomic interactions.")
parser$add_argument(
	"file",
	metavar = 'file',
	type = 'character',
	help = "The SIF file to process."
)
parser$add_argument(
	"-d",
	"--dictionary",
	help = "A tab delimited file containing bioMart fields: hgnc_symbol, chromosome_name."
)
parser$add_argument(
	"-o",
	"--output",
	help = "The name of the file where to write the index."
)
args <- parser$parse_args()

if(is.null(args$output)) {
	idx_fname = gsub(".(txt|sif)", "_chrom.index", basename(args$file))
} else {
	idx_fname = args$output
}

sif <- fread(input = args$file, data.table = FALSE)

genes <- unique(c(unique(sif$V1), unique(sif$V3)))

if(is.null(args$dictionary)) {
	chrom_info = chrom.info.biomart(genes)
} else {
	chrom_info = fread(input = args$dict,
		colClasses=c("character", "character"),
		data.table=FALSE);
}

dict <- gene.chrom.dict(chrom_info)

i1 <- index.chromosome(sif$V1, dict)
i2 <- index.chromosome(sif$V3, dict)
idx <- cbind(i1, i2)

# Write index for intra chromosomic interactions
write.table(
	idx,
	file = idx_fname,
	quote = FALSE,
	sep = "\t",
	row.names = FALSE,
	col.names = FALSE
)
