#!/usr/bin/Rscript

# This script takes a SIF file containting an interaction network
# and makes an index of chromosomes where the involved genes belong.

# Input is of the form:

# gene_1 <tab> M.I. value <tab> gene_2

# The output index is of the form:

# chrom_for_gene_1 <tab> chrom_for_gene_2

source("indexing.R")
library("argparse")
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

idx <- index.chromosome(sif)

idx_fname = gsub(".(txt|sif)", "_chrom.index", basename(args$file))

# Write index for intra chromosomic interactions
write.table(
	idx,
	file = idx_fname,
	quote = FALSE,
	sep = "\t",
	row.names = FALSE,
	col.names = FALSE
)
