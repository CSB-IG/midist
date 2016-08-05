#!/usr/bin/Rscript
# This script takes a SIF file containting an interaction network
# and a chrom.index
#returns Mann-Whitney U test for intrachromosomic and interchromosomic tests for networks
#and a null model list

library("argparse")
library("data.table")
library("sif.subset.by.chrom")
library("midist")

parser <- ArgumentParser(description = "Split a sif file in intra/inter chromosomic interactions.")
parser$add_argument(
	"file",
	type = 'character',
	help = "The SIF file to process."
)

parser$add_argument(
	"index",
	type = 'character',
	help = "The chrom.index for the SIF file."
)

parser$add_argument(
	"--plots",
	action="store_true",
	help = "The chrom.index for the SIF file."
)

parser$add_argument(
	"-o",
	"--output",
	default = "./",
	help = "The directory where the output should go."
)

args <- parser$parse_args()

index <- load.index(args$index)
sif <- fread(args$file, data.table = FALSE)

#clean NAs
k <- known.chrom(index)
sif <- subset(sif, k)
index <- subset(index, k)
rm(k)

# intra and inter subgraphs

intra_test <- subset(sif, is.intra.chromosomic(index))
inter_test <- subset(sif, is.inter.chromosomic(index))

#probability density functions
intra_pdf<-pdf.sif(sif = intra_test)
inter_pdf<-pdf.sif(sif = inter_test)

#probability density functions for unseparated MI distribution
whole_pdf<-pdf.sif(sif = sif) 


#move to output directory
setwd(args$output)

#output prefix
prefix_fname = basename(args$file)

# Plotting


if(args$plots == TRUE){
  #plotting whole pdf for MI distribution
  plot_density.wo(sif, file = paste0("whole_", prefix_fname, "_plot.pdf")) #red 
  plot_density_log.wo(sif, file = paste0("whole_", prefix_fname, "_log_plot.pdf")) #red 
  plot_density_zoom.wo(sif, zoom = 10, file = paste0("whole_", prefix_fname,"_zoom_plot.pdf")) #red 
  #plotting intra and inter pdfs
  plot_density.wo(intra_pdf, inter_pdf, file = paste0(prefix_fname, "_plot.pdf")) #first one red, second blue
  plot_density_log.wo(intra_pdf, inter_pdf, file = paste0(prefix_fname, "_log_plot.pdf")) #first one red, second blue
  plot_density_zoom.wo(intra_pdf, inter_pdf, zoom = 10, file = paste0(prefix_fname,"_zoom_plot.pdf")) #first one red, second blue
}

#compare via U_test

prueba_U<-U_test(pdf1 = intra_pdf, pdf2 = inter_pdf)
#write out
u_out <- unlist(prueba_U)
write.csv(u_out, file = paste0(prefix_fname, ".stats"))

#null model
null_model <- shuffle_repeat(sif, index, n = 100)
write.csv(null_model, file = paste0(prefix_fname, ".nullmodel"))
