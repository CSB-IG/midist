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

prefix = basename(args$file)

index <- load.index(args$index)
sif <- fread(args$file, data.table = FALSE)

# clean NAs
k <- known.chrom(index)
sif <- subset(sif, k)
index <- subset(index, k)
rm(k)

# move to output directory
setwd(args$output)

# Plotting
if(args$plots == TRUE){
	#probability density functions for unseparated MI distribution
	intra_pdf <- sif.density(subset(sif, is.intra.chromosomic(index)))
	inter_pdf <- sif.density(subset(sif, is.inter.chromosomic(index)))
	whole_pdf <- sif.density(sif)

	#plotting whole pdf for MI distribution
	plot_density.wo(whole_pdf, file = paste0("whole_", prefix, "_plot.pdf")) #red
	plot_density_log.wo(whole_pdf, file = paste0("whole_", prefix, "_log_plot.pdf")) #red
	plot_density_zoom.wo(whole_pdf, zoom = 10, file = paste0("whole_", prefix,"_zoom_plot.pdf")) #red
	#plotting intra and inter pdfs
	plot_density.wo(intra_pdf, inter_pdf, file = paste0(prefix, "_plot.pdf")) #first one red, second blue
	plot_density_log.wo(intra_pdf, inter_pdf, file = paste0(prefix, "_log_plot.pdf")) #first one red, second blue
	plot_density_zoom.wo(intra_pdf, inter_pdf, zoom = 10, file = paste0(prefix,"_zoom_plot.pdf")) #first one red, second blue
}

# compare distributions using wilcox test
test.intra.inter <- function(sif, index) {
	# Probability Density Functions (PDF)
	intra_pdf <- sif.density(subset(sif, is.intra.chromosomic(index)))
	inter_pdf <- sif.density(subset(sif, is.inter.chromosomic(index)))

	result <- wilcox.test(
		x = intra_pdf$y,
		y = inter_pdf$y,
		alternative = "two.sided",
		paired = FALSE,
		exact = FALSE
	)
	return(result)
}

# output filenames
result_fname <- paste0(prefix, ".stats")
null_model_fname <- paste0(prefix, ".nullmodel")

result <- test.intra.inter(sif, index)

results_header <- function(wilcox_test_result) {
	value_names = names(wilcox_test_result)
	return(c(
		paste(
			"# ",
			wilcox_test_result$method,
			wilcox_test_result$alternative,
			wilcox_test_result$data.name,
			sep = " "
		),
		paste(
			value_names[1],
			value_names[3],
			value_names[2],
			sep = "\t"
		)
	))
}

results_body <-function(wilcox_test_result) {
	return(paste(
		result$statistic,
		result$p.value,
		result$parameter,
		sep = "\t"
	))
}

write(results_header(result), file = result_fname)
write(results_body(result), file = result_fname, append = TRUE)

## Null model
s <- index
write(results_header(result), file = null_model_fname) # reusing older header because it's just the same
for (i in 1:1000) {
	s$V1 <- sample(s$V1)
	s$V2 <- sample(s$V2)
	result = test.intra.inter(sif, s)
	write(results_body(result), file = null_model_fname, append = TRUE)
}
