library("data.table")
library("argparse")

parser <- ArgumentParser(description = "Plots the MI distribution of two SIF files.")
parser$add_argument(
	"sanos",
	metavar = 'sanos',
	type = 'character',
	help = "The SIF file to process."
)
parser$add_argument(
	"enfermos",
	metavar = 'enfermos',
	type = 'character',
	help = "The SIF file to process."
)
args <- parser$parse_args()

sanos_sif <-fread(input = args$sanos, data.table = FALSE)
enfermos_sif <-fread(input = args$enfermos, data.table = FALSE)

#histogramas
hist_sanos <- hist(sanos_sif$V2, breaks = 200, plot = FALSE)
hist_enfermos <- hist(enfermos_sif$V2, breaks = 200, plot = FALSE)

#TODO: Check whether this is the right way to estimate the distribution.
#density plots
density_sanos<-density(sanos_sif$V2)
density_enfermos<-density(enfermos_sif$V2)

##plot density
plot_fname = paste0(
	gsub(".(txt|sif)", "_blue_", basename(args$sanos)),
	gsub(".(txt|sif)", "_red", basename(args$enfermos)),
	".png"
)
png(plot_fname, width=4, height=4, units="in", res=300)
plot(x = density_enfermos, type = "l", col = "red")
points(x = density_sanos, type = "l", col = "blue")
dev.off()

#compare densities?
ks.test(x = density_sanos$y, y = density_enfermos$y)
