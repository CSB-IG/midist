source("midist.R")
source("indexing.R")

#read sif file
test <- fread(input = "test/test_case.txt", data.table = FALSE)
#setWD to tests
setwd("/labs/chrom_mi_dist/test/tests_scripts")

#index gene interactions into chromosomes
index <- index.chromosome(test)

##write it out
write.table(index,
            file = "test.index",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE
            )
#IF INDEX ALREADY CALCULATED
#index <- load.index(INDEX)


# intra and inter subgraphs

intra_test <- selector(sif, is.intra.chromosomic(index))
inter_test <- selector(sif, is.inter.chromosomic(index))

##write out
write.table(intra_test,
            file = "intra_test.sif",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE
)

write.table(inter_test,
            file = "inter_test.sif",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE
)

#probability density functions
intra_pdf<-pdf.sif(sif = intra_test)
inter_pdf<-pdf.sif(sif = inter_test)

# Plotting
plot_density.wo(intra_pdf, inter_pdf, file = "plot.pdf") #first one red, second blue
plot_density_log.wo(intra_pdf, inter_pdf, file = "log_plot.pdf") #first one red, second blue
plot_density_zoom.wo(intra_pdf, inter_pdf, zoom = 10, file = "log_plot.pdf") #first one red, second blue

#compare via U_test

prueba_U<-U_test(pdf1 = intra_pdf, pdf2 = inter_pdf)
write.table(cbind(prueba_U$statistic, prueba_U$p.value),
            file = "U_test_prueba.txt")
