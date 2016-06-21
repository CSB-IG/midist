library("biomaRt")
#get genes in networks
genes<-unique(c(unique(sanos_sif$V1), unique(sanos_sif$V3), unique(enfermos_sif$V1), unique(enfermos_sif$V3)))
#Timing this code
# ptm <- proc.time()
# genes<-unique(c(unique(sanos_sif$V1), unique(sanos_sif$V3), unique(enfermos_sif$V1), unique(enfermos_sif$V3)))
# proc.time() - ptm

##SLOWER CODE
# ptm <- proc.time()
# genes<-unique(c(sanos_sif$V1, sanos_sif$V3, enfermos_sif$V1, enfermos_sif$V3))
# proc.time() - ptm

#get BiomaRt
listMarts(host="www.ensembl.org")

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

valid_choms = c(1:23, "X", "Y")

genes_in_valid_chroms <- dict_chrom_gene[dict_chrom_gene$chromosome_name%in%valid_choms,]

# not_in_chrom<-genes[!(genes%in%genes_in_valid_chroms)]

