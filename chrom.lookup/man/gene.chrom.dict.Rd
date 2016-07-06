\name{gene.chrom.dict}
\alias{gene.chrom.dict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a dictionary for gene name lookup from a SIF.
}
\description{
Given a vector of 'genes', returns a dictionary that maps each gene name to a chromosome.
}
\usage{
gene.chrom.dict(genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{
A list containing gene names (characters)
}
}
\details{
}
\value{
A 'hash' object containing a map for 'gene' -> 'chromosome'.
}
\references{
}
\author{
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{ \link{chrom.lookup} }
}
\examples{
example_genes <- c("A1BG", "MED27", "ANKRD49", "ARV1", "GEMIN6", "FAM86EP", "FRMPD4")

## The function is currently defined as
function (genes)
{
    mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
        host = "www.ensembl.org")
    genes_in_chrom <- getBM(attributes = c("hgnc_symbol", "chromosome_name"),
        filters = "hgnc_symbol", values = genes, mart = mart)
    valid_chroms = c(1:23, "X", "Y")
    genes_in_valid_chroms <- genes_in_chrom[genes_in_chrom$chromosome_name \%in\%
        valid_chroms, ]
    gene_dict = hash(genes_in_valid_chroms$hgnc_symbol, genes_in_valid_chroms$chromosome_name)
    return(gene_dict)
  }

dict <- gene.chrom.dict(example_genes)
}

\keyword{ chromosome }
\keyword{ gene }
