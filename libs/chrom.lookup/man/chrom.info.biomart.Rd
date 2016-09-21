\name{chrom.info.biomart}
\alias{chrom.info.biomart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get chromosome info for given genes from bioMart.
}
\description{
Given a vector of 'genes', returns a dictionary that maps each gene name to a chromosome.
}
\usage{
chrom.info.biomart(genes)
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
A 'bioMart' object containing 'hgnc_symbol', 'chromosome_name'.
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
    return(genes_in_chrom)
  }

genes_in_chrom <- chrom.info.biomart(example_genes)
}

\keyword{ chromosome }
\keyword{ gene }
