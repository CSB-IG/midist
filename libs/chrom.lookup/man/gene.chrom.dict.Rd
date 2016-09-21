\name{gene.chrom.dict}
\alias{gene.chrom.dict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a dictionary for gene name lookup from a SIF.
}
\description{
Given a bioMart object or data.frame containing the fields hgnc_symbol and chromosome_name.

}
\usage{
gene.chrom.dict(genes_in_chrom)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes_in_chrome}{
bioMart object or data.frame containing the fields hgnc_symbol and chromosome_name.
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
function (genes_in_chrom)
{
    valid_chroms = c(1:23, "X", "Y")
    genes_in_valid_chroms <- genes_in_chrom[genes_in_chrom$chromosome_name \%in\%
        valid_chroms, ]
    gene_dict = hash(genes_in_valid_chroms$hgnc_symbol, genes_in_valid_chroms$chromosome_name)
    return(gene_dict)
  }

dict <- gene.chrom.dict(chrom.info.biomart(example_genes))
}

\keyword{ chromosome }
\keyword{ gene }
