\name{gene.chrom.dict}
\alias{gene.chrom.dict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selects from a chromomose list where the chromosomes are a simple name.
}
\description{
Given data.frame containing the fields hgnc_symbol and chromosome_name (as queried from biomarRt),
it returns a subset where chromosomes are defined as "1" ... "23", "X", "Y".

}
\usage{
gene.chrom.dict(genes_in_chrom)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes_in_chrome}{
data.frame containing the fields hgnc_symbol and chromosome_name (as queried from biomaRt).
}
}
\details{
}
\value{
data.frame containing the fields hgnc_symbol and chromosome_name
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
## The function is currently defined as
valid.chrom.info <- function(chrom_info) {
       valid_chroms = c(1:23, "X", "Y")
       valid_chrom_info <- chrom_info[chrom_info$chromosome_name %in% valid_chroms,]
       return valid_chrom_info
}
}

\keyword{ chromosome }
\keyword{ gene }
