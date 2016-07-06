\name{gene.chrom.dict}
\alias{gene.chrom.dict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a dictionary for gene name lookup from a SIF.
}
\description{
Given a 'SIF', returns a dictionary that maps each gene name to a chromosome.
}
\usage{
gene.chrom.dict(sif)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sif}{
A dataframe containing: GENE1, MI value, GENE2.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A 'hash' object containing a map for 'gene' -> 'chromosome'.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (sif) 
{
    genes <- unique(c(unique(sif$V1), unique(sif$V3)))
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
