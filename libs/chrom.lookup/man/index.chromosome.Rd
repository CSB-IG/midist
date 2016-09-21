\name{index.chromosome}
\alias{index.chromosome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to index the chromosomes for each SIF element.
}
\description{
Given a vector of gene names, and optionally a chrom_dict,
it returns a vector containing the chromosome name
in which the gene is located.
}
\usage{
index.chromosome(genes)
d <- gene.chrom.dict(genes-superset)
index.chromosome(genes, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{
A vector containing gene names.
}
  \item{chrom_dict}{
The dictionary returned from \code{gene.chrom.dict}
}
}
\details{
As this function is used for datasets sets containing millions of gene instances, it is multi threaded.
}
\value{
A vector of chromosome names corresponding to each element on genes.
}
\references{
}
\author{
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (genes, chrom_dict = gene.chrom.dict(genes))
{
	chrom_index <- parallel::mcmapply(
		FUN=gene.chrom,
		genes,
		MoreArgs = list(chrom_dict),
		USE.NAMES=FALSE,
		mc.cores=16
	)
	return(chrom_index)

}
}
