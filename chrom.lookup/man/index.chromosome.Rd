\name{index.chromosome}
\alias{index.chromosome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to index the chromosomes for each SIF element.
}
\description{
Given a SIF and optionally a chrom_dict, it returns a matrix containing the chromosome names.
}
\usage{
index.chromosome(sif, chrom_dict = gene.chrom.dict(sif))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sif}{
A 'data.frame' containing: 'GENE1', 'MI value', 'GENE 2'
}
  \item{chrom_dict}{
The dictionary returned from \code{gene.chrom.dict}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix of the form "chrom1", "chrom2" that contains a row for each interaction on the input SIF.
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
function (sif, chrom_dict = gene.chrom.dict(sif)) 
{
    return(as.matrix(cbind(mcmapply(FUN = gene.chrom, sif$V1, 
        MoreArgs = list(chrom_dict), USE.NAMES = FALSE, mc.cores = 16), 
        mcmapply(FUN = gene.chrom, sif$V3, MoreArgs = list(chrom_dict), 
            USE.NAMES = FALSE, mc.cores = 16)), dimnames = NULL))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
