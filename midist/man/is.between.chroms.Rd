\name{is.between.chroms}
\alias{is.between.chroms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
is.between.chroms(chrom_index, c1, c2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrom_index}{
%%     ~~Describe \code{chrom_index} here~~
}
  \item{c1}{
%%     ~~Describe \code{c1} here~~
}
  \item{c2}{
%%     ~~Describe \code{c2} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
function (chrom_index, c1, c2)
{
    return(known.chrom(chrom_index) & ((chrom_index$V1 == c1 &
        chrom_index$V2 == c2) | (chrom_index$V1 == c2 & chrom_index$V2 ==
        c1)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
