\name{shuffle_interactions_test}
\alias{shuffle_interactions_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
shuffle_interactions_test(sif, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sif}{
%%     ~~Describe \code{sif} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
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
function (sif, index)
{
    dex_shuffle <- sample(as.data.frame(index))
    intranw <- intra(sif = sif, index = dex_shuffle)
    internw <- inter(sif = sif, index = dex_shuffle)
    d.intra <- density(intranw$V2)
    d.inter <- density(internw$V2)
    return(U_test(d.intra, d.inter))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
