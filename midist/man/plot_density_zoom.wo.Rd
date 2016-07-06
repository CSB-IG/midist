\name{plot_density_zoom.wo}
\alias{plot_density_zoom.wo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_density_zoom.wo(..., zoom = 10, file = "plot.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{zoom}{
%%     ~~Describe \code{zoom} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
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
function (..., zoom = 10, file = "plot.pdf")
{
    dots = list(...)
    colores = rainbow(n = length(seq_along(dots)))
    zm = round(length(dots[[1]]$y)/zoom)
    pdf(file = file, paper = "a4r")
    plot(x = dots[[1]]$x[1:zm], y = dots[[1]]$y[1:zm], xlab = "MI",
        ylab = "Frequency", col = colores[1], type = "l", )
    for (i in seq_along(dots)) {
        if (i > 1) {
            points(x = dots[[i]]$x[1:zm], y = dots[[i]]$y[1:zm],
                type = "l", col = colores[i])
        }
    }
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
