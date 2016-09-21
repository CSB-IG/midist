\name{pdf.sif}
\alias{pdf.sif}

\title{
probability density functions}
\description{
Probability density function for mutual information weights of gene regulatory networks}

\usage{
pdf.sif(sif)
}
\arguments{
  \item{sif}{a sif-like data frame with node-interaction-node format
}
}
\value{
A density-class object as returned by density
}


\examples{
function (sif)
{
    return(density(x = sif$V2))
  }
}
\keyword{pdf}
\keyword{pdf.sif}
