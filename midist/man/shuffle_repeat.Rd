\name{shuffle_repeat}
\alias{shuffle_repeat}
\title{Repeated shuffle interaction test
}
\description{
Performs repeated null models for the comparison of intrachromosomal and interchromosomal mutual information distributions in gene regulatory networks or similar 
}
\usage{
shuffle_repeat(sif, index, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{sif}{
a SIF like object, node-interaction-node
}
  \item{index}{
an index obtained by chrom.index 
}
  \item{n}{
number of repeats
}
}

\value{
a list of length n of lists with class "htest" as returned by wilcox.test

}
