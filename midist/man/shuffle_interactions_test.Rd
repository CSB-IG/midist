\name{shuffle_interactions_test}
\alias{shuffle_interactions_test}

\title{
Shuffle interactions test
}
\description{
A null model for the comparison of intrachromosomal and interchromosomal mutual information distributions in gene regulatory networks or similar 
}
\usage{
shuffle_interactions_test(sif, index)
}

\arguments{
  \item{sif}{
a SIF like object, node-interaction-node
}
  \item{index}{
an index obtained by chrom.index 
}
}
\details{
This function shuffles the intra and inter chromosomal labels for each interaction, then performs a Mann-Whitney test for subgraphs obtained with the new labels
}
\value{
A list with class "htest" as returned by wilcox.test

}
