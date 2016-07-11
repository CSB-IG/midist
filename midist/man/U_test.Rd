\name{U_test}
\alias{U_test}

\title{
Mann-Whitney U test
}
\description{
A wrapper for the MAnn-Whitney U test for the comparison of intrachromosomal and interchromosomal mutual information distributions in gene regulatory networks or similar
}
\usage{
U_test(pdf1, pdf2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdf1}{
Objects of class density as produced by density()
}
  \item{pdf2}{
Objects of class density as produced by density()
}
}

\value{
A list with class "htest" as returned by wilcox.test

}
