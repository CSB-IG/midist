\name{sif.subset.by.chrom-package}
\alias{sif.subset.by.chrom-package}
\alias{sif.subset.by.chrom}

\docType{package}

\title{ \packageTitle{sif.subset.by.chrom} }

\description{ \packageDescription{sif.subset.by.chrom} }

\details{

The DESCRIPTION file:
\packageDESCRIPTION{sif.subset.by.chrom}
\packageIndices{sif.subset.by.chrom}

You should have an index made by \code{ \link{chrom.lookup} } for your SIF file.

When you have done this, you can load both the index and the sif to subset it.

}

\keyword{ package }

\seealso{
	\code{\link{is.in.chrom}}
	\code{\link{is.in.chrom.pair}}
	\code{\link{is.in.chroms}}
	\code{\link{is.inter.chromosomic}}
	\code{\link{is.intra.chromosomic}}
	\code{\link{known.chrom}}
	\code{\link{load.index}}
}

\examples{
	#library("data.table")
	#library("sif.subset.by.chrom")
	#sif <- fread("/path/to/sif", data.table = FALSE)
	#sif_idx <- load.index("/path/to/index")
	#intra_chromosomic_sif <- subset(sif, is.intra.chromosomic(sif_idx))
}
