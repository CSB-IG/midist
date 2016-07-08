\name{sif.subsets}
\alias{ is.in.chrom }
\alias{ is.in.chrom.pair }
\alias{ is.in.chroms }
\alias{ is.inter.chromosomic }
\alias{ is.intra.chromosomic }
\alias{ known.chrom }

\title{ Default subsets for \code{\link{sif.subset.by.chrom}} }

\description{
Functions for subsetting a SIF based on chromosome information.
}

\usage{
	is.in.chrom(chrom_index, chrom_name)
	is.in.chrom.pair(chrom_index, c1, c2)
	is.in.chroms(chrom_index, ...)
	is.inter.chromosomic(chrom_index)
	is.intra.chromosomic(chrom_index)
	known.chrom(chrom_index)
}

\arguments{
	\item{chrom_index}{ The object created by \code{\link{load.index}} }
	\item{chrom_name,c1,c2}{ 'string' containing the name of the desired chromosomes. }
	\item{...}{ The names of the desired chromosomes. }
}

\value{
	A 'vector' containing TRUE or FALSE, indicating which interactions described on your SIF fulfill your query.
}

\details{
The index should contain one row for each row of the sif file and two columns,
each column should contain the following info, corresponding to the genes found in your SIF:

	chromosome1 <tab> chromosome2

This object is loaded by the \code{\link{load.index}} function.
}

\examples{
	## The functions are currently defined as:
	# known.chrom <- function(chrom_index) {
	# 	return(!is.na(chrom_index$V1) & !is.na(chrom_index$V2))
	# }

	# is.intra.chromosomic <- function(chrom_index){
	# 	return(
	# 		known.chrom(chrom_index)
	# 		& (chrom_index$V1 == chrom_index$V2)
	# 	)
	# }

	# is.inter.chromosomic <- function(chrom_index){
	# 	return(
	# 		known.chrom(chrom_index)
	# 		& chrom_index$V1 != chrom_index$V2
	# 	)
	# }

	# is.in.chrom <- function(chrom_index, chrom_name) {
	# 	return(
	# 		known.chrom(chrom_index)
	# 		& (chrom_index$V1 == chrom_name | chrom_index$V2 == chrom_name)
	# 	)
	# }

	# is.in.chroms <- function(chrom_index, ...) {
	# 	return(
	# 		known.chrom(chrom_index)
	# 		& (chrom_index$V1 %in% c(...) | chrom_index$V2 %in% c(...))
	# 	)
	# }

	# is.in.chrom.pair <- function(chrom_index, c1, c2) {
	# 	return(
	# 		known.chrom(chrom_index)
	# 		& ((chrom_index$V1 == c1 & chrom_index$V2 == c2)
	# 		| (chrom_index$V1 == c2 & chrom_index$V2 == c1))
	# 	)
	# }
}

\seealso{
	\code{\link{load.index}}
	\code{\link{chrom.lookup-package}}
}

\keyword{SIF}
\keyword{gene}
\keyword{subset}
\keyword{chromosome}
