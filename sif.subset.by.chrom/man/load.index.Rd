\name{load.index}
\alias{load.index}

\title{ Load a chromosome index to subset a SIF }

\descrition{ Given a 'filename' (character), it reads the file and returns a chromosomic_index (data.frame) }

\usage{ load.index("/path/to/index") }

\arguments{
	\item{filename}{ Absolute or relative path to your index file }
}

\value{
	A 'data.frame' containing the chromosome names for your SIF.
}

\details{
The index should contain one row for each row of the sif file and two columns,
each column should contain the following info, corresponding to the genes found in your SIF:

	chromosome1 <tab> chromosome2
}

\examples{
	## The function is currently defined as:
	function(filename){
	  idx = data.table::fread(input = filename, data.table = FALSE)
	  idx$V1 = as.character(idx$V1)
	  idx$V2 = as.character(idx$V2)
	  return(as.data.frame(idx))
	}
}

\keyword{ SIF }
\keyword{ gene }
\keyword{ subset }
\keyword{ chromosome }
