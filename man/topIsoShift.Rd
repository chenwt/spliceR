\name{topIsoShift}
\alias{topIsoShift}
\title{Returns top transcripts in terms of isoform switching}
\description{
Returns top transcripts in terms of isoform switching.
}
\usage{
topIsoShift(spliceRObject, n=10)
}
\arguments{
  \item{spliceRObject}{a \code{SpliceRList} object, that has been successfully analyzed and annotated by \code{\link{spliceR}}.} 
  \item{n}{An integer, giving the number of transcripts to return.} 
}
\details{
This helper function returns the transcripts with the highest delta-isoform fraction (dIF) between samples.
If the data is based on cufflinks (source_id=="cufflinks"), only isoforms flagged significantly changing between samples will be returned.
}

\value{
A dataframe, containing a cast of the GRanges rows of the highest scoring transcripts by dIF.
}

\author{Kristoffer Knudsen, Johannes Waage}

\references{
Knudsen, K, Waage, J (2013).
spliceR, an R-package for splice class and PTC annotation of RNA-seq data
\emph{Journal} xx,xx
}

\examples{
#Load cufflinks example data
cuffDB <- prepareCuffExample()

#Generate SpliceRList from cufflinks data
cuffDB_spliceR <- prepareCuff(cuffDB)

#Reduce dataset size for fast example runtime
cuffDB_spliceR[[1]] <- cuffDB_spliceR[[1]][1:500]

#Run spliceR
mySpliceRList <- spliceR(cuffDB_spliceR, compareTo='preTranscript', filters=c('expressedGenes','geneOK', 'isoOK', 'expressedIso', 'isoClass'))

#Get top dIF transcripts
topIsoShift(mySpliceRList, n=20)
}