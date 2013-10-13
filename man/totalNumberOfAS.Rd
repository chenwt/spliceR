\name{totalNumberOfAS}
\alias{totalNumberOfAS}
\title{Returns total number of alternative splicing events }
\description{
Returns total number of alternative splicing events an SpliceRList.
}
\usage{
totalNumberOfAS(spliceRObject)
}
\arguments{
  \item{spliceRObject}{a \code{SpliceRList} object returned by \code{\link{spliceRPlot}}. 
}}
\details{
This helper function returns number of total number of alternative splicing events. Object must be analyzed by \code{\link{spliceRPlot}} first.
}

\value{
A vector, giving the total number of splicing events for each splice class.
}

\author{Kristoffer Knudsen, Johannes Waage}

\references{
Vitting-Seerup K, Porse BT, Sandelin A, Waage JE. (2013) spliceR: An R package for classification of alternative splicing and prediction of coding potential from RNA-seq data. \emph{PeerJ PrePrints} 1:e80v1
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

#Plot number of exon skipping/inclusion events 
mySpliceRList <- spliceRPlot(mySpliceRList, evaluate="nr_AS", asType="ESI")

totalNumberOfAS(mySpliceRList)
}
