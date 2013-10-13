\name{prepareCuffExample}
\alias{prepareCuffExample}
\title{Prepare the Cufflinks example data}
\description{
Prepare the Cufflinks example data set.
}
\usage{
prepareCuffExample()
}

\details{
This helper function prepares the Cufflinks example dataset, including the example GTF-file.
}

\value{
A CuffSet object.
}

\author{Kristoffer Knudsen, Johannes Waage}

\references{
Vitting-Seerup K, Porse BT, Sandelin A, Waage JE. (2013) spliceR: An R package for classification of alternative splicing and prediction of coding potential from RNA-seq data. \emph{PeerJ PrePrints} 1:e80v1
}


\examples{
#Load cufflinks example data
cuffDB <- prepareCuffExample()
}