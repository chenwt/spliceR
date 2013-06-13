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
Knudsen, K, Waage, J (2013).
spliceR, an R-package for splice class and PTC annotation of RNA-seq data
\emph{Journal} xx,xx
}

\examples{
#Load cufflinks example data
cuffDB <- prepareCuffExample()
}