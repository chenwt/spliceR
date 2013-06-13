\name{conditions}
\alias{conditions}
\title{Returns sample condictions of an SpliceRList or an CuffSet object}
\description{
Returns samples/condictions of an SpliceRList or an CuffSet object.
}
\usage{
conditions(object)
}
\arguments{
  \item{object}{a \code{SpliceRList} object or a \code{CuffSet} object. 
}}
\details{
This helper function returns the "conditions" slot of a  \code{SpliceRList}, or the "sample" slot of a \code{CuffSet}.
}

\value{
A \code{character} vector, giving the samples/conditions.
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

conditions(cuffDB)

#Generate SpliceRList from cufflinks data
cuffDB_spliceR <- prepareCuff(cuffDB)

conditions(cuffDB_spliceR)
}
