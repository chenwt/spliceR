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

\author{Kristoffer Vitting-Seerup, Johannes Waage}

\references{
Vitting-Seerup K , Porse BT, Sandelin A, Waage J. (2014) spliceR: an R package for classification of alternative splicing and prediction of coding potential from RNA-seq data. BMC Bioinformatics 15:81.
}


\examples{
#Load cufflinks example data
cuffDB <- prepareCuffExample()

conditions(cuffDB)

#Generate SpliceRList from cufflinks data
cuffDB_spliceR <- prepareCuff(cuffDB)

conditions(cuffDB_spliceR)
}
