\name{preSpliceRFilter}
\alias{preSpliceRFilter}
\title{Filters on spliceR-lists for reduction of data sets}
\description{
Applies a number of filters on a spliceR object to reduce data set size before running downstream analyses.
}
\usage{
preSpliceRFilter(spliceRobject, filters, expressionCutoff=0)
}
\arguments{
  \item{spliceRobject}{a \code{SpliceRList} object, either created manually from transcript and exon information (see \code{\link{SpliceRList}}), or created by \code{\link{prepareCuff}} from CuffLinks data.}
  \item{filters}{vector, giving the filters that should be applied - any combinations of 'geneOK', 'expressedGenes', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass' and/or 'sigIso'. Works only for data from cufflinks, as a manually generated \code{SpliceRList} does not include these metacolumns.}
  \item{expressionCutoff}{Numeric, giving the expression threshold (often in FPKM) used for the 'expressedGenes' and 'expressedIso' filter. Default value is 0.}
}
\details{
Often, many genes and isoforms are flagged as not "OK" or "LOWDATA" by Cufflinks, indicating low confidence in these. This function is handy for reducing the data size of a Cufflinks data set to reduce running times for downstream analyses.

Note, that \code{preSpliceRFilter} removes trancsripts from the dataset permanently, reducing size, while the filter options of \code{\link{spliceR}} and \code{\link{annotatePTC}} only selects transcripts for analysis, but does not remove any data.
}


\value{
A \code{SpliceRList} with transcripts after filtering.
}

\author{Kristoffer Vitting-Seerup, Johannes Waage}

\references{
Vitting-Seerup K , Porse BT, Sandelin A, Waage J. (2014) spliceR: an R package for classification of alternative splicing and prediction of coding potential from RNA-seq data. BMC Bioinformatics 15:81.
}


\examples{
#Load cufflinks example data
cuffDB <- prepareCuffExample()

#Generate SpliceRList from cufflinks data
cuffDB_spliceR <- prepareCuff(cuffDB)

#Filter 
cuffDB_spliceR_filtered <- preSpliceRFilter(cuffDB_spliceR, filters=c("expressedIso", "isoOK", "expressedGenes", "geneOK"))

}
