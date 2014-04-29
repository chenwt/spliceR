\name{transcripts}
\alias{transcripts}
\alias{exons}
\title{Returns the transcript or exon GRanges from a SpliceRList object}
\description{
Returns the transcript or exon GRanges object from a SpliceRList object.
}
\usage{
transcripts(transcriptData)
exons(transcriptData)
}
\arguments{
  \item{transcriptData}{a \code{SpliceRList} object. 
}}
\details{
These helper functions returns either the "transcript_features" or "exon_features" object of a SpliceRList object.
}

\value{
A \code{GRanges} object. See \code{SpliceRList} for a full description of the contents of a SpliceRList.
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

myTranscripts <- transcripts(cuffDB_spliceR)
myExons <- exons(cuffDB_spliceR)
}
