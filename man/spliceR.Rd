\name{spliceR}
\alias{spliceR}
\title{Splice class detection from assembled RNA-seq data}
\description{
Splice class detection from assembled RNA-seq data.
}
\usage{
spliceR(transcriptData, compareTo, filters, expressionCutoff=0, useProgressBar=T)
}
\arguments{
  \item{transcriptData}{a \code{SpliceRList} object, either created manually from transcript and exon information (see \code{\link{SpliceRList}}), or created by \code{\link{prepareCuff}} from CuffLinks data.}
  \item{compareTo}{a character, either 'preTranscript', for comparison to the hypothetical pre-splicing transcript for each gene, or a character, indicating the reference sample against which to classify splicing events.}
  \item{filters}{vector, giving the filters that should be applied - any combinations of 'geneOK', 'expressedGenes', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass' and/or 'sigIso'. Works only for data from Cufflinks, as a manually generated \code{SpliceRList} does not include these metacolumns.}
  \item{expressionCutoff}{Numeric, giving the expression threshold (often in FPKM) used for the 'expressedGenes' and 'expressedIso' filter. Default value is 0.}
  \item{useProgressBar}{Boolean, indicating whether to use progressbars. For compatibility. Default = TRUE.}

}

\details{
	The following filters are allowed for \code{filters}:
	\code{geneOK} requires Cufflinks to have reported the quantification of the gene as OK. Only works on transcript data from Cufflinks.
	\code{expressedGenes} requires the parent gene to be expressed.
	\code{sigGenes} requires the parent gene to be expressed in at least one sample.
	\code{isoOK} requires cufflinks to have reported the quantification of the isoform as OK. Only works on transcript data from Cufflinks.
	\code{expressedIso} requires the isoform to be expressed in at least one sample.
	\code{isoClass} removed transcripts marked by cufflinks to be either 'possible pre-mRNA fragment', 'Possible polymerase run-on fragment', or 'Repeat'. Only works on transcript data from Cufflinks.
	\code{sigIso} requires cufflinks to have reported the isoform as significant deregulated between samples. Only works on transcript data from Cufflinks.
}

\value{
A \code{SpliceRList}, identical to input SpliceRList transcriptData, with the \code{transcript_features} slot containing the following additional columns:
	\item{spliceR.major}{yes/no, indicating if this isoform is the major isoform expressed of the relevant gene for the reference sample.}
	\item{spliceR.IF1}{Isoform Fraction of total gene expression for sample 1}
	\item{spliceR.IF2}{Isoform Fraction of total gene expression for sample 2}
	\item{spliceR.dIF}{Delta IF (sample 2-sample 1)}
	\item{spliceR.ESI}{Number of exon skipping/inclusion events for this isoform}
	\item{spliceR.MEE}{Number of mutually exclusive exon events for this isoform}
	\item{spliceR.MESI}{Number of mutliple exon skipping/inclusion events for this isoform}
	\item{spliceR.ISI}{Number of intron skipping/retention events for this isoform}
	\item{spliceR.A5}{Number of alternative 5' splice site events for this isoform}
	\item{spliceR.A3}{Number of alternative 3' splice site events for this isoform}
	\item{spliceR.ATSS}{0/1, 1 indicating that this isoform uses an alternative transcription start site}
	\item{spliceR.ATTS}{0/1, 1 indicating that this isoform uses an alternative transcription terminating site}
	\item{spliceR.analyzed}{Yes/no, indicating if this isoform was analyzed(yes), or removed in filtering(no)}
	\item{spliceR.ESI.start}{Genomic start location(s) of ESI elements spliced in/out}
	\item{spliceR.ESI.end}{Genomic end location(s) of ESI elements spliced in/out}
	\item{spliceR.MEE.start}{Genomic start location(s) of MEE elements spliced in/out}
	\item{spliceR.MEE.end}{Genomic end location(s) of MEE elements spliced in/out}
	\item{spliceR.MESI.end}{Genomic end location(s) of MESI elements spliced in/out}
	\item{spliceR.MESI.start}{Genomic start location(s) of MESI elements spliced in/out}
	\item{spliceR.ISI.start}{Genomic start location(s) of ISI elements spliced in/out}
	\item{spliceR.ISI.end}{Genomic end location(s) of ISI elements spliced in/out}
	\item{spliceR.A5.start}{Genomic start location(s) of A5 elements spliced in/out}
	\item{spliceR.A5.end}{Genomic end location(s) of A5 elements spliced in/out}
	\item{spliceR.A3.start}{Genomic start location(s) of A3 elements spliced in/out}
	\item{spliceR.A3.end}{Genomic end location(s) of A3 elements spliced in/out}
	\item{spliceR.ATSS.start}{Genomic start location(s) of ATSS elements spliced in/out}
	\item{spliceR.ATSS.end}{Genomic end location(s) of ATSS elements spliced in/out}
	\item{spliceR.ATTS.start}{Genomic start location(s) of ATTS elements spliced in/out}
	\item{spliceR.ATTS.end}{Genomic end location(s) of ATTS elements spliced in/out}
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
}