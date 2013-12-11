\name{spliceRPlot}
\alias{spliceRPlot}
\title{Plot venn diagrams of alternative splicing events}
\description{Plot venn diagrams of alternative splicing events vs. samples.}
\usage{
spliceRPlot(spliceRobject, evaluate="nr_transcript", asType="All",colors=NULL, alpha=NULL, reset=FALSE, filters=NULL, expressionCutoff=0)
}
\arguments{
  \item{spliceRobject}{A \code{SpliceRList} object, processed and returned by \code{\link{spliceR}}.}
  \item{evaluate}{A character, giving the evaulation criteria (see details).}
  \item{asType}{The alternative splicing type to visualize, either 'ESI','MEE','MESI','ISI','A5', 'A3','ATSS','ATTS' or 'All'. See \code{\link{spliceR}} for a full description of alternative splicing types.}
  \item{colors}{Character, giving plot colors for each condition. Must be same length as number of conditions. If NULL, colors from the ColorBrewer "Dark2" pallette is used.}
  \item{alpha}{A numeric between 0 and 1, giving the transparency of the plot. If NULL, the alpha will be set optimally depending on number of samples.}
  \item{reset}{A boolean, indicating whether to reinitialize the \code{SpliceRList} object for faster replotting.}
  \item{filters}{vector, giving the filters that should be applied - any combinations of 'geneOK', 'expressedGenes', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass' and/or 'sigIso'. Works only for data from cufflinks, as a manually generated \code{SpliceRList} does not include these metacolumns.}
  \item{expressionCutoff}{Numeric, giving the expression threshold (often in FPTKM) used for the 'expressedGenes' and 'expressedIso' filter. Default value is 0.}

}

\details{
	Upon inital usage of spliceRPlot, the \code{SpliceRList} is initiated with internal data, allowing for faster replotting. If the \code{SpliceRList} changes because of filtering or other manipulation, rerun spliceRPlot with \code{reset=T}.
	For the evaulate parameter, the following are valid:
	'nr_transcript','nr_genes', 'nr_transcript_pr_gene', 'nr_AS', 'mean_AS_gene', 'mean_AS_transcript', 'mean_transcript_exp',  'mean_gene_exp'.
  'nr_transcript' outputs number of transcripts, 'nr_AS' outputs number of alternative splicing events, 'mean_as' outputs the average number of AS events per gene, 'mean_transcript_exp' outputs the mean transcript expression and 'mean_gene_exp' output the mean gene expression.
	For a detailed description of filters, see \code{\link{spliceR}}.
}

\value{
A \code{SpliceRList}, contianing additional temporary data for fast subsequent re-plotting.
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
}