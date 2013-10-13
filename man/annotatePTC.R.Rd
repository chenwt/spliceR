\name{annotatePTC}
\alias{annotatePTC}
\title{PTC and NMD-sensitivity detection from assembled RNA-seq data}
\description{
PTC and NMD-sensitivity detection from assembled RNA-seq data.
}
\usage{
annotatePTC(transcriptData, cds, genomeObject, PTCDistance=50, filters, expressionCutoff=0)
}
\arguments{
  \item{transcriptData}{A \code{SpliceRList} object, containing transcript and exon information.} 
  \item{cds}{A \code{CDSSet} object, containing CDS information.} 
  \item{genomeObject}{A \code{BSgenome} object, containing sequence for the relevant genome. Contained in BSGenome objects, downloadable from BioConductor.} 
  \item{PTCDistance}{A numeric giving the premature termination codon-distance: The minimum distance from a STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive.} 
  \item{filters}{Character, giving the filters that should be applied - any combinations of 'geneOK', 'expressedGenes', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass' and/or 'sigIso'. Works only for data from cufflinks, as a manually generated \code{SpliceRList} does not include these metacolumns.}
  \item{expressionCutoff}{Numeric, giving the expression threshold (often in FPKM) used for the 'expressedGenes' and 'expressedIso' filter. Default value is 0.}

  }
\details{
\code{annotatePTC} retrieves sequence data for all exons given in \code{transcriptData}, uses the CDS-information in \code{cds} to scan for the most upstream reading frame, and translates the mRNA, storing information about the first codon in relation to distance from TTS, distance to the final exon-exon junction, etc.
If the STOP distance to the final exon-exon junction is larger than the threshold given in \code{PTCDistance} (and the STOP does not fall in the last exon), the STOP is considered premature and the transcript is marked as NMD (nonsense mediated decay) sensitive. For a review of the PTC and NMD mechanism, see Weischenfeldt et al. 2012. 
}

\value{
A \code{SpliceRList}, with the \code{transcript_features} object containing additional columns:
\code{spliceR.cds_pos}, the nucleotide position (in genomic coordinates) of the start of the most upstream open reading frame.
\code{spliceR.stop_pos}, the nucleotide position (in mRNA coordinates) of the end (first STOP codon position) of the most upstream open reading frame.
\code{spliceR.stop_distance}, the distance in nt (in mRNA coordinates) from the stop position to the end of the transcript.
\code{spliceR.junction_index}, the index of exon in which the stop codon falls, where 0 is the last exon of the transcript, -1 is the second-last, etc; NA, if \code{annotatePTC} was not able to find a ORF.
\code{spliceR.PTC}, a boolean, indicating whether the transcript is (theoretically) susceptible to nonsense mediated decay. annotatePTC sets this value to TRUE if the stop codon falls if any exon other than the last, and the distance to the final downstream exon-exon junction is larger than \code{PTCDistance}.
}

\author{Kristoffer Knudsen, Johannes Waage}

\references{
Vitting-Seerup K, Porse BT, Sandelin A, Waage JE. (2013) spliceR: An R package for classification of alternative splicing and prediction of coding potential from RNA-seq data. \emph{PeerJ PrePrints} 1:e80v1
}

\examples{
\dontrun{

#Rebuild cummeRbund's internal dataset
cuffDB <- readCufflinks(dir=system.file("extdata", package="cummeRbund"), gtf=system.file("extdata/chr1_snippet.gtf", package="cummeRbund"), genome='hg19' ,rebuild=TRUE)

#Generate SpliceRList from cufflinks data
cuffDB_spliceR <- prepareCuff(cuffDB)

# Require BSgenome object, containing genomic sequence
require("BSgenome.Hsapiens.UCSC.hg19",character.only = TRUE)

#Get CDS from UCSC
ucscCDS <- getCDS(selectedGenome="hg19", repoName="UCSC")

#Annotate with PTCs
cuffDB_spliceR_PTC <- annotatePTC(cuffDB_spliceR, cds=ucscCDS, Hsapiens, PTCDistance=50)

}
}
