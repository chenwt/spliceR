\name{annotatePTC}
\alias{annotatePTC}
\title{PTC and NMD-sensitivity detection from assembled RNA-seq data}
\description{
PTC and NMD-sensitivity detection from assembled RNA-seq data.
}
\usage{
annotatePTC(transcriptData, cds, genomeObject, PTCDistance=50)
}
\arguments{
  \item{transcriptData}{A \code{SpliceRList} object, containing transcript and exon information.} 
  \item{cds}{A \code{CDSSet} object, containing CDS information.} 
  \item{genomeObject}{A \code{BSgenome} object, containing sequence for the relevant genome. Contained in BSGenome objects, downloadable from BioConductor.} 
  \item{PTCDistance}{A numeric giving the premature termination codon-distance: The minimum distance from a STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive.} 

  }
\details{
\code{annotatePTC} retrieves sequence data for all exons given in \code{transcriptData}, uses the CDS-information in \code{cds} to scan for the most upstream reading frame, and translates the mRNA, storing information about the first codon in relation to distance from TTS, distance to the final exon-exon junction, etc.
If the STOP distance to the final exon-exon junction is larger than the threshold given in \code{PTCDistance} (and the STOP does not fall in the last exon), the STOP is considered premature and the transcript is marked as NMD (nonsense mediated decay) sensitive. For a review of the PTC and NMD mechanism, see Weischenfeldt et al. 2012. 
}

\value{
A \code{SpliceRList}, with the \code{transcript_features} object containing additional columns:

\code{spliceR.cdsPosGenomic}, the genomic position of the used START codon.

\code{spliceR.stopPosGenomic}, the genomic position of the identified STOP codon,

\code{spliceR.ExonWithStart}, the exon which the used START codon falls within,

\code{spliceR.ExonWithStop}, the exon which the STOP codon falls within.

\code{spliceR.cdsPosTranscipt}, the position relative to transcript start where the used START codon is (measured in nucleotides).

\code{spliceR.stopPosTranscipt}, the position relative to transcript start where the found STOP codon is (measured in nucleotides).

\code{spliceR.stopDistance}, the distance from the found STOP codon to the last exon-exon junction, realative to transcript start.

\code{spliceR.junctionIndex}, the exon number in which the found STOP codon falls when compared to the last exon-exon junction, where 0 is the last exon of the transcript, -1 is the second-last, etc NA, if annotatePTC was not able to find a ORF.

\code{spliceR.PTC}, a boolean, indicating whether the transcript is (theoretically) susceptible to nonsense mediated decay. annotatePTC sets this value to TRUE if the stop codon falls if any exon other than the last, and the distance to the final downstream exon-exon junction is larger than PTCDistance (default 50 nt).

}

\author{Kristoffer Vitting-Seerup, Johannes Waage}

\references{
Vitting-Seerup K , Porse BT, Sandelin A, Waage J. (2014) spliceR: an R package for classification of alternative splicing and prediction of coding potential from RNA-seq data. BMC Bioinformatics 15:81.
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
