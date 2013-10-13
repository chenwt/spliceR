\name{prepareCuff}
\alias{prepareCuff}
\title{Prepare assembled RNA-seq data from Cufflinks for spliceR}
\description{
Prepare assembled RNA-seq data from Cufflinks for spliceR.
}
\usage{
prepareCuff(cuffDB)
}
\arguments{
  \item{cuffDB}{a \code{cuffDB} object, produced by cummeRbund. This object must have been generated with cummeRbund, using the \code{gtf} parameter (see example), for spliceR to extract transcript model and exon information.}
}

\details{
	Use \code{prepareCuff} to prepare a cummeRbund/Cufflinks DB object for use by spliceR (see example).
	Often, it's appropriate to prefilter cufflinks data after running prepareCuff with \code{\link{preSpliceRFilter}} to reduce overhead on downstream analyses.
}

\value{
A \code{SpliceRList} containing a \code{transcript_features} \code{GRanges} object with the following additional metacolumns extracted from the cufflinks DB:
	\item{spliceR.isoform_id}{Cufflinks unique isoform id}
	\item{spliceR.sample_1}{Sample 1 identifier}
	\item{spliceR.sample_2}{Sample 2 identifier}
	\item{spliceR.gene_id}{Cufflinks unique gene id}
	\item{spliceR.CDS_id}{Cufflinks unique CDS id}
	\item{spliceR.gene_short_name}{Cufflinks unique short gene name}
	\item{spliceR.TSS_group_id}{Cufflinks unique TSS id}
	\item{spliceR.class_code}{Cufflinks class code (see cufflinks documentation)}
	\item{spliceR.nearest_ref_id}{Nearest reference id}
	\item{spliceR.length}{Length of the transcript}
	\item{spliceR.gene_status}{Cufflinks gene quantification status}
	\item{spliceR.gene_value_1}{Gene FPKM value for sample 1}
	\item{spliceR.gene_value_2}{Gene FPKM value for sample 2}
	\item{spliceR.gene_log2_fold_change}{Log2 fold change of gene expression (sample2 / sample1)}
	\item{spliceR.gene_p_value}{P-value for differential testing of difference of gene expression between samples}
	\item{spliceR.gene_q_value}{Adjusted p-value for differential testing of difference of gene expression between samples}
	\item{spliceR.gene_significant}{Yes/no; yes if difference of gene expression is significant}
	\item{spliceR.iso_status}{Cufflinks isoform quantification status}
	\item{spliceR.iso_value_1}{Isoform FPKM value for sample 1}
	\item{spliceR.iso_value_2}{Isoform FPKM value for sample 2}
	\item{spliceR.iso_log2_fold_change}{Log2 fold change of isoform expression (sample2 / sample1)}
	\item{spliceR.iso_p_value}{P-value for differential testing of difference of isoform expression between samples}
	\item{spliceR.iso_q_value}{P-value for differential testing of difference of isoform expression between samples}
	\item{spliceR.iso_significant}{Yes/no; yes if difference of isoform expression is significant}
and a \code{exon_features} GRanges object containing exon model information.
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
}