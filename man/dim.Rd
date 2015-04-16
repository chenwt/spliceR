\name{dim}
\alias{dim.SpliceRList}
\alias{length.SpliceRList}
\title{Retrieve the Dimensions of a SpliceRList}
\description{
Retrieve the number of (transcripts) contained in \code{SpliceRList}.
}
\usage{
\method{dim}{SpliceRList}(x)
\method{length}{SpliceRList}(x)
}
\arguments{
  \item{x}{an object of class \code{SpliceRList}.}
}
\details{
As documented in \code{SpliceRList}, a SpliceRList contains two objects; a 'transcript_features' GRanges object, containing transcript information, and a 'exon_features' object, containing exon information. Dim and length currently only gives information about the number of transcripts in a SpliceRList object, i.e. the length() of the 'transcript_features' GRanges object.

}
\value{
Numeric vector of length 1, indicating the number of trancript comparisons in the SpliceRList.
}

\seealso{
  \code{\link[base]{dim}} in the base package.
}
