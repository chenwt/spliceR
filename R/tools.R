    transcripts <- function(transcriptData)
{
    if(!is(transcriptData,"SpliceRList")) stop(class(transcriptData), " is not a SpliceRList")
	return(transcriptData[["transcript_features"]])
}

exons <- function(transcriptData)
{
	if(!is(transcriptData,"SpliceRList")) stop(class(transcriptData), " is not a SpliceRList")
	return(transcriptData[["exon_features"]])
}

conditions <- function(object) 
{
	if(class(object)[1] == 'SpliceRList')     {
		return( object$conditions )
	}
	else if( class(object)[1] == 'CuffSet' )  {
		return( samples(genes(object)) )
	} 
	else {
		stop('The object supplied must be a \'SpliceRList\' or a \'CuffSet\'.')
	}
}

topIsoShift <- function(spliceRObject, n=10) {
	#To avoid R CMD build NOTEs (subset problem)
	spliceR.iso_significant = NULL

	if(!is(spliceRObject,"SpliceRList")) stop(class(spliceRObject), " is not a SpliceRList")

	if(!"spliceR.analyzed"%in%colnames(mcols(spliceRObject[["transcript_features"]])))
		stop("SpliceRList has not yet been analyzed. Run spliceR() first...")

	transcriptData                <- GenomicRanges::as.data.frame( spliceRObject[['transcript_features']] )
	analyzedTranscriptData        <- transcriptData[ which(transcriptData$spliceR.analyzed == 'yes') ,]
	analyzedTranscriptDataSorted  <- analyzedTranscriptData[ sort.list( abs(analyzedTranscriptData$spliceR.dIF), decreasing=TRUE),]
	
	if (spliceRObject[["source_id"]]=="cufflinks") {
	    analyzedTranscriptDataSorted  <- subset(analyzedTranscriptDataSorted, spliceR.iso_significant=="yes")
	}
  
    if(nrow(analyzedTranscriptDataSorted) < n) { # If the user as for more than there is
        message(paste('spliceR only found',nrow(analyzedTranscriptDataSorted),'sigificant isoforms in the dataset'))
        n <- nrow(analyzedTranscriptDataSorted) # set n lower
    }

	return(analyzedTranscriptDataSorted[1:n,])
}

totalNumberOfAS <- function(spliceRObject) {
    if(!is(spliceRObject,"SpliceRList")) stop(class(spliceRObject), " is not a SpliceRList")

    if(!"spliceR.analyzed"%in%colnames(mcols(spliceRObject[["transcript_features"]])))
        stop("SpliceRList has not yet been analyzed. Run spliceR() first...")

    if(is.null(spliceRObject$transcripts_plot))
        stop("SpliceRList has not yet been analyzed using spliceRPlot, run spliceRPlot() first...")

    return(colSums(spliceRObject$transcripts_plot$isoforms[, c("ESI", "MEE", "MESI", "ISI", "A5", "A3", "ATSS", "ATTS", "All")]))
}

preSpliceRFilter <- function(spliceRobject, filters, expressionCutoff=0) {
    # Check class and GRanges
    if (!class(spliceRobject)[1]=="SpliceRList") stop("spliceRobject argument is not of class SpliceRList")
    if ( class(spliceRobject$"transcript_features") != "GRanges" || class(spliceRobject$"exon_features") != "GRanges" ) stop("spliceRobject must have GRanges objects in slots 'transcript_features' and 'exon_features'") 
    
    # Validate required columns in spliceRList
    t_colNames <- colnames(mcols(spliceRobject$"transcript_features"))
    if(!all(c("isoform_id", "sample_1", "sample_2", "gene_id", "iso_value_1", "iso_value_2", "iso_q_value") %in% substr(t_colNames, 9, nchar(t_colNames))))
        stop("Transcript features GRanges not compatible with spliceR - see documentation for more details")
    
    e_colNames <- colnames(mcols(spliceRobject$"exon_features"))
    if(!all(c("isoform_id","gene_id") %in% substr(e_colNames, 9, nchar(e_colNames))))
        stop("Exon features GRanges not compatible with spliceR - see documentation for more details")
    
    dataOrigin <- spliceRobject[["source_id"]]
    
    if(! dataOrigin %in% c('cufflinks', 'granges')  ) {
        stop('The input data was not recogniced, please see ?SpliceRList for more information about the input files')
    }
    
    if(any(!filters %in% c('geneOK','expressedGenes','sigGenes','isoOK','expressedIso','isoClass','sigIso'))) { # if one or more of the supplied filters are not recogniced
        stop('One or more of the supplied filters are not recogniced, please see ?preSpliceRfilter for more information about the filters')
    }
    
    message(length(spliceRobject$transcript_features$spliceR.isoform_id), " entries pre-filtering...")
    
    # Rename the GRange object (without converting it to data.frame)
    temp <- colnames(mcols(spliceRobject$transcript_features))
    temp2 <- substr(temp,9,nchar(temp))
    temp2[which(temp2 %in% c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", "element"))] <- paste('_',temp2[which(temp2 %in% c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", "element"))], sep='') # nessesary because of name duplications
    colnames(mcols(spliceRobject$transcript_features)) <- temp2
    
    isoformsToAnalyzeIndex <- 1:length(spliceRobject$transcript_features$isoform_id)
    
    # Optional filters 
    if('geneOK'         %in% filters) { isoformsToAnalyzeIndex <- .filterOKGenes(           spliceRobject, isoformsToAnalyzeIndex) }
    if('expressedGenes' %in% filters) { isoformsToAnalyzeIndex <- .filterExpressedGenes(    spliceRobject, isoformsToAnalyzeIndex, expressionCutoff) }
    if('sigGenes'       %in% filters) { isoformsToAnalyzeIndex <- .filterSigGenes(          spliceRobject, isoformsToAnalyzeIndex) }
    if('isoOK'          %in% filters) { isoformsToAnalyzeIndex <- .filterOKIso(             spliceRobject, isoformsToAnalyzeIndex) }
    if('expressedIso'   %in% filters) { isoformsToAnalyzeIndex <- .filterExpressedIso(      spliceRobject, isoformsToAnalyzeIndex, expressionCutoff) }
    if('isoClass'       %in% filters) { isoformsToAnalyzeIndex <- .filterIsoClassCode(      spliceRobject, isoformsToAnalyzeIndex) }
    if('sigIso'         %in% filters) { isoformsToAnalyzeIndex <- .filterSigIso(            spliceRobject, isoformsToAnalyzeIndex) }  
    
    # use the index to reduce the number of lines
    spliceRobject$transcript_features <- spliceRobject$transcript_features[isoformsToAnalyzeIndex,]
    # Add original colnames
    colnames(mcols(spliceRobject$transcript_features)) <- temp
    message(length(spliceRobject$transcript_features$spliceR.isoform_id), " entries post-filtering...")
    
    # remove unwanted in the exon GRange as well
    spliceRobject$exon_features <- spliceRobject$exon_features[which(spliceRobject$exon_features$spliceR.isoform_id %in% spliceRobject$transcript_features$spliceR.isoform_id)]
    spliceRobject$transcript_features <- spliceRobject$transcript_features[which(spliceRobject$transcript_features$spliceR.isoform_id %in% spliceRobject$exon_features$spliceR.isoform_id)]
    
    spliceRobject$filter_params <- filters
    
    # return object
    return(spliceRobject)
}

prepareCuffExample <- function()
{
    dir <- tempdir()
    extdata <- system.file("extdata", package="cummeRbund")
    file.copy(file.path(extdata, dir(extdata)), dir)

    cuffDB <- readCufflinks(dir=dir, gtfFile=system.file("extdata/chr1_snippet.gtf", package="cummeRbund"), genome="hg19", rebuild=TRUE)
    cuffDB
}
