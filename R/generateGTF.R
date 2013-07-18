generateGTF <- function(transcriptData, filters=NULL,expressionCutoff=0, scoreMethod='local', filePrefix='spliceR_transcripts', shortDescription="SpliceR Transcripts", longDescription="Transcripts output by SpliceR",  useProgressBar=T)
{
  

    # Check class and GRanges
    if (!class(transcriptData)[1]=="SpliceRList") stop("transcriptData argument is not of class SpliceRList")
    if ( class(transcriptData$"transcript_features") != "GRanges" || class(transcriptData$"exon_features") != "GRanges" ) stop("transcriptData must have GRanges objects in slots 'transcript_features' and 'exon_features'") 

    # Validate required columns in transcriptData
    tColNames <- colnames(mcols(transcriptData$"transcript_features"))
    if(!all(c(
        "isoform_id", "sample_1", "sample_2", "gene_id", "iso_value_1", "iso_value_2", "iso_q_value") %in% substr(tColNames, 9, nchar(tColNames))
        )
    ) stop("Transcript features GRanges not compatible with spliceR - see documentation for more details")

    eColNames <- colnames(mcols(transcriptData$"exon_features"))
    if(!all(c(
        "isoform_id","gene_id") %in% substr(eColNames, 9, nchar(eColNames))
        )
    ) stop("Exon features GRanges not compatible with spliceR - see documentation for more details")

    message("Preparing transcript data...")

    # Check if the filters supplied are OK:
    dataOrigin <- transcriptData[["source_id"]]

    if(dataOrigin == 'cufflinks') { okFilters <- c('none','expressedGenes','geneOK', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass', 'sigIso', 'PTC', 'singleExon') }
    if(dataOrigin == 'granges') { okFilters <- c('none','SingleExon') }
    
  if(is.null(filters)) {
    if(is.null(transcriptData[['filter_params']] )) {
      filters <- 'none'
    } else {
      filters <- transcriptData[['filter_params']]
    }
  }
  
  if('PTC' %in% filters) { # if asked to filter on PTC
    if('spliceR.PTC' %in% colnames(as.data.frame(transcriptData$"transcript_features"[1,]))) { # check whether the spliceR object contain PTC info
      okFilters <- c(okFilters, 'PTC')
    } else {
      stop('spliceR cannot filter on PTC since no PTC info is advailable. PTC information can be obtained through annotatePTC() ')
    }
  }

    if(any(!filters %in% okFilters)) { # if one or more of the supplied filters are not recogniced
        stop('One or more of the supplied filters are not recogniced, please see ?determineAStypes for more information about the filters')
    }
  if(length(scoreMethod) != 1 ) {
    stop('the scoreMethod paramter must be of length 1')
  }
  if(!scoreMethod %in% c('local','global')) {
    stop('the scoreMethod parameter must be one of \'local\' or \'global\' ' )
  }
  
    message(length(unique(transcriptData[["transcript_features"]]$spliceR.isoform_id)), " isoforms pre-filtering...")
  
    message("Converting to internal objects...")

    #Convert to dataframe
    tempDF <- GenomicRanges::as.data.frame(transcriptData[["transcript_features"]], stringsAsFactors=F)
    tempDF <- data.frame(lapply(tempDF, function(x) {if (class(x)=="factor") as.character(x) else (x)}), stringsAsFactors=FALSE)
    colnames(tempDF) <- c(colnames(tempDF)[1:5], substr(colnames(tempDF)[6:ncol(tempDF)],9,nchar(colnames(tempDF)[6:ncol(tempDF)])))
    transcriptData[["transcript_features"]] <- tempDF

    tempDF <- GenomicRanges::as.data.frame(transcriptData[["exon_features"]])
    tempDF <- data.frame(lapply(tempDF, function(x) {if (class(x)=="factor") as.character(x) else (x)}), stringsAsFactors=FALSE)
    colnames(tempDF) <- c(colnames(tempDF)[1:5], substr(colnames(tempDF)[6:ncol(tempDF)],9,nchar(colnames(tempDF)[6:ncol(tempDF)])))
    transcriptData[["exon_features"]] <- tempDF

    rm(tempDF)

    message("Filtering...")

    isoformsToAnalyzeIndex <- 1:nrow(transcriptData[["transcript_features"]])

    # Apply the chosen filters
    if('geneOK'         %in% filters) { isoformsToAnalyzeIndex <- spliceR:::.filterOKGenes(        transcriptData, isoformsToAnalyzeIndex) }
    if('expressedGenes' %in% filters) { isoformsToAnalyzeIndex <- spliceR:::.filterExpressedGenes( transcriptData, isoformsToAnalyzeIndex, expressionCutoff) }
    if('sigGenes'       %in% filters) { isoformsToAnalyzeIndex <- spliceR:::.filterSigGenes(       transcriptData, isoformsToAnalyzeIndex) }
    if('isoOK'          %in% filters) { isoformsToAnalyzeIndex <- spliceR:::.filterOKIso(          transcriptData, isoformsToAnalyzeIndex) }
    if('expressedIso'   %in% filters) { isoformsToAnalyzeIndex <- spliceR:::.filterExpressedIso(   transcriptData, isoformsToAnalyzeIndex, expressionCutoff) }
    if('isoClass'       %in% filters) { isoformsToAnalyzeIndex <- spliceR:::.filterIsoClassCode(   transcriptData, isoformsToAnalyzeIndex) }
    if('sigIso'         %in% filters) { isoformsToAnalyzeIndex <- spliceR:::.filterSigIso(         transcriptData, isoformsToAnalyzeIndex) }  
    if('singleExon'     %in% filters) { isoformsToAnalyzeIndex <- spliceR:::.filterSingleExonIsoAll( transcriptData, isoformsToAnalyzeIndex) }
    if('PTC'            %in% filters) { isoformsToAnalyzeIndex <- spliceR:::.filterPTC(               transcriptData, isoformsToAnalyzeIndex) }
  
    message(length(unique(transcriptData[["transcript_features"]]$isoform_id[isoformsToAnalyzeIndex])), " isoforms post-filtering...")

    # Get isoform indexes

    isoformNamesToAnalyze <- unique(transcriptData[["transcript_features"]]$"isoform_id"[isoformsToAnalyzeIndex])
    nrOfIsoformsToAnalyze <- length(isoformNamesToAnalyze)

    # Splitting the dataset like this is faster than doing it one-at-the-time in the loop
    # extract the features of those isoforms those that should be analyzed
    isoformsToAnalyze <- transcriptData[["exon_features"]][ which( transcriptData[["exon_features"]]$"isoform_id" %in% isoformNamesToAnalyze ) ,]
    # Split these features
    isoformFeaturesSplit <- split( isoformsToAnalyze , f=isoformsToAnalyze$"isoform_id")
  
    ### Calculate gene expressen based on 
    #geneScoreSplit <- split( transcriptData[["transcript_features"]][isoformsToAnalyzeIndex,c('isoform_id','gene_id','')] , f=isoformsToAnalyze$"isoform_id")
  
    # Create dataframe for storing the GTF file (must be done here so i dont have to use rbind)
    numberOfisoforms <- nrow(isoformsToAnalyze)
    gtfFileGlobal <- data.frame(matrix(NA, ncol=9, nrow=numberOfisoforms), stringsAsFactors=F)
    colnames(gtfFileGlobal) <- c('chr','seqSource','feature','start','end','score','strand','frame','group')

    # Get condtion names
    conditionNames <- transcriptData$conditions
  
    # Make data.frame to store scores
    scoresGlobal <- data.frame(matrix(NA, ncol=length(conditionNames)+1, nrow=numberOfisoforms), stringsAsFactors=F)
    colnames(scoresGlobal) <- c('gene', conditionNames)

    message('Generating GTFs...')
    if (useProgressBar) pb <- txtProgressBar(min = 1, max = nrOfIsoformsToAnalyze, style = 3)

    rowCounter <- 1 # the row number form which to insert the localGTF file into the Global GTF file
    for(i in 1:nrOfIsoformsToAnalyze) { # loop over all isoforms to analyze
        # Extract expression/test info  
        transcriptInfo <- transcriptData[["transcript_features"]][ which( transcriptData[["transcript_features"]]$"isoform_id" == isoformNamesToAnalyze[i]  ), ] # just grap the first transcript with the correct name

        myExonInfo <- isoformFeaturesSplit[[isoformNamesToAnalyze[i]]]
        if(is.null(myExonInfo)) {
          warning('Some exons were not found. Please make sure transcriptData[["exon_features"]] contains the nessesary information')
          next # else it will result in errors
        }
        nrOfExons <- nrow(myExonInfo)


        ### Annotate scores
        myLocalScore <- NULL
        for(j in 2:(length(conditionNames)+1)) {
          localScoreIndex <- match( conditionNames[j-1], c(transcriptInfo$sample_1, transcriptInfo$sample_2 ))
          localScore <- c(transcriptInfo$iso_value_1, transcriptInfo$iso_value_2 )[localScoreIndex]
          myLocalScore <- c(myLocalScore, localScore)
          scoresGlobal[rowCounter:(rowCounter + (nrOfExons-1)),j] <- localScore
        }
        scoresGlobal[rowCounter:(rowCounter + (nrOfExons-1)),1] <- transcriptInfo$"gene_id"[1]
    
        # Determine feature type
        if(nrOfExons == 1) {
          featureType <- 'Transcript'
        } else {
          featureType <- 'Exon'
        }

        # Correct chr names from Ensamble
        if(length(grep('chr', myExonInfo$seqnames[1])) == 1) {
          myChr <-  myExonInfo$seqnames[1]
        } else {
          if(myExonInfo$seqnames[1] == 'MT') {
            myChr <- 'chrM'
          } else {
            myChr <- paste('chr',myExonInfo$seqnames[1],sep='')
          }
        }

        # Generate local GTF
        gtfFileLocal <- data.frame( chr=myChr,
                             seqSource='SpliceR',
                             feature=featureType,
                             start=myExonInfo$start,
                             end=myExonInfo$end,
                             score= NA,
                             strand=myExonInfo$strand,
                             frame='.',
                             group=paste(
                               paste('isoform id', transcriptInfo$isoform_id[1], sep=' '),
                               paste('gene id', transcriptInfo$gene_id[1], sep=' '),
                               paste('gene short name', transcriptInfo$gene_short_name[1], sep=' '),
                               paste('nearest refrence id', transcriptInfo$nearest_ref_id[1], sep=' '),
                               paste(paste('FPKM score',conditionNames, myLocalScore, sep=' '),collapse='; '),
                               sep='; '),
                             stringsAsFactors=F, row.names=NULL)

        # Append local GTF to global GTF
        gtfFileGlobal[rowCounter:(rowCounter + (nrOfExons-1)),] <- gtfFileLocal

        # Add to counter so the next itteration will set in correctly
        rowCounter <- rowCounter + nrOfExons

        if (useProgressBar) setTxtProgressBar(pb, i)

    } # End of loop over isoforms
    if (useProgressBar) close(pb)
  
  
    ### Calculate Scores
  if(scoreMethod == 'local') {
    ## log the fpkm values
    scoresGlobal[,2:ncol(scoresGlobal)] <- log2(scoresGlobal[,2:ncol(scoresGlobal)] + 1)
    ## split on gene level
    scoresGlobalSplit <- split(scoresGlobal, f=scoresGlobal$gene)
    
    calculateGeneScores <- function(df) {
      #print(df$gene[1])
      maxForThisGene <- max(df[2:ncol(df)])
      if(maxForThisGene == 0) {
        return(df[,-1])
      } else {
        normalizedScores <- apply(df[,2:ncol(df)], 2, function(x) x* 1000/maxForThisGene)
      }
      
      #print(normalizedScores)
      return(normalizedScores)
    }
    
    ## Normalize FPKM values to geneome browser scores, for each gene.
    #myNormalizedScores <- ldply(scoresGlobalSplit, function(x) calculateGeneScores(x)) # gives an error I cannot expain
    temp <- lapply(scoresGlobalSplit, function(x) calculateGeneScores(x))
    myNormalizedScores <- do.call(rbind, temp)
    
  } else {
    ### Calculate Scores
    maxVal <- quantile( log2(c(transcriptData$transcript_features$iso_value_1[isoformsToAnalyzeIndex],transcriptData$transcript_features$iso_value_2[isoformsToAnalyzeIndex])+1) , 0.99)
    multiplyFactor <- 1000/maxVal
    myNormalizedScores <- scoresGlobal
    myNormalizedScores[,2:ncol(scoresGlobal)] <- log2(scoresGlobal[,2:ncol(scoresGlobal)] +1) * multiplyFactor
  }
    
    ### Qucik fix to correct the error in the Genome browser with score 0 beeing interpred as score 1000
    correctUCSCscoreProblem <- function(x) {
      x[which(x == 0)] <- 1 # which will result in the same color coding
      return( x )
    }
    myNormalizedScores <- apply(myNormalizedScores, 2, function(x) correctUCSCscoreProblem(x) )

    
    ### Generate list of GTF files
    gtfList <- lapply(1:length(conditionNames), function(x) NA)
    for(j in 1:length(conditionNames)) {
        localGTF <- gtfFileGlobal
        localGTF$score <- myNormalizedScores[,j]
    gtfList[[j]] <- localGTF
    }
    names(gtfList) <- conditionNames

    message('Writing GTFs...')

    for(i in 1:length(conditionNames)) {
        # Generate header
        shortDescrip <- paste(shortDescription,'from',conditionNames[i],sep=' ')
        longDescrip <- paste(longDescription,'form',conditionNames[i],sep=' ')
        header <- sprintf('track name="%s" description="%s" useScore=1 visibility=full', shortDescrip, longDescrip)

        # Write header
        myFilename <- paste(filePrefix,"_",conditionNames[i],'.GTF',sep='')
        fileConn<-file(myFilename, open='w') # overwrites the file if it exists - else it creates it
        writeLines(header, fileConn)
        close(fileConn)

        # Write GTF
        write.table(gtfList[[i]], file=myFilename, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE) # append the GTF info to the file
    }
}