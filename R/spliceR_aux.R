######################################################################################################
######################################## Filter functions ############################################
######################################################################################################

# Filter 1: only test genes with OK status (cufflinks specific)
.filterOKGenes <- function(dataList, isoIndex) {
  return( isoIndex[which(  dataList[["transcript_features"]]$gene_status[isoIndex] == 'OK'  )]   )
}

# Filter 2: Pnly test significantly differentially expressed genes (cufflinks specific)
.filterSigGenes <- function(dataList, isoIndex) {
  return( isoIndex[which(  dataList[["transcript_features"]]$gene_significant[isoIndex] == 'yes'  )]   )
}

# Filter 3: Only analyze those with isoform quant status = OK
.filterOKIso <- function(dataList, isoIndex) {
  return( isoIndex[which(  dataList[["transcript_features"]]$iso_status[isoIndex] == 'OK'  )]   )
}

# Filter 4: Only analyze those genes that are expressed (requires gene expression values)
.filterExpressedGenes <- function(dataList, isoIndex, expressionCutoff) {
  return( isoIndex[which(  dataList[["transcript_features"]]$gene_value_1[isoIndex] > expressionCutoff | dataList[["transcript_features"]]$gene_value_2[isoIndex] > expressionCutoff )] )
}

# Filter 4: Get index of those isoforms that are expressed in at least one of the samples
.filterExpressedIso <- function(dataList, isoIndex, expressionCutoff) {
  return( isoIndex[which(  dataList[["transcript_features"]]$iso_value_1[isoIndex] > expressionCutoff | dataList[["transcript_features"]]$iso_value_2[isoIndex] > expressionCutoff )] )# Remove those isoforms that are not expressed
}

.filterIsoClassCode <- function(dataList, isoIndex) {
  classCodeIndex <-  isoIndex[which(  dataList[["transcript_features"]]$class_code[isoIndex] != 'e' & dataList[["transcript_features"]]$class_code[isoIndex] != 'p' & dataList[["transcript_features"]]$class_code[isoIndex] != 'r'  )]
  # Remove those isoforms that are classified as :
  ### Single exon transfrag     (overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment)
  ### Possible polymerase run-on fragment     (within 2Kbases of a reference transcript)
  ### Repeat.     (Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case)
  naIndex <- which(is.na(dataList[["transcript_features"]]$class_code[isoIndex])) # NA's are generated for novel genes, or genes that are not assigned to a known gene.
  return(sort(c(classCodeIndex, naIndex),decreasing=F))
}

# Filter 6: Only analyze significant isoforms (only analyze genes which have two significant differential regulated isoforms (rare))
.filterSigIso <- function(dataList, isoIndex) {
  return( isoIndex[which(  dataList[["transcript_features"]]$iso_significant[isoIndex] == 'yes'  )] )
}

# Filter 7: remove isoforms with only 1 exon
.filterSingleExonIsoAll <- function(dataList,isoformsToAnalyzeIndex) {
  isoformIDsToAnalyze     <- dataList[["transcript_features"]]$"isoform_id"[isoformsToAnalyzeIndex] # get names of those isoforms that I should analyze
  noOfExonsPerIsoform     <- table(dataList[["exon_features"]]$"isoform_id"[ which(dataList[["exon_features"]]$"isoform_id" %in% isoformIDsToAnalyze) ] ) # create table with number of expons for those isoforms (here the number of conditions does not matter since there are no replicates in the feature table)
  isoformIndexesToAnalyze <- isoformsToAnalyzeIndex[ which( dataList[["transcript_features"]]$"isoform_id"[isoformsToAnalyzeIndex] %in% names(noOfExonsPerIsoform[noOfExonsPerIsoform>1]) ) ] # the those indexes that belongs to transcripts with more than one exon
  return(isoformIndexesToAnalyze)
}

# Filter 8: remove genes with less than two isoforms
.filterSingleIsoform <- function(dataList,isoformsToAnalyzeIndex, knownSamples) {
  getUniqGeneIsoformCombinations <- unique(dataList[["transcript_features"]][isoformsToAnalyzeIndex,c('isoform_id','gene_id')]) # extract unique combinations of geneID and isoformID from those indexes that should be analyzed (the uniqe part removes dublicates du to multiple sample comparasons but includes combinations that are only in the index to analyze for some of the comparasons)
  noOfisoformsPerGene     <- table(getUniqGeneIsoformCombinations$gene_id ) 
  isoformIndexesToAnalyze <- isoformsToAnalyzeIndex[ which( dataList[["transcript_features"]]$"gene_id"[isoformsToAnalyzeIndex] %in% names(noOfisoformsPerGene[noOfisoformsPerGene>1]) ) ] # overwrite the expression data with the expression data of those isoform features that have more than 1 exon
  return(isoformIndexesToAnalyze)
}

# Filter 9: remove all those rows that are PTC sensitive
.filterPTC <- function(dataList,isoIndex) {
  return( isoIndex[which( ! dataList[["transcript_features"]]$PTC[isoIndex]) ] ) # NA's are removed
}
########################################################################################################
############################ Functions for determining Major Isoform ###################################
########################################################################################################
.getMajorIsoCuffDB <- function(isoformSubset, referenceSample) {
  # if no refrence condition have been chosen take the isoform with the highest expression in any sample
  if(referenceSample == 'none') {
    maxValue <- max( 
      isoformSubset$iso_value_1, 
      isoformSubset$iso_value_2 
    ) # get the index of the the max value
    maxIsoformIndex <- c(
      which( isoformSubset$iso_value_1 == maxValue), # I cannot use which.max because its two distinct vectors and I need the index
      which( isoformSubset$iso_value_2 == maxValue)  # I cannot use which.max because its two distinct vectors and I need the index
    ) # get the indexes belonging to that value
  } else { # if a refrence condition have been chosen take the isoform with the highest expression in that condition  
    if( length(c(which(isoformSubset$sample_1 == referenceSample),which(isoformSubset$sample_2 == referenceSample))) > 0 ) {
      maxValue <- max( 
        isoformSubset$iso_value_1[which(isoformSubset$sample_1 == referenceSample)], 
        isoformSubset$iso_value_2[which(isoformSubset$sample_2 == referenceSample)] 
      ) # get the index of the the max value
      maxIsoformIndex <- c(
        which( isoformSubset$iso_value_1 == maxValue), # I cannot use which.max because its two distinct vectors
        which( isoformSubset$iso_value_2 == maxValue)  # I cannot use which.max because its two distinct vectors
      )  # get the indexes belonging to the max calue
    } else {  # if the refrence isoform is not in the subset analyzed (probably because the quantification is not ok whereby its filtered out)
      return(NULL)
    } 
  }
  
  # If there are several (different) isoforms with same expression level (that is not as rare as I thought it was) - chose the one with the most exons - if still eqiual chose random.
  if(length(unique(isoformSubset$isoform_id[maxIsoformIndex])) > 1) {
    # chose the longest one
    maxIsoformIndex <- maxIsoformIndex[which.max(isoformSubset$width[maxIsoformIndex])]
    
    # If there are several (different) isoforms with same length
    if(length(unique(isoformSubset$isoform_id[maxIsoformIndex])) > 1) {
      # chose one at random
      maxIsoformIndex <- sample(maxIsoformIndex,size=1)
    }
  }
  return(maxIsoformIndex)
}

###############################################################################################
################################ Core functions ############# #################################
###############################################################################################

.findOverlappingExons <- function(exon1, exon2) {
  if( exon1$start <  exon2$end   & exon1$start >= exon2$start ) { return(TRUE) } # test whether start of exon1 is contained within exon2
  if( exon1$end   <= exon2$end   & exon1$end   >  exon2$start ) { return(TRUE) } # test whether end of exon1 is contained within exon2
  if( exon1$start <  exon2$start & exon1$end   >  exon2$end   ) { return(TRUE) } # test whether exon1 contains exon2
  return(FALSE)
}

.findIdenticalExons <- function(exon1, exon2) {
  if(exon1$start == exon2$start & exon1$end == exon2$end) {
    return(TRUE)
  }
  return(FALSE)
}

### Functions to determin the type of Alternative splicing of overlapping exons in the LEFT side of the exons
.determineLeftOverlappingAStype <- function(exon1, exon2, isStrandEqualToPlus, ExonIndexesAnalyzed, asTypes) {
  localExonList <- list(exon1, exon2)
  if(exon1$start != exon2$start) {
    if(! any( ExonIndexesAnalyzed == c(1,1) )) { 
      ## Annotate AS type
      if(isStrandEqualToPlus) { 
        asTypes$A3 <- asTypes$A3 + 1
        asTypeStart <- 'A3.start' # string to help assign skipped part to the right type
        asTypeEnd <- 'A3.end' # string to help assign skipped part to the right type
      } else {
        asTypes$A5 <- asTypes$A5 + 1 
        asTypeStart <- 'A5.start' # string to help assign skipped part to the right type
        asTypeEnd <- 'A5.end' # string to help assign skipped part to the right type
      }
      ## Annotate spliced out
      minIndex <- which.min( c(exon1$start, exon2$start))
      notMinIndex <- (1:2)[!1:2 %in% minIndex]
      
      if(is.na(asTypes[asTypeStart])) {
        asTypes[asTypeStart] <- paste( localExonList[[minIndex]]$start )
        asTypes[asTypeEnd]   <- paste( localExonList[[notMinIndex]]$start )
      } else {
        asTypes[asTypeStart] <- paste(asTypes[asTypeStart], localExonList[[minIndex]]$start, sep=';')
        asTypes[asTypeEnd]   <- paste(asTypes[asTypeEnd],   localExonList[[notMinIndex]]$start,   sep=';')
      }
      
    }
  }
  return(asTypes)
}



### Functions to determin the type of Alternative splicing of overlapping exons in the RIGHT side of the exons
.determineRightOverlappingAStype <- function(exon1, exon2, isStrandEqualToPlus, ExonIndexesAnalyzed, numberOfExons, asTypes) { # needs exon numbers because it looks from the right side
  localExonList <- list(exon1, exon2)
  if(exon1$end != exon2$end) { # Check whether the coordinates are identical - nessesary since I only filter for completly identical exons
    if(! any( ExonIndexesAnalyzed == numberOfExons ) ) { # Chech wheter both exons are end exons (since the order of ExonIndexesAnalyzed and numberOfExons are the same == can be used)
      ## Annotate AS type
      if(isStrandEqualToPlus) { 
        asTypes$A5 <- asTypes$A5 + 1
        asTypeStart <- 'A5.start' # string to help assign skipped part to the right type
        asTypeEnd <- 'A5.end' # string to help assign skipped part to the right type
      } else {
        asTypes$A3 <- asTypes$A3 + 1
        asTypeStart <- 'A3.start' # string to help assign skipped part to the right type
        asTypeEnd <- 'A3.end' # string to help assign skipped part to the right type
      }
      ## Annotate spliced out
      minIndex <- which.min( c(exon1$end, exon2$end))
      notMinIndex <- (1:2)[!1:2 %in% minIndex]
      
      if(is.na(asTypes[asTypeStart])) {
        asTypes[asTypeStart] <- paste( localExonList[[minIndex]]$end )
        asTypes[asTypeEnd]   <- paste( localExonList[[notMinIndex]]$end )
      } else {
        asTypes[asTypeStart] <- paste(asTypes[asTypeStart], localExonList[[minIndex]]$end, sep=';')
        asTypes[asTypeEnd]   <- paste(asTypes[asTypeEnd],   localExonList[[notMinIndex]]$end,   sep=';')
      }
    }
  }
  return(asTypes)
}


### Function to evaluate counter of non-overlapping exons
.evaluateCounter <- function(count, asTypes, coordinats) {
  if( count == 1) {
    asTypes$ESI <- asTypes$ESI + 1
    
    if(is.na(asTypes$ESI.start)) {
      asTypes$ESI.start <- paste( coordinats$start )
      asTypes$ESI.end   <- paste( coordinats$end )
    } else {
      asTypes$ESI.start <- paste(asTypes$ESI.start, coordinats$start, sep=';')
      asTypes$ESI.end   <- paste(asTypes$ESI.end,   coordinats$end,   sep=';')
    }
    
  } else {
    asTypes$MESI <- asTypes$MESI + 1
    for(i in 1:nrow(coordinats)) {
      if(is.na(asTypes$MESI.start)) {
        asTypes$MESI.start <- paste( coordinats$start[i] )
        asTypes$MESI.end   <- paste( coordinats$end[i] )
      } else {
        asTypes$MESI.start <- paste(asTypes$MESI.start, coordinats$start[i], sep=';')
        asTypes$MESI.end   <- paste(asTypes$MESI.end,   coordinats$end[i],   sep=';')
      }
    }
    
  }
  
  return(asTypes)
}


### Function to determine non overlapping AS types
.determineNonOverlappingAStype <- function(transcript1, transcript2, numberOfExons, index1, index2, exonsToIgnore, asTypes) {
  ### Extract exon info from the skipping in transcript 1
  if(index1[2]-index1[1] >1) { # are there a skipping event for transcript1
    index11 <- seq(index1[1]+1,index1[2]-1) # extract the exon(s) index that should be compared
    
    if(any(index11 %in% exonsToIgnore[[1]])) { # if the exon(s) extracted are among thoes that should be ignored
      removeIndex <- which(index11 %in% exonsToIgnore[[1]])
      index11 <- index11[-removeIndex] # remve those to be ignored
      
      if(length(index11) > 0) { # if there is any exons left
        t1 <- data.frame(start = transcript1[index11,'start'], end = transcript1[index11,'end'] , startExon = (index11 %in% 1), endExon = (index11 %in% (numberOfExons[1])), transcript = 1, stringsAsFactors=F)
      } else {
        t1 <- data.frame()
      }
      
    } else { # if the exons are not in the ignore exons list
      t1 <- data.frame(start = transcript1[index11,'start'], end = transcript1[index11,'end'], startExon = (index11 %in% 1), endExon = (index11 %in% (numberOfExons[1])), transcript = 1, stringsAsFactors=F)
    }    
  } else { # if there is no skipping in this transcript
    t1 <- data.frame()
  }
  
  ### Extract exon info from the skipping in transcript 2
  if(index2[2]-index2[1] >1) { # are there a skipping event for transcript1
    index22 <- seq(index2[1]+1,index2[2]-1) # extract the exon(s) index that should be compared
    
    if(any(index22 %in% exonsToIgnore[[2]])) { # if the exon(s) extracted are among thoes that should be ignored
      removeIndex <- which(index22 %in% exonsToIgnore[[2]])
      index22 <- index22[-removeIndex] # remve those to be ignored
      
      if(length(index22) > 0) { # if there is any exons left
        t2 <- data.frame(start = transcript2[index22,'start'], end = transcript2[index22,'end'], startExon = (index22 %in% 1), endExon = (index22 %in% (numberOfExons[1])), transcript = 1, stringsAsFactors=F)
      } else {
        t2 <- data.frame()
      }
      
    } else { # if the exons are not in the ignore exons list
      t2 <- data.frame(start = transcript2[index22,'start'], end = transcript2[index22,'end'] ,startExon = (index22 %in% 1), endExon = (index22 %in% (numberOfExons[1])), transcript = 1, stringsAsFactors=F)
    }    
  } else { # if there is no skipping in this transcript
    t2 <- data.frame()
  }
  
  # combine and sort transcript (since they are all non-overlapping)
  t3 <- rbind(t1,t2)
  if(nrow(t3) == 0) { # are there any exon left (there migth not be due to the MEE)
    return(asTypes)
  }
  
  t3 <- t3[sort.list(t3$start),]
  
  ### remove all exons that are not within the overlapping body of the two transcripts (including ends - since they are determined elsewhere) 
  cuttOff1 <- 0
  if(any(t3$startExon)) {
    cuttOff1 <- tail(which(t3$startExon),1) +1 # get cutoff for start
  }
  cutOff2 <- nrow(t3)
  if(any(t3$endExon)) {
    cutOff2 <- which(t3$endExon)[1] -1 # get cuttoff for end
  }
  if(cutOff2 >= cuttOff1) { # chech that the cutoffs leves anything (else cuttOff1:cutOff2 whould create a negative sequence)
    t3 <- t3[cuttOff1:cutOff2,]
  } else {
    return(asTypes)
  }
  
  ### determine AS events
  if(nrow(t3) > 1) { # if there is MORE than 1 exon skipping event left to compare after first/last exons have been removed
    # Use the origine to determine the order of the events
    count <- 1
    for(j in 2:nrow(t3)) {
      if( t3$transcript[j] == t3$transcript[(j-1)] ) { # if the skipping occured from the same transcript as the previous
        count <- count +1 # add one to the counter
      } else {
        asTypes <- .evaluateCounter(count, asTypes, t3[(j-count+1):j,c('start','end')])
        count <- 1
      }
    }
    asTypes <- .evaluateCounter(count, asTypes, t3[(j-count+1):j,c('start','end')])
    
  } else if (nrow(t3) == 1) { # if there is only 1 exon skipping left to compare after first/last exons have been removed
    asTypes$ESI <- asTypes$ESI + 1
    
    if(is.na(asTypes$ESI.start)) {
      asTypes$ESI.start <- paste( t3$start )
      asTypes$ESI.end   <- paste( t3$end )
    } else {
      asTypes$ESI.start <- paste(asTypes$ESI.start, t3$start, sep=';')
      asTypes$ESI.end   <- paste(asTypes$ESI.end,   t3$end,   sep=';')
    }
    
    
  }
  return(asTypes)
} # end of determineNonOverlappingAStype function


### Function to create pre-RNA
.getPreRNA <- function(transcript) {
  # We want to extract the the pre mRNA but intron retensions must be excluded since they will give cause to missing classificantions.
  # For this purpose intron retension are defined as an exon, that overlaps at least two other exons, and at least two of the exon
  # the intron retension overlaps does not overlap with one another.
  
  # get unique exons - is better to do this before passing the df to the function
  #myExonInfoUniq <- unique(transcript[,c('start','end')])
  
  # sort them according to start and end coordinats (makes the overlap comparason  much faster)
  myExonInfoUniqSort <- transcript[order(transcript$start,transcript$end ),] 
  
  # Make a pairwise comparason of exons to find overlaps - have been tested thoroughly - works
  overlapIndex <- matrix(nrow=0, ncol=2) # data frame to store the overlapping exon coordinats
  #overlapIndexDF <- data.frame()
  count=2 # counter to avoid creating the same plot twice and from making plots against oneself
  for(i in 1:(nrow(myExonInfoUniqSort)-1)) { # loop over the samples (except the last since that have already been compared with by the second loop)
    for(j in count:nrow(myExonInfoUniqSort)) { # loop over the samples starting from 2, since I dont want to compare an exon with itself
      # Check whether the two exons are overlapping
      #if( .findOverlappingExons(exon1=myExonInfoUniqSort[i,], exon2=myExonInfoUniqSort[j,]) ) {
      if( .C("findOverlappingExons", as.integer(myExonInfoUniqSort[i,c("start","end")]), as.integer(myExonInfoUniqSort[j,c("start","end")]), as.integer(0))[[3]] ) {
        #overlapIndexDF <- rbind(overlapIndexDF, data.frame(i,j))      
        overlapIndex <- rbind(overlapIndex, matrix(c(i,j),nrow=1))      
      }
      ## make sure I skip the rest of the j's when the end coordinat of exon(i) is smaller than the start value of the next j exon (the exons are ordered)
      if(j !=nrow(myExonInfoUniqSort)) { 
        if( myExonInfoUniqSort[i,'end'] < myExonInfoUniqSort[j+1,'start'] ) {
          break
        }
      }
    } # j loop
    count=count+1 # add one to the counter to keep only looking at one half of the matrix of posibilities
  } # i loop
  # now overlapIndex contains the indexes of exons that overlap one another
  
  ### interpret the overlapIndex generated above
  # Find exons represented multiple times (meaning exons that overlaps with more than one other exon)
  #dublicatedExons <- unique(unlist(overlapIndexDF)[duplicated( unlist(overlapIndexDF) )])
  dublicatedExons <- unique(as.integer(overlapIndex)[duplicated( as.integer(overlapIndex) )])
  
  intronRetensions <- NULL
  if(length(dublicatedExons) > 0) { # if there are any exons that overlap more than one other exon
    for(dublicatedExon in dublicatedExons) { # loop over each of the exons
      # extract the indexs of the exons that this exon overlap with
      whichRowsContainsTheDuplicated <- which(apply(overlapIndex,1,function(x) dublicatedExon %in% x))
      RowsOverlappingWith <- overlapIndex[whichRowsContainsTheDuplicated,]
      #overlappingWith <- unlist(RowsOverlappingWith)[unlist(RowsOverlappingWith) != dublicatedExon]
      overlappingWith <- as.integer(RowsOverlappingWith)[unlist(RowsOverlappingWith) != dublicatedExon]
      
      ## make a pairwise comparason of the indexes that this exon overlaps with (nessesary since there migth be more than one)
      count <- 2
      for(i in 1:(length(overlappingWith)-1)) { # loop over the samples (except the last since that have already been compared with by the second loop)
        for(j in count:length(overlappingWith)) { # loop over the samples to compare with      
          # Compare the two exons to find thos that are NOT overlapping
          #if( ! .findOverlappingExons(exon1=myExonInfoUniqSort[overlappingWith[i],], exon2=myExonInfoUniqSort[overlappingWith[j],]) ) {  # it is faster to compare the two exons again thant to look for the indexes in the overlapping table
          if(!.C("findOverlappingExons", as.integer(myExonInfoUniqSort[overlappingWith[i],c("start","end")]), as.integer(myExonInfoUniqSort[overlappingWith[j],c("start","end")]), as.integer(0))[[3]] ) {
            # add it to the intron retension list
            intronRetensions <- c(intronRetensions, dublicatedExon )
          }
          
        }
        count=count+1 # add one to the counter to keep only looking at one half of the matrix of posibilities
      }
    }
  }  
  ### Create the pre-RNA, and make sure to exclude intron retensions
  # If any intron retensions were found exclude them from the preRNA
  if(!is.null(intronRetensions)) {
    myIRange <- IRanges(start=myExonInfoUniqSort$start[-intronRetensions], end=myExonInfoUniqSort$end[-intronRetensions])
  } else {
    myIRange <- IRanges(start=myExonInfoUniqSort$start, end=myExonInfoUniqSort$end)
  }
  end(myIRange) <- end(myIRange) -1 # hack to bypass the fact that IRanges start and end coordinats are 0 based and 1 based respectively, but the coordinats outputted by cufflinks are both 1-based. Does not matter that the 
  myReducedIrange <- reduce(myIRange , min.gapwidth=0 )
  end(myReducedIrange) <- end(myReducedIrange) +1 # add one to the ends again
  myPreRNA <- GenomicRanges::as.data.frame( myReducedIrange )
  
  ## make sure the end have not been cut off (happens somtimes due to intron retensions)
  myPreRNA$start[1] <- min(myExonInfoUniqSort$start)
  myPreRNA$end[nrow(myPreRNA)] <- max(myExonInfoUniqSort$end)
  
  ## Add strand (is removed by the IRange reduce)
  myPreRNA$strand <- transcript$strand[1]
  
  return(myPreRNA)
}

###############################################################################################
############## Core function that loops over the two transcripts comparing them ###############
###############################################################################################
.findOverlap <- function(transcript1, transcript2) {
  numberOfExons1 <- nrow(transcript1) # get number of exons in transcript 1
  numberOfExons2 <- nrow(transcript2) # get number of exons in transcript 2

  ##################################################################################################################
  ### Find overlapping and identical exons (using a while loop is much faster than for loops and IRanges
  # this is the time consuming step (the others take no time)
  exonIndex1 <- 1 # counter for the outer while loop
  exonIndex2 <- 1 # counter for the inner while loop
  overlappingExons <- data.frame()
  identicalExons <- data.frame()
  while(exonIndex1 <= numberOfExons1) {
    while(exonIndex2 <= numberOfExons2) {
      #print(c(exonIndex1,exonIndex2))
      # Test for identical exons (and therfore also overlapping)
      
      ##### C-FUNCTION
      if(.C("findIdenticalExons", as.integer(transcript1[exonIndex1,c('start','end')]), as.integer(transcript2[exonIndex2,c('start','end')]), as.integer(0))[[3]]) {
      #if(.findIdenticalExons(transcript1[exonIndex1,c('start','end')], transcript2[exonIndex2,c('start','end')])) {
        identicalExons <- rbind(identicalExons, c(exonIndex1,exonIndex2))     # if identical add both to identical and overlapping
        overlappingExons <- rbind(overlappingExons, c(exonIndex1,exonIndex2))
        
        if( exonIndex1 ==  numberOfExons1 & exonIndex2 == numberOfExons2) { # make sure I break if the last two are identical (because I only add to the counter if its not the alst two)
          break
        }
        
        if(exonIndex1 < numberOfExons1) {
          exonIndex1 <- exonIndex1 +1 # I only add if it means the counter does not exceed the number of exons in this transcript
        }
        if(exonIndex2 < numberOfExons2) {
          exonIndex2 <- exonIndex2 +1 # I only add if it means the counter does not exceed the number of exons in this transcript
        }
        next
      }
      
#       if(.findIdenticalExons(transcript1[exonIndex1,c('start','end')], transcript2[exonIndex2,c('start','end')])) {
#         identicalExons <- rbind(identicalExons, c(exonIndex1,exonIndex2))     # if identical add both to identical and overlapping
#         overlappingExons <- rbind(overlappingExons, c(exonIndex1,exonIndex2))
#       } else 
      # C-FUNCTION
      # Test for overlapping exons
      if(.C("findOverlappingExons", as.integer(transcript1[exonIndex1,c('start','end')]), as.integer(transcript2[exonIndex2,c('start','end')]), as.integer(0))[[3]]) {
      #if(.findOverlappingExons(transcript1[exonIndex1,c('start','end')], transcript2[exonIndex2,c('start','end')])) {
        overlappingExons <- rbind(overlappingExons, c(exonIndex1,exonIndex2))
      }
      # Determine which of the next exon have the lowest start coordinat so I know which one to go to (to make sure I dont skip exons)
      if(exonIndex1 < numberOfExons1 & exonIndex2 < numberOfExons2) { # nessesary because else there is no coordinat to compare
        if( transcript1[(exonIndex1+1),'start'] < transcript2[(exonIndex2+1),'start']  ) {
          break        
        } else if(transcript1[(exonIndex1+1),'start'] == transcript2[(exonIndex2+1),'start'] ) {
          # in the special case where the start coordinats are identical look at the end coordinats (nessesary because else intron retensions might cause problems)
          if(transcript1[(exonIndex1+1),'end'] < transcript2[(exonIndex2+1),'end'] ) {
            break
          }
        }
      } else if(exonIndex2 == numberOfExons2) { # make sure that exonIndex2 is not made larger than n2 (meaning its not untill the outer loop is done that comparasons will not be made anymore)
        break
      }
      exonIndex2 <- exonIndex2 + 1
    } # end of inner while loop
    exonIndex1 <- exonIndex1 + 1
  } # end of outer while loop
  # add names
  if(nrow(overlappingExons) > 0) {
    colnames(overlappingExons) <- c('isoform1','isoform2')
  }
  if(nrow(identicalExons)> 0) {
    colnames(identicalExons) <- c('isoform1','isoform2')
  }
  
  return(list(overlap=overlappingExons, idenctial=identicalExons))
} # end of find overlap


# transcript1 <- temp[[2]]
# transcript2 <- temp[[1]][which(temp[[1]]$transcript==2),]
# 
# overlappingExons <- test[[1]]
# identicalExons <-  test[[2]]
# #overlappingExons <- overlapList[[i]]
# #identicalExons <- identicalList[[i]]
# 
# exonsToIgnore <- list(NULL,NULL)
# 
# rm(list=c('transcript1','transcript2','overlappingExons','identicalExons','exonsToIgnore'))
# .determineAStypeOverlap(transcript1, transcript2, overlappingExons, identicalExons, exonsToIgnore)

.determineAStypeOverlap <- function(transcript1, transcript2, overlappingExons, identicalExons, exonsToIgnore) {
############################# Inital (nessesary) data analysis  ##############################
  numberOfExons1 <- nrow(transcript1) # get number of exons in transcript 1
  numberOfExons2 <- nrow(transcript2) # get number of exons in transcript 2
  numberOfExons <- c(numberOfExons1, numberOfExons2) # Combine into vecotr
  transcriptList <- list(transcript1,transcript2) # combine transcripts into a list
  ### Determine strand
  strandedness <- unique(c(transcript1$strand, transcript2$strand)) # extract stranness
 
  # Create vector to store AS type findings
  asTypes <- data.frame(ESI=0, MESI=0, ISI=0, A5=0, A3=0, ATSS=0, ATTS=0, ESI.start=NA, ESI.end=NA, MESI.start=NA, MESI.end=NA, ISI.start=NA, ISI.end=NA,
                        A5.start=NA, A5.end=NA, A3.start=NA, A3.end=NA, ATSS.start=NA, ATSS.end=NA, ATTS.start=NA, ATTS.end=NA ,stringsAsFactors=F)
 
 if(length(strandedness)>1) { 
    warning('In pairwise comparison - transcripts are not from the same strand');
    #break
    return(asTypes) }
  isStrandEqualToPlus <- (strandedness == '+') # logic indicating whether the strand is plus
  
   
#   ### check whether the transcripts are completely identical (they should not be)
#   if( numberOfExons1 == numberOfExons2 & nrow(identicalExons) == numberOfExons1 ) { 
#     return(asTypes)
#   } # if the number of exons in each transcrip are the same and all are listed in the identical list
  

  ################################## Test for alternative TSS and TTS  ###############################  
  #Check left side
  if(transcript1$start[1] != transcript2$start[1]) { # if left coordinats are different
    if(isStrandEqualToPlus) {
      asTypes$ATSS <- 1
      asTypeStart <- 'ATSS.start' # string to help assign skipped part to the right type
      asTypeEnd <- 'ATSS.end' # string to help assign skipped part to the right type
    } else {
      asTypes$ATTS <- 1
      asTypeStart <- 'ATTS.start' # string to help assign skipped part to the right type
      asTypeEnd <- 'ATTS.end' # string to help assign skipped part to the right type
    }
    ### Annotate parts spliced out 
    # get indexes
    minIndex <- which.min( c(transcript1$start[1], transcript2$start[1]) )
    notMinIndex <- (1:2)[!1:2 %in% minIndex]
    
    CoordinatsSmallerThanStart <- data.frame(rbind(
      transcriptList[[minIndex]]$start < transcriptList[[notMinIndex]]$start[1],
      transcriptList[[minIndex]]$end < transcriptList[[notMinIndex]]$start[1]
      ))
    
    for(i in 1:ncol(CoordinatsSmallerThanStart)) {
      # if both start and end coordinats are smaller
      if(all(CoordinatsSmallerThanStart[,i])) {
        
        if(is.na(asTypes[asTypeStart])) {
          asTypes[asTypeStart] <- paste( transcriptList[[minIndex]]$start[i] )
          asTypes[asTypeEnd]   <- paste( transcriptList[[minIndex]]$end[i] )
        } else {
          asTypes[asTypeStart] <- paste(asTypes[asTypeStart], transcriptList[[minIndex]]$start[i], sep=';')
          asTypes[asTypeEnd]   <- paste(asTypes[asTypeEnd],   transcriptList[[minIndex]]$end[i],   sep=';')
        }
              
      } else if( any(CoordinatsSmallerThanStart[,i]) ) { # if only start coordinat is smaller
        
        if(is.na(asTypes[asTypeStart])) {
          asTypes[asTypeStart] <- paste( transcriptList[[minIndex]]$start[i] )
          asTypes[asTypeEnd]   <- paste( transcriptList[[notMinIndex]]$start[1] )
        } else {
          asTypes[asTypeStart] <- paste(asTypes[asTypeStart], transcriptList[[minIndex]]$start[i],    sep=';')
          asTypes[asTypeEnd]   <- paste(asTypes[asTypeEnd],   transcriptList[[notMinIndex]]$start[1], sep=';')
        }
        
      } else {
        break # since there are no coordinats smaller thant the start
      }
    }
  }
  
  ### Check rigth side
  if(tail(transcript1$end,1) != tail(transcript2$end,1)) {
    if(isStrandEqualToPlus) {
      asTypes$ATTS <- 1
      asTypeStart <- 'ATTS.start' # string to help assign skipped part to the right type
      asTypeEnd <- 'ATTS.end' # string to help assign skipped part to the right type
    } else {
      asTypes$ATSS <- 1
      asTypeStart <- 'ATSS.start' # string to help assign skipped part to the right type
      asTypeEnd <- 'ATSS.end' # string to help assign skipped part to the right type
    }

    ### Annotate parts spliced out 
    # get indexes
    maxIndex <- which.max( c(tail(transcript1$end,1),  tail(transcript2$end,1)) )
    notMaxIndex <- (1:2)[!1:2 %in% maxIndex]
    
    CoordinatsLargerThanEnd <- data.frame(rbind(
      transcriptList[[maxIndex]]$start > tail(transcriptList[[notMaxIndex]]$end,1),
      transcriptList[[maxIndex]]$end > tail(transcriptList[[notMaxIndex]]$end,1)
    ))
    
    for(i in ncol(CoordinatsLargerThanEnd):1) { # loop over them backwards since it is here the action is
      # if both start and end coordinats are larger
      if(all(CoordinatsLargerThanEnd[,i])) {
        
        if(is.na(asTypes[asTypeStart])) {
          asTypes[asTypeStart] <- paste( transcriptList[[maxIndex]]$start[i] )
          asTypes[asTypeEnd]   <- paste( transcriptList[[maxIndex]]$end[i] )
        } else {
          asTypes[asTypeStart] <- paste(transcriptList[[maxIndex]]$start[i], asTypes[asTypeStart] , sep=';') # switched because i loop over then backwards
          asTypes[asTypeEnd]   <- paste(transcriptList[[maxIndex]]$end[i], asTypes[asTypeEnd],   sep=';') # switched because i loop over then backwards
        }
                
      } else if( any(CoordinatsLargerThanEnd[,i]) ) { # if only the end coordinat is larger
        
        if(is.na(asTypes[asTypeStart])) {
          asTypes[asTypeStart] <- paste( tail(transcriptList[[notMaxIndex]]$end,1) )
          asTypes[asTypeEnd]   <- paste( transcriptList[[maxIndex]]$end[i] )
        } else {
          asTypes[asTypeStart] <- paste( tail(transcriptList[[notMaxIndex]]$end,1), asTypes[asTypeStart], sep=';') # switched because i loop over then backwards
          asTypes[asTypeEnd]   <- paste( transcriptList[[maxIndex]]$end[i],         asTypes[asTypeEnd],   sep=';') # switched because i loop over then backwards
        }
        
      } else {
        break # since there are no coordinats smaller thant the start
      }            
    }
  }
  

  ### Vector of the same length as the number of overlapping exons containing logicis indicating whether an Intron retion have occured 
  duplicatedExons <-(duplicated(overlappingExons$isoform1) | duplicated(overlappingExons$isoform1, fromLast=T )) | (duplicated(overlappingExons$isoform2) | duplicated(overlappingExons$isoform2, fromLast=T )) # Logic indicating whether I find one exon in one trancript overlapping two or more exons in the other transcript)

  # Make an index of the number of exons skipped based on index of overlapping exons - this index contains the number of skipping events between the overlapping exons (0 equals no skipping event)
  exonSkippingIndex <- data.frame(transcript1=diff(c(0,overlappingExons$isoform1,numberOfExons1+1)), transcript2=diff(c(0,overlappingExons$isoform2,numberOfExons2+1)) ) -1
  numberOfSkippingComparasons <- nrow(exonSkippingIndex)


  #################################### Analyze overlapping exons  ###############################
  t1 <- Sys.time()
  
  c1 <- 1 # counter to use in looping
  while(c1 <= nrow(overlappingExons)) {
    if(!duplicatedExons[c1]) { # if the exons analyzed are not part of a intron retension
      # Single exon comparason
      if( any(apply(identicalExons,MARGIN=1,function(x) all(overlappingExons[c1,] %in% x ))) ) { # dont do anything if the exons are identical (meaning found in the identical list)
        c1 <- c1 + 1 # add one to the counter to go to the next level  
        next 
      } else {
        # analyze left and rigth
        #asTypes[,c('A5','A3')] <- asTypes[,c('A5','A3')] + .C("determineLeftOverlappingAStype", as.integer(transcript1[overlappingExons[c1,1],c(2,3)]), as.integer(transcript2[overlappingExons[c1,2],c(2,3)]), as.integer(isStrandEqualToPlus), as.integer(overlappingExons[c1,]), as.integer(c(0,1)))[[5]]
        #asTypes[,c('A5','A3')] <- asTypes[,c('A5','A3')] + .C("determineRightOverlappingAStype", as.integer(transcript1[overlappingExons[c1,1],c(2,3)]), as.integer(transcript2[overlappingExons[c1,2],c(2,3)]), as.integer(isStrandEqualToPlus), as.integer(overlappingExons[c1,]), as.integer(numberOfExons), as.integer(c(0,1)))[[6]]
        asTypes[,c('A5','A3','A5.start','A5.end','A3.start','A3.end')] <- .determineLeftOverlappingAStype(  transcript1[overlappingExons[c1,1],] , transcript2[overlappingExons[c1,2],], isStrandEqualToPlus , overlappingExons[c1,], asTypes[,c('A5','A3','A5.start','A5.end','A3.start','A3.end')]  ) 
        asTypes[,c('A5','A3','A5.start','A5.end','A3.start','A3.end')] <- .determineRightOverlappingAStype(  transcript1[overlappingExons[c1,1],] , transcript2[overlappingExons[c1,2],], isStrandEqualToPlus , overlappingExons[c1,], numberOfExons, asTypes[,c('A5','A3','A5.start','A5.end','A3.start','A3.end')]  )
        c1 <- c1 + 1 # add one to the counter to go to the next level  
      }    
    } else {
      #### intron retension
      # annotate ISI
      asTypes$ISI <- asTypes$ISI + 1
      
      # Get the exon indexe(s) involved in the intron retension (with respect to the transcript not the overlapping list)
      numberOfReplicates1 <- overlappingExons[which(overlappingExons$isoform1 == overlappingExons$isoform1[c1]),'isoform2'] # get the exon(s) from transcript1 for the exon currently under investigation
      numberOfReplicates2 <- overlappingExons[which(overlappingExons$isoform2 == overlappingExons$isoform2[c1]),'isoform1'] # get the  exon(s) from transcript2 for the exon currently under investigation
      numberOfReplicatesList <- list(numberOfReplicates1, numberOfReplicates2)
      # these vectors indicates which exons is the intron retension (the one with length > 1) and the one with two exons
      # therefore these vectors needs to be changed so that numberOfReplicates1 is used with transcript 2 and vice versa
      # e.g. if transcript 1 have a intron retension the overlappingExons look like
      # isoform1   isoform2
      #   1            1
      #   1            2
      # Then the numberOfReplicates1 will have length 2, whereas numberOfReplicates2 will have length 1
      # Since i want isoform1 exon 1 to be compared to both isoform2 exon1 and exon2 I will have to use numberOfReplicates1 with transcript2 and vice versa
      
      # test for 3' overhang in the intron retension
      if( all( c(numberOfReplicates1[1], numberOfReplicates2[1]) != c(1,1) )) { # only test if non of the exons are first exon
        if( transcript1[numberOfReplicates2[1],'start'] != transcript2[numberOfReplicates1[1],'start'] ) {
          minIndex <- which.min( c(transcript1[numberOfReplicates2[1],'start'], transcript2[numberOfReplicates1[1],'start']) )
          notMinIndex <- (1:2)[!1:2 %in% minIndex]
          
          if(is.na(asTypes$ISI.start)) {
            asTypes$ISI.start <- paste( transcriptList[[minIndex]][numberOfReplicatesList[[notMinIndex]][1],'start'] )
            asTypes$ISI.end   <- paste( transcriptList[[notMinIndex]][numberOfReplicatesList[[minIndex]][1],'start'] )
          } else {
            asTypes$ISI.start <- paste(asTypes$ISI.start, transcriptList[[minIndex]][numberOfReplicatesList[[notMinIndex]][1],'start'], sep=';')
            asTypes$ISI.end   <- paste(asTypes$ISI.end,   transcriptList[[notMinIndex]][numberOfReplicatesList[[minIndex]][1],'start'], sep=';')
          }
        }
      }
      
      
      # anotate spliced out
      transcriptWithIntronRetension <- which.max( sapply(numberOfReplicatesList, length ) )
      transcriptWithOUTintronRetension <- (1:2)[!1:2 %in% transcriptWithIntronRetension]
      
      for(i in 1:( max(sapply(numberOfReplicatesList, length )) -1) ) {
        
        if(is.na(asTypes$ISI.start)) {
          asTypes$ISI.start <- paste( transcriptList[[transcriptWithOUTintronRetension]]$end[ numberOfReplicatesList[[transcriptWithIntronRetension]][i]] )
          asTypes$ISI.end   <- paste( transcriptList[[transcriptWithOUTintronRetension]]$start[ numberOfReplicatesList[[transcriptWithIntronRetension]][i] +1 ] )
        } else {
          asTypes$ISI.start <- paste(asTypes$ISI.start, transcriptList[[transcriptWithOUTintronRetension]]$end[ numberOfReplicatesList[[transcriptWithIntronRetension]][i]], sep=';')
          asTypes$ISI.end   <- paste(asTypes$ISI.end,   transcriptList[[transcriptWithOUTintronRetension]]$start[ numberOfReplicatesList[[transcriptWithIntronRetension]][i] +1 ],   sep=';')
        }
      }
      
      # test for 5' overhang in the intron retension
      if( all( c(tail(numberOfReplicates2,1), tail(numberOfReplicates1,1)) != numberOfExons) ) { # only test if non of the exons are last exon
        if( transcript1[tail(numberOfReplicates2,1),'end'] != transcript2[tail(numberOfReplicates1,1),'end'] ) {
          minIndex <- which.min( c(transcript1[tail(numberOfReplicates2,1),'end'], transcript2[tail(numberOfReplicates1,1),'end']) )
          notMinIndex <- (1:2)[!1:2 %in% minIndex]
          
          # no need to test whether it is NA - since the intron retension have already been anotated
          asTypes$ISI.start <- paste(asTypes$ISI.start, transcriptList[[minIndex]][tail(numberOfReplicatesList[[notMinIndex]],1),'end'], sep=';')
          asTypes$ISI.end   <- paste(asTypes$ISI.end,   transcriptList[[notMinIndex]][tail(numberOfReplicatesList[[minIndex]],1),'end'], sep=';')
        }
      }
      
      
      c1 <- c1 + max(length(numberOfReplicates1), length(numberOfReplicates2)) # add the number of exons overlapped by the intron retionsion so that in next itteration these are skipped (this is the reason a while loop is used)
    }
  } # end of overlapping exons while loop

  ################################ Analyze non-overlapping exons  ###############################
  ## A data.frame to help retrive the exons corresponding to the skipping indexes
  overlappingExons2 <- rbind(c(0,0), overlappingExons, (numberOfExons+1))
  colnames(overlappingExons2) <- c('isoform1','isoform2')

  for(i in 1:numberOfSkippingComparasons) {
    if(sum(exonSkippingIndex[i,]) == 0) { # if no exons were skipped
      next
    } 
    else if(sum(exonSkippingIndex[i,]) == 1) {
      if(i == 1 | i == numberOfSkippingComparasons) {
        next # the ATSS/ATTS have already been found - they cannot be ommitted because they might contain more than alternative TSS/TTS
      } else {
        ### Check whether the exon skipping is part of the exons to ignore
        # check whith of the transcripts the single scipping occured in
        skippingInTranscript <- which(exonSkippingIndex[i,] == 1)
        # Extract the exon number of the exon to be skipped
        localExonSkippingIndex <- overlappingExons2[i:(i+1),skippingInTranscript]
        LocalExonToSkip <- seq(localExonSkippingIndex[1]+1,localExonSkippingIndex[2]-1) # extract the exon(s) index that should be compared
        # Check whether the exon to be skipped is part of those that should be ignored
        if(LocalExonToSkip %in% exonsToIgnore[[skippingInTranscript]]) {
          next
        } else {
          asTypes$ESI <- asTypes$ESI + 1
          
          if(is.na(asTypes$ESI.start)) {
            asTypes$ESI.start <- paste( transcriptList[[skippingInTranscript]]$start[LocalExonToSkip] )
            asTypes$ESI.end   <- paste( transcriptList[[skippingInTranscript]]$end[LocalExonToSkip] )
          } else {
            asTypes$ESI.start <- paste(asTypes$ESI.start, transcriptList[[skippingInTranscript]]$start[LocalExonToSkip],sep=';')
            asTypes$ESI.end   <- paste(asTypes$ESI.end,   transcriptList[[skippingInTranscript]]$end[LocalExonToSkip],  sep=';')
          } 
        }
      } # end of if not start or end
    } else if(sum(exonSkippingIndex[i,]) > 1) { # if the number of exons skipped is larger than 1 (nessesary to evaluate since intron retension cause negative numbers in the index)
      if(i == 1 | i == numberOfSkippingComparasons) { # if the skipping is in the first or last step
        # Make sure it is not just one of the trasnscripts having many exons in the alternative TSS or TTS
        if(exonSkippingIndex$transcript1[i] > 0 & exonSkippingIndex$transcript2[i] > 0) {
          asTypes[,c('ESI','MESI','ESI.start','ESI.end','MESI.start','MESI.end')]  <- .determineNonOverlappingAStype(transcript1, transcript2,numberOfExons , overlappingExons2$isoform1[i:(i+1)], overlappingExons2$isoform2[i:(i+1)], exonsToIgnore, asTypes[,c('ESI','MESI','ESI.start','ESI.end','MESI.start','MESI.end')] )
        } else {
          next # continue since ATSS/ATTS have already been found
        }
      } else { # if its not the first or last step
        asTypes[,c('ESI','MESI','ESI.start','ESI.end','MESI.start','MESI.end')]  <- .determineNonOverlappingAStype(transcript1, transcript2, numberOfExons,  overlappingExons2$isoform1[i:(i+1)], overlappingExons2$isoform2[i:(i+1)], exonsToIgnore, asTypes[,c('ESI','MESI','ESI.start','ESI.end','MESI.start','MESI.end')] )
      } 
    }
  } # end of non-overlapping loop
  ################################ Finish up################################
  
  return(asTypes)
}