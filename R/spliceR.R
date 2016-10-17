spliceR <- function(transcriptData, compareTo, filters, expressionCutoff=0, useProgressBar=T)
{
    startTime <- Sys.time()
    
    # Check class and GRanges
    if (!class(transcriptData)[1]=="SpliceRList") stop("transcriptData argument is not of class SpliceRList")
    if ( class(transcriptData$"transcript_features") != "GRanges" || class(transcriptData$"exon_features") != "GRanges" ) stop("transcriptData must have GRanges objects in slots 'transcript_features' and 'exon_features'") 
    
    # Validate required columns in spliceRList
    t_colNames <- colnames(mcols(transcriptData$"transcript_features"))
    if(!all(c(
        "isoform_id", "sample_1", "sample_2", "gene_id", "iso_value_1", "iso_value_2", "iso_q_value") %in% substr(t_colNames, 9, nchar(t_colNames))
    )
    ) stop("Transcript features GRanges not compatible with spliceR - see documentation for more details")
    
    e_colNames <- colnames(mcols(transcriptData$"exon_features"))
    if(!all(c(
        "isoform_id","gene_id") %in% substr(e_colNames, 9, nchar(e_colNames))
    )
    ) stop("Exon features GRanges not compatible with spliceR - see documentation for more details")
    
    # check correct user input
    conditionNames <- transcriptData[['conditions']]
    
    if(!compareTo %in% c('preTranscript', conditionNames)) {
        stop(paste('Error in determining compareTo, must be one of: \'preTranscript\' \'', paste(conditionNames, collapse="', '"), "\'. See ?spliceR for more information", sep='') )
    }
    
    dataOrigin <- transcriptData[["source_id"]]
    
    if(! dataOrigin %in% c('cufflinks', 'granges')  ) {
        stop('The input data was not recogniced, please see ?SpliceRList for more information about the input files')
    }
    
    # Check if the filters supplied are OK:
    if(dataOrigin == 'cufflinks') { okFilters <- c('none','expressedGenes','geneOK', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass', 'sigIso', 'singleExon') }
    if(dataOrigin == 'granges') { okFilters <- c('none', 'SingleExon') } # ok since it forces the user to acknowledge that no filters are used
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
    
    
    message("Preparing transcript data...")
    # Create placeholder rows
    transcriptData$"transcript_features"$"spliceR.major"=NA
    
    transcriptData$"transcript_features"$"spliceR.IF1"=NA
    transcriptData$"transcript_features"$"spliceR.IF2"=NA
    transcriptData$"transcript_features"$"spliceR.dIF"=NA
    
    transcriptData$"transcript_features"$"spliceR.ESI"=NA
    transcriptData$"transcript_features"$"spliceR.MEE"=NA
    transcriptData$"transcript_features"$"spliceR.MESI"=NA
    transcriptData$"transcript_features"$"spliceR.ISI"=NA
    transcriptData$"transcript_features"$"spliceR.A5"=NA
    transcriptData$"transcript_features"$"spliceR.A3"=NA
    transcriptData$"transcript_features"$"spliceR.ATSS"=NA
    transcriptData$"transcript_features"$"spliceR.ATTS"=NA
    
    transcriptData$"transcript_features"$"spliceR.analyzed"='no'
    
    transcriptData$"transcript_features"$"spliceR.ESI.start"=NA
    transcriptData$"transcript_features"$"spliceR.ESI.end"=NA
    transcriptData$"transcript_features"$"spliceR.MEE.start"=NA
    transcriptData$"transcript_features"$"spliceR.MEE.end"=NA
    transcriptData$"transcript_features"$"spliceR.MESI.start"=NA
    transcriptData$"transcript_features"$"spliceR.MESI.end"=NA
    transcriptData$"transcript_features"$"spliceR.ISI.start"=NA
    transcriptData$"transcript_features"$"spliceR.ISI.end"=NA
    transcriptData$"transcript_features"$"spliceR.A5.start"=NA
    transcriptData$"transcript_features"$"spliceR.A5.end"=NA
    transcriptData$"transcript_features"$"spliceR.A3.start"=NA
    transcriptData$"transcript_features"$"spliceR.A3.end"=NA
    transcriptData$"transcript_features"$"spliceR.ATSS.start"=NA
    transcriptData$"transcript_features"$"spliceR.ATSS.end"=NA
    transcriptData$"transcript_features"$"spliceR.ATTS.start"=NA
    transcriptData$"transcript_features"$"spliceR.ATTS.end"=NA
    
    #Create backup spliceRList with GRanges before converting to dataframes for output
    originalTranscriptData <- transcriptData
    
    message("Converting to internal objects...")
    
    
    #Convert Granges to dataframe
    tempDF <- GenomicRanges::as.data.frame(transcriptData[["transcript_features"]])
    tempDF <- data.frame(lapply(tempDF, function(x) {if (class(x)=="factor") as.character(x) else (x)}), stringsAsFactors=FALSE) # remove factors
    colnames(tempDF) <- c(colnames(tempDF)[1:5], substr(colnames(tempDF)[6:ncol(tempDF)],9,nchar(colnames(tempDF)[6:ncol(tempDF)])))
    #remove columns not needed here
    transcriptData[["transcript_features"]] <- tempDF
    
    tempDF <- GenomicRanges::as.data.frame(transcriptData[["exon_features"]])
    tempDF <- data.frame(lapply(tempDF, function(x) {if (class(x)=="factor") as.character(x) else (x)}), stringsAsFactors=FALSE) # remove factors
    colnames(tempDF) <- c(colnames(tempDF)[1:5], substr(colnames(tempDF)[6:ncol(tempDF)],9,nchar(colnames(tempDF)[6:ncol(tempDF)])))
    tempDF <- tempDF[,c("start", "end", "strand", "isoform_id")]
    transcriptData[["exon_features"]] <- tempDF
    
    rm(tempDF)
    
    message(length(unique(transcriptData[["transcript_features"]]$isoform_id)), " isoforms pre-filtering...")
    
    message("Filtering...")
    
    
    ### Filter transcript info
    isoformsToAnalyzeIndex <- 1:nrow(transcriptData[["transcript_features"]])
    
    # Optional filters 
    if('geneOK'         %in% filters) { isoformsToAnalyzeIndex <- .filterOKGenes(           transcriptData, isoformsToAnalyzeIndex) }
    if('expressedGenes' %in% filters) { isoformsToAnalyzeIndex <- .filterExpressedGenes(    transcriptData, isoformsToAnalyzeIndex, expressionCutoff) }
    if('sigGenes'       %in% filters) { isoformsToAnalyzeIndex <- .filterSigGenes(          transcriptData, isoformsToAnalyzeIndex) }
    if('isoOK'          %in% filters) { isoformsToAnalyzeIndex <- .filterOKIso(             transcriptData, isoformsToAnalyzeIndex) }
    if('expressedIso'   %in% filters) { isoformsToAnalyzeIndex <- .filterExpressedIso(      transcriptData, isoformsToAnalyzeIndex, expressionCutoff) }
    if('isoClass'       %in% filters) { isoformsToAnalyzeIndex <- .filterIsoClassCode(      transcriptData, isoformsToAnalyzeIndex) }
    if('sigIso'         %in% filters) { isoformsToAnalyzeIndex <- .filterSigIso(            transcriptData, isoformsToAnalyzeIndex) }  
    if('singleExon'     %in% filters) { isoformsToAnalyzeIndex <- .filterSingleExonIsoAll(  transcriptData, isoformsToAnalyzeIndex) }
    if('PTC'            %in% filters) { isoformsToAnalyzeIndex <- .filterPTC(               transcriptData, isoformsToAnalyzeIndex) }
    # Mandatory filters
    isoformsToAnalyzeIndex <- .filterSingleIsoform(transcriptData, isoformsToAnalyzeIndex, conditionNames) # Remove genes with only one isoform left
    
    
    message(length(unique(transcriptData[["transcript_features"]]$isoform_id[isoformsToAnalyzeIndex])), " isoforms post-filtering...")
    
    ### Extract unique gene name
    geneIDs <- unique(transcriptData[["transcript_features"]]$gene_id[isoformsToAnalyzeIndex])
    numberOfGenes <- length(geneIDs)
    
    ### Extract gene names of all genes that ought to be analyzed
    geneIdsToAnalyze  <- transcriptData[["transcript_features"]]$gene_id[isoformsToAnalyzeIndex]
    
    message("Preparing exons...")
    
    ### Split exon info
    # It is fasters to split on isoform id than on genID when scaling up, probably because no exoninfo not used is extracted.
    isoformIDs <- unique(transcriptData[["transcript_features"]]$isoform_id[isoformsToAnalyzeIndex])
    temp <- transcriptData[["exon_features"]][which(transcriptData[["exon_features"]]$isoform_id %in% isoformIDs),]
    isoformFeaturesSplit <- split(temp, f=temp$"isoform_id")
    rm(temp)
    
    ## Determine whether major is chosen or not
    if(compareTo == 'preTranscript') { 
        # A logic indicating whether preTranscrip comparason or major is chosen
        major <- FALSE
    } else {
        ### Find major isoform if that option is toggeled
        major <- TRUE
    }
    
    
    message('Analyzing transcripts...')
    
    # Create statusbar (this statement also automaticlly prints the statusbar)
    if (useProgressBar) pb <- txtProgressBar(min = 1, max = numberOfGenes, style = 3)
    
    for(geneIndex in 1:numberOfGenes) {
        #################### Extract indexs of isoforms belonging to the genes ####################
        ### extract information about the gene
        isoformsToAnalyzeWithinGeneIndexGlobal <- isoformsToAnalyzeIndex[which(geneIdsToAnalyze == geneIDs[geneIndex])] # get the global indexes for the gene analyzed now    #Indexing moved outside of loop 
        isoformsToAnalyze <- transcriptData[["transcript_features"]][isoformsToAnalyzeWithinGeneIndexGlobal,] # extract info about the gene analyzed now    #OBS, IMPROVE SPEED
        
        # extract unique isoforms indexes
        
        uniqIsoformNames <- unique(isoformsToAnalyze$isoform_id)
        uniqueIsoformsIndex <- match(uniqIsoformNames, isoformsToAnalyze$isoform_id)
        
        ### annotate which have been analyzed
        isoformsToAnalyze$analyzed <- 'yes'
        isoformsToAnalyze$MEE <- 0 # these are not nessesarely overwritten else
        
        #0.08
        
        
        #### Determine whether major or preTranscript
        ## Extract all exons to analyze belonging to this gene
        exonList <- isoformFeaturesSplit[ uniqIsoformNames ] # this list is used several times
        ### Create preTranscript
        
        ################################ SLOW ######################
        # !!! Improved !!! KVS
        allUniqueExons <- unique(do.call(rbind,exonList)[,c('start','end','strand')])
        if(nrow(allUniqueExons) > 1) {
            exonInfoPreTranscript <- .getPreRNA( allUniqueExons ) # only pass unique coordinates to the function
        } else {
            next # for the special occation where a gene with multiple identical transcripts with only one exon
        }
        
        ################################ SLOW ######################
        
        ### see whether the minimum requirements for MEE are there (to speed up the calculations)
        areMEEposible <- all( length(which(sapply(exonList, nrow) >= 3)) >= 2 , nrow(exonInfoPreTranscript) >=4)
        
        
        if(major) { ## Check if major is toggeled
            # ###Determine which isoform is major
            # if(dataOrigin == 'cufflinks') {
            # maxIsoformIndex <- .getMajorIsoCuffDB(isoformsToAnalyze, compareTo)
            # } else {
            maxIsoformIndex <- .getMajorIsoCuffDB(isoformsToAnalyze, compareTo)
            # }
            if(length(maxIsoformIndex) == 0) { next } # since it means that no transcript from refrence sample is expressed in for this gene 
            
            # make sure all rows with the transcripts are annotated included in the maxIsoformIndex (nessesary when having more than two samples since else there will be rows with the isoform, but not containing the refrence sample)
            maxIsoformIndex <- which(isoformsToAnalyze$isoform_id == isoformsToAnalyze$isoform_id[maxIsoformIndex[1]])
            
            # extract exon info from the major isoform
            majorExonInfo <-  exonList[[ isoformsToAnalyze$isoform_id[maxIsoformIndex[1]] ]]
            
            # annotate which is major and which is not
            isoformsToAnalyze[,'major'] <- 'no'
            isoformsToAnalyze[maxIsoformIndex,'major'] <- 'yes'
        }
        #############################################################
        ####################### Compare isoforms ####################
        #############################################################
        
        ################# Check for MEE #################
        ### Create an empty list with the exons to ignore for each of the isoforms
        exonsToIgnoreList <- lapply(1:length(uniqueIsoformsIndex),function(x) NULL)
        
        if(areMEEposible) { # if not major Generate list to store the overlapping info from
            ### Create a dataframe to indicate whether the exons are included or not - used to find mutually exclusive exons
            exonIncluded <- data.frame(matrix(0, ncol=length(uniqueIsoformsIndex), nrow=nrow(exonInfoPreTranscript)))
            colnames(exonIncluded) <- uniqIsoformNames
            
            if(!major) { 
                overlapList <- list(NULL) # it is faster to store the overlaps in a list than to redo them - even for small transcript
                identicalList <- list(NULL) # it is faster to store the overlaps in a list than to redo them - even for small transcript
            }
            
            # Loop over genes and extract info of which exons are expressed in which transcripts
            for(i in 1:length(uniqueIsoformsIndex)) {
                ## extract exon features of minor isoform
                isoformExonInfo <-  exonList[[ uniqIsoformNames[i] ]]
                
                overlapIdenticalList <- .findOverlap(exonInfoPreTranscript,isoformExonInfo)
                
                if(!major) { 
                    ## save the overlap tables so I dont need to create them again (they are not used again if major is choseccn)
                    overlapList[[i]] <- overlapIdenticalList$overlap
                    identicalList[[i]] <- overlapIdenticalList$idenctial
                }
                
                exonIncluded[overlapIdenticalList$overlap$isoform1,i] <- overlapIdenticalList$overlap$isoform1
                # Add collum with expressed info to the data.frame
            } # end of loop over unique isoforms
            
            # Change to binary table
            exonIncludedTF <- apply((exonIncluded > 0),2,function(x) as.integer(x))
            
            # Determine MEE based on the exons included table pair else they are ends
            startExons <- apply(exonIncludedTF,2,function(x) match(1, x))
            endExons <- nrow(exonIncludedTF) + 1 - apply(exonIncludedTF,2,function(x) match(1, rev(x)))
            
            for(i in 1:(nrow(exonIncludedTF)-1)) {
                if( sum(exonIncludedTF[i,]) == 1 &  sum(exonIncludedTF[i+1,]) == 1 ) {
                    expressedIn1 <- which(as.logical(exonIncludedTF[i,])) # get the transcript index
                    expressedIn2 <- which(as.logical(exonIncludedTF[i+1,])) # get the transcript index
                    if( expressedIn1 != expressedIn2 ) {
                        if(i != startExons[expressedIn1] & i != endExons[expressedIn1] & i+1 != startExons[expressedIn2] & i+1 != endExons[expressedIn2]) {
                            # annotate the number of MEE
                            MEEindexes1 <- which(isoformsToAnalyze$isoform_id == isoformsToAnalyze$isoform_id[uniqueIsoformsIndex[expressedIn1]])
                            MEEindexes2 <- which(isoformsToAnalyze$isoform_id == isoformsToAnalyze$isoform_id[uniqueIsoformsIndex[expressedIn2]])
                            isoformsToAnalyze$MEE[c(MEEindexes1,MEEindexes2)] <- isoformsToAnalyze$MEE[c(MEEindexes1,MEEindexes2)]  + 1
                            
                            # add the exons to the ignore list
                            exonsToIgnoreList[[expressedIn2]] <-  c(exonsToIgnoreList[[expressedIn2]], exonIncluded[i,expressedIn1]) #exonsToIgnoreList[[expressedIn2]] is used instead of expressedIn1 since I want the exon not expressed in this isoform
                            exonsToIgnoreList[[expressedIn1]] <-  c(exonsToIgnoreList[[expressedIn1]], exonIncluded[i+1,expressedIn2])
                            
                            # annotate the positions of MEE
                            if(is.na(isoformsToAnalyze$MEE.start[MEEindexes1[1]])) {
                                isoformsToAnalyze$MEE.start[MEEindexes1] <- paste(exonInfoPreTranscript$start[i])
                                isoformsToAnalyze$MEE.end[MEEindexes1] <- paste(exonInfoPreTranscript$end[i])
                            } else {
                                isoformsToAnalyze$MEE.start[MEEindexes1] <- paste(isoformsToAnalyze$MEE.start[MEEindexes1],exonInfoPreTranscript$start[i],sep=';')
                                isoformsToAnalyze$MEE.end[MEEindexes1] <- paste(isoformsToAnalyze$MEE.end[MEEindexes1],exonInfoPreTranscript$end[i],sep=';')
                            }
                            if(is.na(isoformsToAnalyze$MEE.end[MEEindexes2[1]])) {
                                isoformsToAnalyze$MEE.start[MEEindexes2] <- paste(exonInfoPreTranscript$start[i+1])
                                isoformsToAnalyze$MEE.end[MEEindexes2] <- paste(exonInfoPreTranscript$end[i+1])
                            } else {
                                isoformsToAnalyze$MEE.start[MEEindexes2] <- paste(isoformsToAnalyze$MEE.start[MEEindexes2],exonInfoPreTranscript$start[i+1],sep=';')
                                isoformsToAnalyze$MEE.end[MEEindexes2] <- paste(isoformsToAnalyze$MEE.end[MEEindexes2],exonInfoPreTranscript$end[i+1],sep=';')
                            } 
                        }
                    }
                }
            }
            
            ## Make sure every exon is only reppresented once (a special case where an exon is envolved in two MEE events)
            exonsToIgnoreList <- lapply(exonsToIgnoreList, function(x) unique(x))
        } # end of is MEE posible
        
        ### Determine which of the indexes is major (only nessesary if any exons ought to be skipped)
        if(major) {
            majorIndex <- which(uniqIsoformNames %in% isoformsToAnalyze$isoform_id[maxIsoformIndex[1]]) # get the index of major
        }
        
        ######## Classify the rest of the AS types ########
        # loop over the unique isoforms to make the AS type comparason
        for(i in 1:length(uniqueIsoformsIndex)) {
            
            if(major) { # if I compare to major
                if(uniqueIsoformsIndex[i] %in% maxIsoformIndex) {next} # if this isoform is the major
            }
            
            ## extract exon features of minor isoform
            isoformExonInfo <-  exonList[[ uniqIsoformNames[i] ]]
            # Get indexes for all rows containing this transcript so they can all be annotated (many samples == many rows)
            thisIsoformIndex <- which(isoformsToAnalyze$isoform_id == isoformsToAnalyze$isoform_id[uniqueIsoformsIndex[i]])
            
            ### Determine AS classification and overlap
            if(!major) { # if pre-transcript
                ## Determine which exons to ignore
                exonsToIgnore <- list(exonsToIgnoreList[[i]], NULL) # NULL since I know no exons is skipped in pre-transcripted, switched since the skipping is going to occure in the OPPOSITE transcript
                ### annotate all instances of this isoform with the Alternative splicing found
                if(areMEEposible) {
                    isoformsToAnalyze[thisIsoformIndex, c(
                        'ESI','MESI','ISI','A5','A3','ATSS','ATTS',
                        'ESI.start', 'ESI.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end',
                        'A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end'
                    )]  <- .determineAStypeOverlap(exonInfoPreTranscript,isoformExonInfo,overlapList[[i]],identicalList[[i]], exonsToIgnore)
                } else {
                    overlapListLocal <- .findOverlap(exonInfoPreTranscript,isoformExonInfo)
                    isoformsToAnalyze[thisIsoformIndex, c(
                        'ESI','MESI','ISI','A5','A3','ATSS','ATTS',
                        'ESI.start', 'ESI.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end',
                        'A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end'
                    )]  <- .determineAStypeOverlap(exonInfoPreTranscript,isoformExonInfo,overlapListLocal[[1]],overlapListLocal[[2]], exonsToIgnore)
                }
            } else { # if major
                ### Determine which exons to ignore
                exonsToIgnore <- list(exonsToIgnoreList[[i]], exonsToIgnoreList[[majorIndex]]) # switched since the skipping is going to occure in the OPPOSITE transcript
                
                ### Determine overlapping exons between major and minor
                temp <- .findOverlap(majorExonInfo,isoformExonInfo)
                
                ### annotate all instances of this isoform with the Alternative splicing found
                isoformsToAnalyze[thisIsoformIndex, c(
                    'ESI','MESI','ISI','A5','A3','ATSS','ATTS',
                    'ESI.start', 'ESI.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end',
                    'A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end'
                )]  <- .determineAStypeOverlap(majorExonInfo,isoformExonInfo,temp[[1]],temp[[2]], exonsToIgnore)
            }     
        } # end of loop over isoforms
        
        ### Anotate IF and dIF values
        # get total expression of all the isoforms to analyze for EACH condition (is different than the expression of the gene - because i exclude some transcripts)
        totalIsoformExpression <- NULL
        for(conditionName in conditionNames) {
            maxOFthisIsoform <- sum( 
                isoformsToAnalyze$iso_value_1[which(isoformsToAnalyze$sample_1 == conditionName)], 
                isoformsToAnalyze$iso_value_2[which(isoformsToAnalyze$sample_2 == conditionName)] 
            )
            totalIsoformExpression <- c(totalIsoformExpression , maxOFthisIsoform)
        }
        
        # Extract info about which collums are the ones that contain the wanted info
        sampleCol1 <- which(colnames(transcriptData[["transcript_features"]])=="sample_1")   	#spliceR.sample_1
        isoValCol1 <- which(colnames(transcriptData[["transcript_features"]])=="iso_value_1")	#spliceR.iso_value_1
        IFvalCol1  <- which(colnames(transcriptData[["transcript_features"]])=="IF1")        	#spliceR.IF1
        
        # print(cat(colnames(transcriptData[["transcript_features"]])))
        # print(cat(isoValCol1))
        # print(cat(PSIvalCol1))
        
        
        # if(dataOrigin == 'cufflinks') { # in this way we can controle the different indexes for different input files
        #   # sampleCol1 <- 7                                                                             	#spliceR.sample_1
        #   # isoValCol1 <- 24                                                                            	#spliceR.iso_value_1
        #   # PSIvalCol1 <- 31                                                                            	#spliceR.PSI1
        #   sampleCol1 <- which(colnames(transcriptData[["transcript_features"]]))=="spliceR.sample_1")   	#spliceR.sample_1
        #   isoValCol1 <- which(colnames(transcriptData[["transcript_features"]]))=="spliceR.iso_value_1")	#spliceR.iso_value_1
        #   PSIvalCol1 <- which(colnames(transcriptData[["transcript_features"]]))=="spliceR.PSI1")       	#spliceR.PSI1
        # }
        # if(dataOrigin == 'granges') { # in this way we can controle the different indexes for different input files
        #   sampleCol1 <- NULL
        #   isoValCol1 <- NULL
        #   PSIvalCol1 <- NULL
        # }
        
        # 	    # loop over all indexes to analyze to calculte IF values
        # 	    for(isoformIndex in 1:nrow(isoformsToAnalyze)) {
        # 	      # if isoform is major annotate it
        # 	      if(major) {
        # 	        if(isoformIndex %in% maxIsoformIndex) { next }
        # 	      } 
        #         
        #         for(i in 0:1) { #loop over indexes so i can calculate IF for both the sample_1 and sample_2 collumns
        # 	        # Get condition name
        # 	        myCondition <- isoformsToAnalyze[isoformIndex,(sampleCol1+i)] # get condition name (which is in collumn 2 and 3)
        # 	          # get total expression of isoforms within that gene
        # 	        totalExpValue <- totalIsoformExpression[conditionNames %in% myCondition]
        # 	        if(totalExpValue == 0) {
        # 	          isoformsToAnalyze[isoformIndex,(IFvalCol1+i)] <- 0 # else I would devide by zero
        # 	        } else {
        # 	          # calculate IF value
        # 	          isoformsToAnalyze[isoformIndex,(IFvalCol1+i)] <- round( isoformsToAnalyze[isoformIndex,(isoValCol1+i)] / totalExpValue * 100 ,digits = 2)
        # 	        }
        # 	      }
        # 	    }
        # 	    # annotate dIF
        # 	    isoformsToAnalyze$dIF <- isoformsToAnalyze$IF2 - isoformsToAnalyze$IF1
        
        # write local data to global dataframe (so everything is stored and can be returned)    
        # this is faster than replacing the full dataset and also faster (and more readiable) than using c(26:40)
        transcriptData[["transcript_features"]][isoformsToAnalyzeWithinGeneIndexGlobal,c("major","IF1","IF2","dIF","ESI","MEE","MESI","ISI","A5","A3","ATSS","ATTS","analyzed",'ESI.start', 'ESI.end','MEE.start','MEE.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end','A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end')] = isoformsToAnalyze[,c("major","IF1","IF2","dIF","ESI","MEE","MESI","ISI","A5","A3","ATSS","ATTS","analyzed",'ESI.start', 'ESI.end','MEE.start','MEE.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end','A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end')] #slow index - THE RATE LIMITING STEP
        #paste(difftime(Sys.time(),t10,u='sec'),'Time to write to global file',sep=' ')
        
        ### Update progressbar
        if (useProgressBar) setTxtProgressBar(pb, geneIndex)
    } # belongs to loop over genes
    
    ### Annotate IF values
    transcriptData$transcript_features$IF1[isoformsToAnalyzeIndex] <- round(transcriptData$transcript_features$iso_value_1[isoformsToAnalyzeIndex] / transcriptData$transcript_features$gene_value_1[isoformsToAnalyzeIndex] * 100, digits=4)
    transcriptData$transcript_features$IF2[isoformsToAnalyzeIndex] <- round(transcriptData$transcript_features$iso_value_2[isoformsToAnalyzeIndex] / transcriptData$transcript_features$gene_value_2[isoformsToAnalyzeIndex] * 100, digits=4)
    transcriptData$transcript_features$dIF[isoformsToAnalyzeIndex] <- transcriptData$transcript_features$IF2[isoformsToAnalyzeIndex] - transcriptData$transcript_features$IF1[isoformsToAnalyzeIndex]
    
    #close progress bar
    if (useProgressBar) close(pb)
    
    message('Preparing output...')
    
    ori_col_names <- colnames(mcols(originalTranscriptData[[1]]))
    ori_col_names_no_spliceR <- substr(ori_col_names, 9, nchar(ori_col_names))
    for (i in 1:length(ori_col_names))
    {
        mcols(originalTranscriptData[[1]])[ori_col_names[i]] <- transcriptData[[1]][,ori_col_names_no_spliceR[i]]
    }
    
    #Add filters to spliceR object
    originalTranscriptData[['filter_params']] <- filters
    
    endTime <- Sys.time()
    message('Done in ', format(difftime(endTime,startTime), digits=2))
    
    # return data list to give back all annotation
    return(originalTranscriptData)
    
    #return GRanges
}