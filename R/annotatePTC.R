annotatePTC <- function(transcriptData, cds, genomeObject, PTCDistance=50) 
{
    if (!class(transcriptData)[1]=="SpliceRList") stop("transcriptData argument is not of class SpliceRList")
    if (!class(genomeObject)[1]=="BSgenome") stop("genomeObject argument is not of class BSgenome")
    if (!class(cds)[1]=="CDSSet") stop("cds argument is not of class CDSSet")
    
    message("Initializing....", sep="")
    
    # Remove novel gene isoforms having "*" as strand
    # transcriptData[["transcript_features"]] <- transcriptData[["transcript_features"]][strand(transcriptData[["transcript_features"]])!="*"]
    
    message("Scanning ", length(unique(transcriptData[["transcript_features"]]$"spliceR.isoform_id")) , " transcript models for ORFs....", sep="")
    
    #Get genomic sequence
    message("Fetching exon sequence....", sep="")
    
    #Establish exon features full copy
    exon_features_full <- transcriptData[["exon_features"]]
    
    
    #Get only + and - strands
    #transcriptData[["exon_features"]] <- transcriptData[["exon_features"]][match(strand(transcriptData[["exon_features"]]),c("+", "-"), nomatch=0)>0]
    transcriptData[["exon_features"]] <- transcriptData[["exon_features"]][match(as.vector(strand(transcriptData[["exon_features"]])),c("+", "-"), nomatch = 0) > 0]
    # sort (because some annotation database (such as Biocondocturs annotation packages) give exons in oposite order on '-' strand, which will cause problems with our transcript sequence assembly)
    transcriptData[["exon_features"]] <- sort(transcriptData[["exon_features"]])
    
    #cuffDB_spliceR[["exon_features"]] <- cuffDB_spliceR[["exon_features"]][ match(strand(cuffDB_spliceR[["exon_features"]]),c("+", "-"))>0]
    
    # Correct chr names from Ensamble
    if(length(grep('chr',as.character(seqnames(transcriptData[["exon_features"]]))))== 0)
    {
        # Correct all chromosomes
        seqlevels(transcriptData[["exon_features"]]) <- unique(paste("chr",seqnames(transcriptData[["exon_features"]]), sep=""))
        seqlevels(transcriptData[["transcript_features"]]) <- unique(paste("chr",seqnames(transcriptData[["transcript_features"]]), sep=""))
        # Correct Mitochondria
        seqlevels(transcriptData[["exon_features"]])<- sub('chrMT','chrM',seqlevels(transcriptData[["exon_features"]]))
        seqlevels(transcriptData[["transcript_features"]])<- sub('chrMT','chrM',seqlevels(transcriptData[["transcript_features"]]))   
    }
    
    exonSeq  <- getSeq(genomeObject, transcriptData[["exon_features"]])
    
    unFactor <- function(df)
    {
        classes <- unlist(lapply(df, class))
        for (i in 1:ncol(df))
        {
            if (classes[i] == "factor") 
                df[,i] <- as.character(df[,i])
        }
        df
    }
    transcriptDF     <- unFactor(GenomicRanges::as.data.frame(transcriptData[["transcript_features"]]))
    exonDF 			<- unFactor(GenomicRanges::as.data.frame(transcriptData[["exon_features"]]))
    exonDF$"seq"	<- as.character(exonSeq)
    
    #Split exons
    exonDFSplit 	<- split(exonDF, exonDF$"spliceR.isoform_id", drop=T)
    
    # Variables
    cdsPosGenomic	<- rep(-1, length(transcriptData[["transcript_features"]]))	# ATG position in genomic coordinats
    cdsPosTranscipt <- rep(-1, length(transcriptData[["transcript_features"]])) # ATG position in comparison to transcript length
    cdsPosExon      <- rep(-1, length(transcriptData[["transcript_features"]])) # ATG is in exon
    stopPosGenomic	<- rep(NA, length(transcriptData[["transcript_features"]]))	# Stop codon position (in aa) in genomic coordinats
    stopPosTranscipt<- rep(NA, length(transcriptData[["transcript_features"]]))	# Stop codon position (in aa) compared to transcript length
    stopPosExon     <- rep(NA, length(transcriptData[["transcript_features"]]))  # Stop codon position (in aa) compared to transcript length
    stopDistance 	<- rep(NA, length(transcriptData[["transcript_features"]]))	# Distance from first stop to downstream junction
    junctionIndex	<- rep(NA, length(transcriptData[["transcript_features"]]))	# Downstream junction index. 0 = last, -1 = second last...
    noExons			<- rep(NA, length(transcriptData[["transcript_features"]]))	# No of exons of transcript
    
    #Make sure only to loop over unique transcript ID's
    uniqueIsoformIds <- unique(transcriptData[[1]]$"spliceR.isoform_id")
    
    # Message
    message("Analyzing....")
    
    #Initialize progress bar
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    
    #Main loop
    for(ii in 1:length(uniqueIsoformIds))
    {
        setTxtProgressBar(pb, floor(ii/length(uniqueIsoformIds)*100))
        
        i 							<- which(transcriptDF$"spliceR.isoform_id"==uniqueIsoformIds[ii])
        # currentCuff	 			<- transcriptDF[i,]
        currentCuff					<- transcriptDF[transcriptDF$"spliceR.isoform_id"==uniqueIsoformIds[ii],][1,]
        currentExonRows		 		<- exonDFSplit[[currentCuff$"spliceR.isoform_id"]]
        currentExonSeq				<- currentExonRows$"seq"
        
        
        
        if (is.null(currentExonRows)) next	#For transcript models without exons (could be removed because of "*" strand)
        
        noExons	[i] <- nrow(currentExonRows)
        
        # Next, if one exon gene
        if (nrow(currentExonRows		) == 1) { junctionIndex[i] <- -1;next }
        
        # To avoid CMD Build NOTEs
        cdsStart = cdsEnd = NULL
        
        # Identify the most upstream annotated start codon
        if (currentExonRows[1,"strand"] == "+")
        {
            # Find compatible annotated start codons that are within the transcript - pre-filtering
            qualifiedCDS <- cds[cds$"chrom"==currentExonRows[1, "seqnames"]     &
                                    cds$"strand"==currentExonRows[1, "strand"]      & 
                                    cds$"cdsStart">= min(currentExonRows$"start"-1) & 
                                    cds$"cdsStart"<= max(currentExonRows$"end"+1),]
            
            # Find compatible annotated start codons that are within the exons of the transcript - final-filtering
            if (!is.null(qualifiedCDS)) 
            {
                for (j in 1:nrow(currentExonRows))
                {
                    subset(qualifiedCDS,
                           chrom		== currentExonRows[j, "seqnames"] &
                               cdsStart 	>= currentExonRows[j,"start"]-1 & 
                               cdsStart 	<= currentExonRows[j,"end"]+1 &
                               strand		== currentExonRows[j, "strand"]
                    ) -> result
                    if (nrow(result))
                    {
                        cdsPosGenomic[i] <- min(result$"cdsStart") +1 # KVS corrected +1 to get to 1 based system
                        break
                    }
                }
            }
            
        } else if(currentExonRows[1,"strand"] == "-")
        {
            # Find compatible annotated start codons that are within the transcript - pre-filtering
            qualifiedCDS <- cds[cds$"chrom"==currentExonRows[1, "seqnames"] &
                                    cds$"strand"==currentExonRows[1, "strand"]  & 
                                    cds$"cdsEnd">= min(currentExonRows$"start"-1) & 
                                    cds$"cdsEnd"<= max(currentExonRows$"end"+1),]
            
            # Find compatible annotated start codons that are within the exons of the transcript - final-filtering
            if (!is.null(qualifiedCDS))
            {
                for (j in nrow(currentExonRows):1)
                {
                    subset(qualifiedCDS,
                           chrom		== currentExonRows		[j, "seqnames"] &
                               cdsEnd 		>= currentExonRows		[j,"start"]-1 & 
                               cdsEnd	 	<= currentExonRows		[j,"end"]+1 &
                               strand		== currentExonRows		[j, "strand"]
                    ) -> result
                    if (nrow(result))
                    {
                        cdsPosGenomic[i] <- max(result$"cdsEnd")
                        break
                    }
                }
            }
        } # end of finding the most upstream annotated start codon
        
        # Identify stop codon in translated sequence
        if (cdsPosGenomic[i][1] != -1)
        {
            # Cumulative exon sizes
            exonCumsum         <- cumsum(c(0,(currentExonRows$"end"-currentExonRows$"start"+1)))
            exonRevCumsum    <- cumsum(c(0,rev(currentExonRows$"end"-currentExonRows$"start"+1)))
            
            if (currentExonRows[1,"strand"] == "+")
            {
                
                #CDS exon
                cdsExon		 	<- max(which(cdsPosGenomic[i][1] >= (currentExonRows$"start")-1)) # get exon number
                cdsPosExon[i]   <- cdsExon
                cdsOffset		<- cdsPosGenomic[i][1]+1 - currentExonRows[cdsExon,"start"] + exonCumsum[cdsExon] # KVS corr
                # save the position of the start codon
                cdsPosTranscipt[i] <- cdsOffset
                
                # extract DNA sequence
                dnaSequencePre	<- paste(currentExonSeq	, collapse="")
                dnaSequence 	<- DNAString(substr(dnaSequencePre ,cdsOffset, nchar(dnaSequencePre)))
                
                # jump to next if it is a corrupted DNA stand
                if (length(grep('[^ATCG]', dnaSequence, perl=T))>0) next
                
                # Translate sequence
                suppressWarnings(proteinSequence 	<- as.character(translate(dnaSequence)))
                # find stop codons
                stopPosTranscipt[i]					<- regexpr("\\*", proteinSequence)[1]
                if (stopPosTranscipt[i][1] == -1) stopPosTranscipt[i]	<- nchar(proteinSequence) # if no found set it to end
                stopPosTranscipt[i]					<- (stopPosTranscipt[i]*3)+cdsOffset-1 # translate info to nt
                if(!any(stopPosTranscipt[i][1]<=exonCumsum))    {stopPosTranscipt[i] <- max(exonCumsum)}
                
                stopIsLocatedInExon <- max(which(stopPosTranscipt[i][1] > exonCumsum))
                stopPosExon[i] <- stopIsLocatedInExon
                distanceToLocalExonstart <- stopPosTranscipt[i][1] - exonCumsum[stopIsLocatedInExon]
                stopPosGenomic[i] <- currentExonRows$start[stopIsLocatedInExon] + distanceToLocalExonstart -2 # -3 since the *3 to get the stopPosTranscipt brings us to the end of the stop codon
                
                #stopDistance[i]						<- exonCumsum[min(which((stopPos[i])<=exonCumsum))]-(stopPos[i])
                stopDistance[i]						<- exonCumsum[length(exonCumsum)-1]-(stopPosTranscipt[i][1])
                junctionIndex[i]					<- (min(which((stopPosTranscipt[i][1])<=exonCumsum)))-length(exonCumsum)
                
            } else if (currentExonRows[1,"strand"] == "-")
            {
                cdsExon							 	<- min(which(cdsPosGenomic[i][1] <= (currentExonRows$"end")+1))
                cdsPosExon[i]                       <- (nrow(currentExonRows):1)[ cdsExon ] # KVS corr - the index needs to be turned around (so we count in minus direction) since we compare to exon structure (annotated as if on plus strand)
                cdsOffset							<- currentExonRows[cdsExon,"end"] -  cdsPosGenomic[i][1]+1 + exonRevCumsum[length(exonRevCumsum)-cdsExon] # KVS corr
                
                # save the position of the start codon
                cdsPosTranscipt[i] <- cdsOffset
                
                # extract DNA sequence
                dnaSequencePre						<- paste(rev(currentExonSeq	), collapse="")
                dnaSequence 						<- DNAString(substr(dnaSequencePre,cdsOffset, nchar(dnaSequencePre))) # KVS corr
                
                # jump to next if it is a corrupted DNA stand
                if (length(grep("N", dnaSequence, perl=T))>0) next
                
                suppressWarnings(proteinSequence 	<- as.character(translate(dnaSequence)))
                stopPosTranscipt[i] 							<- regexpr("\\*", proteinSequence)[1]
                if (stopPosTranscipt[i][1] == -1) stopPosTranscipt[i]	<- nchar(proteinSequence)
                stopPosTranscipt[i]							<- (stopPosTranscipt[i]*3)+cdsOffset-1 # -1 because both positions included # KVS corr
                if(!any(stopPosTranscipt[i][1]<=exonCumsum))	{stopPosTranscipt[i] <- max(exonRevCumsum)}
                
                stopIsLocatedInExon <- max(which(stopPosTranscipt[i][1] > exonRevCumsum))
                stopPosExon[i] <- stopIsLocatedInExon
                
                cdsEndExonIndex     <- max(which(stopPosTranscipt[i][1]   >=  exonRevCumsum ))
                
                distanceToLocalExonstart <- stopPosTranscipt[i][1] - exonRevCumsum[stopIsLocatedInExon]
                stopPosGenomic[i] <- currentExonRows$end[ (nrow(currentExonRows):1)[stopIsLocatedInExon] ] - distanceToLocalExonstart +1 # KVS Corr
                
                # stopDistance[i]						<- exonRevCumsum[min(which((stopPos[i])<=exonRevCumsum))]-(stopPos[i])
                stopDistance[i]						<- exonRevCumsum[length(exonRevCumsum)-1]-(stopPosTranscipt[i][1]) # KVS corr
                junctionIndex[i]					<- (min(which((stopPosTranscipt[i][1])<=exonRevCumsum)))-length(exonRevCumsum)	
            }
        }	#End of any compatible cds found
    }	#End of iteration over unique transcripts
    
    message("\nProcessed ", sum(cdsPosGenomic[!duplicated(transcriptData[[1]]$"spliceR.isoform_id")]!=-1, na.rm=T) , " putative ORFs...", sep="")
    
    transcriptData[["transcript_features"]]$"spliceR.cdsPosGenomic"     <- cdsPosGenomic
    transcriptData[["transcript_features"]]$"spliceR.stopPosGenomic"    <- stopPosGenomic
    transcriptData[["transcript_features"]]$"spliceR.ExonWithStart"     <- cdsPosExon
    transcriptData[["transcript_features"]]$"spliceR.ExonWithStop"      <- stopPosExon
    transcriptData[["transcript_features"]]$"spliceR.cdsPosTranscipt" 	<- cdsPosTranscipt
    transcriptData[["transcript_features"]]$"spliceR.stopPosTranscipt" 	<- stopPosTranscipt
    transcriptData[["transcript_features"]]$"spliceR.stopDistance" 		<- stopDistance
    transcriptData[["transcript_features"]]$"spliceR.junctionIndex" 	<- junctionIndex
    transcriptData[["transcript_features"]]$"spliceR.PTC"           	<-(stopDistance>=PTCDistance) & junctionIndex != 0
    
    #restore full copy of exon features
    
    #Establish exon features full copy
    transcriptData[["exon_features"]] <- exon_features_full
    
    return(transcriptData)
}