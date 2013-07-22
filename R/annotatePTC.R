annotatePTC <- function(transcriptData, cds, genomeObject, PTCDistance=50, filters, expressionCutoff=0)
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


	#Get only + and minus strand
	transcriptData[["exon_features"]] <- transcriptData[["exon_features"]][match(strand(transcriptData[["exon_features"]]),c("+", "-"), nomatch=0)>0]
	

	#cuffDB_spliceR[["exon_features"]] <- cuffDB_spliceR[["exon_features"]][ match(strand(cuffDB_spliceR[["exon_features"]]),c("+", "-"))>0]

	# Correct chr names from Ensamble
	if(length(grep('chr',as.character(seqnames(transcriptData[["exon_features"]]))))== 0)
	{
		seqlevels(transcriptData[["exon_features"]]) <- unique(paste("chr",seqnames(transcriptData[["exon_features"]]), sep=""))
		seqlevels(transcriptData[["transcript_features"]]) <- unique(paste("chr",seqnames(transcriptData[["transcript_features"]]), sep=""))
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
	transcriptDF 	<- unFactor(GenomicRanges::as.data.frame(transcriptData[["transcript_features"]]))
	exonDF 			<- unFactor(GenomicRanges::as.data.frame(transcriptData[["exon_features"]]))
	exonDF$"seq"	<- as.character(exonSeq)

	#Split exons
	exonDFSplit 	<- split(exonDF, exonDF$"spliceR.isoform_id")

	# Variables
	cdsPos 			<- rep(-1, length(transcriptData[["transcript_features"]]))	# ATG position
	stopPos 		<- rep(NA, length(transcriptData[["transcript_features"]]))	# Stop codon position (in aa)
	stopDistance 	<- rep(NA, length(transcriptData[["transcript_features"]]))	# Distance from first stop to downstream junction
	junctionIndex	<- rep(NA, length(transcriptData[["transcript_features"]]))	# Downstream junction index. 0 = last, -1 = second last...
	noExons			<- rep(NA, length(transcriptData[["transcript_features"]]))	# No of exons of transcript

	#Initialize progress bar
	pb <- txtProgressBar(min = 0, max = 100, style = 3)


	#Make sure only to loop over unique transcript ID's

	uniqueIsoformIds <- unique(transcriptData[[1]]$"spliceR.isoform_id")

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

		if (currentExonRows[1,"strand"] == "+")
		{
			qualifiedCDS <- cds[cds$"chrom"==currentExonRows[1, "seqnames"]&cds$"strand"==currentExonRows[1, "strand"] & cds$"cdsStart">= min(currentExonRows$"start"-1) & cds$"cdsStart"<= max(currentExonRows$"end"+1),]

			if (!is.null(qualifiedCDS))
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
					cdsPos[i] <- min(result$"cdsStart")
					break
				}
			}
		} else if(currentExonRows[1,"strand"] == "+")
		{
			qualifiedCDS <- cds[cds$"chrom"==currentExonRows[1, "seqnames"]&cds$"strand"==currentExonRows[1, "strand"] & cds$"cdsEnd">= min(currentExonRows$"start"-1) & cds$"cdsEnd"<= max(currentExonRows$"end"+1),]
			if (!is.null(qualifiedCDS))
			{
				for (j in nrow(currentExonRows		):1)
				{
					subset(qualifiedCDS,
						chrom		== currentExonRows		[j, "seqnames"] &
						cdsEnd 		>= currentExonRows		[j,"start"]-1 & 
						cdsEnd	 	<= currentExonRows		[j,"end"]+1 &
						strand		== currentExonRows		[j, "strand"]
					) -> result
					if (nrow(result))
					{
						cdsPos[i] <- max(result$"cdsEnd")
						break
					}
				}
			}
		}	

		if (cdsPos[i][1] != -1)
		{
			exonCumsum	 	<- cumsum(c(0,(currentExonRows$"end"-currentExonRows$"start"+1)))
			exonRevCumsum	<- cumsum(c(0,rev(currentExonRows$"end"-currentExonRows$"start"+1)))
			if (currentExonRows[1,"strand"] == "+")
			{
				# Cumulative exon sizes
				#CDS exon
				cdsExon		 	<- max(which(cdsPos[i][1] >= (currentExonRows$"start")-1))
				cdsOffset		<- cdsPos[i][1]+1 - currentExonRows[cdsExon,"start"] + exonCumsum[cdsExon]+1

				dnaSequencePre	<- paste(currentExonSeq	, collapse="")
				dnaSequence 	<- DNAString(substr(dnaSequencePre ,cdsOffset, nchar(dnaSequencePre)))

				# Break if any position has N 
				if (length(grep("N", dnaSequence, perl=T))>0) next

				# Translate sequence
				suppressWarnings(proteinSequence 	<- as.character(translate(dnaSequence)))
				stopPos[i] 							<- regexpr("\\*", proteinSequence)[1]
				if (stopPos[i][1] == -1) stopPos[i]	<- nchar(proteinSequence)

				if(!max(stopPos[i][1]<=exonCumsum))	{stopPos[i] <- max(exonCumsum)}
				stopPos[i]							<- (stopPos[i]*3)+cdsOffset-1

				#stopDistance[i]						<- exonCumsum[min(which((stopPos[i])<=exonCumsum))]-(stopPos[i])
				stopDistance[i]						<- exonCumsum[length(exonCumsum)-1]-(stopPos[i][1])
				junctionIndex[i]					<- (min(which((stopPos[i][1])<=exonCumsum)))-length(exonCumsum)
			} else
			{
				cdsExon							 	<- min(which(cdsPos[i][1] <= (currentExonRows$"end")+1))
				cdsOffset							<- currentExonRows[cdsExon,"end"] -  cdsPos[i][1]+1 + exonRevCumsum[length(exonCumsum)-cdsExon] +1

				dnaSequencePre						<- paste(rev(currentExonSeq	), collapse="")
				dnaSequence 						<- DNAString(substr(dnaSequencePre,cdsOffset-1, nchar(dnaSequencePre)))

				# Break if any position has N 
				if (length(grep("N", dnaSequence, perl=T))>0) next

	 			suppressWarnings(proteinSequence 	<- as.character(translate(dnaSequence)))
				stopPos[i] 							<- regexpr("\\*", proteinSequence)[1]
				if (stopPos[i][1] == -1) stopPos[i]	<- nchar(proteinSequence)
				stopPos[i]							<- (stopPos[i]*3)+cdsOffset-1
				if(!max(stopPos[i][1]<=exonCumsum))	{stopPos[i] <- max(exonRevCumsum)}

				# stopDistance[i]						<- exonRevCumsum[min(which((stopPos[i])<=exonRevCumsum))]-(stopPos[i])
				stopDistance[i]						<- exonCumsum[length(exonCumsum)-1]-(stopPos[i][1])
				junctionIndex[i]					<- (min(which((stopPos[i][1])<=exonRevCumsum)))-length(exonRevCumsum)	
			}
		}	#End of cds positive
	}	#End of iteration

	message("\nProcessed ", sum(cdsPos[!duplicated(transcriptData[[1]]$"spliceR.isoform_id")]!=-1, na.rm=T) , " putative ORFs...", sep="")

	transcriptData[["transcript_features"]]$"spliceR.cds_pos"       		              	<-cdsPos 		
	transcriptData[["transcript_features"]]$"spliceR.stop_pos"      		              	<-stopPos		
	transcriptData[["transcript_features"]]$"spliceR.stop_distance" 		<-stopDistance	
	transcriptData[["transcript_features"]]$"spliceR.junction_index"	<-junctionIndex
	transcriptData[["transcript_features"]]$"spliceR.PTC"           				<-(stopDistance>=PTCDistance) & junctionIndex != 0

	return(transcriptData)
}