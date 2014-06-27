prepareCuff <- function(cuffDB, fixCufflinksAnnotationProblem=TRUE, removeNonCanonical=TRUE)
{
    if(!is.logical(fixCufflinksAnnotationProblem)) { stop("The fixCufflinksAnnotationProblem parameter must be set to either TRUE or FALSE, indicating whether to try to correct the Cufflinks annotation problem") }
    ### Get isoform and gene info
	message("Reading cuffDB, isoforms...")

	# Get gene and isoform pointers
	cuffGenes 			<- genes(cuffDB)
	cuffIsoforms 		<- isoforms(cuffDB)
	
	### Get isoform annotation
	isoformAnnotation 			<- data.frame(annotation(cuffIsoforms),stringsAsFactors=F)
	# Unique + removeal of colums is nessesary of the new cummeRbund devel (where these colums are included making duplication of rows)
    isoformAnnotation          	<- unique( isoformAnnotation[,-which( colnames(isoformAnnotation) %in% c('start','end','width','exon_number'))])
    
    # Get gene diff analysis
	geneDiffanalysis 			<- data.frame(diffData(cuffGenes),stringsAsFactors=F)[,-8]
	colnames(geneDiffanalysis) 	<- c('gene_id','sample_1', 'sample_2' ,unlist(lapply(colnames(geneDiffanalysis[,-1:-3]), function(x) paste("gene_",x,sep="")))) # add gene to the colnames so they can be destinquished from the isoform diff data
	
    # Merge data
    isoformData 				<- data.frame(merge(isoformAnnotation, geneDiffanalysis, by='gene_id'),stringsAsFactors=F)
	
	# Get isoform diff analysis and merge with the rest
	isoformDiffanalysis 			<- data.frame(diffData(cuffIsoforms),stringsAsFactors=F)[,-8]
	colnames(isoformDiffanalysis) 	<- c('isoform_id','sample_1', 'sample_2' ,unlist(lapply(colnames(isoformDiffanalysis[,-1:-3]), function(x) paste("iso_",x,sep="")))) # add gene to the colnames so they can be destinquished from the gene diff data
	isoformData 					<- data.frame(merge(isoformData, isoformDiffanalysis, by=c('isoform_id','sample_1','sample_2')), stringsAsFactors=F) #CORR JW 
    
    ### Get exon info
    message("Reading cuffDB, exons...")
    isoformFeatureQuery		<-paste("SELECT y.* FROM features y JOIN genes x on y.gene_id = x.gene_id ",sep="")
    isoformFeatures			<-data.frame(dbGetQuery(cuffDB@DB,isoformFeatureQuery),stringsAsFactors=F)
    # Another faster way of doing it would be to use the colums removed from isoformAnnotation, but this approach is not backwards compatible
    # isoformFeatures <- isoformAnnotation[,c("seqnames", "start", "end", "strand", "isoform_id", "gene_id")]
	
	### Remove unknown chromosomes
    if( removeNonCanonical ) {
    	isoformData <- isoformData[grep('^[1-9mc]', isoformData$locus, ignore.case=T, perl=T),] # only that those that starts with either a number or c or m (meaning random ect are removed - they cause problems in annotatePTC)
    	isoformFeatures <- isoformFeatures[grep('^[1-9mc]', isoformFeatures$seqnames, ignore.case=T, perl=T),] # only that those that starts with either a number or c or m (meaning random ect are removed - they cause problems in annotatePTC)
    }
	
	### Make sure none of the data.frames only contain isoforms not found in the other data.frame
	isoformData <- isoformData[which(isoformData$isoform_id %in% isoformFeatures$isoform_id),]
	isoformFeatures <- isoformFeatures[which(isoformFeatures$isoform_id %in% isoformData$isoform_id),]
    
	### Fix to correct for Cufflinks annotation problem where cufflinks assignes transcripts from several annotated genes to 1 cuffgene (this solution does not correct for gene expression or differential analysis)
    if(fixCufflinksAnnotationProblem) {
        message("Analyzing cufflinks annotation problem...")
        
        geneName <- unique(isoformData[,c('gene_id','gene_short_name')])
        geneNameSplit <- split(geneName$gene_short_name , f=geneName$gene_id)
        # remove all unique
        geneNameSplit <- geneNameSplit[which(sapply(geneNameSplit, length) > 1)]
        
        if(length(geneNameSplit) > 0 ) { # if there are any problems
            message("Fixing cufflinks annotation problem...")
            #get indexes of those affected
            geneNameIndexesData     <- which(isoformData$gene_id %in% names(geneNameSplit))
            geneNameIndexesFeatures <- which(isoformFeatures$gene_id %in% names(geneNameSplit))
            # combine names of cuffgenes and 
            isoformData$gene_id[geneNameIndexesData]            <- paste(isoformData$gene_id[geneNameIndexesData]        , isoformData$gene_short_name[geneNameIndexesData],     sep=':')
            isoformFeatures$gene_id[geneNameIndexesFeatures]    <- paste(isoformFeatures$gene_id[geneNameIndexesFeatures], isoformFeatures$gene_name[geneNameIndexesFeatures],   sep=':')
            
            ## Correct gene expression levels and differntial analysis
            problematicGenes <- isoformData[geneNameIndexesData, c('isoform_id','gene_id','sample_1','sample_2','gene_value_1','gene_value_2','gene_log2_fold_change','gene_p_value','gene_q_value','gene_significant','iso_value_1','iso_value_2')]
            problematicGenesSplit <- split(problematicGenes, f=problematicGenes[,c('gene_id','sample_1','sample_2')], drop=T )
            
            correctGeneExpression <- function(df) {
                df$gene_value_1 <- sum(df$iso_value_1)
                df$gene_value_2 <- sum(df$iso_value_2)
                df$gene_log2_fold_change <- log2( (df$gene_value_2[2]+0.0001)/(df$gene_value_1[1]+0.0001) )
                df$gene_p_value <- 1
                df$gene_q_value <- 1
                df$gene_significant <- 'no'
                return(df)
            }
            correctedGenes <- ldply(problematicGenesSplit, correctGeneExpression)
            
            # sort so genes end up being in correct order for overwriting
            correctedGenes <- correctedGenes[order(correctedGenes$isoform_id, correctedGenes$gene_id, correctedGenes$sample_1, correctedGenes$sample_2),-1] # -1 removes the index created by ldply
            # overwrite problematic genes
            isoformData[geneNameIndexesData, c('gene_id','sample_1','sample_2','gene_value_1','gene_value_2','gene_log2_fold_change','gene_p_value','gene_q_value','gene_significant','iso_value_1','iso_value_2')] <- correctedGenes[,-1] # -1 removes the id that ldply uses
                        
            message( paste("Cufflinks annotation problem was fixed for", length(geneNameSplit), "Cuff_genes", sep=' ') )
        } else {
            message( "No instances of a Cufflinks annotation problem found - no changes were made" )
        }
    } # end of fix cufflinks annotatopn problem
	message("Creating spliceRList...")
	    
	### Extract conditions info from cufflinks
	conditions <- samples(genes(cuffDB))

	# Create GRanges for transcript features
	isoformDataSplit 	<- strsplit(isoformData$"locus", split=":")
	isoformChr 			<- unlist(lapply(isoformDataSplit, FUN=function(x){x[1]}))
	posSplit			<- strsplit(unlist(lapply(isoformDataSplit, FUN=function(x){x[2]})), split="-")
	isoformStart		<- as.integer(lapply(posSplit, function(x){x[1]}))
	isoformEnd			<- as.integer(lapply(posSplit, function(x){x[2]}))

	transcriptFeatures <- GRanges(
		seqnames=isoformChr,
		ranges=IRanges(
			start=isoformStart,
			end=isoformEnd),
		spliceR=isoformData[,-10]		#add spliceR to metadata colnames
	)

	# Create GRanges for exon features
	exonFeatures <- GRanges(
		seqnames=isoformFeatures$"seqnames",
		strand=isoformFeatures$"strand",
		ranges=IRanges(
			start=isoformFeatures$"start",
			end=isoformFeatures$"end"),
		spliceR=isoformFeatures[,c("isoform_id", "gene_id")]		#add spliceR to metadata colnames
	)

	if (length(exonFeatures)==0) stop("No exon information extracted from GTF")

	# Return SpliceRList
	return(
	    new("SpliceRList",
			list(
				transcript_features=transcriptFeatures,
				exon_features=exonFeatures,
				assembly_id=runInfo(cuffDB)[5,2],
				source_id="cufflinks",
				conditions=conditions,
				transcripts_plot=NULL,
				filter_params=NULL
			)		
		)
	)
}