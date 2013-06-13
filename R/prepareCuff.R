prepareCuff <- function(cuffDB)
{

	message("Reading cuffDB, isoforms...")

	# Get gene and isoform pointers
	cuffGenes 			<- genes(cuffDB)
	cuffIsoforms 		<- isoforms(cuffDB)
	
	# Get isoform annotation and gene diff analysis; merge
	isoformAnnotation 			<- data.frame(annotation(cuffIsoforms),stringsAsFactors=F)[,-10]
	geneDiffanalysis 			<- data.frame(diffData(cuffGenes),stringsAsFactors=F)[,-8]
	colnames(geneDiffanalysis) 	<- c('gene_id','sample_1', 'sample_2' ,unlist(lapply(colnames(geneDiffanalysis[,-1:-3]), function(x) paste("gene_",x,sep="")))) # add gene to the colnames so they can be destinquished from the isoform diff data
	isoformData 				<- data.frame(merge(isoformAnnotation, geneDiffanalysis, by='gene_id'),stringsAsFactors=F)
	
	# Get isoform diff analysis and merge with the rest
	isoformDiffanalysis 			<- data.frame(diffData(cuffIsoforms),stringsAsFactors=F)[,-8]
	colnames(isoformDiffanalysis) 	<- c('isoform_id','sample_1', 'sample_2' ,unlist(lapply(colnames(isoformDiffanalysis[,-1:-3]), function(x) paste("iso_",x,sep="")))) # add gene to the colnames so they can be destinquished from the gene diff data
	isoformData 					<- data.frame(merge(isoformData, isoformDiffanalysis, by=c('isoform_id','sample_1','sample_2')), stringsAsFactors=F) #CORR JW
	
	message("Reading cuffDB, exons...")
	# Get exon info
	isoformFeatureQuery		<-paste("SELECT y.* FROM features y JOIN genes x on y.gene_id = x.gene_id ",sep="")
	isoformFeatures			<-data.frame(dbGetQuery(cuffDB@DB,isoformFeatureQuery),stringsAsFactors=F)
	
	# Extract conditions info from cufflinks
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
	isoformFeatures <- isoformFeatures[,c("seqnames", "start", "end", "strand", "isoform_id", "gene_id")]
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