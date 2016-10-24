# DOCUMENTATION NEEDS TO BE UPDATED

getCDS <- function(selectedGenome, repoName)
{
	# Check classes and params
	if(length(repoName)>1)
		stop("getCDS: Please supply only one repository")

	if(!any(selectedGenome%in%c(
		"hg19", "hg38", "mm9", "mm10"))
	) stop("getCDS: supported genomes are currently: hg19, hg38, mm9, mm10")

	if(!any(repoName%in%c(
		"ensemble", "UCSC", "refseq", "GENCODE"))
	) stop("getCDS: Supported repositories are currently: Ensemble, UCSC, Refseq, GENCODE (hg19/hg38 only)")

	message("Retrieving CDS tables for ", repoName ,"...", sep="")

	repoName <- c("ensGene", "knownGene", "refGene", "wgEncodeGencodeV19")[which(c(
		"ensemble", "UCSC", "refseq", "GENCODE") %in% repoName)]

	# Get CDS
	session      		<- browserSession("UCSC",url="http://genome-euro.ucsc.edu/cgi-bin/")
	genome(session)   	<- selectedGenome
	query         		<- ucscTableQuery(session, repoName)
	
	# Encode table name fix
	if (repoName == "wgEncodeGencodeV19") repoName = "wgEncodeGencodeCompV19"

	tableName(query)  	<- repoName
	cdsTable      		<- getTable(query)

	# Repository specific filtering of transcripts (remove transcripts without CDS)
	if (repoName=="ensGene") cdsTable[cdsTable$"cdsStartStat"!="none",] -> cdsTable
	if (repoName=="knownGene") cdsTable[cdsTable$"cdsStart" != cdsTable$"cdsEnd",] -> cdsTable
	if (repoName=="refGene") cdsTable[cdsTable$"cdsStart" != cdsTable$"cdsEnd",] -> cdsTable

	cdsTable      <- cdsTable[,c("chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "name")]
	message("Retrieved ", nrow(cdsTable) ," records...", sep="")
	flush.console()
	return(new("CDSSet",cdsTable))
}



# UNIT TESTS

# hg19_ens  <- getCDS("hg19","ensemble")
# hg19_ref  <- getCDS("hg19","refseq")
# hg19_ucsc <- getCDS("hg19","UCSC")
# hg19_gen  <- getCDS("hg19","GENCODE")

# hg38_ens  <- getCDS("hg38","ensemble")	#ERROR, track name fault
# hg38_ref  <- getCDS("hg38","refseq")
# hg38_ucsc <- getCDS("hg38","UCSC")
# hg38_gen  <- getCDS("hg38","GENCODE")

# mm9_ens  <- getCDS("mm9","ensemble")
# mm9_ref  <- getCDS("mm9","refseq")
# mm9_ucsc <- getCDS("mm9","UCSC")

# mm10_ens  <- getCDS("mm10","ensemble")
# mm10_ref  <- getCDS("mm10","refseq")
# mm10_ucsc <- getCDS("mm10","UCSC")
