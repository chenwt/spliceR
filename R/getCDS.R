getCDS <- function(selectedGenome, repoName)
{
	# Check classes and params
	if(length(repoName)>1)
		stop("getCDS: Please supply only one repository")

	if(!any(selectedGenome%in%c(
		"hg19", "mm9"))
	) stop("getCDS: supported genomes are currently: hg19, mm9")

	if(!any(repoName%in%c(
		"ensemble", "UCSC", "refseq"))
	) stop("getCDS: Supported repositories are currently: Ensemble, UCSC, Refseq")

	message("Retrieving CDS tables for ", repoName ,"...", sep="")

	repoName <- c("ensGene", "knownGene", "refGene")[which(c(
		"ensemble", "UCSC", "refseq", "gencode")%in%repoName)]

	# Get CDS
	session      		<- browserSession("UCSC")
	genome(session)   	<- selectedGenome
	query         		<- ucscTableQuery(session, repoName)
	
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
