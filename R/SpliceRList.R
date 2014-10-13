SpliceRList <- function(transcript_features, exon_features, assembly_id, source_id, conditions,transcripts_plot=NULL,filter_params=NULL) 
{
	x <- new(
		"SpliceRList",
		list(
			transcript_features=transcript_features,
			exon_features=exon_features,
			assembly_id=assembly_id,
			source_id=source_id,
			conditions=conditions,
			transcripts_plot=transcripts_plot,
			filter_params=filter_params
		)
	)
	return(x)
}
