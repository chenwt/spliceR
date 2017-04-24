CDSSet <- function(cds) 
{
	x <- new(
		"CDSSet",
		as.data.frame(cds)
	)
	return(x)
}
