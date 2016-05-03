require(methods)

setClass("SpliceRList",
	representation("list")
)

setClass("CDSSet",
	representation("data.frame")
)


dim.SpliceRList <- function(x) {
	length(x[[1]])
}

length.SpliceRList <- function(x) {
	length(x[[1]])
}

