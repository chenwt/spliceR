################ Use the score table to extract the numbers that can be passed to venn diagram functions #####################
.getIntersectVennData <- function(scoreDataFrame, conditionsToAnalyze , evaluate='nr_transcript', asType='All', expressionCutoff=0) {
  # notes:
  # created with inspiration from the venn.diagram() function. 
  # This function have been thoroughly tested and gives the same results as venn.diagram()
  
  # The idea is that I create a vector with the size of the regions analyzed. These can the be used to overwrite the 
  # numbers in a dummy venn diagram. By calculateing the region sizes myself (instead of just using venn.diagram with 
  # a list) I can weight each of the inputs for example by expression level or number of AS events (or the combination)
  
  # This allows me to create highly costumized Venn plots. E.g number of AS evenst / nr_Transcript
    
  ### Manipulate the score dataframe to the use I need
  expressed <- scoreDataFrame[,c(3:(length(conditionsToAnalyze)+2))] > expressionCutoff

  ########################### Function of how to ecaluate the numbers ##########################
  if(evaluate == 'nr_transcript') { # function to return number of transcripts
    evaluateFunction <- function(expressedTranscripts, scoreDataFrame, asType, index) {
      return(sum(expressedTranscripts))
    }
  } else if( evaluate == 'nr_genes') {
    evaluateFunction <- function(expressedTranscripts, scoreDataFrame, asType, index) { 
      return(length(unique(scoreDataFrame$gene_id[which(expressedTranscripts)])))
    }
  } else if(evaluate == 'nr_AS') { # function to return number of AS events
    evaluateFunction <- function(expressedTranscripts, scoreDataFrame, asType, index) {
      return(sum(scoreDataFrame[,asType][expressedTranscripts]))
    }
  } else if(evaluate == 'exp_mass') { # function to return the fpkm mass
    evaluateFunction <- function(expressedTranscripts, scoreDataFrame, asType, index) {
      if(length(index) == 1) {
        return(sum(          log2(scoreDataFrame[,(index+2)]+1)[expressedTranscripts]  ))
      } else {
        return(sum( rowMeans(log2(scoreDataFrame[,(index+2)]+1)) [expressedTranscripts] ))
      }
    }
  }
  
  # Function to evaluate whether at a transcript is expressed in some samples but not in other samples
  evaluateExpressedIN <- function(myRow, expressed_in, not_expressed_in) {
    
    expressedIn <- all(myRow[expressed_in]) # check that the isoform is expressed in these conditions
    
    if(!is.null(not_expressed_in)) {
      notExpressed <- !any(myRow[not_expressed_in]) # check that the isoform is NOT expressed in these conditions
    } else {
      notExpressed <- NULL
    }
    return(all(c(expressedIn,notExpressed))) # combine the two statements into one
  }
  
  ############################# Get numbers for Venn diagrams #############################
  ### pairwise Venn diagram (2)
  
  if(length(conditionsToAnalyze)==2) { 
    A   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1),   not_expressed_in=c(2))    ), scoreDataFrame, asType , c(1)   )
    B   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2),   not_expressed_in=c(1))    ), scoreDataFrame, asType , c(2)   )
    n12 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2), not_expressed_in=c(NULL)) ), scoreDataFrame, asType , c(1,2) )
        
    vennData <- c( A,B ,n12 ) # same order as draw.pairwise.venn() wants the input
  }
  
  ### triple Venn diagram (3)
  else if(length(conditionsToAnalyze)==3) { 
    A   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1),   not_expressed_in=c(2,3)) ), scoreDataFrame, asType , c(1)   )
    B   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2),   not_expressed_in=c(1,3)) ), scoreDataFrame, asType , c(2)   )
    C   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(3),   not_expressed_in=c(1,2)) ), scoreDataFrame, asType , c(3)   )
    
    n12 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2), not_expressed_in=c(3))   ), scoreDataFrame, asType , c(1,2) )
    n23 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,3), not_expressed_in=c(1))   ), scoreDataFrame, asType , c(2,3) )
    n13 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,3), not_expressed_in=c(2))   ), scoreDataFrame, asType , c(1,3) )
    
    n123 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,3), not_expressed_in=c(NULL)) ), scoreDataFrame, asType , c(1,2,3) )
    
    vennData <- c(A ,B ,C ,n12 ,n23 ,n13 ,n123 ) # same order as draw.triple.venn() wants the input
    # Rearage to the order that the numbers are actually plotted
    vennData <- vennData[c(7,2,1,3,5,6,4)]
   
  }
  
  ### quad venn diagram (4)
  else if(length(conditionsToAnalyze)==4) {
    A   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1),   not_expressed_in=c(2,3,4))    ), scoreDataFrame, asType , c(1)   )
    B   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2),   not_expressed_in=c(1,3,4))    ), scoreDataFrame, asType , c(2)   )
    C   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(3),   not_expressed_in=c(1,2,4))    ), scoreDataFrame, asType , c(3)   )
    D   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(4),   not_expressed_in=c(1,2,3))    ), scoreDataFrame, asType , c(4)   )
    
    n12 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2), not_expressed_in=c(3,4)) ), scoreDataFrame, asType , c(1,2) )
    n13 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,3), not_expressed_in=c(2,4)) ), scoreDataFrame, asType , c(1,3) )
    n14 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,4), not_expressed_in=c(2,3)) ), scoreDataFrame, asType , c(1,4) )
    n23 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,3), not_expressed_in=c(1,4)) ), scoreDataFrame, asType , c(2,3) )
    n24 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,4), not_expressed_in=c(1,3)) ), scoreDataFrame, asType , c(2,4) )
    n34 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(3,4), not_expressed_in=c(1,2)) ), scoreDataFrame, asType , c(3,4) )
    
    n123 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,3), not_expressed_in=c(4)) ), scoreDataFrame, asType , c(1,2,3) )
    n124 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,4), not_expressed_in=c(3)) ), scoreDataFrame, asType , c(1,2,4) )
    n134 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,3,4), not_expressed_in=c(2)) ), scoreDataFrame, asType , c(1,3,4) )
    n234 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,3,4), not_expressed_in=c(1)) ), scoreDataFrame, asType , c(2,3,4) )
    
    n1234 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,3,4), not_expressed_in=c(NULL)) ), scoreDataFrame, asType , c(1,2,3,4) )
    
    vennData <- c(A ,B ,C ,D ,n12 ,n13 ,n14 ,n23 ,n24 ,n34 ,n123 ,n124 ,n134 ,n234 ,n1234 ) # same order as draw.quad.venn() wants the input
    # Rearage to the order that the numbers are actually plotted
    vennData <- vennData[c(3,10,4,6,13,15,14,9,1,7,12,11,8,2,5)]
  }
  
  ### quintuple venn (5)
  else if(length(conditionsToAnalyze)==5) { 
    A     <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1),   not_expressed_in=c(2,3,4,5))    ), scoreDataFrame, asType , c(1)   )
    B     <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2),   not_expressed_in=c(1,3,4,5))    ), scoreDataFrame, asType , c(2)   )
    C     <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(3),   not_expressed_in=c(1,2,4,5))    ), scoreDataFrame, asType , c(3)   )
    D     <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(4),   not_expressed_in=c(1,2,3,5))    ), scoreDataFrame, asType , c(4)   )
    E     <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(5),   not_expressed_in=c(1,2,3,4))    ), scoreDataFrame, asType , c(5)   )
    
    n12   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2), not_expressed_in=c(3,4,5)) ), scoreDataFrame, asType , c(1,2) )
    n13   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,3), not_expressed_in=c(2,4,5)) ), scoreDataFrame, asType , c(1,3) )
    n14   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,4), not_expressed_in=c(2,3,5)) ), scoreDataFrame, asType , c(1,4) )
    n15   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,5), not_expressed_in=c(2,3,4)) ), scoreDataFrame, asType , c(1,5) )
    n23   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,3), not_expressed_in=c(1,4,5)) ), scoreDataFrame, asType , c(2,3) )
    n24   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,4), not_expressed_in=c(1,3,5)) ), scoreDataFrame, asType , c(2,4) )
    n25   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,5), not_expressed_in=c(1,3,4)) ), scoreDataFrame, asType , c(2,5) )
    n34   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(3,4), not_expressed_in=c(1,2,5)) ), scoreDataFrame, asType , c(3,4) )
    n35   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(3,5), not_expressed_in=c(1,2,4)) ), scoreDataFrame, asType , c(3,5) )
    n45   <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(4,5), not_expressed_in=c(1,2,3)) ), scoreDataFrame, asType , c(4,5) )
    
    n123  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,3), not_expressed_in=c(4,5)) ), scoreDataFrame, asType , c(1,2,3) )
    n124  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,4), not_expressed_in=c(3,5)) ), scoreDataFrame, asType , c(1,2,4) )
    n125  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,5), not_expressed_in=c(3,4)) ), scoreDataFrame, asType , c(1,2,5) )
    n134  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,3,4), not_expressed_in=c(2,5)) ), scoreDataFrame, asType , c(1,3,4) )
    n135  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,3,5), not_expressed_in=c(2,4)) ), scoreDataFrame, asType , c(1,3,5) )
    n145  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,4,5), not_expressed_in=c(2,3)) ), scoreDataFrame, asType , c(1,4,5) )
    n234  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,3,4), not_expressed_in=c(1,5)) ), scoreDataFrame, asType , c(2,3,4) )
    n235  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,3,5), not_expressed_in=c(1,4)) ), scoreDataFrame, asType , c(2,3,5) )
    n245  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,4,5), not_expressed_in=c(1,3)) ), scoreDataFrame, asType , c(2,4,5) )
    n345  <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(3,4,5), not_expressed_in=c(1,2)) ), scoreDataFrame, asType , c(3,4,5) )
    
    n1234 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,3,4), not_expressed_in=c(5)) ), scoreDataFrame, asType , c(1,2,3,4) )
    n1235 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,3,5), not_expressed_in=c(4)) ), scoreDataFrame, asType , c(1,2,3,5) )
    n1245 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,4,5), not_expressed_in=c(3)) ), scoreDataFrame, asType , c(1,2,4,5) )
    n1345 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,3,4,5), not_expressed_in=c(2)) ), scoreDataFrame, asType , c(1,3,4,5) )
    n2345 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(2,3,4,5), not_expressed_in=c(1)) ), scoreDataFrame, asType , c(2,3,4,5) )
    
    n12345 <- evaluateFunction( apply(expressed, 1, function(x) evaluateExpressedIN(x, expressed_in=c(1,2,3,4,5), not_expressed_in=c(NULL)) ), scoreDataFrame, asType , c(1,2,3,4,5) )
    
    vennData <- c(A, B, C, D, E, n12, n13, n14, n15, n23, n24, n25, n34, n35, n45, n123, n124, n125, n134, 
                  n135, n145, n234, n235, n245, n345, n1234, n1235, n1245, n1345, n2345, n12345) # same order as draw.quintuple.venn() wants the input
    # Rearage to the order that the numbers are actually plotted
    vennData <- vennData[c(1,2,3,4,5,14,9,8,6,12,10,7,13,11,15,25,20,21,17,18,23,16,19,22,24,30,29,28,27,26,31 )]
  }
  return(vennData)
}

### Function to combine it all
spliceRPlot <- function(spliceRobject, evaluate='nr_transcript', asType='All', colors=NULL, alpha=NULL, reset=FALSE, filters=NULL, expressionCutoff=0) {
  # Check class and GRanges
  if (!class(spliceRobject)[1]=="SpliceRList") stop("spliceRobject argument is not of class SpliceRList")
  if ( class(spliceRobject$"transcript_features") != "GRanges" || class(spliceRobject$"exon_features") != "GRanges" ) stop("spliceRobject must have GRanges objects in slots 'transcript_features' and 'exon_features'") 
  
  # Validate required columns in spliceRList
  t_colNames <- colnames(mcols(spliceRobject$"transcript_features"))
  if(!all(c(
    "isoform_id", "sample_1", "sample_2", "gene_id", "iso_value_1", "iso_value_2", "iso_q_value") %in% substr(t_colNames, 9, nchar(t_colNames))
  )
  ) stop("Transcript features GRanges not compatible with spliceR - see documentation for more details")
  
  # Vailidate that the spliceR object have the AS data
  AScols <- c("spliceR.major","spliceR.IF1","spliceR.IF2","spliceR.dIF","spliceR.ESI","spliceR.MEE","spliceR.MESI","spliceR.ISI","spliceR.A5","spliceR.A3","spliceR.ATSS","spliceR.ATTS","spliceR.analyzed","spliceR.ESI.start","spliceR.ESI.end","spliceR.MEE.start","spliceR.MEE.end","spliceR.MESI.start","spliceR.MESI.end","spliceR.ISI.start","spliceR.ISI.end","spliceR.A5.start","spliceR.A5.end","spliceR.A3.start","spliceR.A3.end","spliceR.ATSS.start","spliceR.ATSS.end","spliceR.ATTS.start","spliceR.ATTS.end")
  if( !all( AScols %in% t_colNames)) {
    stop("SpliceRList has not yet been analyzed. Run spliceR() first...")
  }
  e_colNames <- colnames(mcols(spliceRobject$"exon_features"))
  if(!all(c(
    "isoform_id","gene_id") %in% substr(e_colNames, 9, nchar(e_colNames))
  )
  ) stop("Exon features GRanges not compatible with spliceR - see documentation for more details")
    
  conditionsToAnalyze <- spliceRobject[['conditions']]
  
  ### Test correctness of input
  if(class(evaluate) != "character") {
    stop('Evaluate must be a text string descriping the plot type. Advailable plots are \'nr_transcript\', \'nr_AS\', \'mean_AS\', \'mean_transcript_exp\' and \'weighted_mean_transcript_exp\'')
  }
  if(class(asType) != "character") {
    stop('asType must be a text string descriping the type of alternative splicing to be analyzed. Advailable AS types are \'ESI\', \'MEE\', \'MESI\', \'ISI\', \'A5\', \'A3\', \'ATSS\', \'ATTS\' and \'All\' (sum of All the others) ')
  }
  if(! evaluate %in% c('nr_transcript', 'nr_genes' , 'nr_AS', 'mean_AS_gene', 'mean_AS_transcript', 'mean_transcript_exp', 'mean_gene_exp', 'nr_transcript_pr_gene')) {
    stop('Plot type not recogniced. Advailable plots are \'nr_transcript\',\'nr_genes\' \'nr_transcript_pr_gene\', \'nr_AS\', \'mean_AS_gene\', \'mean_AS_transcript\', \'mean_transcript_exp\',  \'mean_gene_exp\' ')
  }
  if(length(evaluate) != 1) {
    stop('Evaluate can one be ONE of the following \'nr_transcript\',\'nr_genes\', \'nr_transcript_pr_gene\' , \'nr_AS\', \'mean_AS_gene\', \'mean_AS_transcript\', \'mean_transcript_exp\',  \'mean_gene_exp\' ')
  } 
  if( any( !is.null(asType), evaluate=='nr_AS', evaluate=='AS_mass') ) {
    if(! asType %in% c('ESI','MEE','MESI','ISI','A5', 'A3','ATSS','ATTS','All')) {
      stop('AS type type not recogniced. Advailable AS types are \'ESI\', \'MEE\', \'MESI\', \'ISI\', \'A5\', \'A3\', \'ATSS\', \'ATTS\' and \'All\' (sum of All the others)')
    }
  }
  if(!length(conditionsToAnalyze) %in% 2:5) {
    stop(paste('spliceRplot only supports 2,3,4 or 5 conditions. SpliceRplot found',length(conditionsToAnalyze),'conditions',sep=" "))
  }
  if(!is.null(colors)) {
    if(length(conditionsToAnalyze) != length(colors)) {
      stop(paste('Splicerplot needs one color per condition. There are', length(conditionsToAnalyze), 'conditions, but', length(colors), 'colors were supplied', sep=' '))
    }
  }
  if(!is.null(alpha)) {
    if(length(alpha) != 1) {
      stop('The \'alpha\' value must be a vector with length 1')
    }
    if(alpha > 1 | alpha < 0) {
      stop('The \'alpha\' value must be a number between 0 and 1')
    }
  }
  ## check filters
  dataOrigin <- spliceRobject[["source_id"]]
  
  if(! dataOrigin %in% c('cufflinks', 'granges')  ) {
    stop('The input data was not recogniced, please see ?SpliceRList for more information about the input files')
  }
  
  # Check if the filters supplied are OK:
  if(dataOrigin == 'cufflinks') { okFilters <- c('none','expressedGenes','geneOK', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass', 'sigIso', 'singleExon') }
  if(dataOrigin == 'granges') { okFilters <- c('none', 'singleExon') } # ok since it forces the user to acknowledge that no filters are used
  if('PTC' %in% filters) { # if asked to filter on PTC
    if('spliceR.PTC' %in% colnames(as.data.frame(spliceRobject$"transcript_features"[1,]))) { # check whether the spliceR object contain PTC info
      okFilters <- c(okFilters, 'PTC')
    } else {
      stop('spliceR cannot filter on PTC since no PTC info is advailable. PTC information can be obtained through annotatePTC() ')
    }
  }
  if(any(!filters %in% okFilters)) { # if one or more of the supplied filters are not recogniced
    stop('One or more of the supplied filters are not recogniced, please see ?spliceRPlot for more information about the filters')
  }
  
  # check whether the splice data have already been added to the splicer object - if not add it
  if(is.null(spliceRobject[['transcripts_plot']])) {
    reset <- TRUE
  }
  
  if(reset) { 
    message('Initializing (this will only happen once)...')
    
    if(is.null(filters)) { # check wheter filters are supplied
      if(is.null(spliceRobject[['filter_params']])) { # check whether any filters are stored in the spliceR object
        stop('spliceRplot needs data generated by spliceR(), please use spliceR() on the data before using spliceRplot')
      }
      filters <- spliceRobject[['filter_params']]
    }
    
    # Extract transcript_features from spliceRobject
    transcript_features <- GenomicRanges::as.data.frame(spliceRobject[["transcript_features"]])
    transcript_features <- data.frame(lapply(transcript_features, function(x) {if (class(x)=="factor") as.character(x) else (x)}), stringsAsFactors=FALSE)
    colnames(transcript_features) <- c(colnames(transcript_features)[1:5], substr(colnames(transcript_features)[6:ncol(transcript_features)],9,nchar(colnames(transcript_features)[6:ncol(transcript_features)])))
    
    
    isoformsToAnalyzeIndex <- 1:nrow(transcript_features)
    ##################################### Apply the chosen filters #####################################
    # Optional filters 
    if('geneOK'         %in% filters) { isoformsToAnalyzeIndex <- .filterOKGenes(        list("transcript_features"=transcript_features), isoformsToAnalyzeIndex) }
    if('expressedGenes' %in% filters) { isoformsToAnalyzeIndex <- .filterExpressedGenes( list("transcript_features"=transcript_features), isoformsToAnalyzeIndex, expressionCutoff) }
    if('sigGenes'       %in% filters) { isoformsToAnalyzeIndex <- .filterSigGenes(       list("transcript_features"=transcript_features), isoformsToAnalyzeIndex) }
    if('isoOK'          %in% filters) { isoformsToAnalyzeIndex <- .filterOKIso(          list("transcript_features"=transcript_features), isoformsToAnalyzeIndex) }
    if('expressedIso'   %in% filters) { isoformsToAnalyzeIndex <- .filterExpressedIso(   list("transcript_features"=transcript_features), isoformsToAnalyzeIndex, expressionCutoff) }
    if('isoClass'       %in% filters) { isoformsToAnalyzeIndex <- .filterIsoClassCode(   list("transcript_features"=transcript_features), isoformsToAnalyzeIndex) }
    if('sigIso'         %in% filters) { isoformsToAnalyzeIndex <- .filterSigIso(         list("transcript_features"=transcript_features), isoformsToAnalyzeIndex) }  
    if('singleExon'     %in% filters) { 
      exonDF <- GenomicRanges::as.data.frame(spliceRobject[["exon_features"]])
      exonDF <- data.frame(lapply(exonDF, function(x) {if (class(x)=="factor") as.character(x) else (x)}), stringsAsFactors=FALSE) # remove factors
      colnames(exonDF) <- c(colnames(exonDF)[1:5], substr(colnames(exonDF)[6:ncol(exonDF)],9,nchar(colnames(exonDF)[6:ncol(exonDF)])))
      
      isoformsToAnalyzeIndex <- .filterSingleExonIsoAll( list("transcript_features"=transcript_features, "exon_features"=exonDF), isoformsToAnalyzeIndex) 
      rm(exonDF)
    }
    if('PTC'            %in% filters) { isoformsToAnalyzeIndex <- .filterPTC(               list('transcript_features'=transcript_features), isoformsToAnalyzeIndex) }
    
    analyzedIsoformData <- transcript_features[isoformsToAnalyzeIndex,]
    
    ### function remodle the expression df by using apply
    myExtractData <- function(dfForThisIso) {
      myResults <- replicate(length(conditionsToAnalyze), 0)
      
      for(index2 in 1:length(conditionsToAnalyze)) {
        indexOfCondition <- match(conditionsToAnalyze[index2], c(dfForThisIso[,'sample_1'], dfForThisIso[,'sample_2'])) 
        if(is.na(indexOfCondition)) {next}
        myResults[index2] <- c(dfForThisIso[,'iso_value_1'], dfForThisIso[,'iso_value_2'])[indexOfCondition]
        
      }
      return(myResults)
    }
    
    ### Extract a data.frame with transcript in rows and expression values from each condidition in seperat collumns
    splittedIsoform <- split(analyzedIsoformData[,c('sample_1','sample_2','iso_value_1','iso_value_2')], f=analyzedIsoformData$isoform_id)
    myIsoformScore <- ldply(splittedIsoform, function(x) myExtractData(x))
    colnames(myIsoformScore) <- c('isoform_id',conditionsToAnalyze)
    
    ### merge with info from the known table
    myASisoformScores <- merge(myIsoformScore, unique(analyzedIsoformData[,c('isoform_id','gene_id','ESI','MEE','MESI','ISI','A5', 'A3','ATSS','ATTS')]),by ='isoform_id')
    myASisoformScores$All <- apply(myASisoformScores[,c('ESI','MEE','MESI','ISI','A5', 'A3','ATSS','ATTS')],1,sum)
    
    ### Replace NA with 0
    for(i in (ncol(myASisoformScores)-8):ncol(myASisoformScores) ) {
      myASisoformScores[ which(is.na(myASisoformScores[,i])) ,i] <- 0  
    }
    
    ### Rearange collumns
    myASisoformScores <- myASisoformScores[,c('isoform_id','gene_id',conditionsToAnalyze,'ESI','MEE','MESI','ISI','A5', 'A3','ATSS','ATTS','All')]
    
    ### Calculate gene expression
    temp <- split(myASisoformScores, f=myASisoformScores$gene_id)
    # function to use on the splitted data
    getGeneExpression <- function(df) {
        if(nrow(df) == 1 ) {
            return(df)
        } else {
            for(i in 3:(3+length(conditionsToAnalyze)-1)) {
                df[,i]      <- sum(df[,i])
            }
            return(df)
        }
    }
    myGeneExpression <- ldply(temp, getGeneExpression)[,-1]
    
    spliceRobject[['transcripts_plot']] <- list(isoforms=myASisoformScores, genes=myGeneExpression)
  } # end of is.null(cummerSpliceRanalyzed[['plotData']])
  else {
    myASisoformScores <- spliceRobject$transcripts_plot$isoforms
    myGeneExpression <- spliceRobject$transcripts_plot$genes
  }
  

  ### Generate numbers to overwrite those in the dummy venn diagram and headers
  mainText <- 'spliceR venn diagram'
  
  if(evaluate == 'nr_transcript') {
    overwriteValues <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='nr_transcript', asType=asType, expressionCutoff=expressionCutoff)
    subText         <- 'Number of transcipts'
  } else if(evaluate == 'nr_genes' ) {
    overwriteValues <- .getIntersectVennData(myGeneExpression, conditionsToAnalyze , evaluate='nr_genes', asType=asType, expressionCutoff=expressionCutoff)
    subText         <- 'Number of genes'
  } else if(evaluate == 'nr_transcript_pr_gene' ) {
    nrTranscripts   <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='nr_transcript', asType=asType, expressionCutoff=expressionCutoff)
    nrGenes         <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='nr_genes', asType=asType, expressionCutoff=expressionCutoff)
    overwriteValues <- round(nrTranscripts/nrGenes,digits=2)
    subText         <- 'Number of transcript per gene'
  } else if(evaluate == 'nr_AS') {
    overwriteValues <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='nr_AS', asType=asType, expressionCutoff=expressionCutoff)
    subText         <- paste('Number of',asType,'AS events',sep=' ')
  } else if(evaluate == 'mean_AS_gene') {
    nrAS            <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='nr_AS', asType=asType, expressionCutoff=expressionCutoff)
    nrTranscripts   <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='nr_genes', asType=asType, expressionCutoff=expressionCutoff)
    overwriteValues <- round(nrAS/nrTranscripts,digits=2)
    subText         <- paste('Average number of',asType,'AS events per gene',sep=' ')
  } else if(evaluate == 'mean_AS_transcript') {
    nrAS            <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='nr_AS', asType=asType, expressionCutoff=expressionCutoff)
    nrTranscripts   <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='nr_transcript', asType=asType, expressionCutoff=expressionCutoff)
    overwriteValues <- round(nrAS/nrTranscripts,digits=2)
    subText         <- paste('Average number of',asType,'AS events per transcript',sep=' ')
  } else if(evaluate == 'mean_transcript_exp') {
    expressionMass  <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='exp_mass', asType=asType, expressionCutoff=expressionCutoff)
    nrTranscripts   <- .getIntersectVennData(myASisoformScores, conditionsToAnalyze , evaluate='nr_transcript', asType=asType, expressionCutoff=expressionCutoff)
    overwriteValues <- round(expressionMass/nrTranscripts,digits=2)
    subText         <- 'Mean transcript expression'
  } else if(evaluate == 'mean_gene_exp') {
    expressionMass  <- .getIntersectVennData(myGeneExpression, conditionsToAnalyze , evaluate='exp_mass', asType=asType, expressionCutoff=expressionCutoff)
    nrGenes         <- .getIntersectVennData(myGeneExpression, conditionsToAnalyze , evaluate='nr_genes', asType=asType, expressionCutoff=expressionCutoff)
    overwriteValues <- round(expressionMass/nrGenes,digits=2)
    subText <- 'Mean gene expression'
  }
  
  # Correct for deviding with zero
  overwriteValues[which(is.na(overwriteValues))] <- 0
 
  ### Create a dummy venn diagram
  if(length(conditionsToAnalyze) == 2) {
    if(is.null(colors)) {
      vennColors <- brewer.pal(n=3,name='Dark2')[1:2]
    } else {
      vennColors <- colors
    }
    if(is.null(alpha)) {
      myAlpha <- 0.7
    } else {
      myAlpha <- alpha
    }
  } else {
    if(is.null(colors)) {
      vennColors <- brewer.pal(n=length(conditionsToAnalyze),name='Dark2')
    } else {
      vennColors <- colors
    }
    if(is.null(alpha)) {
      myAlpha <- 0.4
    } else {
      myAlpha <- alpha
    }
  }
  
  myVennList <- list()
  for(i in 1:length(conditionsToAnalyze)) {
    myVennList[[conditionsToAnalyze[i]]] <- 0
  }
  myVenn <- venn.diagram(myVennList, category.names=conditionsToAnalyze , euler.d=FALSE, scale=FALSE, 
                              col='transparent', alpha=myAlpha, fill=vennColors, filename=NULL,
                              main=mainText, sub=subText, main.cex = 1.5)
  
  ### Modify the numbers in the venn diagram
  # get indexes of which elements in the gList (myVenn) is text elements (have a label)
  textIndexes <- which(!sapply( sapply(myVenn,'[[',"label") , is.null ))
  # extract info about how many values are in the venn diagram
  numberOfnumbersInVennDiagram <- c(3,7,15,31)[length(conditionsToAnalyze)-1]
  # loop over the values in the venn diagram and change them
  for(i in 1:numberOfnumbersInVennDiagram) {
    `[[`(`[[`( myVenn , textIndexes[i] ), "label") <- overwriteValues[i]
  }
  
  ### Draw the venn diagram
  grid.newpage()
  grid.draw(myVenn)
  
  
  return(spliceRobject) # since I add the plot data.frame and it should be stored
}