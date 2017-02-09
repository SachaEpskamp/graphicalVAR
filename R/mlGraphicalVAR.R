# Multi-level like graphical VAR:
mlGraphicalVAR <- function(
  data,
  vars,
  beepvar,
  dayvar,
  idvar,
  scale = TRUE,
  centerWithin = TRUE,
  glasso_gamma = 0.5, # Gamma used in glasso in qgraph
  verbose = TRUE,
  subjectNetworks = TRUE,
  ... # Args sent to graphicalVAR
){
  if (missing(idvar)) stop("'idvar' must be assigned")
  
  # Prep data:
  dataPrepped <- tsData(data,vars=vars,beepvar=beepvar,dayvar=dayvar,idvar=idvar,scale=scale,centerWithin=centerWithin)
  
  if (verbose){
    message("Estimating fixed networks")
  }
  
  # Fixed effects:
  ResFixed <- graphicalVAR(dataPrepped, ...)
  
  # Between-subjects:
  if (verbose){
    message("Estimating between-subjects network")
  }
  meansData <- dataPrepped$data_means
  meansData <- meansData[,names(meansData) != idvar]
  meansData <- meansData[rowMeans(is.na(meansData))!=1,]
  ResBetween <- qgraph::EBICglasso(cov(meansData),nrow(meansData),glasso_gamma,returnAllResults = TRUE)
  
  # Computing model per person:
 
  
  IDs <- unique(dataPrepped$data[[idvar]])
  idResults <- list()
  if (subjectNetworks){
    if (verbose){
      message("Estimating subject-specific networks")
      pb <- txtProgressBar(max=length(IDs),style=3)
    }
    for (i in seq_along(IDs)){
      capture.output({idResults[[i]] <- try(suppressWarnings(graphicalVAR(dataPrepped$data[dataPrepped$data[[idvar]] == IDs[i],],
                                                          vars=vars,
                                                          beepvar=beepvar,
                                                          dayvar=dayvar,
                                                          idvar=idvar,
                                                          scale = scale,
                                                          centerWithin = centerWithin,...,verbose = FALSE)))})
      if (verbose){
        setTxtProgressBar(pb,i)
      }
      if (is(idResults[[i]], "try-error")){
        idResults[[i]] <- list()
      }
    }   
    if (verbose){
      close(pb)
    }
  } else {
    idResults <- lapply(seq_along(IDs),function(x)list())
  }
  

  
  # Aggregate results:
  Results <- list(fixedPCC = ResFixed$PCC, 
                  fixedPDC = ResFixed$PDC,
                  fixedResults = ResFixed,
                  betweenNet = ResBetween$optnet,
                  betweenResults = ResBetween,
                  ids = IDs,
                  subjectPCC = lapply(idResults, '[[', 'PCC'),
                  subjecResults = idResults)
  class(Results) <- "mlGraphicalVAR"
  return(Results)
}