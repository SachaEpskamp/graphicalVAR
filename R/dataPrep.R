# Data preperation function:
tsData <- function(data,
                     vars,
                     beepvar,
                     dayvar,
                     idvar,
                     scale = TRUE,
                     centerWithin = TRUE,
                   deleteMissings = TRUE){
  . <- NULL
  
  data <- as.data.frame(data)
  
  # Add subject:
  if (missing(idvar)){
    idvar <- "ID"
    data[[idvar]] <- 1
  }
  
  # Add day:
  if (missing(dayvar)){
    dayvar <- "DAY"
    data[[dayvar]] <- 1
  }

  # Add beepvar:
  if (missing(beepvar)){
    beepvar <- "BEEP"
    data <- data %>% dplyr::group_by_(dayvar,idvar) %>% 
      dplyr::mutate_(BEEP = ~seq_len(n()))
  }
  
  # Vars:
  if (missing(vars)){
    vars <- names(data[!names(data)%in%c(idvar,dayvar,beepvar)])
  }
  
  # Center and scale data:
  data[,vars] <- scale(data[,vars], TRUE, scale)
  
  # Obtain person specific means:
  MeansData <- data %>% dplyr::group_by_(idvar) %>% dplyr::summarize_each_(funs(mean(.,na.rm=TRUE)),vars = vars)
  
  # Within-person center:
  if (centerWithin){
    # Only if N > 1 (very minimal floating point error can lead to different layout to older version otherwise)
    if (length(unique(data[[idvar]])) > 1){
      data <- data %>% dplyr::group_by_(idvar) %>% dplyr::mutate_each_(funs(scale(.,center=TRUE,scale=FALSE)),vars = vars)          
    }
  }

  # From mlVAR: Augment data:
  # Augment the data
  augData <- data
  
  # Add missing rows for missing beeps
  beepsPerDay <-  eval(substitute(dplyr::summarize_(data %>% group_by_(idvar,dayvar), 
                                                    first = ~ min(beepvar,na.rm=TRUE),
                                                    last = ~ max(beepvar,na.rm=TRUE)), 
                                  list(beepvar = as.name(beepvar))))
  
  # all beeps:
  allBeeps <- expand.grid(unique(data[[idvar]]),unique(data[[dayvar]]),seq(min(data[[beepvar]],na.rm=TRUE),max(data[[beepvar]],na.rm=TRUE))) 
  names(allBeeps) <- c(idvar,dayvar,beepvar)
  
  # Left join the beeps per day:
  allBeeps <- eval(substitute({
    allBeeps %>% dplyr::left_join(beepsPerDay, by = c(idvar,dayvar)) %>% 
      dplyr::group_by_(idvar,dayvar) %>% dplyr::filter_(~BEEP >= first, ~BEEP <= last)%>%
      dplyr::arrange_(idvar,dayvar,beepvar)
  },  list(BEEP = as.name(beepvar))))
  
  
  # Enter NA's:
  augData <- augData %>% dplyr::right_join(allBeeps, by = c(idvar,dayvar,beepvar))
  
  # Obtain data_c (slice away first row per day/subject):
  data_c <- augData %>% dplyr::group_by_(idvar,dayvar) %>% dplyr::slice(-1)
  
  # Obtain data_l (slice away last row per day/subject):
  data_l <- augData %>% dplyr::group_by_(idvar,dayvar) %>% dplyr::slice(-n())

  # # Remove rows with missings:
  if (deleteMissings){
    isNA <- rowSums(is.na(data_c)) > 0 | rowSums(is.na(data_l)) > 0
    data_c <- data_c[!isNA,]
    data_l <- data_l[!isNA,]
  }

  # Return datasets:
  Results <- list(
    data = augData,
    data_c = data_c[,vars],
    data_l = cbind(1,data_l[,vars]),
    data_means = MeansData,
    vars=vars,
    idvar=idvar,
    dayvar=dayvar,
    beepvar=beepvar
  )
  
  return(Results)
}