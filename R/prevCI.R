#' Prevalence Confidence Intervals
#'
#' @export

prevCI <- function(data,
                   column,
                   measure = "ClopPear",
                   group = NULL,
                   conf = 0.95){

  # Data
  if(!inherits(data, c("tbl", "tbl_df", "data.frame"))){
    stop("Data must be of classes 'tbl', 'tbl_df', or 'data.frame'")
  }

  # Column
  if(!is.character(column)){
    stop("Column must be a character. For example 'flea_presence'")
  }

  if(!is.numeric(data %>% pull(column))){
    stop("Input parasite presence values must be numerical (0 or 1)")
  }

  if(isTRUE(!any(rodent %>% pull(column) %in% c(0,1)))){
    stop("Input parasite presence values must only be 0 or 1")
  }

  # Confidence
  if(!is.numeric(conf)){
    stop("Input value must be numeric")
  }

  # Group
  if(!is.null(group)){
    if(!is.character(group)){
      stop("Group must be a character. For example 'sex' ")
    }
  }

  #If any values are NA
  if(sum(is.na(data[[column]])) > 0){
    message("NAs present in parasite presence column. These values were removed for calculations.")
    data <- data[!is.na(data[[column]]),]
  }

  # If length of measure > 1
  if(length(measure) > 1){
    list <- lapply(1:length(measure), function(i){
      switch(measure[i],
             ClopPear = cpCI(data = data,
                             column = column,
                             conf = conf,
                             group = group))
    })
  }

  # If length measure = 1
  else{
    # Run the switch
    return(switch(measure,
                  ClopPear = cpCI(data = data,
                                  column = column,
                                  conf = conf,
                                  group = group)))
  }
}
