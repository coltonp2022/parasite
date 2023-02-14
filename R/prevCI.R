#' Prevalence Confidence Intervals
#'
#' @export

prevCI <- function(data,
                   column,
                   measure = c("ClopPear", "Blaker", "Stern"),
                   group = NULL,
                   conf = 0.95,
                   tolerance = 1e-05){

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

  if(isTRUE(!any(data %>% pull(column) %in% c(0,1)))){
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

  # If any values in group are NA
  if(!is.null(group)){
    if(sum(is.na(data[[group]])) > 0){
      message("NAs present in grouping column. These rows were removed for calculations.")
      data <- data[!is.na(data[[group]]),]
    }
  }

  #If any values are NA
  if(sum(is.na(data[[column]])) > 0){
    message("NAs present in parasite presence column. These rows were removed for calculations.")
    data <- data[!is.na(data[[column]]),]
  }

  # If no parasitized individuals
  if(sum(data %>% pull(column)) == 0){
    stop("No parasitized individuals. All values for presence are 0.")
  }

  # Print confidence interval
  message(paste0((conf * 100), "% Confidence Intervals"))

  # If length of measure > 1
  if(length(measure) > 1){
    list <- lapply(1:length(measure), function(i){
      switch(measure[i],
             ClopPear = cpCI(data = data,
                             column = column,
                             conf = conf,
                             group = group),
             Blaker = blaker(data = data,
                             column = column,
                             conf = conf,
                             group = group,
                             tolerance = tolerance),
             Stern = stern(data = data,
                           column = column,
                           group = group,
                           conf = conf))
    })
    names(list) <- measure
    list <- lapply(list, function(x) x %>% mutate_if(is.numeric, round, 3))
    return(list)
  }

  # If length measure = 1
  else{
    # Run the switch
    out <- switch(measure,
                  ClopPear = cpCI(data = data,
                                  column = column,
                                  conf = conf,
                                  group = group),
           Blaker = blaker(data = data,
                           column = column,
                           conf = conf,
                           group = group,
                           tolerance = tolerance),
           Stern = stern(data = data,
                         column = column,
                         group = group,
                         conf = conf))
    # Now round the output
    return(out %>% mutate_if(is.numeric, round, 3))

  }
}
rodent[nrow(rodent) + 1,] <- NA

n <- (prevCI(rodent, "flea", group = "fire"))

