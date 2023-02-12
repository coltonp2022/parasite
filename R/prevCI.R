#' Prevalence Confidence Intervals
#'
#' @export

prevCI <- function(data,
                   column,
                   measure = "ClopPear",
                   group = NULL,
                   conf = 0.95){

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
