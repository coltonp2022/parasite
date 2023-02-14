#' Aggregation Measures
#'
#' @export

Agg_meas <- function(data,
                     column,
                     measure = c("k", "Hoover", "PoulinD", "QPcrowd",
                                 "Lloydcrowd"),
                     group = NULL,
                     conf = 0.95,
                     r = 2000,
                     limit = 2000){

  # Data
  if(!is.data.frame(data)){
    stop("Data must be of classes 'tbl', 'tbl_df' or 'data.frame'")
  }

  # Column
  if(!is.character(column)){
    stop("Input must be a character. Ex. 'num_parasite'.")
  }

  if(!is.numeric(data %>% pull(column))){
    stop("Input parasite intensities must be numerical")
  }

  # Measure
  if(!is.character(measure)){
    stop("Input measure must be a character. Ex. 'Hoover'")
  }

  # Group
  if(!is.null(group)){
    if(!is.character(group)){
      stop("Input group must be a character. Ex. 'Sex'")
    }
  }

  # Conf
  if(!is.numeric(conf)){
    stop("Confidence must be a numeric value")
  }

  # Reps
  if(!is.numeric(r)){
    stop("r must be a numeric value")
  }

  # Limit
  if(!is.numeric(limit)){
    stop("Limit must be a numeric value")
  }

  #If any values are NA
  if(sum(is.na(data[[column]])) > 0){
    message("NAs present in parasite intensity column. These values were removed for calculations.")
    data <- data[!is.na(data[[column]]),]
  }

  # Print a message regarding confidence
  message(paste0((conf * 100), "% ", "Confidence Intervals"))

  # If k
  if("k" %in% measure & !any(measure %in% c("Hoover", "PoulinD", "QPcrowd",
                                            "Lloydcrowd"))){
    message("Confidence Intervals for k estimated by MLE")
  }

  # If k and any other
  if("k" %in% measure & any(measure %in% c("Hoover", "PoulinD", "QPcrowd",
                                           "Lloydcrowd"))){
    message("Confidence Intervals for k estimated by MLE, others by BCa Bootstrap")
  }

  # If any other
  if(!"k" %in% measure & any(measure %in% c("Hoover", "PoulinD", "QPcrowd",
                                            "Lloydcrowd"))){
    message("Confidence Intervals estimated using BCa bootstrap")
  }

  # If more than one measure
  if(length(measure) > 1){
    # Set up the output list
    list <- lapply(1:length(measure), function(i){
      # Run the switch
      switch(measure[i],
             k = k_est(data = data,
                       column = column,
                       conf = conf,
                       limit = limit,
                       group = group),
             QPcrowd = QPcrowding(data = data,
                                  column = column,
                                  group = group,
                                  r = r,
                                  conf = conf),
             Hoover = H_index(data = data,
                              column = column,
                              group = group,
                              conf = conf,
                              r = r),
             PoulinD = inDISC(data = data,
                              column = column,
                              group = group,
                              conf = conf,
                              r = r),
             Lloydcrowd = Lloydcrowding(data = data,
                                        column = column,
                                        r = r,
                                        conf = conf,
                                        group = group))
    })
    # Now set the names of the list
    names(list) <- measure
    invisible(list)
  } else{
    # Run the switch
    return(switch(measure,
           k = k_est(data = data,
                     column = column,
                     conf = conf,
                     limit = limit,
                     group = group),
           QPcrowd = QPcrowding(data = data,
                                column = column,
                                group = group,
                                r = r,
                                conf = conf),
           Hoover = H_index(data = data,
                            column = column,
                            group = group,
                            conf = conf,
                            r = r),
           PoulinD = inDISC(data = data,
                            column = column,
                            group = group,
                            conf = conf,
                            r = r),
           Lloydcrowd = Lloydcrowding(data = data,
                                      column = column,
                                      r = r,
                                      conf = conf,
                                      group = group)))
  }
}
