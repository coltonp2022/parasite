#' Aggregation Measures
#'

Agg_meas <- function(data,
                     column,
                     measure = c("k", "Hoover", "PoulinD", "QPcrowd",
                                 "Lloydcrowd"),
                     group = NULL,
                     conf.int = FALSE,
                     conf = 0.95,
                     r = 2000,
                     limit = 2000){

  # Print a message regarding confidence
  if(isTRUE(conf.int)){
    print(paste0((conf * 100), "% ", "Confidence Intervals"))
  }

  # If k
  if("k" %in% measure & !any(measure %in% c("Hoover", "PoulinD", "QPcrowd",
                                            "Lloydcrowd"))){
    print("Confidence Intervals for k estimated by MLE")
  }

  # If k and any other
  if("k" %in% measure & any(measure %in% c("Hoover", "PoulinD", "QPcrowd",
                                           "Lloydcrowd"))){
    print("Confidence Intervals for k estimated by MLE, others by BCa Bootstrap")
  }

  # If any other
  if(!"k" %in% measure & any(measure %in% c("Hoover", "PoulinD", "QPcrowd",
                                            "Lloydcrowd"))){
    print("Confidence Intervals estimated using BCa bootstrap")
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
                       limit = limit),
             QPcrowd = QPcrowding(data = data,
                                  column = column,
                                  group = group,
                                  r = r,
                                  conf = conf),
             Hoover = H_index(data = data,
                              column = column,
                              group = group,
                              conf.int = conf.int,
                              conf = conf,
                              r = r),
             PoulinD = inDISC(data = data,
                              column = column,
                              group = group,
                              conf.int = conf.int,
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
                     limit = limit),
           QPcrowd = QPcrowding(data = data,
                                column = column,
                                group = group,
                                r = r,
                                conf = conf),
           Hoover = H_index(data = data,
                            column = column,
                            group = group,
                            conf.int = conf.int,
                            conf = conf,
                            r = r),
           PoulinD = inDISC(data = data,
                            column = column,
                            group = group,
                            conf.int = conf.int,
                            conf = conf,
                            r = r),
           Lloydcrowd = Lloydcrowding(data = data,
                                      column = column,
                                      r = r,
                                      conf = conf,
                                      group = group)))
  }
}
