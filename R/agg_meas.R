#' Parasite Aggregation Measures
#'
#' @description `Agg_meas()` calculates multiple different metrics quantifying the aggregation of parasites within a host population.
#'
#' @param data A data frame `('tbl', 'tbl_df', or 'data.frame')` consisting of at least one column of parasite intensity values for each host.
#' @param column Character, indicating the column consisting of parasite intensity values
#' @param measure Character, indicating which measures to be calculated.
#' @param group Character, indicating the column that data should be grouped by
#' @param conf Numerical, indicating coverage of confidence intervals for these measures
#' @param r Numerical, number of replicates for bootstrap resampling
#' @param limit Numerical, passed into theta.ml() as maximum number of iterations for MLE.
#'
#' @return If only a single measure is chosen, a data frame consisting of the estimate and confidence interval bound for the measure is returned. If multiple measures are chosen, returns a named_list of data frames structured as mentioned.
#'
#' @examples
#' data(sex)
#'
#' # Single measure no group
#' Agg_meas(sex, "intensity", measure = "PoulinD")
#'
#' # Multiple Measures no group
#' (Agg_meas(sex, "intensity", measure = c("Hoover", "PoulinD")))
#'
#' # Adding a group
#' Agg_meas(sex, "intensity", measure = "PoulinD", group = "sex")
#'
#' @references
#' Guyatt, H.L. and Bundy, D.A.P., 1991. Estimating prevalence of community morbidity due to intestinal helminths: prevalence of infection as an indicator of the prevalence of disease. Transactions of the Royal Society of Tropical Medicine and Hygiene 85:778-782.
#' Lloyd, M. 1967. Mean crowding. Journal of Animal Ecology 36: 1–30.
#' Poulin, R., 1993. The disparity between observed and uniform distributions: a new look at parasite aggregation. International journal for parasitology 23:937-944.
#' Reiczigel, J., Lang, Z., Rózsa, L. and Tóthmérész, B., 2005. Properties of crowding indices and statistical tools to analyze parasite crowding data. Journal of Parasitology 91:245-252.
#' Shaw, D.J. and Dobson, A.P., 1995. Patterns of macroparasite abundance and aggregation in wildlife populations: a quantitative review. Parasitology 111:S111-S133.
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
