#' Prevalence Confidence Intervals
#'
#' @description `prevCI()` calculates 4 different types of confidence intervals for prevalence of parasites.
#'
#' @param data A data frame `('tbl', 'tbl_df', or 'data.frame')` that consists of at least one column of parasite presence datas in binomial form (0 or 1)
#' @param column Character, indicating the column of parasite presence values
#' @param method Character, indicating which which method or methods should be used to calculate CI for prevalence.
#' @param group Character, indicating which column to group data by.
#' @param conf Numerical, indicating coverage of the confidence interval
#' @param tolerance Numerical, Value passed into `blakerci()` for tolerance of Blaker confidence intervals.
#' @param correct Logical, indicating whether or not to apply continuity correction to Wilson score intervals.
#'
#' @return If only a single method is used, returns a data frame consisting of naive estimate, lower bound, and upper bound of prevalence. If multiple methods are specified, a list of data frames structured as mentioned is returned.
#'
#' @examples
#' data(sex)
#'
#' # Only a single measure
#' prevCI(sex, "presence", method = "ClopPear")
#'
#' # Multiple measures
#' prevCI(sex, "presence", method = c("Blaker", "Stern"))
#'
#' # Add a group & get all measures
#' (prevCI(sex, "presence", group = "sex"))
#'
#' @references
#' Blaker, H., 2000. Confidence curves and improved exact confidence intervals for discrete distributions. Canadian Journal of Statistics 28: 783-798.
#' Clopper, C.J. and Pearson, E.S., 1934. The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika 26: 404-413.
#' Sterne, T.E., 1954. Some remarks on confidence or fiducial limits. Biometrika 41: 275-278.
#' Wilson, E.B., 1927. Probable inference, the law of succession, and statistical inference. Journal of the American Statistical Association 22:209-212.
#'
#'
#' @export

prevCI <- function(data,
                   column,
                   method = c("ClopPear", "Blaker", "Sterne", "Wilson"),
                   group = NULL,
                   conf = 0.95,
                   tolerance = 1e-05,
                   correct = FALSE){

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

  # If any values in group are NA
  if(!is.null(group)){
    if(sum(is.na(data[[group]])) > 0){
      print("NAs present in grouping column. These rows were removed for calculations.")
      data <- data[!is.na(data[[group]]),]
    }
  }

  #If any values are NA
  if(sum(is.na(data[[column]])) > 0){
    print("NAs present in parasite presence column. These rows were removed for calculations.")
    data <- data[!is.na(data[[column]]),]
  }

  if(isTRUE(!all(data %>% pull(column) %in% c(0,1)))){
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

  # Correct
  if(any(method %in% "Wilson") & isFALSE(correct)){
    message("Wilson Score Interval without continuity correction")
  }

  if(any(method %in% "Wilson") & isTRUE(correct)){
    message("Wilson Score Interval with continuity correction")
  }

  # Print confidence interval
  message(paste0((conf * 100), "% Confidence Intervals"))

  # If length of method > 1
  if(length(method) > 1){
    list <- lapply(1:length(method), function(i){
      switch(method[i],
             ClopPear = cpCI(data = data,
                             column = column,
                             conf = conf,
                             group = group),
             Blaker = blaker(data = data,
                             column = column,
                             conf = conf,
                             group = group,
                             tolerance = tolerance),
             Sterne = sterne(data = data,
                           column = column,
                           group = group,
                           conf = conf),
             Wilson = Wilson(data = data,
                             column = column,
                             group = group,
                             conf = conf,
                             correct = correct))
    })
    names(list) <- method
    list <- lapply(list, function(x) x %>% mutate_if(is.numeric, round, 3))
    return(list)
  }

  # If length method = 1
  else{
    # Run the switch
    out <- switch(method,
                  ClopPear = cpCI(data = data,
                                  column = column,
                                  conf = conf,
                                  group = group),
           Blaker = blaker(data = data,
                           column = column,
                           conf = conf,
                           group = group,
                           tolerance = tolerance),
           Stern = sterne(data = data,
                         column = column,
                         group = group,
                         conf = conf),
           Wilson = Wilson(data = data,
                           column = column,
                           group = group,
                           conf = conf,
                           correct = correct))
    # Now round the output
    out <- out %>% mutate_if(is.numeric, round, 3)
    return(out)
  }
}
