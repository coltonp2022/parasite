#' Chi-Squared Tests for Data Frames
#'
#' Takes input data in the form of a data frame, manipulates that data, and then runs the chisq.test() function for specified groups.
#'


X2test <- function(data,
                   column,
                   group,
                   simulate.p.value = FALSE,
                   B = 2000){
  # Data
  if(!is.data.frame(data)){
    stop("Data must be input in data frame form")
  }

  # Column
  if(!is.character(column)){
    stop("Column must be a character value (i.e. 'yourcolumnnamehere')")
  }

  # Group
  if(!is.character(group)){
    stop("Column must be a character value (i.e. 'yourgroupnamehere')")
  }

  # Now start functions
  data <- data %>%
    dplyr::select(.data[[group]], .data[[column]]) %>%
    dplyr::count(.data[[group]], .data[[column]]) %>%
    tidyr::spread(., key = group, value = "n") %>%
    as.matrix() %>%
    chisq.test(.)

  # Reset the value for the dataframe name
  data[["data.name"]] <- deparse(substitute(rodent))

  # Print
  print(data)
}



