#' Chi-Squared Tests for Data Frames
#'
#' Takes input data in the form of a data frame, manipulates that data, and then runs the chisq.test() function for specified groups.
#'
#'
#'
#' @import


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

  # Now restructure data and run chisq.test()
  data <- data %>%
    dplyr::select(.data[[group]], .data[[column]]) %>% # Select needed columns
    dplyr::count(.data[[group]], .data[[column]]) %>% # Use the count function to get n
    tidyr::spread(., key = group, value = "n") %>% # Spread the data to wide
    as.matrix() %>% # Make this data a matrix
    chisq.test(.) # Run the chisq.test() function

  # Reset the value for the dataframe name
  data[["data.name"]] <- deparse(substitute(rodent))

  # Print
  print(data)
}



