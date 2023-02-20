#' Chi-Squared Tests for Data Frames
#'
#' @description Takes input data in the form of a data frame, manipulates that data, and then runs the chisq.test() function for specified groups.
#'
#' @param data An input data frame consisting of at least binomial (0 or 1) parasitism data and a grouping variable column.
#' @param column Character, indicating the column of parasite presence data.
#' @param group Character, indicating the column name of the grouping variable column
#' @param simulate.p.value Logical, indicating whether to compute p-values by Monte Carlo simulation. Passed to chisq.test().
#' @param B Numerical, specifying the number of replicates used for the Monte Carlo test
#'
#' @return Returns the same data as the chisq.test() function. See documentation for descriptions.
#'
#' @export



X2test <- function(data,
                   column,
                   group,
                   simulate.p.value = FALSE,
                   B = 2000){
  # Data
  if(!inherits(data, c("tbl_df", "tbl", "data.frame"))){
    stop("Data must be of classes 'tbl_df', 'tbl', or 'data.frame'")
  }

  # Column
  if(!is.character(column)){
    stop("Column must be a character value (i.e. 'yourcolumnnamehere')")
  }

  if(!is.numeric(data %>% pull(column))){
    stop("Input parasite presence values must be numerical (0 or 1)")
  }

  # Group
  if(!is.character(group)){
    stop("Column must be a character value (i.e. 'yourgroupnamehere')")
  }

  # Now restructure data and run chisq.test()
  test <- data %>%
    dplyr::select(.data[[group]], .data[[column]]) %>% # Select needed columns
    dplyr::count(.data[[group]], .data[[column]]) %>% # Use the count function to get n
    tidyr::spread(., key = group, value = "n") %>% # Spread the data to wide
    as.matrix() %>% # Make this data a matrix
    chisq.test(., simulate.p.value = simulate.p.value, B = B) # Run the chisq.test() function

  # Reset the value for the dataframe name
  test[["data.name"]] <- deparse(substitute(data))

  # Print
  print(test)
}



