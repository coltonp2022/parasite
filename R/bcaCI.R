#' Calculate bca bootstrap confidence intervals
#'
#' This function allows you to calculate bootstrap confidence intervals using the bca method.
#' @param data A data frame with parasite data
#' @param column Character, indicating the column of parasite intensity data where values are >= 0
#' @param conf Numerical, indicating the confidence level desired
#' @param r Numerical, indicated the number of replicates for the bootstrapping function
#' @param group Character, allows for grouping the data together and calculating intervals for each group
#'
#'
#' @return Returns an object in the form of a data frame that includes mean and confidence intervals. If group = T, a column with your grouping variables will also be returned within this data frame.
#'
#' @export

bcaCI <- function(data, column, conf = 0.95, r = 2000, group = NULL){

  # Data
  if(!is.data.frame(data)){
    stop("Input must be a data frame")
  }

  # Column
  if(!is.character(column)){
    stop("Column must be a character. For example 'flea_intensity' ")
  }

  # Confidence
  if(!is.numeric(conf)){
    stop("Confidence must be numeric.")
  }

  # Run the bootstrap in a pipe
  df1 <- data %>%
    pull(.data[[column]]) %>% # Pull the number of fleas column
    boot(data = ., # Use that column
         statistic = function(x, i) mean(x[i]), # Statistic is the mean
         R = r) %>% # Bootstrap 2000 times
    boot.ci(boot.out = ., # Use that bootstrap sample for confidence intervals
            type = "bca",
            conf = conf) # Use the bias corrected and accelerated bootstrap

  # Now get all the estimates and make a dataframe
  final_df <- data.frame(
    Mean = df1$t0,
    Lower = df1$bca[4],
    Upper = df1$bca[5]
  )

  # Now return that final df
  return(final_df)
}
