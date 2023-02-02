#' Calculate bca bootstrap confidence intervals
#'
#' This function allows you to calculate bootstrap confidence intervals using the bca method.
#' @param data A data frame with parasite data
#' @param column Character, indicating the column of parasite intensity data where values are >= 0
#' @param measure Character, indicating either abundance or intensity to be used for calculation
#' @param conf Numerical, indicating the confidence level desired
#' @param r Numerical, indicated the number of replicates for the bootstrapping function
#' @param group Character, allows for grouping the data together and calculating intervals for each group
#' @param print Logical, allows for printing of data frame when creating an object from the data frame
#'
#'
#' @return Returns an object in the form of a data frame that includes mean, and confidence intervals.
#'
#' @export

bcaCI <- function(data,
                  column,
                  measure,
                  conf = 0.95,
                  r = 2000,
                  group = NULL){

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

  # Abundance or Intensity
  if(measure == "abun"){
    data <- data
  }
  else if(measure == "int"){
    data <- data %>% filter(.data[[column]] > 0)
  }
  else{
    stop("Input must be either 'abun' or 'int'")
  }

  # Run without grouping
  if(is.null(group)){
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

  else {
    # Make the input df easier to work with
    df <- data %>%
      dplyr::select(.data[[column]], .data[[group]])

    # Now create a looping function
    final_df <- do.call(rbind, lapply(1:nrow(unique(df[2])), function(i){
      # Run the bootstrap in a pipe
      df1 <- df %>%
        dplyr::filter(.data[[group]] == as.character(unique(df[2])[i,])) %>%
        pull(.data[[column]]) %>% # Pull the number of fleas column
        boot(data = ., # Use that column
             statistic = function(x, i) mean(x[i]), # Statistic is the mean
             R = r) %>% # Bootstrap 2000 times
        boot.ci(boot.out = ., # Use that bootstrap sample for confidence intervals
                type = "bca",
                conf = conf) # Use the bias corrected and accelerated bootstrap

      # Now get all the estimates and make a dataframe
      df2 <- data.frame(
        Group = unique(df[2])[i,],
        Mean = df1$t0,
        Lower = df1$bca[4],
        Upper = df1$bca[5]
      )

      # Now return this df
      return(df2)
    }))

    # Now return that final df
    return(final_df)
  }
}



