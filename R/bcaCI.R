#' Calculate bca bootstrap confidence intervals
#'
#' @description This function allows you to calculate bootstrap confidence intervals using the bca method.
#'
#' @param data A data frame `('tbl', 'tbl_df', or 'data.frame')` consisting of at least one column of parasite intensity data
#' @param column Character, indicating the column of parasite intensity data where values are >= 0
#' @param measure Character, indicating either abundance ("abun") or intensity ("int") to be used for calculation
#' @param conf Numerical, indicating the confidence level desired
#' @param r Numerical, indicated the number of replicates for the bootstrapping function
#' @param group Character, allows for grouping the data together and calculating intervals for each group
#'
#'
#' @return Returns an object in the form of a data frame that includes mean and confidence intervals. If data is grouped using the group argument, a column with your grouping variables will also be returned within this data frame.
#'
#' @examples
#' data(sex)
#'
#' # Intensity
#' bcaCI(sex, "intensity", measure = "int")
#'
#' # Abundance
#' bcaCI(sex, "intensity", measure = "abun")
#'
#' # Add a group in
#' bcaCI(sex, "intensity", group = "sex")
#'
#' @export

bcaCI <- function(data,
                  column,
                  group = NULL,
                  measure = "int",
                  conf = 0.95,
                  r = 2000){

  # Data
  if(!inherits(data, c("tbl", "tbl_df", "data.frame"))){
    stop("Input must be of classes 'tbl', 'tbl_df', or 'data.frame'")
  }

  # Column
  if(!is.character(column)){
    stop("Column must be a character. For example 'flea_intensity' ")
  }

  if(!is.numeric(data %>% pull(column))){
    stop("Input parasite intensities must be numerical")
  }

  # Confidence
  if(!is.numeric(conf)){
    stop("Confidence must be numeric.")
  }

  # Abundance or Intensity
  if(measure == "int"){
    data <- data %>% dplyr::filter(.data[[column]] > 0)
  }

  # Run without grouping
  if(is.null(group)){
    # Run the bootstrap in a pipe
    df1 <- data %>%
      dplyr::pull(.data[[column]]) %>% # Pull the number of fleas column
      boot::boot(data = ., # Use that column
           statistic = function(x, i) mean(x[i]), # Statistic is the mean
           R = r) %>% # Bootstrap 2000 times
      boot::boot.ci(boot.out = ., # Use that bootstrap sample for confidence intervals
              type = "bca",
              conf = conf) # Use the bias corrected and accelerated bootstrap

    # Now get all the estimates and make a dataframe
    final_df <- data.frame(
      Measure = df1$t0,
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
        dplyr::pull(.data[[column]]) %>% # Pull the number of parasites column
        boot::boot(data = ., # Use that column
             statistic = function(x, i) mean(x[i]), # Statistic is the mean
             R = r) %>% # Bootstrap r times
        boot::boot.ci(boot.out = ., # Use that bootstrap sample for confidence intervals
                type = "bca", # Use the bias corrected and accelerated bootstrap
                conf = conf) # Set the confidence level

      # Now get all the estimates and make a dataframe
      df2 <- data.frame(
        Group = unique(df[2])[i,],
        Measure = df1$t0,
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

