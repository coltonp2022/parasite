#' Calculate Clopper-Pearson confidence intervals
#'
#' This function allows you to calculate Clopper-Pearson confidence intervals.
#' @param data A data frame with parasite data
#' @param column Character, indicating the column consisting of binomial presence data in 1,0 form
#' @param conf Numerical, indicating the confidence level desired
#' @param group Character, allows for grouping the data together and calculating intervals for each group
#'
#' @return Returns an object in the form of a data frame that includes naive prevalence, and confidence intervals. If data is grouped using the group argument, this data frame also includes a column with your original groups.
#'
#'

cpCI <- function(data, column, conf = 0.95, group = NULL){

  # Data
  if(!inherits(data, c("tbl", "tbl_df", "data.frame"))){
    stop("Input must be of classes 'tbl', 'tbl_df', or 'data.frame'")
  }

  # Column
  if(!is.character(column)){
    stop("Column must be a character. For example 'flea_presence'")
  }

  if(!is.numeric(data %>% pull(column))){
    stop("Input parasite presence values must be numerical (0 or 1)")
  }

  if(isTRUE(!any(rodent %>% pull(column) %in% c(0,1)))){
    stop("Input values must only be 0 or 1")
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

  if(is.null(group)){
    # Take the input data and summarize it
    df <- data %>%
      dplyr::summarize(tot_parasite = sum(.data[[column]]), # Get a total number of parasites
                n = dplyr::n(), # Get a sample size
                naive_prev = tot_parasite / n) # Calculate a Naive Prevalence

    # Now get your input parameters
    alpha <- 1 - conf # Alpha
    n <- df$n # Sample size
    n1 <- df$tot_parasite # Number of "successes"
    f1 <- qf(1 - alpha/2, (2 * n1), 2*(n - n1 + 1), lower.tail = F) # Actual lower coverage
    f2 <- qf(alpha/2, 2*(n1 + 1), 2*(n - n1), lower.tail = F) # Actual upper coverage

    # Now get the confidence intervals
    pLow <- (1 + ((n - n1 + 1) / (n1 * f1))) ^ (-1)
    pUpp <- (1 + ((n - n1)) / ((n1 + 1) * f2)) ^ (-1)

    # Now make a final df
    final_df <- data.frame(
      Naive_Prev = df$naive_prev, # Mean prevalence
      Lower = pLow, # Lower bound
      Upper = pUpp # Upper bound
    )

    # Return
    return(final_df)
  }

  else{
    # Take the input data and summarize it
    df <- data %>%
      dplyr::group_by(.data[[group]]) %>%
      dplyr::summarize(tot_parasite = sum(.data[[column]]), # Get a total number of parasites
                n = dplyr::n(), # Get a sample size
                naive_prev = tot_parasite / n) # Calculate a Naive Prevalence

    # Loop through each group
    final_df <- do.call(rbind, lapply(1:nrow(df), function(i){

      # Subset the dataframe
      df1 <- df %>%
        dplyr::filter(.data[[group]] == as.character(unique(df[group])[i,]))

      # Now get your input parameters
      alpha <- 1 - conf # Alpha
      n <- df1$n # Sample size
      n1 <- df1$tot_parasite # Number of "successes"
      f1 <- qf(1 - alpha/2, (2 * n1), 2*(n - n1 + 1), lower.tail = F) # Actual lower coverage
      f2 <- qf(alpha/2, 2*(n1 + 1), 2*(n - n1), lower.tail = F) # Actual upper coverage

      # Now get the confidence intervals
      pLow <- (1 + ((n - n1 + 1) / (n1 * f1))) ^ (-1)
      pUpp <- (1 + ((n - n1)) / ((n1 + 1) * f2)) ^ (-1)

      # Now make a final df
      df2 <- data.frame(
        Group = unique(data[[group]][i]), # Grouping variable
        Naive_Prev = df1$naive_prev, # Mean Prevalence
        Lower = pLow, # Upper
        Upper = pUpp # Lower
      )
      colnames(df2)[1] <- group
      return(df2)
    }))

    # Return
    return(final_df)
  }
}
