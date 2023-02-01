#' This function allows you to calculate Clopper-Pearson confidence intervals.
#' @param data A data frame with parasite data
#' @param column Character, indicating the column consisting of binomial presence data in 1,0 form
#' @param conf Numerical, indicating the confidence level desired
#' @param group Character, allows for grouping the data together and calculating intervals for each group
#' @param print Logical, allows for printing of data frame when creating an object from the data frame
#'
#' @return Returns an object in the form of a data frame that includes naive prevalence, and confidence intervals. If group = T, this data frame also includes

#' @export

cpCI <- function(data, column, conf, group = NULL, print = FALSE){

  # Data type
  if(!is.data.frame(data)){
    stop("Input must be a data frame")
  }

  # Column
  if(!is.character(column)){
    stop("Column must be a character. For example 'flea_presence'")
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
      summarize(tot_parasite = sum(.data[[column]]), # Get a total number of parasites
                n = n(), # Get a sample size
                naive_prev = tot_parasite / n) # Calculate a Naive Prevalence

    # Now create the Confidence intervals
    df1 <- exactci(
      df$tot_parasite, # Total parasites
      df$n, # Total Sample Size
      conf.level = conf # Confidence
    )

    # Now get the CI and Naive Prevalence into a vector
    vector <- c(df$naive_prev, df1$conf.int)

    # Now make a final df
    final_df <- data.frame(
      Naive_Prev = vector[1],
      Lower = vector[2],
      Upper = vector[3]
    )

    if(isTRUE(print)){
      print(final_df)
    }

    # Return
    return(final_df)
  }

  else{
    # Take the input data and summarize it
    df <- data %>%
      group_by(.data[[group]]) %>%
      summarize(tot_parasite = sum(.data[[column]]), # Get a total number of parasites
                n = n(), # Get a sample size
                naive_prev = tot_parasite / n) # Calculate a Naive Prevalence

    # Loop through each group
    final_df <- do.call(rbind, lapply(1:nrow(df), function(i){

      # Subset the dataframe
      df1 <- df %>%
        filter(.data[[group]] == as.character(unique(df[group])[i,]))

      # Now create the Confidence intervals
      df2 <- exactci(
        df1$tot_parasite, # Total parasites
        df1$n, # Total Sample Size
        conf.level = conf # Confidence
      )

      # Now get the CI and Naive Prevalence into a vector
      vector <- c(df1$naive_prev,
                  df2$conf.int)

      # Now make a final df
      df3 <- data.frame(
        Group = as.character(unique(df[group])[i,]),
        Naive_Prev = vector[1],
        Lower = vector[2],
        Upper = vector[3]
      )

      return(df3)
    }))

    if(isTRUE(print)){
      print(final_df)
    }

    # Return
    return(final_df)
  }
}






