#' Sterne Confidence Intervals
#'
#' @description Calculates confidence intervals as described by Stern (1954)
#'
#' @param data A data frame `('tbl', 'tbl_df', or 'data.frame')` that consists of at least one column of parasite presence datas in binomial form (0 or 1)
#' @param column Character, indicating the column of parasite presence values
#' @param group Character, indicating which column to group data by.
#' @param alternative Character, indicating what type of test to run in binom.exact(). Is not important for the output.
#' @param method Character, indicating the 'minlike' method within binom.exact() to get Stern confidence intervals
#' @param conf Numerical, indicating coverage of the confidence interval
#'
#' @noRd

sterne <- function(data,
                  column,
                  group = NULL,
                  alternative = "t",
                  method = "minlike",
                  conf = 0.95){
  # If group is null
  if(is.null(group)){

    # Take the input data and summarize it
    df <- data %>%
      dplyr::summarize(tot_parasite = sum(.data[[column]]), # Get a total number of parasites
                       n = dplyr::n(), # Get a sample size
                       naive_prev = tot_parasite / n) # Calculate a Naive Prevalence

    # Run the binom.exact function
    df1 <- exactci::binom.exact(df$tot_parasite,
                                df$n,
                                p = 0.5,
                                alternative = alternative,
                                tsmethod = method,
                                conf.level = conf)

    # Now reformat the data into a DF
    out <- data.frame(
      Naive_Prev = df$naive_prev,
      Lower = df1[["conf.int"]][1],
      Upper = df1[["conf.int"]][2],
      N = df$n
    )
  }
  else{

    # Take the input data and summarize it
    df <- data %>%
      group_by(.data[[group]]) %>%
      dplyr::summarize(tot_parasite = sum(.data[[column]]), # Get a total number of parasites
                       n = dplyr::n(), # Get a sample size
                       naive_prev = tot_parasite / n) # Calculate a Naive Prevalence

    # Loop through each unique group
    out <- do.call(rbind, lapply(1:length(unique(data[[group]])), function(i){

      # Subset the original df
      df1 <- df %>%
        filter(.data[[group]] == unique(data[[group]])[i])

      # Run the binom.exact function
      df2 <- exactci::binom.exact(df1$tot_parasite,
                                  df1$n,
                                  p = 0.5,
                                  alternative = alternative,
                                  tsmethod = method,
                                  conf.level = conf)

      # Now reformat the data into a DF
      df3 <- data.frame(
        Group = unique(data[[group]])[i],
        Naive_Prev = df1$naive_prev,
        Lower = df2[["conf.int"]][1],
        Upper = df2[["conf.int"]][2],
        N = df1$n
      )
      return(df3)
    }))
  }
  return(out)
}

