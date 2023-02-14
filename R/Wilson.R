#' Wilson Score Interval
#'
#' @param data A data frame `('tbl', 'tbl_df', or 'data.frame')` that consists of at least one column of parasite presence datas in binomial form (0 or 1)
#' @param column Character, indicating the column of parasite presence values
#' @param group Character, indicating which column to group data by.
#' @param alternative Character, indicating the type of test passed to prop.test(). Is not important for confidence interval construction.
#' @param conf Numerical, indicating coverage of the confidence interval
#' @param correct Logical, indicating whether or not to apply continuity correction to Wilson score intervals.

Wilson <- function(data,
                  column,
                  group = NULL,
                  alternative = "t",
                  conf = 0.95,
                  correct = FALSE){
  # If group is null
  if(is.null(group)){

    # Take the input data and summarize it
    df <- data %>%
      dplyr::summarize(tot_parasite = sum(.data[[column]]), # Get a total number of parasites
                       n = dplyr::n(), # Get a sample size
                       naive_prev = tot_parasite / n) # Calculate a Naive Prevalence

    # Run the binom.exact function
    df1 <- prop.test(df$tot_parasite,
                     df$n,
                     alternative = alternative,
                     conf.level = conf,
                     correct = correct)

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

      # Run the Wilson test
      df2 <- prop.test(df1$tot_parasite,
                       df1$n,
                       alternative = alternative,
                       conf.level = conf,
                       correct = correct)

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

