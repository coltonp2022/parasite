#' Blaker Confidence Intervals
#'
#'

blaker <- function(data,
                   column,
                   group = NULL,
                   conf = 0.95,
                   tolerance = 1e-05){

  # If group is null
  if(is.null(group)){

    # Take the input data and summarize it
    df <- data %>%
      dplyr::summarize(tot_parasite = sum(.data[[column]]), # Get a total number of parasites
                       n = dplyr::n(), # Get a sample size
                       naive_prev = tot_parasite / n) # Calculate a Naive Prevalence

    # Run the blakerci function
    df1 <- PropCIs::blakerci(df$tot_parasite,
                             df$n,
                             conf.level = conf,
                             tolerance = tolerance)

    # Now reformat the data into a DF
    out <- data.frame(
      Naive_Prev = df$naive_prev,
      Lower = df1[[1]][1],
      Upper = df1[[1]][2],
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
        dplyr::filter(.data[[group]] == unique(data[[group]])[i])

      # Run the blakerci function
      df2 <- PropCIs::blakerci(df1$tot_parasite,
                               df1$n,
                               conf.level = conf,
                               tolerance = tolerance)

      # Now reformat the data into a DF
      df3 <- data.frame(
        Group = unique(data[[group]])[i],
        Naive_Prev = df1$naive_prev,
        Lower = df2[[1]][1],
        Upper = df2[[1]][2],
        N = df1$n
      )
      return(df3)
    }))
  }
  return(out)
}

