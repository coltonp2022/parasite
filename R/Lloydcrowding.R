#' Calculate Lloyd's Mean Crowding and BCa Confidence Intervals
#'
#' Calculates mean crowding as defined by Lloyd 1967 from parasite intensity data.
#'
#' @param data A data frame (`"tbl", "tbl_df", "data.frame"`) consisting of at least one column of parasite intensity
#' @param column Character, indicating the name of the column of parasite intensity
#' @param r Numeric, number of replicates to use within the BCa bootstrap
#' @param conf Numeric, confidence level for the BCa bootstrap.
#' @param group Character, indicating the name of the column that contains different groups within your data. Ex. `"sex", "fire", "age"`.
#'
#' @return Returns a data frame object with three columns: mean crowding, lower and upper limits of the confidence interval.
#'
#' @noRd


Lloydcrowding <- function(data,
                     column,
                     group = NULL,
                     conf = 0.95,
                     r = 5000){
  # Data
  if(!inherits(data, c("tbl", "tbl_df", "data.frame"))){
    stop("Data must be class 'tbl', 'tbl_df', or 'data.frame'")
  }

  # Column
  if(!is.character(column)){
    stop("Column name must be a character i.e. 'num_parasites'")
  }

  if(!is.numeric(data %>% pull(column))){
    stop("Input parasite intensities must be numerical")
  }

  # Group
  if(!is.null(group)){
    if(!is.character(group)){
      stop("Group name must be a character i.e. 'sex'")
    }
  }

  if(is.null(group)){

    # Now lets calculate a bca bootstrap for crowding
    df1 <- data %>%
      dplyr::pull(.data[[column]]) %>%
      boot::boot(data = ., # Use that crowding df
               statistic = function(x, i){
                (mean(x[i])) + (var(x[i]) / mean(x[i])) - 1 # Statistic is the Lloyds mean crowding
               },
               R = r) %>% # Bootstrap 2000 times
      boot::boot.ci(boot.out = ., # Use that bootstrap sample for confidence intervals
                    type = "bca",
                    conf = conf) # Use the bias corrected and accelerated bootstrap

    # Now get all the estimates and make a dataframe
    final_df <- data.frame(
      Mean = df1$t0,
      Lower = df1$bca[4],
      Upper = df1$bca[5]
    )
  } else{
    # Make the input df easier to work with
    df <- data %>%
      dplyr::select(.data[[column]], .data[[group]])

    # Now create a looping function
    final_df <- do.call(rbind, lapply(1:nrow(unique(df[2])), function(i){

      # Filter the data to the group
      df2 <- data %>%
        dplyr::filter(.data[[group]] == as.character(unique(df[2])[i,])) %>%
        dplyr::pull(.data[[column]]) %>%
        boot::boot(data = ., # Use that column
                   statistic = function(x, i){
                     (mean(x[i])) + (var(x[i]) / mean(x[i])) - 1 # Statistic is the Lloyd's mean crowding
                   },
                   R = r) %>% # Bootstrap r times
        boot::boot.ci(boot.out = ., # Use that bootstrap sample for confidence intervals
                      type = "bca", # Use the bias corrected and accelerated bootstrap
                      conf = conf) # Set the confidence level

      # Now get all the estimates and make a dataframe
      df3 <- data.frame(
        column = unique(df[2])[i,],
        Measure = df2$t0,
        Lower = df2$bca[4],
        Upper = df2$bca[5]
      )
      colnames(df3)[1] <- paste0(group)

      # Now return this df
      return(df3)
    }))
  }
  return(final_df)
}
