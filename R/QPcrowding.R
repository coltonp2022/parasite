#' Calculate Mean Crowding and BCa Confidence Intervals
#'
#' @description Calculates parasite crowding as defined by Reiczigel et al. 2005 from normal host parasite intensity data. This is the same crowding index as what is used in QPweb.
#'
#' @param data A data frame (`"tbl", "tbl_df", "data.frame"`) consisting of at least one column of parasite intensity values
#' @param column Character, indicating the name of the column of parasite intensity values
#' @param r Numeric, Number of replicates to use within the BCa bootstrap
#' @param conf Numeric, Confidence level for the BCa bootstrap.
#' @param group Character, indicating the name of the column that contains different groups within your data. Ex. `"sex"`, `"fire"`, `"age"`.
#'
#' @return Returns a data frame object with three columns: mean crowding, lower and upper limits of the confidence interval.
#'
#' @export


QPcrowding <- function(data,
                     column,
                     r = 5000,
                     conf = 0.95,
                     group = NULL){
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

  # r
  if(!is.numeric(r)){
    stop("Value for r must be numerical")
  }

  # Confidence
  if(!is.numeric(conf)){
    stop("Value for confidence must be numerical")
  }

  # Group
  if(!is.null(group)){
    if(!is.character(group)){
      stop("Group name must be a character i.e. 'sex'")
    }
  }

  # Create a crowding data set
  df <- do.call(c, lapply(1:length(data[[column]]), function(i){
    d <- rep(data[[column]][i], times = data[[column]][i])
    return(d)
  }))

  if(is.null(group)){
    # Now lets calculate a bca bootstrap for crowding
    df1 <- boot::boot(data = df, # Use that crowding df
                      statistic = function(x, i) mean(x[i]), # Statistic is the mean
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
      dat <- data %>% filter(.data[[group]] == as.character(unique(df[2])[i,]))

      # Create the crowding data set
      df1 <- do.call(c, lapply(1:length(dat[[column]]), function(j){
        d <- rep(dat[[column]][j], times = dat[[column]][j])
        return(d)
      }))

      # Run the bootstrap in a pipe
      df2 <- boot::boot(data = df1, # Use that column
                   statistic = function(x, i) mean(x[i]), # Statistic is the mean
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
