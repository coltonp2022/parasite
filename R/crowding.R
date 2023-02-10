#' Calculate Mean Crowding and Confidence Intervals
#'
#'

crowding <- function(data,
                     column,
                     r = 2000,
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

  # Make the input a dataframe
  if(inherits(data, c("tbl", "tbl_df"))){
    data <- as.data.frame(data)
  }

  # Create a crowding data set
  df <- do.call(c, lapply(1:length(data[,column]), function(i){
    d <- rep(data[i, column], times = data[i, column])
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
      df1 <- do.call(c, lapply(1:length(dat[,column]), function(j){
        d <- rep(dat[j, column], times = dat[j, column])
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
        Group = unique(df[2])[i,],
        Measure = df2$t0,
        Lower = df2$bca[4],
        Upper = df2$bca[5]
      )

      # Now return this df
      return(df3)
    }))
  }
  return(final_df)
}


