#' Estimating the Negative Binomial Exponent, k
#'
#' @description Calculates both estimates and confidence intervals for the negative binomial exponent, k.
#'
#' @param data A data frame `('tbl', 'tbl_df', or 'data.frame')` consisting of at least one column of parasite intensity data for hosts
#' @param column `Character`, indicating the column of parasite intensity values
#' @param group `Character`, indicating the column that the data should be grouped by
#' @param conf `Numerical`, indicating the coverage of the confidence interval.
#' @param limit `Numerical`, maximum number of iterations for MLE
#'
#' @noRd

k_est <- function(data,
                  column,
                  conf = 0.975,
                  limit = 2000,
                  group = NULL){

  if(is.null(group)){
    # Create a dataset
    df1 <- data %>%
      dplyr::pull(.data[[column]]) # Pull the number of parasites column

    # Get point estimate of k
    k_point <- MASS::theta.ml(df1, mean(df1), length(df1), limit = 2000)

    # Get the zscore needed
    z <- qnorm(conf)

    # Now calculate an upper and lower
    lower <- k_point - (z * attributes(k_point)[[1]])
    upper <- k_point + (z * attributes(k_point)[[1]])
    if(lower < 0){
      message("Lower bound estimate for k overlapping 0")
      lower = 0
    }

    # Now put this all in a DF
    k <- data.frame(
      Estimate = k_point,
      Lower = lower,
      Upper = upper
    )
  }
  else{
    k <- do.call(rbind, lapply(1:length(unique(data[[group]])), function(i){
      # Create a dataset
      df1 <- data %>%
        dplyr::filter(.data[[group]] == unique(data[[group]])[i]) %>% # Subset
        dplyr::pull(.data[[column]]) # Pull the number of parasites column

      # Get point estimate of k
      k_point <- MASS::theta.ml(df1, mean(df1), length(df1), limit = 2000)

      # Get the zscore needed
      z <- qnorm(conf)

      # Now calculate an upper and lower
      lower <- k_point - (z * attributes(k_point)[[1]])
      upper <- k_point + (z * attributes(k_point)[[1]])
      if(lower < 0){
        message(paste0("Lower bound k estimate for '",
                       group, " = ",
                       unique(data[[group]])[i],
        "' overlapping 0"))
        lower = 0
      }

      # Now put this all in a DF
      out <- data.frame(
        Group = unique(data[[group]])[i],
        Estimate = k_point,
        Lower = lower,
        Upper = upper
      )
      colnames(out)[1] <- group
      return(out)
    }))
  }
  return(k)
}
