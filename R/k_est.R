#' Estimating the Negative Binomial Exponent, k
#'
#' @noRd

k_est <- function(data,
                  column,
                  conf = 0.975,
                  limit = 2000){

    # Create a dataset
    df1 <- data %>%
      pull(.data[[column]]) # Pull the number of parasites column

    # Get point estimate of k
    k_point <- theta.ml(df1, mean(df1), length(df1), limit = 2000)

    # Get the zscore needed
    z <- qnorm(conf)

    # Now calculate an upper and lower
    lower <- k_point - (z * attributes(k)[[1]])
    upper <- k_point + (z * attributes(k)[[1]])

    # Now put this all in a DF
    k <- data.frame(
      Estimate = k_point,
      Lower = lower,
      Upper = upper
    )
    return(k)
}





