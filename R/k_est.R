#' Estimating the Negative Binomial Exponent, k
#'
#'

k_est <- function(data,
                  column,
                  conf.int = NULL,
                  conf = 0.975,
                  rep = NULL,
                  limit = 2000){

  if(isTRUE(conf.int)){
    # Create a dataset
    df1 <- data %>%
      count(.data[[column]]) %>% # Pull the number of parasites column
      as.data.frame()

    # Restructure the data set
    df2 <- rep(df1[,1], df1[,2])

    # Get point estimate of k
    k_point <- theta.ml(df2, mean(df2), length(df2), limit = 2000)

    # Now run the bootstrap
    fun <- function(){
      A <- sample(df2, length(df2), replace = T) # Sample
      k <- theta.ml(A, mean(A), limit = limit) # Get k estimate
    }

    # Now get the needed values
    k_boot <- replicate(rep, fun())

    # Now choose between regular bootstrap CI and bca bootstrap CI
    if(shapiro.test(k_boot)$p.value > 0.05){
      quant <- quantile(k_boot, c((1 - conf), conf))
    } else{

      # Now for a bca bootstrap
      # Set cutoffs
      cutoffs <- c((1 - 0.975), 0.975)

      # Compute constants for bca bootstrap
      z0 <- qnorm(mean(k_boot <= k_point))
      zu <- qnorm(cutoffs)

      # Calculate acceleration
      # Set an empty repetition
      I <- rep(NA, length(df2))

      # Start the jackknife
      for(i in 1:length(df2)){
        # Remove the ith data point
        dfnew <- df2[-i]
        # Estimate k
        k_jack <- theta.ml(dfnew, # New data frame
                           mean(dfnew), # Mean of that data
                           length(dfnew), # Sample size
                           limit = 5000) # Iteration limit
        I[i] <- (length(df2) - 1) * (k_point - k_jack)
      }
      # Now Estimate acceleration
      a_hat <- (sum(I^3) / sum(I^2)^1.5) / 6

      # Now adjust the quantiles
      u_adjusted <- pnorm(z0 + (z0 + zu) / (1 - a_hat * (z0 + zu)))

      # Now get the quantiles
      quant <- quantile(k_boot, u_adjusted)

      # Now create a data frame
      dat <- data.frame(cbind(k_point[1], unname(quant[1]), unname(quant[2])))
      colnames(dat) <- c("k_est", paste0(cutoffs[1] * 100, "%"), paste0(cutoffs[2] * 100, "%"))

      # Return the data
      invisible(dat)

      # Now print a message about the interval
      cat(paste0("BCa method used\n\n", "k = ", round(k_point, 3), "\n\n",
                 as.character(paste0(conf *100, "%")), " Confidence Interval: ",
                 "(", unname(round(quant[1], 3)), " - ",
                 unname(round(quant[2], 3)), ")"))
    }
  } else{
    # Set up the parasite count data
    df1 <- data %>%
      count(.data[[column]]) %>% # Count the data
      as.data.frame() # Make it a normal df

    # Restructure the data
    df2 <- rep(df1[,1], df1[,2])

    # Get a mean
    m <- mean(df2)

    # Get MLE of k
    kmle <- theta.ml(df2, m, limit = limit)
    return(kmle)
  }
}








