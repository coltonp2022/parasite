#' Hoover Index
#'
#' @description Calculates the Hoover index to quantify parasite aggregation.
#'
#' @param data A data frame `('tbl', 'tbl_df', or 'data.frame')` consisting of at least one column of parasite intensity data for hosts
#' @param column `Character`, indicating the column of parasite intensity values
#' @param group `Character`, indicating the column that the data should be grouped by
#' @param conf.int `Logical`, indicating whether or not to construct a confidence interval for the value.
#' @param conf `Numerical`, indicating the coverage of the confidence interval.
#' @param r `Numerical`, number of bootstrap replicates for confidence interval construction
#'
#' @noRd
#'


H_index <- function(data,
                    column,
                    group = NULL,
                    conf.int = FALSE,
                    conf = 0.95,
                    r = 2000){

  if(is.null(group)){
      # Calculate the hoover index using the bca bootstrap
      dat <- data %>% # Start with the data
        dplyr::pull(column) %>% # Pull the column into a vector
        boot::boot(., # Data from pipe
                   statistic = function(x, i){ # Hoover index function
                     do.call(sum, lapply(1:length(x[i]), function(j){
                       return(abs(x[i][j] - mean(x[i])))
                     })) / (2 * length(x[i]) * mean(x[i]))
                   },
                   R = r) %>% # Run r replicates
        boot::boot.ci(., # Output bootstrap
                      conf = conf, # Confidence level
                      type = "bca") # Use the bca bootstrap

      # Now get a full data frame
      hoov <- data.frame(
        Mean = dat$t0, # Mean value
        Lower = dat$bca[4], # Lower CI
        Upper = dat$bca[5] # Upper CI
      )
  } else{

      # Loop through multiple dfs
      hoov <- do.call(rbind, lapply(1:length(unique(data[[group]])), function(i){

        # Get a new data frame
        dat <- data %>%
          dplyr::filter(.data[[group]] == unique(data[[group]])[i])

        # Calculate the hoover index using the bca bootstrap
        dat1 <- dat %>%
          dplyr::pull(column) %>%
          boot::boot(.,
                     statistic = function(x, j){ # Hoover index function
                       do.call(sum, lapply(1:length(x[j]), function(k){
                         return(abs(x[j][k] - mean(x[j])))
                       })) / (2 * length(x[j]) * mean(x[j]))
                     },
                     R = r) %>%
          boot::boot.ci(.,
                        conf = conf,
                        type = "bca")

        # Now get a full data frame
        dat2 <- data.frame(
          Group = unique(data[[group]])[i],
          Mean = dat1$t0,
          Lower = dat1$bca[4],
          Upper = dat1$bca[5]
        )
        colnames(dat2)[1] <- group
        return(dat2)
      }))
  }
  return(hoov)
}
