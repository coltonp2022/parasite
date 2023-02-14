#' Hoover Index
#'
#' @description Calculates the Hoover index to quantify parasite aggregation.
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
        pull(column) %>% # Pull the column into a vector
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
          filter(.data[[group]] == unique(data[[group]])[i])

        # Calculate the hoover index using the bca bootstrap
        dat1 <- dat %>%
          pull(column) %>%
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
