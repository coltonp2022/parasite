#' Hoover Index
#'
#' @description Calculates the Hoover index to quantify parasite aggregation.
#'
#' @export


H_index <- function(data,
                    column,
                    group = NULL,
                    conf.int = FALSE,
                    conf = 0.95,
                    r = 2000){

  if(is.null(group)){
    if(!isTRUE(conf.int)){
      # Calculate the hoover index
      hoov <- do.call(sum, lapply(1:length(data[[column]]), function(i){
        return(abs(data[[column]][i] - mean(data[[column]])))
      })) / (2 * length(data[[column]]) * mean(data[[column]]))
    } else{
      # Calculate the hoover index using the bca bootstrap
      hoov <- data %>%
        pull(column) %>%
        boot::boot(.,
                   statistic = function(x, i){
                     do.call(sum, lapply(1:length(x[i]), function(j){
                       return(abs(x[i][j] - mean(x[i])))
                     })) / (2 * length(x[i]) * mean(x[i]))
                   },
                   R = r) %>%
        boot::boot.ci(.,
                      conf = conf,
                      type = "bca")
    }
  } else{
    if(!isTRUE(conf.int)){

      # Loop this through multiple dfs
      hoov <- do.call(rbind, lapply(1:length(unique(data[[group]])), function(i){

        # Get the data frames for each unique group
        dat <- data %>%
          filter(.data[[group]] == unique(data[[group]][i]))

        # Calculate the hoover index
        dat1 <- do.call(sum, lapply(1:length(dat[[column]]), function(j){
          return(abs(dat[[column]][j] - mean(dat[[column]])))
        })) / (2 * length(dat[[column]]) * mean(dat[[column]]))

        # Now get a full data frame
        dat2 <- data.frame(
          Group = unique(data[[group]][i]),
          Mean = dat1$t0,
          Lower = dat1$bca[4],
          Upper = dat1$bca[5]
        )
        colnames(dat2)[1] <- group
        return(dat2)
      }))


    } else{
      # Loop through multiple dfs
      hoov <- do.call(rbind, lapply(1:length(unique(data[[group]])), function(i){

        # Get a new data frame
        dat <- data %>%
          filter(.data[[group]] == unique(data[[group]][i]))

        # Calculate the hoover index using the bca bootstrap
        dat1 <- dat %>%
          pull(column) %>%
          boot::boot(.,
                     statistic = function(x, j){
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
          Group = unique(data[[group]][i]),
          Mean = dat1$t0,
          Lower = dat1$bca[4],
          Upper = dat1$bca[5]
        )
        colnames(dat2)[1] <- group
        return(dat2)
      }))
    }
  }
  return(hoov)
}
