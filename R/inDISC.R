#' Poulin's Index of Discrepancy
#'
#' @description Calculates the index of discrepancy as described by Poulin (1993)
#'
#' @param data A data frame `('tbl', 'tbl_df', or 'data.frame')` consisting of at least one column of parasite intensity data for hosts
#' @param column `Character`, indicating the column of parasite intensity values
#' @param group `Character`, indicating the column that the data should be grouped by
#' @param conf `Numerical`, indicating the coverage of the confidence interval.
#' @param r `Numerical`, number of bootstrap replicates for confidence interval construction
#'
#' @noRd


inDISC <- function(data,
                  column,
                  group = NULL,
                  conf = 0.95,
                  r = 2000){

  # If group is NULL
  if(is.null(group)){

      # Get the data and run the bootstrap
      out <- data %>%
        dplyr::pull(.data[[column]]) %>%
        boot::boot(.,
             statistic = function(x, i){
               # Get the index
               do.call(sum, lapply(1:length(x[i]), function(j, k){
                 d <- abs(x[i][j] - x[i][k])
               })) / (2 * (length(x[i]) ^ 2) * mean(x[i]))
             },
             R = r) %>%
        boot::boot.ci(.,
                conf = conf,
                type = "bca")

      # Format the index
      index <- data.frame(
        Mean = out$t0,
        Lower = out$bca[4],
        Upper = out$bca[5]
      )
  }
  # If group is specified
  else {
      # Get values for each group
      index <- do.call(rbind, lapply(1:length(unique(data[[group]])), function(g){

        # Run the bootstrap
        out <- data %>%
          dplyr::filter(.data[[group]] == unique(data[[group]])[g]) %>%
          dplyr::pull(.data[[column]]) %>%
          boot(.,
               statistic = function(x, i){
                 # Get the index
                 do.call(sum, lapply(1:length(x[i]), function(j, k){
                   d <- abs(x[i][j] - x[i][k])
                 })) / (2 * (length(x[i]) ^ 2) * mean(x[i]))
               },
               R = 2000) %>%
          boot.ci(.,
                  conf = 0.95,
                  type = "bca")

        # Format the index
        i <- data.frame(
          Group = unique(data[[group]])[g],
          Mean = out$t0,
          Lower = out$bca[4],
          Upper = out$bca[5]
        )
        colnames(i)[1] <- group
        return(i)
      }))
  }
  return(index)
}
