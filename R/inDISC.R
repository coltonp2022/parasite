#' Poulin's Index of Discrepancy
#'
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
        pull(.data[[column]]) %>%
        boot(.,
             statistic = function(x, i){
               # Get the index
               do.call(sum, lapply(1:length(x[i]), function(j, k){
                 d <- abs(x[i][j] - x[i][k])
               })) / (2 * (length(x[i]) ^ 2) * mean(x[i]))
             },
             R = r) %>%
        boot.ci(.,
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
        out <- rodent %>%
          dplyr::filter(.data[[group]] == unique(data[[group]])[g]) %>%
          pull(.data[[column]]) %>%
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
