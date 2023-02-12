#' Poulin's Index of Discrepancy
#'
#'
#' @noRd

inDISC <- function(data,
                  column,
                  group = NULL,
                  conf.int = FALSE,
                  conf = 0.95,
                  r = 2000){

  # If group is NULL
  if(is.null(group)){
    # If confidence interval = FALSE
    if(!isTRUE(conf.int)){
      # Bring in the data and restructure it
      dat <- data %>%
        dplyr::pull(.data[[column]])

      # Now calculate the top
      top <- do.call(sum, lapply(1:length(dat), function(i, j){
        d <- abs(dat[i] - dat[j])
      }))

      # Now calculate the bottom
      bottom <- (2 * length(dat) ^2) * mean(dat)

      # Now calculate the Index of discrepancy
      index <- top / bottom
      return(index)
    }

    # If confidence interval = TRUE
    else{

      # Get the data and run the bootstrap
      index <- data %>%
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
    }
  }
  # If group is specified
  else {
    # If confidence interval = FALSE
    if(!isTRUE(conf.int)){

      # Create a new df
      index <- do.call(rbind, lapply(1:length(unique(data[[group]])), function(g){

        # Bring in the data and restructure it
        dat <- data %>%
          dplyr::filter(.data[[group]] == unique(data[[group]][g])) %>%
          dplyr::pull(.data[[column]])

        # Now calculate the top
        top <- do.call(sum, lapply(1:length(dat), function(i, j){
          d <- abs(dat[i] - dat[j])
        }))

        # Now calculate the bottom
        bottom <- (2 * length(dat) ^2) * mean(dat)

        # Now calculate the Index of discrepancy
        i <- top / bottom
        return(i)
      }))
    }

    # If confidence interval = TRUE
    else{

      # Get values for each group
      index <- do.call(rbind, lapply(1:length(unique(data[[group]])), function(g){

        # Run the bootstrap
        i <- rodent %>%
          dplyr::filter(sex == unique(data[[group]])[g]) %>%
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
        return(i)
      }))
    }
  }
 return(index)
}




