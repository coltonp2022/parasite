#' Bootstrapped 2-sample t-test for mean abundance or intensity
#'
#' @description Calculates a bootstrap two-sample t-test to test for differences in mean abundance or intensity
#'
#' @param data A data frame `('tbl', 'tbl_df', or 'data.frame')` consisting of at least one column of parasite intensity data
#' @param column Character, indicating the column of parasite intensity data where values are >= 0
#' @param group Character, indicating the column to group the data by for testing. Must be only two unique values for group.
#' @param alternative Character, indicating the type of test to be conducted. Can be "two.sided", "greater", or "less".
#' @param r Numerical, indicating the number of bootstrap replicates.
#'
#' @return Returns a printed statement indicating the alternative hypothesis and p-value. Returns the value for the p-value if assigning to an object.
#'
#' @examples
#' data(sex)
#'
#' # Two sided test
#' boot2samp_t_test(sex, "intensity", group = "sex")
#'
#' # Greater than
#' boot2samp_t_test(sex, "intensity", group = "sex", alternative = "greater")
#'
#' # Less than
#' boot2samp_t_test(sex, "intensity", group = "sex", alternative = "less")
#'
#' @references
#' Tibshirani, R.J. and Efron, B., 1993. An introduction to the bootstrap. Monographs on statistics and applied probability 57.
#'
#' @export

boot2samp_t_test <- function(data,
                             column,
                             group,
                             alternative = "two.sided",
                             r = 2000){
  # Data
  if(!inherits(data, c("tbl_df", "tbl", "data.frame"))){
    stop("Data must be class 'data.frame', 'tbl', or 'tbl_df'")
  }

  # Column
  if(!is.character(column)){
    stop("Input must be a character (i.e. 'intensity'")
  }

  if(!is.numeric(data %>% pull(column))){
    stop("Input parasite intensities must be numerical")
  }

  # Group
  if(!is.character(group)){
    stop("Input must be a character (i.e. 'sex'")
  }

  # Alternative
  if(!alternative %in% c("two.sided", "greater", "less")){
    stop("Alternative hypothesis must be specified")
  }

  # Get the unique categories of the group
  cat <- data %>%
    dplyr::select(.data[[group]]) %>%
    unique(.) %>%
    pull(1)

  # Get the first data frame
  df1 <- data %>%
    dplyr::filter(.data[[group]] == cat[1]) %>%
    pull(.data[[column]])

  # Get the second data frame
  df2 <- data %>%
    dplyr::filter(.data[[group]] == cat[2]) %>%
    pull(.data[[column]])

  # Get a t estimate
  t.est <- t.test(df1, df2, var.equal = F)$stat

  # Now create bootstrapping function
  boot_fun <- function(){
      A <- sample(df1, nrow(data), replace = T)
      mean_A <- mean(A)
      B <- sample(df2, nrow(data), replace = T)
      mean_B <- mean(B)
      test <- t.test(A, B, var.equal = F)$stat
      dat <- rbind(mean_A, mean_B, test) %>% as.data.frame()
      return(dat)
  }

  # Now replicate that function
  t.vect <- data.frame(do.call(rbind, replicate(r, boot_fun())))
  colnames(t.vect) <- c("MeanA",
                        "MeanB",
                        "Stat")

  # Now for two.sided
  if(alternative == "two.sided"){
   p_value <- (1 - mean(abs(t.vect$Stat) > abs(t.est)))
   message(paste0("Alternative Hypothesis: ",
                 cat[1], " != ", cat[2]))

   if(p_value < 0.001){
     cat("P-value < 0.001 \n")
   } else cat("P-value = ", p_value, "\n")

   cat(paste0("'", cat[1], "'", " Estimated Mean = ", round(mean(t.vect$MeanA), 3),
                "\n",
                "'", cat[2], "'", " Estimated Mean = ", round(mean(t.vect$MeanB), 3),
              "\n",
              "Estimated Difference in Means = ", round(mean(t.vect$MeanA) - mean(t.vect$MeanB), 3)))
  }

  # Greater
  if(alternative == "greater"){
   p_value <- (1 - mean(t.vect$Stat > t.est))
   message(paste0("Alternative Hypothesis: ",
                  cat[1], " > ", cat[2]))

   if(p_value < 0.001){
     cat("P-value < 0.001 \n")
   } else cat("P-value = ", p_value, "\n")

   cat(paste0("'", cat[1], "'", " Estimated Mean = ", round(mean(t.vect$MeanA), 3),
              "\n",
              "'", cat[2], "'", " Estimated Mean = ", round(mean(t.vect$MeanB), 3),
              "\n",
              "Estimated Difference in Means = ", round(mean(t.vect$MeanA) - mean(t.vect$MeanB), 3)))
  }

  # Less
  if(alternative == "less"){
    p_value <- (1 - mean(t.vect$Stat < t.est))
    message(paste0("Alternative Hypothesis: ",
                   cat[1], " < ", cat[2]))

    if(p_value < 0.001){
      cat("P-value < 0.001 \n")
    } else cat("P-value = ", p_value, "\n")

    cat(paste0("'", cat[1], "'", " Estimated Mean = ", round(mean(t.vect$MeanA), 3),
               "\n",
               "'", cat[2], "'", " Estimated Mean = ", round(mean(t.vect$MeanB), 3),
               "\n",
               "Estimated Difference in Means = ", round(mean(t.vect$MeanA) - mean(t.vect$MeanB), 3)))
  }

  invisible(p_value)
}
