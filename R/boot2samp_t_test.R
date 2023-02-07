#' Bootstrapped 2-sample t-test for means
#'
#' Takes normal input data and calculates a bootstrap two-sample t-test from that data to test for differences in means or medians
#'
#'

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
      B <- sample(df2, nrow(data), replace = T)
      test <- t.test(A, B, var.equal = F)
      test$stat
  }

  # Now replicate that function
  t.vect <- replicate(r, boot_fun())

  # Now for two.sided
  if(alternative == "two.sided"){
   p_value <- round((1 - mean(abs(t.vect) > abs(t.est))), 3)
   cat(paste0("Hypothesis: ",
                 cat[1], " != ", cat[2],
                 "\n\nP-value = ", p_value))
  }

  # Greater
  if(alternative == "greater"){
   p_value <- round((1 - mean(t.vect > t.est)), 3)
   cat(paste0("Hypothesis: ",
                 cat[1], " > ", cat[2],
                 "\n\nP-value = ", p_value))
  }

  # Less
  if(alternative == "less"){
    p_value <- round((1 - mean(t.vect < t.est)), 3)
    cat(paste("Hypothesis: ",
                  cat[1], " < ", cat[2],
                  "\n\nP-value = ", p_value))
  }

  invisible(p_value)
}

