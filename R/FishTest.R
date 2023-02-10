#' Calculating Fisher's Exact Test
#'
#' This function takes normal input data and calculates p-values from Fisher's Exact Test on your specified groups. If looking for greater or less hypotheses, then make sure you set your factor levels correctly. This function depends on the fisher.test() function from base R. See that function for some extra information.
#'
#' @param data A data frame consisting of at least one column of parasite presence data in binomial form (0 or 1) as well as a column of groups to be tested against each other.
#' @param column Character, indicating the column of binomial parasite presence values.
#' @param group Character, indicating the column of groups wanting to be compared.
#' @param alternative Character, indicating the type of test to be conducted within the fisher.test() function. Options are "two.sided", "greater", "less".
#' @param simulate.p.value Logical, see fisher.test() function for description of this option.
#' @param B Numerical, indicating the number of replicates used in the Monte Carlo test
#' @param conf.int Logical, indicating whether or not to calculate the odds ratio confidence interval. Only works for data with 2 groups. Data with more groups will not have this option.
#' @param conf.level Numerical, indicating the level of confidence to calculate the confidence interval. Only works for data with 2 groups.
#'
#' @return Returns the list object described within the fisher.test() function help page.
#'
#' @export
#'

# Function

FishTest <- function(data,
                     column,
                     group,
                     alternative = "two.sided",
                     simulate.p.value = FALSE,
                     B = 2000,
                     conf.int = TRUE,
                     conf.level = 0.95){


  # Data
  if(!is.data.frame(data)){
    stop("Data must be a data frame")
  }

  # Column
  if(!is.character(column)){
    stop("Column must be in character form (i.e. 'flea')")
  }

  if(!is.numeric(data %>% pull(column))){
    stop("Input parasite presence values must be numerical (0 or 1)")
  }

  # Group
  if(!is.character(group)){
    stop("Group must be in character form (i.e. 'fire')")
  }

  # Alternative
  if(!alternative %in% c("two.sided", "greater", "less")){
    stop("Alternative must be either 'two.sided', 'greater', or 'less'")
  }


  # Run the function
  df <- data %>%
    dplyr::group_by(.data[[group]]) %>% # Group the data by the fire
    dplyr::count(.data[[column]]) %>% # Get a count of fleas
    tidyr::spread(., key = column, value = "n") %>% # Spread the data to diff format
    tibble::column_to_rownames(group) %>% # Send a column to rownames
    fisher.test(.,
                alternative = alternative,
                conf.int = conf.int,
                conf.level = conf.level,
                simulate.p.value = simulate.p.value,
                 B = B) # Run the fisher exact test

  # Reset the value for the dataframe name
  df[["data.name"]] <- deparse(substitute(data))

  # View the df at the end
  print(df)

}
