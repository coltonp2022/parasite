% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot2samp_t_test.R
\name{boot2samp_t_test}
\alias{boot2samp_t_test}
\title{Bootstrapped 2-sample t-test for mean abundance or intensity}
\usage{
boot2samp_t_test(data, column, group, alternative = "two.sided", r = 2000)
}
\arguments{
\item{data}{A data frame \verb{('tbl', 'tbl_df', or 'data.frame')} consisting of at least one column of parasite intensity data}

\item{column}{Character, indicating the column of parasite intensity data where values are >= 0}

\item{group}{Character, indicating the column to group the data by for testing. Must be only two unique values for group.}

\item{alternative}{Character, indicating the type of test to be conducted. Can be "two.sided", "greater", or "less".}

\item{r}{Numerical, indicating the number of bootstrap replicates.}
}
\value{
Returns a printed statement indicating the alternative hypothesis and p-value. Returns the value for the p-value if assigning to an object.
}
\description{
Calculates a bootstrap two-sample t-test to test for differences in mean abundance or intensity
}
\examples{
data(sex)

# Two sided test
boot2samp_t_test(sex, "intensity", group = "sex")

# Greater than
boot2samp_t_test(sex, "intensity", group = "sex", alternative = "greater")

# Less than
boot2samp_t_test(sex, "intensity", group = "sex", alternative = "less")

}
