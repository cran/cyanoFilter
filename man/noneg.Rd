% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noneg.R
\name{noneg}
\alias{noneg}
\title{Removes negative values from the expression matrix}
\usage{
noneg(x)
}
\arguments{
\item{x}{is the flowframe whose expression matrix contains negative values}
}
\value{
flowframe with non-negative values in its expression matrix
}
\description{
Removes negative values from the expression matrix
}
\examples{
flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
              mustWork = TRUE)
flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
                               transformation = FALSE, emptyValue = FALSE,
                               dataset = 1) #FCS file contains only one data object
flowfile_nona <- cyanoFilter::nona(x = flowfile)
noneg(x = flowfile_nona)


}
