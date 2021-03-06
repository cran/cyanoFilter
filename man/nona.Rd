% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nona.R
\name{nona}
\alias{nona}
\title{Removes NA values from the expression matrix of a flow cytometer file.}
\usage{
nona(x)
}
\arguments{
\item{x}{flowframe with expression matrix containing NAs.}
}
\value{
flowframe with expression matrix rid of NAs.
}
\description{
Removes NA values from the expression matrix of a flow cytometer file.
}
\examples{
flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
              mustWork = TRUE)
flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
                               transformation = FALSE, emptyValue = FALSE,
                               dataset = 1) #FCS file contains only one data object
nona(x = flowfile)


}
