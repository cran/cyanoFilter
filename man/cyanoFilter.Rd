% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cyanoFilter.R
\docType{package}
\name{cyanoFilter}
\alias{cyanoFilter}
\alias{cyanoFilter-package}
\title{cyanoFilter: A package to identify and/or assign indicators to BS4, BS5 cyanobacteria cells contained in water sample.}
\description{
The package provides two categories of functions:
\emph{metafile} preprocessing functions and \emph{fcsfile} processing functions.
}
\section{metafile preprocessing functions}{

          This set of functions (\code{\link{goodfcs}} and \code{\link{retain}}) helps to identify the appropriate fcs file to read.
}

\section{fcsfile processing functions}{

         These functions (\code{\link{nona}} and \code{\link{noneg}}, \code{\link{noneg}}, \code{\link{celldebris_nc}}, \code{\link{celldebris_emclustering}})
         works on the fcs file to identify the cell populations contained in the sample that generated this file.
}

