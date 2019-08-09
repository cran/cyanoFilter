#' Removes negative values from the expression matrix
#'
#' @param x is the flowframe whose expression matrix contains negative values
#' @return flowframe with non-negative values in its expression matrix
#'
#' @examples
#' flowfile_path <- system.file("extdata", "text.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' noneg(x = flowfile_nona)
#'
#'
#'@importFrom methods new
#'@export noneg

noneg <- function(x) {
    dtest <- !apply(flowCore::exprs(x), 1, function(row) any(row <= 0))
    exx <- flowCore::exprs(x)[dtest == T, ]
    paraa <- x@parameters
    describe <- x@description
    return(flowCore::flowFrame(exprs = exx, parameters = paraa, description = describe))
}
