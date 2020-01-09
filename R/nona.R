#' Removes NA values from the expression matrix of a flow cytometer file.
#'
#' @param x flowframe with expression matrix containing NAs.
#' @return flowframe with expression matrix rid of NAs.
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' nona(x = flowfile)
#'
#'
#' @importFrom methods new
#' @export nona

nona <- function(x) {
    dtest <- !apply(flowCore::exprs(x), 1, function(row) any(is.na(row) | is.nan(row)))
    exx <- flowCore::exprs(x)[dtest == T, ]
    paraa <- x@parameters
    describe <- x@description
    return(flowCore::flowFrame(exprs = exx, parameters = paraa, description = describe))
}
