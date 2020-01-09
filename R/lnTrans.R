#' log transforms the expression matrix of a flowframe
#'
#' @param x flowframe to be transformed
#' @param notToTransform columns not to be transformed
#' @return \strong{flowframe} with log transformed expression matrix
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#'
#'
#' @importFrom methods new
#' @export lnTrans

lnTrans <- function(x, notToTransform = c("SSC.W", "TIME")) {
    exx <- cbind(log(flowCore::exprs(x)[, which(!(flowCore::colnames(x) %in% notToTransform))]),
                 flowCore::exprs(x)[, which(flowCore::colnames(x) %in% notToTransform)])
    colnames(exx) <- flowCore::colnames(x)
    paraa <- x@parameters
    describe <- x@description
    return(flowCore::flowFrame(exprs = exx, parameters = paraa, description = describe))
}
