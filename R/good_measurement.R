#' indicates if measurement from a flowfile is good or bad.
#'
#' @param metafile associated metafile to the supplied fcsfile. This is a csv file containig computed stats from the flow cytometer.
#' @param col_cpml column name or column number in metafile containing cell per microlitre measurements.
#' @param mxd_cellpML maximal accepted cell per microlitre. Flowfiles with larger cell per microlitre are termed bad.
#'                    Defaults to 1000.
#' @param mnd_cellpML minimum accepted cell per microlitre. Flowfiles with lesser cell per microlitre are termed bad.
#'                    Defaults to 50.
#'
#' @return character vector with length same as the number of rows in the metafile whose entries are
#'          \strong{good} for good files and \strong{bad} for bad files.
#'
#' @description This function examines the column containig \eqn{cells/\mu L}  and determins if the measurement can be used for further analysis or not based on a supplied range.
#'
#' @details Most flow cytometer makers will always inform clients within which range can measurements from the machine be trusted. The machines normally stores the amount of
#'          \eqn{cells/\mu L} it counted in a sample. Too large value could mean possible doublets and too low value could mean too little cells.
#'
#' @examples
#'  metadata <- system.file("extdata", "2019-03-25_Rstarted.csv", package = "cyanoFilter",
#'               mustWork = TRUE)
#'  metafile <- read.csv(metadata, skip = 7, stringsAsFactors = FALSE,
#'                      check.names = TRUE, encoding = "UTF-8")
#'  metafile <- metafile[, 1:65] #first 65 columns contains useful information
#'  #extract the part of the Sample.ID that corresponds to BS4 or BS5
#'  metafile$Sample.ID2 <- stringr::str_extract(metafile$Sample.ID, "BS*[4-5]")
#'  #clean up the Cells.muL column
#'  names(metafile)[which(stringr::str_detect(names(metafile), "Cells."))] <- "CellspML"
#'  goodfcs(metafile = metafile, col_cpml = "CellspML", mxd_cellpML = 1000, mnd_cellpML = 50)
#'
#' @export goodfcs

goodfcs <- function(metafile, col_cpml = "CellspML", mxd_cellpML = 1000, mnd_cellpML = 50) {
    if (!is.null(metafile) & !is.null(mxd_cellpML & !is.null(mnd_cellpML))) {

        goodfile <- ifelse((metafile[, col_cpml] < mxd_cellpML & metafile[, col_cpml] > mnd_cellpML), "good", "bad")

    } else stop("At least metafile is empty")
    return(goodfile)
}
