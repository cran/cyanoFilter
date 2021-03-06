#' Decides if a file should be retiained or removed based on its status.
#'
#' @param meta_files dataframe from meta file that has been preprocessed by the \code{\link{goodfcs}} function.
#' @param make_decision decision to be made should more than one \eqn{cells/\mu L} be good.
#' @param Status column name in meta_files containing status obtained from the \code{\link{goodfcs}} function.
#' @param CellspML column name in meta_files containing \eqn{cells/\mu L} measurements.
#'
#' @return a character vector with entries "Retain" for a file to be retained or "No!" for a file to be discarded.
#'
#' @description Function to determine what files to retain and finally read from the flow cytometer FCS file.
#'
#' @details It is typically not known in advance which dilution level would result in the desired \eqn{cells/\mu L}, therefore
#'          the samples are ran through the flow cytometer at two or more dilution levels. Out of these, one has to decide which
#'          to retain and finally use for further analysis. This function and \code{\link{goodfcs}} are to help you decide that.
#'          If more than one of the dilution levels are judged good, the option \emph{make_decision = "maxi"} will give "Retain" to the
#'          row with the maximum \eqn{cells/\mu L} while the opposite occurs for \emph{make_decision = "mini"}.
#'          \emph{make_decision = "unique"} i there is only one measurement for that particular sample, while \emph{make_decision = "maxi"}
#'          and \emph{make_decision = "mini"} should be used for files with more than one measurement for the sample in question.
#'
#' @seealso \code{\link{goodfcs}}
#'
#' @examples
#'
#'  metadata <- system.file("extdata", "2019-03-25_Rstarted.csv", package = "cyanoFilter",
#'               mustWork = TRUE)
#'  metafile <- read.csv(metadata, skip = 7, stringsAsFactors = FALSE,
#'                       check.names = TRUE, encoding = "UTF-8")
#'  metafile <- metafile[, 1:65] #first 65 columns contain useful information
#'  #extract the part of the Sample.ID that corresponds to BS4 or BS5
#'  metafile$Sample.ID2 <- stringr::str_extract(metafile$Sample.ID, "BS*[4-5]")
#'  #clean up the Cells.muL column
#'  names(metafile)[which(stringr::str_detect(names(metafile), "Cells."))] <- "CellspML"
#'  metafile$Status <- cyanoFilter::goodfcs(metafile = metafile, col_cpml = "CellspML",
#'                             mxd_cellpML = 1000, mnd_cellpML = 50)
#'  metafile$Retained <- NULL
#'  # first 3 rows contain BS4 measurements at 3 dilution levels
#'  metafile$Retained[1:3] <- cyanoFilter::retain(meta_files = metafile[1:3,], make_decision = "maxi",
#'                    Status = "Status", CellspML = "CellspML")
#'  # last 3 rows contain BS5 measurements at 3 dilution levels as well
#'  metafile$Retained[4:6] <- cyanoFilter::retain(meta_files = metafile[4:6,], make_decision = "maxi",
#'                    Status = "Status", CellspML = "CellspML")
#'
#'
#' @export retain
#'
retain <- function(meta_files, make_decision = c("maxi", "mini", "unique"), Status = "Status", CellspML = "CellspML") {
    are_all_good <- sum(meta_files[Status] == "good")

    decision <- rep(NA, nrow(meta_files))
    decision[which(meta_files[Status] != "good")] <- "No!"
    good_pos <- which(meta_files[Status] == "good")

    if (make_decision == "mini") {

        if (are_all_good == 1) {

            decision <- ifelse(meta_files[Status] == "good", "Retain", "No!")

        } else if (are_all_good >= 2) {

            mini_pos <- which(meta_files[good_pos, CellspML] == min(meta_files[good_pos, CellspML]))
            nmini_pos <- which(meta_files[good_pos, CellspML] != min(meta_files[good_pos, CellspML]))

            decision[good_pos[mini_pos]] <- "Retain"

            decision[good_pos[nmini_pos]] <- "No!"


        } else {

            decision <- "All Files are bad!"

        }

    } else if (make_decision == "maxi") {

        if (are_all_good == 1) {

            decision <- ifelse(meta_files[Status] == "good", "Retain", "No!")

        } else if (are_all_good >= 2) {

            maxi_pos <- which(meta_files[good_pos, CellspML] == max(meta_files[good_pos, CellspML]))
            nmaxi_pos <- which(meta_files[good_pos, CellspML] != max(meta_files[good_pos, CellspML]))

            decision[good_pos[maxi_pos]] <- "Retain"

            decision[good_pos[nmaxi_pos]] <- "No!"


        } else {

            decision <- "All Files are bad!"

        }

    } else if (make_decision == "unique") {

                decision[good_pos] <- "Retain"

    } else stop("Supply make_decision")


    return(decision)
}
