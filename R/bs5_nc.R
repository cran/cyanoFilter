#' gates out or assign indicators to BS5 cyanobacteria cells.
#'
#' @param bs4bs5 flowframe with debris removed.
#' @param p1 first flowcytometer channel that can be used to separate BS5 cells from the rest, e.g. "RED.B.HLin".
#' @param p2 second flowcytometer channel that can be used to separate BS5 cells from the rest, e.g. "YEL.B.HLin"
#' @param others row numbers for non-debris events. This is provided by the debris_nc or debris_inc function.
#' @param  to_retain should potential candidates be retained or further gating be applied to filter out only certain BS5 cells.
#'  @return list containing; \itemize{
#' \item \strong{bs5_reduced -} flowframe containing only BS5s
#' \item \strong{others_nk -} unidentified particle positions
#' \item \strong{bs5_pos -} BS5 positions
#' \item \strong{others_nk2 -} other unidentified particle positions
#' }
#'
#' @description This function takes in a flowframe with debris removed and identifies BS5 populations in the provided frame.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and \code{\link[flowDensity]{deGate}} functions in the \emph{flowDensity} package to
#'          identify peaks and identify cut-off points between these peaks. This function is not designed to be called in isolation, if called
#'          in isolation an error will be returned. It is preferably called on the results from \code{\link{debris_nc}} or \code{\link{debris_inc}} function. A graph with horizontal
#'          and vertical lines used in separating the populations is returned and if \emph{to_retain="potential"}, all BS5 points are coloured red while a
#'          circle made of dashed lines is drawn around BS5 points if \emph{to_retain="refined"}.
#'
#'
#' @import Rdpack
#' @importFrom stats sd
#' @importFrom utils capture.output
#' @export bs5_nc


bs5_nc <- function(bs4bs5, p1, p2, others, to_retain = "potential") {

    yel.bhlin_peaks <- flowDensity::getPeaks(bs4bs5, p2)
    red.bhlin_peaks <- flowDensity::getPeaks(bs4bs5, p1)
    msg <- capture.output(flowDensity::deGate(bs4bs5, p1, all.cuts = T))[1]
    msg_yel <- capture.output(flowDensity::deGate(bs4bs5, p2, all.cuts = T))[1]

    ## checking for invasion on the YEL.B.HLin axis
    if (stringr::str_detect(msg_yel, "Only one peak") == T) {

        yel_cuts <- flowDensity::deGate(bs4bs5, p2, use.percentile = T, percentile = 0.03)
        bs5_pot <- bs4bs5[which(bs4bs5@exprs[, p2] > yel_cuts), ]
        others_pot <- others[which(bs4bs5@exprs[, p2] > yel_cuts)]
        others_nk <- others[which(!(bs4bs5@exprs[, p2] > yel_cuts))]
        ptt <- "1"

    } else if (length(yel.bhlin_peaks$Peaks) == 2) {

        mgroup1 <- yel.bhlin_peaks$Peaks[1] - 0
        if (mgroup1 > 2.9) {
            yel_cuts <- flowDensity::deGate(bs4bs5, p2, use.percentile = T, percentile = 0.03)
            bs5_pot <- bs4bs5[which(bs4bs5@exprs[, p2] > yel_cuts), ]
            others_pot <- others[which(bs4bs5@exprs[, p2] > yel_cuts)]
            others_nk <- others[which(!(bs4bs5@exprs[, p2] > yel_cuts))]
            ptt <- "2a"

        } else {

            yel_cuts <- flowDensity::deGate(bs4bs5, p2, all.cuts = T)
            bs5_pot <- bs4bs5[which(bs4bs5@exprs[, p2] > yel_cuts), ]
            others_pot <- others[which(bs4bs5@exprs[, p2] > yel_cuts)]
            others_nk <- others[which(!(bs4bs5@exprs[, p2] > yel_cuts))]
            ptt <- "2b"

        }
    } else {

        yel_cuts <- min(yel.bhlin_peaks$Peaks) + 0.5 * sd(bs4bs5@exprs[, p2])
        bs5_pot <- bs4bs5[which(bs4bs5@exprs[, p2] > yel_cuts), ]
        others_pot <- others[which(bs4bs5@exprs[, p2] > yel_cuts)]
        others_nk <- others[which(!(bs4bs5@exprs[, p2] > yel_cuts))]
        ptt <- "3"

    }

    abline(h = yel_cuts, lty = 2, col = 2)
    # points(bs5_pot@exprs[, c(p1, p2)], pch = '.', col = 2) text(0, 2.5, paste('BS5', ptt, sep = ':'), col = 2)

    if (to_retain == "refined") {
        bs5s <- flowDensity::flowDensity(bs5_pot, channels = c(p1, p2), position = c(F, F),
                                         ellip.gate = T, use.upper = c(T, T),
                                         upper = c(T, T))
        # plotting BS5
        points(bs5s@filter, type = "l", col = 2, lwd = 2, lty = 4)
        text(mean(bs5s@filter[, 1]), mean(bs5s@filter[, 2]), "BS5", col = 2)
        bs5_reduced <- bs5_pot[which(!is.na(bs5s@flow.frame@exprs[, 1])), ]
        bs5_pos <- others_pot[which(!is.na(bs5s@flow.frame@exprs[, 1]))]
        others_nk2 <- others_pot[which(is.na(bs5s@flow.frame@exprs[, 1]))]

    } else if (to_retain == "potential") {

        text(mean(bs5_pot@exprs[, p1]), mean(bs5_pot@exprs[, p2]), "BS5", col = "red4")
        bs5_reduced <- bs5_pot
        bs5_pos <- others_pot
        others_nk2 <- NA

    } else stop("supply to_retain")

    # yel_peaks <- flowDensity::getPeaks(bs5_pot, p2) if((length(yel_peaks$Peaks) == 2)) { if(abs(diff(yel_peaks$Peaks)) > sd(bs5_pot@exprs[, p2])) { bs5s <-
    # flowDensity::flowDensity(bs5_pot, channels = c(p1, p2), position = c(T, T), ellip.gate = T, use.upper = c(T, T), upper = c(F, F)) ptt <- '1a' } else { bs5s
    # <- flowDensity::flowDensity(bs5_pot, channels = c(p1, p2), position = c(T, T), ellip.gate = T, use.upper = c(T, T), upper = c(F, F)) ptt <- '1b' } } else {
    # bs5s <- flowDensity::flowDensity(bs5_pot, channels = c(p1, p2), position = c(T, T), ellip.gate = T, use.upper = c(T, T), upper = c(F, F)) ptt <- '2' }
    # #plotting BS5 points(bs5s@filter, type = 'l', col = 2, lwd = 2, lty = 4) text(mean(bs5s@filter[,1]), mean(bs5s@filter[, 2]), paste('BS5', ptt2, ptt, sep =
    # ':'), col = 2)

    return(list(bs5_reduced = bs5_reduced, others_nk = others_nk, others_nk2 = others_nk2, bs5_pos = bs5_pos))
}
