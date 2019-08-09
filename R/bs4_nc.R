#' gates out or assign indicators to BS4 cyanobacteria cells.
#'
#' @param bs4bs5 flowframe with debris removed.
#' @param p1 first flowcytometer channel that can be used to separate BS4 cells from the rest, e.g. "RED.B.HLin".
#' @param p2 second flowcytometer channel that can be used to separate BS4 cells from the rest, e.g. "YEL.B.HLin"
#' @param others row numbers for non-debris events. This is provided by the \code{\link{debris_nc}} or \code{\link{debris_inc}} function.
#' @param to_retain should potential candidates be retained or further gating be applied to filter out only certain BS4 cells.
#' @return list containing; \itemize{
#' \item \strong{bs4_reduced -} flowframe containing only BS4s
#' \item \strong{others_nk -} unidentified particle positions
#' \item \strong{bs4_pos -} BS4 positions
#' \item \strong{others_nk2 -} other unidentified particle positions
#' }
#'
#' @description This function takes in a flowframe with debris removed and identifies BS4 populations in the provided frame.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and \code{\link[flowDensity]{deGate}} functions in the \emph{flowDensity} package to
#'          identify peaks and identify cut-off points between these peaks. This function is not designed to be called in isolation, if called
#'          in isolation an error will be returned. It is preferably called on the results from \code{\link{debris_nc}} or \code{\link{debris_inc}} function. A graph with horizontal
#'          and vertical lines used in separating the populations is returned and if \emph{to_retain = "potential"}, all BS4 points are coloured red while a
#'          circle made of dashed lines is drawn around BS4 points if \emph{to_retain = "refined"}.
#'
#' @seealso \code{\link{bs5_nc}}
#'
#'
#' @importFrom utils capture.output
#' @export bs4_nc

bs4_nc <- function(bs4bs5, p1, p2, others, to_retain = c("refined", "potential")) {

    yel.bhlin_peaks <- flowDensity::getPeaks(bs4bs5, p2)
    red.bhlin_peaks <- flowDensity::getPeaks(bs4bs5, p1)
    msg <- capture.output(flowDensity::deGate(bs4bs5, p2, all.cuts = T))[1]

    if (length(yel.bhlin_peaks$Peaks) == 1) {

        yel_cuts <- flowDensity::deGate(bs4bs5, p2, use.percentile = T, percentile = 0.97)
        bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] < yel_cuts), ]
        others_pot <- others[which(bs4bs5@exprs[, p2] < yel_cuts)]
        others_nk <- others[which(!(bs4bs5@exprs[, p2] < yel_cuts))]
        ptt <- "1"

    } else if (length(yel.bhlin_peaks$Peaks) == 2) {

        mgroup1 <- yel.bhlin_peaks$Peaks[2] - 0
        if (mgroup1 < 2.9) {
            yel_cuts <- flowDensity::deGate(bs4bs5, p2, use.percentile = T, percentile = 0.99)
            bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] < yel_cuts), ]
            others_pot <- others[which(bs4bs5@exprs[, p2] < yel_cuts)]
            others_nk <- others[which(!(bs4bs5@exprs[, p2] < yel_cuts))]
            ptt <- "2a"

        } else {

            yel_cuts <- flowDensity::deGate(bs4bs5, p2, all.cuts = T)
            # BS4 is to the left of the minimum red_cuts
            bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] < yel_cuts), ]
            others_pot <- others[which(bs4bs5@exprs[, p2] < yel_cuts)]
            others_nk <- others[which(!(bs4bs5@exprs[, p2] < yel_cuts))]
            ptt <- "2b"

        }
    } else {

        yel_cuts <- flowDensity::deGate(bs4bs5, p2, use.percentile = T, percentile = 0.05)
        # min(yel.bhlin_peaks$Peaks) + 0.50*sd(bs4bs5@exprs[, p2])
        bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] > yel_cuts), ]
        others_pot <- others[which(bs4bs5@exprs[, p2] > yel_cuts)]
        others_nk <- others[which(!(bs4bs5@exprs[, p2] > yel_cuts))]
        ptt <- "3"

    }
    #
    abline(h = yel_cuts, lty = 2, col = 2)
    # points(bs4_pot@exprs[, c(p1, p2)], pch = '.', col = 2) text(0, 2.5, paste('BS4', ptt, sep = ':'), col = 2)

    if (to_retain == "refined") {

        bs4s <- flowDensity::flowDensity(bs4_pot, channels = c(p1, p2), position = c(F, F), use.upper = c(T, F),
                                         upper = c(T, NA),
                                         use.percentile = c(F, T),
                                         percentile = c(NA, 0.7),
                                         ellip.gate = T)
        points(bs4s@filter, type = "l", col = 2, lwd = 2, lty = 4)
        text(mean(bs4s@filter[, 1]), mean(bs4s@filter[, 2]), paste("BS4", ptt, sep = "-"), col = 2)
        # reduced flowframe for BS4
        bs4_reduced <- bs4_pot[which(!is.na(bs4s@flow.frame@exprs[, 1])), ]
        # positions of BS4s and other unidentified particles
        bs4_pos <- others_pot[which(!is.na(bs4s@flow.frame@exprs[, 1]))]
        others_nk2 <- others_pot[which(is.na(bs4s@flow.frame@exprs[, 1]))]

    } else if (to_retain == "potential") {

        text(mean(bs4_pot@exprs[, p1]), mean(bs4_pot@exprs[, p2]), "BS4", col = "red4")
        # reduced flowframe for BS4
        bs4_reduced <- bs4_pot
        # positions of BS4s and other unidentified particles
        others_nk <- others_nk
        bs4_pos <- others_pot
        others_nk2 <- NA

    } else stop("supply to_retain")

    return(list(bs4_reduced = bs4_reduced, others_nk = others_nk, bs4_pos = bs4_pos, others_nk2 = others_nk2))

}
