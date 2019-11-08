#' gates out or assign indicators to debris particle.
#'
#' @param flowframe flowframe with debris, BS4, BS5 and other cells.
#' @param p1 first flowcytometer channel that can be used to separate debris from the rest, e.g. "RED.B.HLin".
#' @param p2 second flowcytometer channel that can be used to separate debris from the rest, e.g. "YEL.B.HLin"
#'
#' @return list containing; \itemize{
#' \item \strong{syn - flowframe containing non-debris particles}
#' \item \strong{deb_pos - position of particles that are debris}
#' \item \strong{syn_pos - position of particles that are not debris}
#' }
#'
#' @description The function takes in a flowframe and identifies debris contained in the provided flowframe. It is specially designed for flowframe contaning both debris,
#'              BS4, BS5 and possibly other invading populations.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and \code{\link[flowDensity]{deGate}} functions in the flowDensity package to
#'          identify peaks between peaks and identify cut-off points between these peaks. A plot of both channels supplied with horizontal line separating
#'          debris from other cell populations is also returned.
#'
#' @seealso \code{\link{debris_nc}}
#'
#' @examples
#' \donttest{
#' flowfile_path <- system.file("extdata", "text.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flow.frame = flowfile_logtrans, Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' debris_inc(flowframe = flowfile, p1 = "RED.B.HLin", p2 = "YEL.B.HLin")
#' }
#'
#' @importFrom stats sd
#' @export debris_inc


debris_inc <- function(flowframe, p1, p2) {
    # plotting
    flowDensity::plotDens(flowframe, c(p1, p2), main = flowCore::identifier(flowframe), frame.plot = F)
    # debris gating
    p1_peaks <- flowDensity::getPeaks(flowframe, p1)
    if (length(p1_peaks$Peaks) == 2) {
        # all is well
        deb_cut <- flowDensity::deGate(flowframe, p1, bimodal = T)
        if (p1_peaks$Peaks[1] - (1 * sd(flowframe@exprs[, p1])) > 0 & deb_cut >= p1_peaks$Peaks[1]) {
            deb_cut <- flowDensity::deGate(flowframe, p1, use.upper = T, upper = F)
            ptt <- "1a"
        } else {

            deb_cut <- deb_cut
            ptt <- "1b"

        }

    } else if (length(p1_peaks$Peaks) == 3) {
        # invader present
        deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)
        if (p1_peaks$Peaks[1] - (1 * sd(flowframe@exprs[, p1])) <= 0) {

            deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)[1]
            ptt <- "2a"

        } else if (p1_peaks$Peaks[1] - (1 * sd(flowframe@exprs[, p1])) > 0) {

            deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)[1]

            flowframe2 <- flowframe[which(flowCore::exprs(flowframe)[, p1] > deb_cut), ]
            deb_cut <- ifelse(length(flowDensity::getPeaks(flowframe2, p1, tinypeak.removal = 1/5)$Peaks) < 2, flowDensity::deGate(flowframe, p1, use.upper = T,
                upper = F), deb_cut)
            ptt <- "2b"

        } else {

            deb_cut <- flowDensity::deGate(flowframe, p1, bimodal = F)
            ptt <- "2c"

        }


    } else if (length(p1_peaks$Peaks) < 2) {
        # not much debris, hence one peak
        deb_cut <- flowDensity::deGate(flowframe, p1, sd.threshold = T, n.sd = 1)
        if (deb_cut < 2) {

            deb_cut <- flowDensity::deGate(flowframe, p1, after.peak = T)
            ptt <- "3a"

        } else if (p1_peaks$Peaks[1] - (2 * sd(flowframe@exprs[, p1])) > 0 & deb_cut >= p1_peaks$Peaks[1]) {

            deb_cut <- flowDensity::deGate(flowframe, p1, use.upper = T, upper = F)
            ptt <- "3b"

        } else {

            deb_cut <- deb_cut
            ptt <- "3c"
        }
    } else {

        deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)[2]  #p1_peaks$Peaks[2] - (0.5 * sd(flowframe@exprs[, p1]))
        ptt <- "4"
    }

    abline(v = deb_cut, lty = 2, col = 2)
    # abline(v = flowDensity::deGate(flowframe, p1, all.cuts = T), lty = 2, col = 1)

    bs4bs5 <- flowframe[which(flowCore::exprs(flowframe)[, p1] > deb_cut), ]
    other_pos <- which(flowCore::exprs(flowframe)[, p1] > deb_cut)
    deb_pos <- which(flowCore::exprs(flowframe)[, p1] <= deb_cut)

    text(x = mean(flowframe@exprs[which(flowCore::exprs(flowframe)[, p1] <= deb_cut), p1]),
         y = mean(flowframe@exprs[which(flowCore::exprs(flowframe)[, p2] <=
        deb_cut), p2]), "Deb", col = 2)

    return(list(syn = bs4bs5, deb_pos = deb_pos, syn_pos = other_pos))
}
