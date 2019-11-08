#' gates out or assign indicators to debris particle.
#'
#' @param flowframe flowframe with debris, BS4, BS5 and other cells.
#' @param p1 first flowcytometer channel that can be used to separate debris
#'           from the rest, e.g. "RED.B.HLin".
#' @param p2 second flowcytometer channel that can be used to separate debris
#'           from the rest, e.g. "YEL.B.HLin"
#'
#' @return list containing; \itemize{
#' \item \strong{syn - flowframe containing non-debris particles}
#' \item \strong{deb_pos - position of particles that are debris}
#' \item \strong{syn_pos - position of particles that are not debris}
#' }
#'
#' @description The function takes in a flowframe and identifies debris contained in the
#'              provided flowframe.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and
#'          \code{\link[flowDensity]{deGate}} functions in the flowDensity package to
#'          identify peaks between peaks and identify cut-off points between these peaks.
#'          A plot of both channels supplied with horizontal line separating
#'          debris from other cell populations is also returned.
#'
#' @seealso \code{\link{debris_inc}}
#'
#'
#'
#' @export debris_nc

debris_nc <- function(flowframe, p1, p2) {

    # debris gating
    p1_peaks <- flowDensity::getPeaks(flowframe, p1)

    # plotting
    flowDensity::plotDens(flowframe, c(p1, p2), main = paste(flowCore::identifier(flowframe),
                                                             length(p1_peaks$Peaks), sep = "-"),
                          frame.plot = F)

    if (length(p1_peaks$Peaks) == 2) {
        # all is well
        deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)[1]

    } else if (length(p1_peaks$Peaks) > 2) {
        # invader present deb_cut <- flowDensity::deGate(flowframe, p1, bimodal = F)#p1_peaks$Peaks[2]
        deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)[2]

    } else if (length(p1_peaks$Peaks) < 2) {
        # not much debris, hence one peak
        deb_cut <- flowDensity::deGate(flowframe, p1, sd.threshold = T, n.sd = 1)
        if (deb_cut < 2) {
            deb_cut <- flowDensity::deGate(flowframe, p1, after.peak = T)
        } else if (deb_cut >= 2 & deb_cut >= p1_peaks$Peaks[1]) {
            deb_cut <- flowDensity::deGate(flowframe, p1, use.upper = T, upper = F)
        } else {
            deb_cut <- flowDensity::deGate(flowframe, p1, sd.threshold = T, n.sd = 1)
        }
    }

    abline(v = deb_cut, lty = 2, col = 2)

    bs4bs5 <- flowframe[which(flowCore::exprs(flowframe)[, p1] > deb_cut), ]
    other_pos <- which(flowCore::exprs(flowframe)[, p1] > deb_cut)
    deb_pos <- which(flowCore::exprs(flowframe)[, p1] <= deb_cut)

    text(x = mean(flowframe@exprs[which(flowCore::exprs(flowframe)[, p1] <= deb_cut), p1]),
         y = mean(flowframe@exprs[which(flowCore::exprs(flowframe)[, p2] <=
        deb_cut), p2]), "Deb", col = 2)

    return(list(syn = bs4bs5, deb_pos = deb_pos, syn_all_pos = other_pos))
}
