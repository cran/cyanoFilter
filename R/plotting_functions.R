#' plots the expression matrix of a flowframe. Note that, it takes some time to display the plot.
#'
#' @param flowfile flowframe to be plotted
#' @param notToPlot column in expression matrix not to be plotted
#'
#' @examples \donttest{
#' flowfile_path <- system.file("extdata", "text.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' pair_plot(x = flowfile_logtrans,
#'           notToPlot = c("TIME", "FSC.HLin", "RED.R.HLin", "NIR.R.HLin"))
#'
#' }
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline points panel.smooth pairs smoothScatter text
#' @export pair_plot

pair_plot <- function(flowfile, notToPlot = c("TIME")) {
    toplot <- setdiff(flowCore::colnames(flowfile), notToPlot)
    col.palette <- colorRampPalette(c("white", "blue", "cyan", "green", "orange", "red"),
                                    space = "Lab")
    pairs(flowCore::exprs(flowfile)[, toplot], pch = ".",
          panel = function(...) smoothScatter(..., nrpoints = 0,
                  colramp = col.palette,
                  add = TRUE), gap = 0.2, main = flowCore::identifier(flowfile))
}


#' plots the expression matrix of a flowframe analysed with celldebris_emclustering.
#'
#' @param flowfile flowframe to be plotted
#' @param channels channels used in gating
#' @param mus matrix of means obtained from celldebris_emclustering
#' @param tau vector of cluster weights obtained from celldebris_emclustering
#' @param classifier cells will be assigned to a cluster if belongs to that cluster
#'                    by at least this probability. Only for plotting purposes.
#'
#' @importFrom grDevices chull densCols
#' @importFrom graphics plot polygon
#' @export cluster_plot


cluster_plot <- function(flowfile, channels,
                         mus = NULL, tau = NULL,
                         classifier) {

  if(is.null(mus) | is.null(tau)) stop("supply the matrix of mean and
                                        vector of percentages obtained
                                       from the celldebris_emclustering function")

    ddata <- flowCore::exprs(flowfile)

    ddata2 <- ddata[, stringr::str_detect(colnames(ddata), "Prob") ]

    color_code <- apply(ddata2, 1, function(x) {

        rest <- which(x >= classifier & x == max(x))
        # maximum = not sure in other words NA
        frest <- ifelse(length(rest) == 0, 0, rest)

        return(frest)

    })

    #plotting
    color_seq <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "seq", ]
    panel.xy <- function(x, y, ...) {

      for(i in unique(color_code)) {
        if(i != 0) {

          cols.pal <- RColorBrewer::brewer.pal(n = 9, row.names(color_seq)[i])
          col.palette <- colorRampPalette(c("white", cols.pal))
          pal <- densCols(x[color_code == i], y[color_code == i],
                          colramp =  col.palette)
          points(x[color_code == i], y[color_code == i],
                 col = pal, ...)

          #drawign cluster boundaries
          if(i == which(tau == max(tau)) ) {
          plot_data <- cbind(x[color_code == i], y[color_code == i])
          convhull <- chull(plot_data)
          polygon(plot_data[convhull, ], lty = 4, lwd = 2,
                  border = "red"
                  ) }


        } else {

          cols.pal <- RColorBrewer::brewer.pal(n = 9, "Greys")
          col.palette <- colorRampPalette(c("white", cols.pal))
          pal <- densCols(x[color_code == i], y[color_code == i],
                          colramp =  col.palette)
          points(x[color_code == i], y[color_code == i],
                 col = pal, ...)

        }

      }

    }

    pairs(ddata[, channels], pch = ".",
          panel = panel.xy, gap = 0.2,
          main = flowCore::identifier(flowfile)
      )

    msg <- paste("points are assigned to a cluster if they have at least",
                 classifier, "of belonging to that cluster.
                 The boundary of the largest cluster is drawn but the
                 colors differentiate each cluster", sep = " ")
    message(msg)


}
