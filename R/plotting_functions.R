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
    col.palette <- colorRampPalette(c("white", "blue", "cyan", "green", "orange", "red"), space = "Lab")
    pairs(flowCore::exprs(flowfile)[, toplot], pch = ".",
          panel = function(...) smoothScatter(..., nrpoints = 0,
                  colramp = col.palette,
                  add = TRUE), gap = 0.2, main = flowCore::identifier(flowfile))
}


#' plots the expression matrix of a flowframe analysed with celldebris_emclustering.
#'
#' @param flowfile flowframe to be plotted
#' @param channel1 column in expression matrix not to be plotted
#' @param channel2 column in expression matrix not to be plotted
#' @param mus matrix of means obtained from celldebris_emclustering
#' @param tau vector of cluster weights obtained from celldebris_emclustering
#'
#' @importFrom grDevices chull densCols
#' @importFrom graphics plot polygon
#' @export cluster_plot


cluster_plot <- function(flowfile, channel1 = "RED.B.HLin", channel2 = "YEL.B.HLin", mus = NULL, tau = NULL) {

  if(is.null(mus) | is.null(tau)) stop("supply the matrix of mean and vector of percentages obtained
                                       from the celldebris_emclustering function")

    ddata <- flowCore::exprs(flowfile)

    ddata2 <- ddata[, stringr::str_detect(colnames(ddata), "Prob") ]

    color_code <- apply(ddata2, 1, function(x) {

        rest <- which(x >= 0.7 & x == max(x))
        frest <- ifelse(length(rest) == 0, ncol(ddata2)+1, rest) # maximum = not sure in other words NA

        return(frest)

    })

    #plotting
    cols.pal <- RColorBrewer::brewer.pal(n = length(unique(color_code)), "Dark2")
    plot(ddata[, c(channel1, channel2)], pch = ".", main = flowCore::identifier(flowfile),
         frame.plot = F, type = "n"
        )

    for(i in 1:length(unique(color_code))) {

      plotdata <- ddata[color_code==unique(color_code)[i], c(channel1, channel2)]
      pal <- colorRampPalette(c(cols.pal[i], "gray88", "green", "yellow", "orange", cols.pal[i]))
      col <- densCols(plotdata, colramp =  pal)

      points(plotdata, pch = ".", col =  col)
      #points to write cluster names
      x <- ifelse(class(try(mus[channel1, unique(color_code)[i]], silent = TRUE)) == "numeric",
                  mus[channel1, unique(color_code)[i]],
                  mean(plotdata[, channel1]) )

      y <- ifelse(class(try(mus[channel2, unique(color_code)[i]], silent = TRUE)) == "numeric",
                  mus[channel2, unique(color_code)[i]],
                  mean(plotdata[, channel2]) )
      #cluster names
      text(x = x,
           y = y,
           labels = paste("C", unique(color_code)[i], sep = ""),
           col = 2, offset = 0.0, cex = 0.8
           )
      #convexhull
        if(unique(color_code)[i] == which(tau == max(tau)) ) {
          convhull <- chull(plotdata[ ,
                                      c(channel1, channel2)])
          polygon(plotdata[convhull, ], lty = 4, lwd = 2, border = "red")
        }



    }

    message("points are assigned to a cluster if they have at least 70% probability of belonging to that
            cluster. The boundary of the largest cluster is drawn but the outermost colors represent the
            boundaries of each cluster")


}
