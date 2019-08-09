#' Removes or assign indicators to margin events.
#'
#' @param flow.frame Flowframe containing margin events to be filtered out
#' @param Channel The channel on which margin events are. Defaults to SSC.W (side scatter width)
#' @param type The method to be used in gating out the margin cells. Can either be 'manual' where
#' user supplies a cut off point on the channel, 1 = not margin 0 = margin
#' @param cut sould not be NULL if type = 'manual'
#' @param y_toplot channel on y-axis of plot with \emph{Channel} used to gate out margin events
#'
#' @return list containing; \itemize{
#' \item \strong{reducedflowframe -} flowframe without margin events
#' \item \strong{fullflowframe -} flowframe with an Margin.Indicator added as an extra column added to the expression matrix
#' to indicate which particles are margin events. 1 = not margin event, 0 = margin event
#' \item \strong{N_margin -} number of margin events recorded
#' \item \strong{N_cell -} numner of non-margin events
#' \item \strong{N_particle -} is the number of particles in total, i.e. N_cell + N_margin
#' }
#'
#' @description The function identifies margin events, i.e. cells that are too large for the flow cytometer to measure.
#'
#' @details Users can either supply a cut-off point along the channel describing particle width or allow the function to estimate the cut-off point using the
#'          \code{\link[flowDensity]{deGate}} function from the \emph{flowDensity} package. A plot of channel against "FSC.HLin" is provided with a vertical line showing the
#'          cut-off point separating margin events from other cells.
#'
#' @examples
#' flowfile_path <- system.file("extdata", "text.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' cellmargin(flow.frame = flowfile_logtrans, Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#'
#'
#' @importFrom methods new
#' @export cellmargin

cellmargin <- function(flow.frame, Channel = "SSC.W", type = c("manual", "estimate"), cut = NULL, y_toplot = "FSC,HLin") {

    dvarMetadata <- flow.frame@parameters@varMetadata
    ddimnames <- flow.frame@parameters@dimLabels
    describe <- flow.frame@description

    if (type == "manual" & !is.null(cut)) {
        margin.ind <- ifelse(flowCore::exprs(flow.frame)[, Channel] <= cut, T, F)
        n_margin <- sum(margin.ind == F)

        # plotting
        flowDensity::plotDens(flow.frame, c(Channel, y_toplot), main = flowCore::identifier(flow.frame))
        abline(v = cut, lwd = 1, lty = 4, col = 2)

        # constructing full flow frame with both margin and non-margin events
        exx <- as.matrix(cbind(flowCore::exprs(flow.frame), as.numeric(margin.ind)))  #1=not margin, 0=margin
        colnames(exx)[ncol(exx)] <- "Margin.Indicator"

        # constructing the annotated data frame for the parameter
        ddata <- data.frame(rbind(flow.frame@parameters@data, c("Margin.Indicator", "Margin Indicator", 1, 0, 1)))
        paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata, dimLabels = ddimnames)
        row.names(ddata) <- c(row.names(flow.frame@parameters@data), paste0("$", "P", ncol(exx)))

        fflowframe <- flowCore::flowFrame(exprs = exx, parameters = paraa, description = describe)

        # constructing flowframe with only the non margin events
        exx2 <- exx[exx[, "Margin.Indicator"] == 1, ]
        rflowframe <- flowCore::flowFrame(exprs = exx2, parameters = paraa, description = describe)
    } else if (type == "estimate") {

        infl_point1 <- flowDensity::deGate(flow.frame, Channel, bimodal = T) #upper boundary
        infl_point2 <- flowDensity::deGate(flow.frame, Channel, use.upper = T, upper = F) #lower boundary

        # plotting
        flowDensity::plotDens(flow.frame, c(Channel, y_toplot), main = flowCore::identifier(flow.frame))
        abline(v = infl_point1, lwd = 1, lty = 4, col = 2)
        abline(v = infl_point2, lwd = 1, lty = 4, col = 2)

        margin.ind <- ifelse(flowCore::exprs(flow.frame)[, Channel] > infl_point2 &
                             flowCore::exprs(flow.frame)[, Channel] < infl_point1, T, F)  #FALSE = Margin Event
        n_margin <- sum(margin.ind == F)  #part of output

        # constructing full flow frame with both margin and non-margin events
        exx <- as.matrix(cbind(flowCore::exprs(flow.frame), as.numeric(margin.ind)))  #1=not margin, 0=margin
        colnames(exx)[ncol(exx)] <- "Margin.Indicator"  #1 = not amrgin, 0 = margin

        # constructing the annotated data frame for the parameter
        ddata <- data.frame(rbind(flow.frame@parameters@data, c("Margin.Indicator", "Margin Indicator", 1, 0, 1)))
        paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata, dimLabels = ddimnames)
        row.names(ddata) <- c(row.names(flow.frame@parameters@data), paste0("$", "P", ncol(exx)))

        fflowframe <- flowCore::flowFrame(exprs = exx, parameters = paraa, description = describe)

        # constructing flowframe with only the non margin events
        exx2 <- exx[exx[, "Margin.Indicator"] == 1, ]
        rflowframe <- flowCore::flowFrame(exprs = exx2, parameters = paraa, description = describe)
    } else stop("Error: check your inputs")

    return(list(reducedflowframe = rflowframe, fullflowframe = fflowframe, N_margin = n_margin, N_nonmargin = flowCore::nrow(rflowframe),
                N_particle = flowCore::nrow(fflowframe)))
}
