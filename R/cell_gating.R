#' gates out or assign indicators to Synechococcus cyanobacteria cells.
#'
#' @param flowframe flowframe with debris and Synechococcus cells.
#' @param channel1 first flowcytometer channel that can be used to separate cyanobacteria
#'                 cells from the rest, e.g. "RED.B.HLin".
#' @param channel2 second flowcytometer channel that can be used to separate cyanobacteria
#'                 cells from the rest, e.g. "YEL.B.HLin"
#' @param interest a string indicating poistion of population of interest to be gated,
#'                 can be "bottom-right", "top-right" or "both-right".
#' @param to_retain should potential candidates be retained or further
#'                  gating be applied to filter out only certain cyano cells.
#' @return list containing; \itemize{
#' \item \strong{fullframe -} full flowframe with indicator for debris, BS4/BS5 or both.
#' \item \strong{reducedframe -} flowframe with onlySynechococcus cyanobacteria.
#' \item \item \strong{Cell_count -} a vector containing number of Synechococcus cyanobacteria
#'                                   cells. Might be a single value or vector of two values
#'                                   depending on interest.
#' \item \strong{Debris_count -} number of debris particles.
#' }
#'
#' @description This is a top-level function that calls other functions to identify cell population of interest.
#'
#' @details The indicators assigned to the "BS4BS5.Indicator" column in the full flowframe depends
#'          on the \emph{interest} supplied. For \emph{interest="bottom-right"} or
#'          \emph{interest="top-right"}; 0 = Debris, 1 = BS4/BS5, 2 = not-identified while for
#'          \emph{interest="Both"}, 0 = Debris, 1 = Syn-1, 2 = Syn-2, 3 = not-identified.
#'          This function calls the \code{\link{debris_nc}} or \code{\link{debris_inc}} function to
#'          identify debris and afterwards call the \code{\link{bs4_nc}} and/or \code{\link{bs5_nc}}
#'          depending on the interest supplied.
#'
#' @seealso \code{\link{celldebris_emclustering}}
#'
#' @examples
#' \donttest{
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flow.frame = flowfile_logtrans, Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' celldebris_nc(flowframe = cells_nonmargin$reducedframe, channel1 = "RED.B.HLin",
#'               channel2 = "YEL.B.HLin",
#'               interest = "bottom-right", to_retain = "refined")
#' }
#'
#' @import Biobase
#' @export celldebris_nc


celldebris_nc <- function(flowframe, channel1 = "RED.B.HLin", channel2 = "YEL.B.HLin",
                          interest = c("bottom-right", "top-right", "both-right"),
                          to_retain = c("refined", "potential")) {


    if (interest == "both-right") {

        ddata <- data.frame(rbind(flowframe@parameters@data, c("Indicator",
                                                               "Indicator", 1, 0, 1)))
        row.names(ddata) <- c(row.names(flowframe@parameters@data),
                              paste("$P", length(row.names(flowframe@parameters@data))+1,
                                    sep = ""))

    } else if (interest == "bottom-right" | interest == "top-right") {

        ddata <- data.frame(rbind(flowframe@parameters@data, c("Indicator",
                                                               "Indicator", 1, 0, 1)))
        row.names(ddata) <- c(row.names(flowframe@parameters@data),
                              paste("$P", length(row.names(flowframe@parameters@data))+1,
                                    sep = ""))

    } else stop("incorrect option supplied, check interest")

    dvarMetadata <- flowframe@parameters@varMetadata
    ddimnames <- flowframe@parameters@dimLabels
    paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata,
                                         dimLabels = ddimnames)
    describe <- flowframe@description
    BS4BS5_ind <- rep(NA, flowCore::nrow(flowframe))

    # gating debris
    if (interest == "both-right") {

        debris <- debris_inc(flowframe = flowframe, p1 = channel1, p2 = channel2)

    } else debris <- debris_nc(flowframe = flowframe, p1 = channel1, p2 = channel2)

    # BS4s and BS5s
    bs4bs5 <- debris$syn
    # BS4BS5 positions in the expression matrix
    bs4bs5_pos <- debris$syn_all_pos

    # Fill debris positions with 0
    BS4BS5_ind[debris$deb_pos] <- 0

    if (interest == "bottom-right") {

        bs4s <- bs4_nc(bs4bs5, p1 = channel1, p2 = channel2, others = bs4bs5_pos,
                       to_retain = to_retain)

        ### indicators for each cell type
        BS4BS5_ind[bs4s$syn_pos] <- 1  #BS4
        BS4BS5_ind[bs4s$others_nk] <- 2  #others not identified
        BS4BS5_ind[bs4s$others_nk2] <- 2  #others not identified

        # number of BS4
        n_cyano <- sum(BS4BS5_ind == 1)

    } else if (interest == "top-right") {

        bs5s <- bs5_nc(bs4bs5, p1 = channel1, p2 = channel2, others = bs4bs5_pos,
                       to_retain = to_retain)

        ### indicators for each cell type
        BS4BS5_ind[bs5s$syn_pos] <- 1  #BS5
        BS4BS5_ind[bs5s$others_nk] <- 2  #others not identified
        BS4BS5_ind[bs5s$others_nk2] <- 2  #others not identified

        # number of BS5
        n_cyano <- sum(BS4BS5_ind == 1)

    } else if (interest == "both-right") {

        # BS4
        bs4s <- bs4_nc(bs4bs5, p1 = channel1, p2 = channel2, others = bs4bs5_pos,
                       to_retain = to_retain)
        # BS5
        bs5s <- bs5_nc(bs4bs5, p1 = channel1, p2 = channel2, others = bs4bs5_pos,
                       to_retain = to_retain)

        ### indicators for each cell type
        BS4BS5_ind[bs4s$syn_pos] <- 1  #BS4
        BS4BS5_ind[bs5s$syn_pos] <- 2  #BS5
        BS4BS5_ind[is.na(BS4BS5_ind)] <- 3  #others not identified

        # number of BS4, number of BS5
        n_cyano <- c(BS4 = sum(BS4BS5_ind == 1), BS5 = sum(BS4BS5_ind == 2))


    } else stop("Supply interest")

    ### Full flowframe Forming a new expression matrix for the full flowframe with indicator added for BS4 or BS5
    nexp_mat <- as.matrix(cbind(flowCore::exprs(flowframe), as.numeric(BS4BS5_ind)))
    # giving a name to the newly added column to the expression matrix
    colnames(nexp_mat)[length(colnames(nexp_mat))] <- "Indicator"
    # full flow frame with indicator for particly type
    fflowframe <- flowCore::flowFrame(exprs = nexp_mat, parameters = paraa,
                                      description = describe)

    ### reduced flowframe
    if (interest == "both-right") {

        rframe <- fflowframe[fflowframe@exprs[, 13] %in% c(1, 2), ]

    } else {

        rframe <- fflowframe[fflowframe@exprs[, 13] == 1, ]
    }

    # number of debris
    n_debris <- sum(BS4BS5_ind == 0)

    return(list(fullframe = fflowframe, reducedframe = rframe, Cell_count = n_cyano,
                Debris_Count = n_debris))

}
