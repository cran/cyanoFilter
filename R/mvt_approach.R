#' identifies BS4, BS5 and Debris in a flowfile using an EM style algorithm.
#'
#' @param flowfile flowframe to be clustered.
#' @param channels channels to use for the clustering
#' @param mu pre-specified mean matrix for the clusters. Number of rows should equal ncluster and number
#'           of columns should equal length(channels). Defaults to NULL and will be computed from the data internally if left as NULL.
#' @param sigma pre-specified list of variance-covariace matrix for the clusters. Each element of the list should contain a square matrix
#'        of variance-covariance matrix with length equal ncluster. Defaults to NULL and will be computed from the data internally if left as NULL.
#' @param  ncluster number of cluster desired.
#' @param  min.itera minimum number of EM iterations.
#'
#' @return list containing; \itemize{
#' \item \strong{percentages -} percentage of cells in each cluster
#' \item \strong{mus -} matrix of mean vectors for each cluster
#' \item \strong{sigmas -} list of variance-covariance matrix for each cluster
#' \item \strong{result -} flowframe with probabilities of each cluster added as columns to the expression matrix of the flowfile
#' }
#'
#' @description separates BS4, BS5 and Debris population in a flowfile using an EM style algorithm.
#'              Algorithm starts with \emph{ncluster} number of clusters and automatically reduces
#'              this number if need be.
#'
#' @details The function using EM algorithm involving mixtures of multivariate normals to separate the entire cell-population provided into cluster. The \code{\link{mvnorm}}
#'          function is used to compute the densities and only the probabilites of each point belonging to a cluster are returned as additional columns to the
#'          expression matrix of \emph{result}.
#'
#' @seealso \code{\link{celldebris_nc}}
#'
#' @examples
#'\donttest{
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
#'
#' emapproach <- cyano_emclustering(flowfile = cells_nonmargin, channels = c("RED.B.HLin",
#'                    "YEL.B.HLin", "FSC.HLin", "RED.R.HLin"),
#'                    ncluster = 5, min.itera = 20)
#' }
#'
#'
#' @importFrom stats var cov quantile runif
#' @export celldebris_emclustering




celldebris_emclustering <- function(flowfile, channels, mu = NULL, sigma = NULL, ncluster = 5,
                       min.itera = 20) {

      data <- flowCore::exprs(flowfile)[, channels]


      if(is.null(mu)) {
        #generate a uniform data
        rands <- matrix(NA, nrow = length(channels), ncol = ncluster)
        rand_mu <- apply(rands, 2, function(x){runif(length(channels), 0, 1)})
        boundaries <- apply(data, 2, quantile, probs = c(0.05, 0.95) )
        mu <- boundaries[1, ] + rand_mu *(boundaries[2, ] - boundaries[1, ])
        #mu <- matrix(, nrow = ncluster, ncol = length(channels))

      } else {
          mu <- mu
      }

      if(is.null(sigma)) {
          sigma <- vector("list", length = ncluster)
          sigma <- lapply(sigma, function(i){ cov(data)} )
      } else {
        sigma <- sigma
      }

      # probability of each cell to belong to each cluster
      lambda <- matrix(NA, nrow = nrow(data), ncol = ncluster)
      tau <- rep(1/ncluster, ncluster) # wheight of each cluster

      row.names(mu) <- channels

    i <- 0

    while(TRUE) {
      i <- i + 1
      # to check progress
      tau.old <- tau
      mu.old <- mu
      # compute the probability of each cell for each cluster (E step)
      for (p in 1:ncluster) {
        lambda[,p] <- tau[p] * mvnorm(data, mu[,p], sigma[[p]])
      }
      if(anyNA(lambda)) { # one cluster has size 0, repeat algorithm
        #clusters.pres <- clusters.pres[!is.na(lambda[1,clusters.pres])]
        return( celldebris_emclustering(flowfile, ncluster = ncluster - 1))
    }

    # each cell should have probability 1
    lambda <- lambda/rowSums(lambda)

    # probabilities/sizes of the different clusters
    tau <- colMeans(lambda)

    # adapt the mean and the std for the different clusters (M - step)
    for(p in 1:ncluster) {
        mu[,p] <- colSums(lambda[, p] * data) / sum(lambda[, p])
        sigma[[p]] <- t(lambda[,p] * sweep(data, 2, mu[,p])) %*% sweep(data, 2, mu[, p]) /
                    sum(lambda[, p])
    }
    # check whether local maxima has been achieved
    rel.diff.mu <- max((abs(mu - mu.old) / mu))
    rel.diff.tau <- max((abs(tau - tau.old) / tau))

    if (i > min.itera & rel.diff.mu < 1e-3 & rel.diff.tau < 1e-3) break

  }

    ddata <- data.frame()
    for(i in 1:ncol(lambda)) {

      ddata <- data.frame(rbind(ddata, data.frame(name = paste("Cluster_Prob", i, sep = "_"),
                                                  desc = paste("Cluster_Prob", i, sep = "_"), range = 1, minRange = range(lambda[, i])[1],
                                                  maxRange = range(lambda[, i])[2]))
      )
    }

    ddata <- rbind(flowfile@parameters@data, ddata)

    dvarMetadata <- flowfile@parameters@varMetadata
    ddimnames <- flowfile@parameters@dimLabels
    paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata, dimLabels = ddimnames)
    describe <- flowfile@description
    row.names(ddata) <- c(row.names(flowfile@parameters@data),
                          paste("$P", length(row.names(flowfile@parameters@data))+ 1:ncluster, sep = ""))

    ### Full flowframe Forming a new expression matrix for the full flowframe with indicator added for BS4 or BS5
    nexp_mat <- as.matrix(cbind(flowCore::exprs(flowfile), lambda))
    # giving a name to the newly added column to the expression matrix
    colnames(nexp_mat) <- ddata$name

    # full flow frame with indicator for particly type
    fflowframe <- methods::new("flowFrame", exprs = nexp_mat, parameters = paraa, description = describe)
    cluster_plot(fflowframe, channel1 = channels[1], channel2 = channels[2], mus = mu, tau = tau)

    return(list(percentages = tau, mus = mu, sigmas = sigma,
                  result = fflowframe))
}
