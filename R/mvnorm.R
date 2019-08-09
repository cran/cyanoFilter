#' multivariate normal density
#'
#' @param x matrix to compute density on
#' @param mu mean vector
#' @param sigma variance covariance matrix
#' @return vector of density
#'
#' @export mvnorm

# mnrom computes the probability of an event
mvnorm <- function(x, mu, sigma){
  x_center = sweep(x, 2, mu) # subtract mean, normalizing
  sig_inv <- sigma * NaN
  try(sig_inv <- solve(sigma), silent = T)
  dens <- exp(-0.5 * rowSums((x_center %*% sig_inv) * x_center))/
    sqrt(abs(det(2 * pi * sigma)))
  return (dens)
}
