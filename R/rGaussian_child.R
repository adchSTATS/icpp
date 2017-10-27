#' Generate clusters disperced with Gaussian density
#'
#' Generate clusters disperced by a normal distribution around a center point.
#'
#' @param parents An object of class \code{\link{ppp}}. Determining the possition of the centre points
#' @param sd A single non-negative numeric or a vector consisting of two non-negative numerics.
#' May also be a matrix such that squaring each entry yields the covariance matrix.
#' @param mean A non-negative numeric. The mean number of points in a single cluster.
#' @param var A non-negative numeric. The variance of the number of points in a cluster. See details.
#' @param ncores Positive integer specifying the number of cores to be used.
#' @details Number of point are Poisson distributed if \code{var} is \code{NULL} or identical to \code{mean}.
#' It is negativel binomially distributed if \code{var>mean} and binomially distributed if \code{var<mean}.
#' @return A point pattern of class \code{\link{ppp}} in the same window as \code{parents}.
#' @import spatstat parallel
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnbinom rpois
#' @export
rGaussian_child <- function(parents, sd, mean, var = NULL, ncores = 1) {
  verifyclass(parents, "ppp")
  pp <- do.call(
    rbind,
    mclapply(as.data.frame(t(as.data.frame(parents))), function(x) {
      if(is.null(var) || mean == var) {
        np <- rpois(1, mean)
      } else if(var > mean) {
        np <- rnbinom(1, size = mean ^ 2 / (var - mean), prob = mean / var)
      } else if(var < mean) {
        prob <- 1 - var/mean
        np <- rbinom(1, size = round(mean / prob), prob = prob)
      }
      if(np == 0) return()
      mvrnorm(n = np, mu = x, Sigma = sd ^ 2 * diag(1, 2, 2))
    }, mc.cores = ncores)
  )
  as.ppp(pp, W = parents$window)
}
