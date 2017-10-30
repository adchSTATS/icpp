#' 2D, uncorrelated, equal variance, Gaussian density function
#'
#' Function for calculating the density of a uncorrelated 2 dimensional normal distribution with equal variance.
#'
#' @param x vector of quantiles.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @return A vector of density values.
#' @importFrom stats dnorm
my2ddnorm <- function(x, mean = 0, sd = 1) {
  const <- sqrt(2 * pi * sd ^ 2)
  dnorm(x, mean = mean, sd = sd) / const
}
