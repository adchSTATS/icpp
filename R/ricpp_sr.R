#' Simulate iterated cluster point prcocesses with independent imigration
#'
#' Simulate iterated cluster point prcocesses with independent imigration under same reproduction.
#'
#' @param intens Non-negative numeric. The desired intensity for each generation.
#' @param clsiz Numeric between 0 and 1. The expected cluster size.
#' Must be between 0 and 1, due to the assumption of same reproduction.
#' @param sd Non-negative numeric. The offspring density standard deviation.
#' @param noise String describing the type of noise to be added. \code{"Det"}, \code{"Pois"}, or \code{"Per"}.
#' @param win Window of class \code{\link{owin}}.
#' @param kernel_sd Non-negative numeric. Kernel standard deviation.
#' Default is \code{NULL} in which case \code{noise} must be \code{"Pois"}.
#' @param clsiz_var Non-negative numeric. Variance of the number of points in a cluster. See \code{\link{rGaussian_child}}.
#' @param n Positive integer. The number of iterations.
#' @param tol Positive real number.
#' @param ncores Positive integer specifying the number of cores to be used.
#' See the argument \code{mc.cores} in \code{\link{mclapply}}.
#' @import spatstat parallel
#' @export
ricpp_sr <- function(intens, clsiz, sd, noise, win = owin(), kernel_sd = NULL, clsiz_var = NULL, n = NULL, tol = .Machine$double.eps, ncores = 1){
  if (intens < 0) stop("The constant intensity must be a non-negative number.")
  if (sd < 0) stop("The standard deviation for the offspring density must be non-negative.")
  if (clsiz < 0 || clsiz > 1) stop("The expected cluster size (clsiz) must be between 0 and 1.")
  if (!(noise %in% c("det", "pois", "per"))) stop("Type of noise process not recognized. Must be 'Det', 'Pois', or 'Per'")
  if (noise %in% c("det", "per") && is.null(kernel_sd)) stop("Since noise is either 'det' or 'per', 'kernel_sd must be specified.")

  # Determine n
  noise_int <- intens * (1 - clsiz)
  if(is.null(n)){
    n <- 0
    exp_int <- 0
    while (!isTRUE(all.equal(exp_int, intens, tolerance = tol))) {
      n <- n + 1
      exp_int <- noise_int * sum(rep(clsiz, n + 1) ^ (0:n))
    }
  }
  cat(n, "\n")

  # Do the simulations
  noise_pp <- mclapply(0:n, function(i) {
    win_expanded <- expand.owin(win, distance = 4 * sqrt(i) * sd)
    switch(noise,
           pois = {
             out <- rpoispp(noise_int, win = win_expanded)
             cat(i)
             out
           },
           det = {
             catch <- "no"
             while (catch != "ppp") {
               out <- simulate.dppm(dppGauss(lambda = noise_int, alpha = kernel_sd, d = 2), W = win_expanded)
               catch <- class(out)
             }
             cat(i)
             out
           },
           per = {
             out <- rPPP(noise_int, kernel_sd, alpha = 0.5, win = win_expanded)
             cat(i)
             out
           })
  }, mc.cores = ncores)
  current <- noise_pp[[n + 1]]
  for (i in 0:(n - 1)) {
    current <- rGaussian_child(current, sd = sd, mean = clsiz, var = clsiz_var, ncores = ncores)
    noise_pp_tmp <- noise_pp[[n - i]]
    current <- superimpose(current, noise_pp_tmp, W = noise_pp_tmp$window)
    cat("Iteration:", i + 1, "\n",
        "Number of points:", npoints(current), "\n")
  }
  current[win]
}
