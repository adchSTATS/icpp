#' Simulate a weighted-permanental point pattern
#'
#' Simulate an alpha-weighted-permanental point pattern when alpha is half an integer. This may be regarded as a Cox process.
#'
#' @param intens Non-negative numeric. Desired intensity.
#' @param kernel_sd Non-negative numeric. Kernel standard deviation.
#' @param alpha Half integer. Default is \code{1/2}. The weight of the weighted permanental point process.
#' @param win Window of class \code{\link{owin}}. Default is \code{owin()}.
#' @return A pixel image of class \code{\link{im}}.
#' @import spatstat
#' @importFrom RandomFields RFsimulate RMgauss
#' @export
rPPP <- function (intens, kernel_sd, alpha = 0.5, win = owin()) {
  alpha2 <- 2 * alpha
  if(alpha2 != round(alpha2)) stop("Alpha must be a half integer.")
  rfmodel <- RMgauss(var = intens / alpha2, scale = kernel_sd)
  w <- as.mask(w = win)
  z_squared_sum <- 0
  for (i in 1:alpha2) {
    spc <- RandomFields::RFoptions()$general$spConform
    if(spc) RandomFields::RFoptions(spConform=FALSE)
    z <- RFsimulate(rfmodel, w$xcol, w$yrow, grid = TRUE)
    if(spc) RandomFields::RFoptions(spConform=TRUE)
    z_squared_sum <- z_squared_sum + z^2
  }
  Lambda <- as.im(matrix(z_squared_sum, nrow = w$dim[1], ncol = w$dim[2], byrow = TRUE), W = w)
  rpoispp(Lambda)[win]
}
