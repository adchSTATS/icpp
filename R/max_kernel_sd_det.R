#' Maximum kernel standard deviation for the Gaussian determinantal point process
#'
#' The Gaussian determinantal point process is specified by a zero-mean Guassian density function, more generally refered to as the kernel function.
#'
#' @param int The intensity of the point pattern.
#' @param check If \code{FALSE} a rounded value is returned. This ensures that the value will pass a check of whether the value is appropriate.
#' If thrue a non-rounded value is returned, which can be used to ensure that the condition is satisfied.
#' @details The variance determines the 'interaction radius' between points.
#' Large values of this variance means more repulsion for the weighted determinantal point process and more clustering for the weighted permanental point process.
#' It has a maximum for the determinantal point process that depends on the intensity, in order to ensure existence of the process.
#' The more points the less regularity is obtainable.
#' @return The maximal value of the kernel standard deviation, under the determinantal point process.
#' @export
max_kernel_sd_det <- function(int, check = FALSE) {
  if (check == TRUE) {
    1 / sqrt(pi * int)
  } else {
    floor(1 / sqrt(pi * int) * 1000000) / 1000000
  }
}
