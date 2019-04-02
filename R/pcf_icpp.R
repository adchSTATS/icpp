#' PCF for iterated cluster point process
#'
#' Closure function (ie. returns a function) that:
#' Calculats the pcf for the n'th generation iterated Neyman-Scott type process with or without noise.
#' It is asusmed that nothing except for the intensity depends on the generation number.
#' The offspring density is assumed to be Gaussian.
#' The retention probability in 0, meaning that points do not survive to the next generation.
#'
#' @param n Positive integer. The number of iterations.
#' @param init_int Non-negative numeric. Intensity of the initial generation.
#' @param clsiz Non-negative numeric. The expected cluster size.
#' @param c Numeric. Constant related to the over-/underdispersion of number of point in a cluster.
#' Default is 1, corresponding to the number of points in a cluster follow a poisson distriution.
#' See article on iterated cluster point processes.
#' @param sd Non-negative number. Offspring density standard deviation.
#' @param initial Vector of strings or a single string. Decribing the type of the initial process.
#' Can be either "Pois" (default), "det", or "per".
#' @param kernel_sd_det_init Non-negative numeric. Kernel standard deviation.
#' (interaction 'radius' for the initial process).
#' Ignored if \code{initial} is "Pois".
#' @param kernel_sd_per_init Non-negative numeric.
#' Kernel standard deviation. (interaction 'radius' for the initial process).
#' Ignored if \code{initial} is "Pois".
#' @param weight_det_init A positive integer. The alpha weight for an alpha weighted DPP.
#' Default is 1 corresponding to the DPP (most repulsive gaussian DPP).
#' Ignored if \code{initial} is "Pois".
#' @param weight_per_init A positive half integer. The alpha weight for an alpha weighted PPP.
#' Default is 1/2 corresponding to the most clustered gaussian PPP which is a Cox process..
#' Ignored if \code{initial} is "Pois".
#' @param noise Vector of strings or a single string. Decribing the type of the noise process.
#' Can be either "non" (default) "det", "pois", or "per".
#' @param noise_int Non-negative numeric or the string \code{"sr"}.
#' If \code{"sr"} same reproduction system is assumed and \code{noise_int = init_int * (1 - clsiz)}.
#' Ignored if \code{noise} is "non".
#' @param kernel_sd_det_noise Ignored if \code{noise} is "non" or Pois".
#' @param kernel_sd_per_noise Ignored if \code{noise} is "non" or Pois".
#' @param weight_det_noise A positive integer. The alpha weight for an alpha weighted DPP.
#' Default is 1 corresponding to the DPP (most repulsive gaussian DPP).
#' Ignored if \code{noise} is "non" or Pois".
#' @param weight_per_noise A positive half integer. The alpha weight for an alpha weighted PPP.
#' Default is 1/2 corresponding to the most clustered gaussian PPP which is a Cox process.
#' Ignored if \code{noise} is "non" or Pois".
#' @details The function is a closure function.
#' The benefit of making a closure function here is that I can run the output function multiple times
#' without running initial checks and preliminary calculations.
#' We want this here because we have to compile the code for maybe 1000 values of r.
#'
#' The initial processes and noise processes may be either Poisson, weighted determinantal or weighted permanetal point process.
#' @return A function that takes inter-point distance (commonly referred to as r) as the only argument.
#' @export
pcf_icpp <- function(n,
                     init_int, clsiz, c = 1,
                     sd,
                     initial = "Pois",
                     kernel_sd_det_init = NULL,
                     kernel_sd_per_init = NULL,
                     weight_det_init = 1,
                     weight_per_init = 1/2,
                     noise = "non",
                     noise_int = 0,
                     kernel_sd_det_noise = NULL,
                     kernel_sd_per_noise = NULL,
                     weight_det_noise = 1,
                     weight_per_noise = 1/2) {
  # return proper warnings concerning generationindex
  if (n != round(n)) stop("The generation number (n) should be an integer.")

  # return proper warnings concerning intensities
  if (init_int < 0) stop("The initial intensity (init_int) should be a non-negative number.")
  if (clsiz < 0) stop("The expected cluster size (clsiz) should be a non-negative number.")
  if (sd < 0) stop("The off spring density standard deviatoin (sd) should be a non-negative number.")
  if (noise_int != "sr" && numeric(noise_int) && noise_int < 0) {
    stop("The initial intensity (init_int) should be a non-negative number or sr for same reproduction.")
  }
  if (noise_int == "sr") {
    if (clsiz > 1) stop("noise intensity was set to sr (same reproduction) and clsiz must therfore be less than 1.")
    noise_int <- init_int * (1 - clsiz)
  }

  # return proper warnings concerning the initial process
  initial <- unique(initial)
  if (!all(initial %in% c("det", "Pois", "per"))) stop("At least 1 element of initial is not recognized. Should be either det, Pois, or per.")
  if ("det" %in% initial && init_int != 0) {
    if (is.null(kernel_sd_det_init)) {
      kernel_sd_det_init <- max_kernel_sd_det(init_int)
    } else {
      if (kernel_sd_det_init < 0) stop("kernel_sd_det_init should be a non-negative number.")
      if (kernel_sd_det_init > max_kernel_sd_det(init_int, check = TRUE)) {
        stop("kernel_sd_det_init is too large to guarantee existence of the underlying DPP.")
      }
    }
    if (weight_det_init != round(weight_det_init) || weight_det_init < 0) {
      stop("weight_det_init should be a non-negative integer")
    }
  }
  if ("per" %in% initial && init_int != 0) {
    if (is.null(kernel_sd_per_init)) {
      kernel_sd_per_init <- max_kernel_sd_det(init_int)
    } else if (kernel_sd_per_init < 0) {
      stop("kernel_sd_per_init should be a non-negative number.")
    }
    if (weight_per_init < 0) stop("weight_per_init should be a non-negative numeric")
  }

  # return proper warnings concerning the noise processes
  noise <- unique(noise)
  if (length(noise) == 1 && noise == "non" && noise_int != 0) {
    warning("noise is set to non, but the intensity for the noise process is not 0. No noise will be added.")
  }
  if (!all(noise %in% c("non", "det", "Pois", "per"))) stop("At least 1 element of noise is not recognized. Should be either non, det, Pois, or per.")
  if (any(c("det", "Pois", "per") %in% noise) && noise_int == 0) {
    warning("The intensity of the noise processes is 0 and hence no noise will be added.")
  }
  if ("det" %in% noise && noise_int != 0) {
    if (is.null(kernel_sd_det_noise)) {
      kernel_sd_det_noise <- max_kernel_sd_det(noise_int)
    } else {
      if (kernel_sd_det_noise < 0) stop("kernel_sd_det_noise should be a non-negative number.")
      if (kernel_sd_det_noise > max_kernel_sd_det(noise_int, check = TRUE)) {
        stop("kernel_sd_det_noise is too large to guarantee existence of the underlying DPP.")
      }
    }
    if (weight_det_noise != round(weight_det_noise) || weight_det_noise < 0) {
      stop("weight_det_noise should be a non-negative integer")
    }
  }
  if ("per" %in% noise && noise_int != 0) {
    if (is.null(kernel_sd_per_noise)) {
      kernel_sd_per_noise <- max_kernel_sd_det(noise_int)
    } else if (kernel_sd_per_noise < 0) {
      stop("kernel_sd_per_noise should be a non-negative number.")
    }
    if (weight_per_noise < 0) stop("weight_per_noise should be a non-negative numeric")
  }

  if (noise == "non") noise_int <- 0
  if(noise_int == 0) noise <- "non"

  # Determine the combinations to return
  combs_tmp <- expand.grid(initial, noise)
  combs <- paste(combs_tmp[, 1], combs_tmp[, 2], sep = "_")

  # Common preliminary calculations
  out <- numeric(length(combs))
  names(out) <- combs
  index <- 0:(n - 1)
  clsiz_power <- clsiz ^ (0:n)
  sd_pois <- sqrt(2 * (index + 1) * sd ^ 2)
  int_gen_i <- init_int * clsiz_power + noise_int * cumsum(c(0, clsiz^index))

  # Initial specific calculations
  switch(initial,
         "det" = {
           a <- -pi * kernel_sd_det_init ^ 2 / (2 * weight_det_init)
           sd_init <- sqrt(kernel_sd_det_init^2 / 4 + 2 * n * sd^2)
         }, "per" = {
           a <- pi * kernel_sd_per_init ^ 2 / (2 * weight_per_init)
           sd_init <- sqrt(kernel_sd_per_init^2 / 4 + 2 * n * sd^2)
         }, {
           a <- 0
           sd_init <- 1
         })

  # Noise specific calculations
  switch (noise,
          "det" = {
            b <- -pi * kernel_sd_det_noise ^ 2 / (2 * weight_det_noise)
            sd_noise <- sqrt(kernel_sd_det_noise ^ 2 / 4 + 2 * index * sd ^ 2)
          }, "per" = {
            b <- pi * kernel_sd_per_noise ^ 2 / (2 * weight_per_noise)
            sd_noise <- sqrt(kernel_sd_per_noise ^ 2 / 4 + 2 * index * sd ^ 2)
          }, {
            b <- 0
            sd_noise <- 1
          })


  function(r) {
    # Common term
    density_convolution_mat <- sapply(sd_pois, FUN = function(sd) my2ddnorm(r, sd = sd))
    Pois_Pois <- as.vector(1 + (density_convolution_mat %*% (rev(int_gen_i[-(n+1)]) * clsiz_power[-1]^2)) / int_gen_i[n+1]^2)

    # Term for initial generation Initial
    density_convolution <- my2ddnorm(r, sd = sd_init)
    init_term <- (int_gen_i[1] / int_gen_i[n+1])^2 * a * density_convolution * clsiz_power[n+1]^2

    # Terms for noise
    density_convolution_mat <- sapply(sd_noise, FUN = function(sd) my2ddnorm(r, sd = sd))
    noise_term <- (noise_int / int_gen_i[n+1])^2 * b * as.vector(density_convolution_mat %*% clsiz_power[-(n+1)]^2)


    # Final output calculations
    Pois_Pois + init_term + noise_term
  }
}
