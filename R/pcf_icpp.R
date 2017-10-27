#' PCF for iterated cluster point process
#'
#' Closure function (ie. returns a function) that:
#' Calculats the pcf for the n'th generation iterated Neyman-Scott type process with or without noise.
#' It is asusmed that nothing except for the intensity depends on the generation number.
#' The offspring density is assumed to be Gaussian.
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

  # Determine the combinations to return
  combs_tmp <- expand.grid(initial, noise)
  combs <- paste(combs_tmp[, 1], combs_tmp[, 2], sep = "_")

  # Common preliminary calculations
  out <- numeric(length(combs))
  names(out) <- combs
  index <- 0:(n - 1)
  clsiz_power <- clsiz ^ index
  sd_pois <- sqrt(2 * (n - index) * sd ^ 2)

  # Initial specific calculations
  a_det <- -pi * kernel_sd_det_init ^ 2 / weight_det_init
  a_per <- pi * kernel_sd_per_init ^ 2 / weight_per_init
  sd_det_per_init <- sqrt(c(kernel_sd_det_init, kernel_sd_per_init) ^ 2 / 4 + 2 * n * sd ^ 2)
  sd_det_init <- sd_det_per_init[1]
  sd_per_init <- sd_det_per_init[2]

  # Noise specific calculations
  no_noise <- "non" %in% noise || noise_int == 0
  with_noise <- any(noise %in% c("det", "Pois", "per")) && noise_int != 0
  if (with_noise) {
    b_det <- -pi * kernel_sd_det_noise ^ 2 / weight_det_noise
    b_per <- pi * kernel_sd_per_noise ^ 2 / weight_per_noise
    sd_det_noise <- sqrt(kernel_sd_det_noise ^ 2 / 4 + 2 * (n - index[-n]) * sd ^ 2)
    sd_per_noise <- sqrt(kernel_sd_per_noise ^ 2 / 4 + 2 * (n - index[-n]) * sd ^ 2)
    int_gen_i <- c(init_int, sapply(index, function(n) {
      if (n == 0) {
        sum_term <- 0
      } else {
        i <- 1:n
        sum_term <- sum(clsiz ^ (n + 1 - i))
      }
      noise_int * (1 + sum_term) + init_int * clsiz ^ (n + 1)
    }))
    clsiz_squared <- clsiz ^ 2
  }

  function(r) {
    # Define constants that are common for the case with and without noise
    common_conv <- my2ddnorm(r, sd = sd_pois)
    if ("det" %in% initial) {
      init_det_conv <- a_det * my2ddnorm(r, sd = sd_det_init)
    }
    if ("per" %in% initial) {
      init_per_conv <- a_per * my2ddnorm(r, sd = sd_per_init)
    }

    # Calculations when there is no noise
    if (no_noise) {
      Pois_non = 1 + c / init_int * sum(common_conv / clsiz_power)
    }

    # Calculations when there is noise
    if (with_noise) {
      Pois_Pois <- 1 + c / int_gen_i[n + 1] ^ 2 * sum(int_gen_i[-(n + 1)] * clsiz_squared ^ (n - index) * common_conv)

      # Calculations specific for choice of initial process
      if (any(c("det", "per") %in% initial)) {
        init_const <- init_int ^ 2 / int_gen_i[n + 1] ^ 2 *clsiz_squared ^ n
        if ("det" %in% initial) {
          det_init <- init_const * init_det_conv
        }
        if ("per" %in% initial) {
          per_init <- init_const * init_per_conv
        }
      }

      # Calculations specific for choice of noise process
      if (any(c("det", "per") %in% noise)) {
        noise_const <- (noise_int / int_gen_i[n + 1]) ^ 2
        if("det" %in% noise) {
          det_noise <- sum(clsiz_squared ^ (n - index[-n]) * b_det * my2ddnorm(r, sd = sd_det_noise))
          det_noise <- det_noise + b_det * my2ddnorm(r, sd = kernel_sd_det_noise ^ 2 / 4)
          det_noise <- noise_const * det_noise
        }
        if("per" %in% noise) {
          per_noise <- sum(clsiz_squared ^ (n - index[-n]) * b_per * my2ddnorm(r, sd = sd_per_noise))
          per_noise <- per_noise + b_per * my2ddnorm(r, sd = kernel_sd_per_noise ^ 2 / 4)
          per_noise <- noise_const * per_noise
        }
      }
    }

    # Final output calculations
    sapply(combs, function(x) {
      switch(x,
             Pois_non = Pois_non,
             det_non = Pois_non + init_det_conv,
             per_non = Pois_non + init_per_conv,
             Pois_Pois = Pois_Pois,
             det_Pois = Pois_Pois + det_init,
             per_Pois = Pois_Pois + per_init,
             Pois_det = Pois_Pois + det_noise,
             det_det = Pois_Pois + det_init + det_noise,
             per_det = Pois_Pois + per_init + det_noise,
             Pois_per = Pois_Pois + per_noise,
             det_per = Pois_Pois + det_init + per_noise,
             per_per = Pois_Pois + per_init + per_noise)
    })
  }
}
