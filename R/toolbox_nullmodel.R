#' Generate metapopulation
#' @param patch_num Number of patches
#' @param kernel_type which dispersal kernel to use
#' @param charactersitic_distance characteristic dispersal distance
#' @param location_distribution how patches are spatially distributed
#' @param self_regulation whether self-colonization exists
#' @param dispersal_rate area dependence of dispersal
#' @param extinction_rate area dependence of extinction
#' @param area_type distribution of patch area/quality
#' @param area_para parameter for area/quality distribution
#' @param ... other arguments passed to methods
#' @return a list of metapopulation parameters
#' @export
generate_metapopulation <- function(patch_num,
                                    kernel_type,
                                    charactersitic_distance,
                                    location_distribution,
                                    self_regulation,
                                    dispersal_rate,
                                    extinction_rate,
                                    area_type,
                                    area_para,
                                    ...) {
  dispersal_kernel <- generate_dispersal_kernel(
    kernel_type,
    charactersitic_distance
  )
  patch_area <- generate_patch_area(
    patch_num,
    area_type,
    area_para
  )
  patch_location <- generate_patch_location(
    location_distribution,
    patch_num
  )

  distance <- dist(patch_location, diag = T, upper = T) %>%
    as.matrix()

  A <- outer(patch_area^extinction_rate, patch_area^dispersal_rate)

  M <- calculate_landscape_matrix(
    distance = distance,
    patch_area = patch_area,
    dispersal_kernel = dispersal_kernel,
    dispersal_rate = dispersal_rate,
    extinction_rate = extinction_rate,
    self_regulation = self_regulation
  )

  D <- as.matrix(dispersal_kernel(distance))

  list(
    distance = distance,
    distance_matrix = D,
    patch_area = patch_area,
    area_matrix = A,
    landscape_matrix = M,
    patch_location = patch_location
  )
}


#' Generate dispersal kernels
#' @param kernel_type which kernel to use
#' @param para parameter in the kernel
#' @param ... other arguments passed to methods
#' @return a function of dispersal kernel
#' @export
generate_dispersal_kernel <- function(kernel_type, para, ...) {
  dispersal_linear <- function(para) {
    function(d) {
      d <- abs(d)
      ifelse(d < para, 1 - d / para, 0)
    }
  }

  dispersal_negexp <- function(para) {
    function(d) {
      exp(-para * abs(d))
    }
  }

  dispersal_gaussian <- function(para) {
    function(d) {
      exp(-abs(d^2) / (2 * para^2))
    }
  }

  dispersal_nonmono <- function(para, a = 2, b = 2) {
    function(d) {
      gamma(a + b) / (gamma(a) * gamma(b)) * (d / para)^(a - 1) * (sqrt(2) / para - d / para)^(b - 1)
    }
  }

  if (kernel_type == "linear") {
    return(dispersal_linear(para))
  } else if (kernel_type == "exponential") {
    return(dispersal_negexp(para))
  } else if (kernel_type == "gaussian") {
    return(dispersal_gaussian(para))
  } else if (kernel_type == "nonmono") {
    return(dispersal_nonmono(para, ...))
  }
}

#' Generate patch areas/qualities
#' @param patch_num number of patches
#' @param area_type distribution of patch area/quality
#' @param area_para parameter for area/quality distribution
#' @return a vector of patch areas/qualities
#' @export
generate_patch_area <- function(patch_num, area_type, area_para) {
  patch_area <- if (area_type == "uniform") {
    runif(patch_num, 0, area_para)
  } else if (area_type == "normal") {
    rnorm(patch_num, 1, area_para)
  } else if (area_type == "fixed") {
    rep(area_para, patch_num)
  }

  patch_area
}

#' Generate spatial distributions of patches
#' @param patch_type which spatial distribution to use
#' @param patch_num number of patches
#' @param maxRange parameter for spatial interactions
#' @param ncluster number of spatial clusters
#' @return a vector of patch areas/qualities
#' @export
generate_patch_location <- function(patch_type,
                                    patch_num = 100,
                                    maxRange = 1,
                                    ncluster = 5) {
  # Poisson
  if (patch_type == "poisson") {
    patch_location <- matrix(runif(patch_num * 2, min = 0, max = maxRange), ncol = 2)
  }
  # Poisson inhomogeneous
  if (patch_type == "poisson_inhomogeneous") {
    intenfun <- function(x, y) {
      1000 * (x^2 + y)
    }
    patch_location <- rpoispp(intenfun,
      lmax = 1000,
      win = owin(
        c(0, maxRange),
        c(0, maxRange)
      )
    ) %>%
      {
        cbind(.$x, .$y)
      } %>%
      {
        .[sample(1:nrow(.), patch_num), ]
      }
  }
  # Regular
  if (patch_type == "regular") {
    by_range <- floor(sqrt(patch_num))
    patch_location <- expand.grid(
      V1 = seq(0, maxRange, length.out = by_range),
      V2 = seq(0, maxRange, length.out = by_range)
    ) %>%
      as.matrix()
  }
  # homogeneous cluster
  if (patch_type == "cluster_homogeneous") {
    patch_location <- rCauchy(30, 0.01, patch_num / ncluster) %>%
      {
        cbind(.$x, .$y)
      } %>%
      {
        .[sample(1:nrow(.), patch_num), ]
      }
  }
  # inhomogeneous cluster
  if (patch_type == "cluster_inhomogeneous") {
    ff <- function(x, y) {
      exp(2 - 3 * abs(x))
    }
    Z <- as.im(ff, W = owin())
    Y <- rCauchy(100, 0.01, Z)
    patch_location <- Y %>%
      {
        cbind(.$x, .$y)
      } %>%
      {
        .[sample(1:nrow(.), patch_num), ]
      }
  }
  # modular
  if (patch_type == "modular") {
    by_range <- maxRange / ncluster
    patch_location <- 1:ncluster %>%
      map(~ matrix(runif(2 * patch_num / ncluster,
        min = (. - 1) * by_range,
        max = . * by_range
      ), ncol = 2)) %>%
      {
        do.call("rbind", .)
      }
  }
  patch_location
}

# get_subpatch_mc <- function(patch_num,
#                             dispersal_kernel,
#                             area_type,
#                             location_distribution,
#                             self_regulation,
#                             extinction_rate,
#                             omega,
#                             metapopulation,
#                             subpatch_num_max,
#                             Nrand_num,
#                             ...) {
#   if (missing(metapopulation)) {
#     metapopulation <- generate_metapopulation(
#       patch_num,
#       dispersal_kernel,
#       area_type,
#       location_distribution,
#       self_regulation,
#       extinction_rate,
#       omega
#     )
#   }
#
#   patch_area <- metapopulation$patch_area
#   d <- metapopulation$d
#   patch_location <- metapopulation$patch_location
#   capacity_true <- metapopulation$capacity_true
#   mc_true <- metapopulation$mc_true
#
#   dispersal_kernel <- dispersal_negexp(dispersal_kernel)
#
#   mc_random_subnet <- function(subpatch_num) {
#     subpatch_label <- sample(1:patch_num, subpatch_num)
#     d_sub <- d[subpatch_label, subpatch_label]
#     a <- patch_area[subpatch_label]
#
#     calculate_metapopulation_capacity(
#       d_sub,
#       a,
#       dispersal_kernel,
#       self = self_regulation,
#       ex = extinction_rate,
#       omega = omega,
#       patch_mc = T
#     )$patch_mc
#   }
#
#   if (missing(subpatch_num_max)) subpatch_num_max <- 20
#   if (missing(Nrand_num)) Nrand_num <- 20
#   expand_grid(
#     subpatch_num = c(2:subpatch_num_max),
#     Nrand = 1:Nrand_num
#   ) %>%
#     rowwise() %>%
#     mutate(mc = sd(mc_random_subnet(subpatch_num))) %>%
#     ungroup() %>%
#     mutate(mc_true = mc_true)
# }
