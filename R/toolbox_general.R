#' Calculate the landscape/connectivity matrix from given metapopulation parameters
#' @param distance a matrix containing pairwise distances between patches
#' @param patch_area a vector containing the patch areas/qualities
#' @param dispersal_kernel which dispersal kernel to use
#' @param dispersal_rate area dependence of dispersal
#' @param extinction_rate area dependence of extinction
#' @param self_regulation whether self-colonization exists
#' @return landscape/connectivity matrix
#' @export
calculate_landscape_matrix <- function(distance,
                                       patch_area,
                                       dispersal_kernel,
                                       dispersal_rate = 1,
                                       extinction_rate = 0.5,
                                       self_regulation = TRUE) {
  m <- as.matrix(dispersal_kernel(distance))
  diag(m) <- ifelse(self_regulation, 1, 0)
  m <- outer(patch_area^extinction_rate, patch_area^dispersal_rate) * m

  m
}

#' Calculate metapopulation capacity and patch importance
#' @param m landscape/connectivity matrix
#' @param patch_mc whether to compute patch importance
#' @return indicators of metapopulation persistence
#' @export
calculate_metapopulation_capacity <- function(m,
                                              patch_mc = FALSE) {
  e <- eigen(m, only.values = !patch_mc)
  if (patch_mc) {
    cap <- list(capacity = Re(e$values[1]), patch_mc = Re(e$vectors[, 1])^2)
  } else {
    cap <- Re(e$values[1])
  }
  cap
}
