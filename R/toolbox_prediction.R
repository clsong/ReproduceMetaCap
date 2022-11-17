#' Calculate the mean of off-diagonal elements in a matrix
#' @param m a matrix
#' @return mean of off-diagonal elements
#' @export
estimate_offdigonal_mean <- function(m) {
  (sum(m) - sum(diag(m))) / (nrow(m)^2 - nrow(m))
}

#' Calculate the variance of off-diagonal elements in a matrix
#' @param m a matrix
#' @return variance of off-diagonal elements
#' @export
estimate_offdigonal_variance <- function(m) {
  diag(m) <- NA
  var(c(m), na.rm = T)
}

#' Estimate the ratio of dependency
#' @param A0 area matrix (A_i^e A_j^omega)
#' @param self_regulation whether self-colonization exists
#' @return ratio of dependency
#' @export
estimate_ratio_area <- function(A0, self_regulation) {
  A0_mean <- estimate_offdigonal_mean(A0)
  A0_var <- estimate_offdigonal_variance(A0)
  diag_elements <- mean(diag(A0))

  calculate_metapopulation_capacity(A0) / estimate_asymptotic_capacity(A0, nrow(A0))
}

#' Analytic-based inference of metapopulation capacity
#' @param m landscape/connectivity matrix
#' @param total_num total number of patches in a metapopulation
#' @return predicted metapopulation capacity
#' @export
estimate_asymptotic_capacity <- function(m, total_num) {
  m_mean <- estimate_offdigonal_mean(m)
  m_var <- estimate_offdigonal_variance(m)
  diag_element <- mean(diag(m))

  total_num * m_mean + diag_element + m_var / m_mean
}

#' Regression-based inference of metapopulation capacity
#' @param df a tibble of sub-metapopulations and their capacities
#' @param method regression method
#' @return predicted metapopulation capacity
#' @export
predict_capacity <- function(df, method) {
  if (method == "capacity_linear") {
    capacity <- linear_reg() %>%
      fit(capacity_empirical ~ subpatch_num, data = df) %>%
      predict(
        new_data = data.frame(subpatch_num = patch_num)
      ) %>%
      pull(.pred)
  } else if (method == "capacity_nonlinear") {
    capacity <- gen_additive_mod() %>%
      set_engine("mgcv") %>%
      set_mode("regression") %>%
      fit(capacity_empirical ~ s(subpatch_num), data = df) %>%
      predict(
        new_data = data.frame(subpatch_num = patch_num)
      ) %>%
      pull(.pred)
  } else if (method == "ratio_based") {
    capacity <- df %>%
      filter(subpatch_num > .8 * max(subpatch_num)) %>%
      pull(capacity_asymptotic) %>%
      mean()
  }
  capacity
}
