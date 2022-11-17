#' Calculate the metapopulation capacity of sub-metapopulations
#' @param M landscape/connectivity matrix
#' @param subpatch_num_min minimum number of samples
#' @param subpatch_num_max maximum number of samples
#' @param Nrand_num number of replicates with a given number of sampled patches
#' @return metapopulation capacity for sub-metapopulations
#' @export
sample_sub_metapopulation <- function(M,
                                      subpatch_num_min = 2,
                                      subpatch_num_max = patch_num / 2,
                                      Nrand_num = 20) {
  samplable_patches <- sample(1:nrow(M), nrow(M))

  expand_grid(
    subpatch_num = c(subpatch_num_min:subpatch_num_max),
    Nrand = 1:Nrand_num
  ) %>%
    mutate(sampled_patch = map(
      subpatch_num,
      ~ sample(samplable_patches, .)
    )) %>%
    mutate(m = map(sampled_patch, ~ M[., .])) %>%
    mutate(capacity_empirical = map_dbl(
      m,
      calculate_metapopulation_capacity
    )) %>%
    mutate(capacity_asymptotic = map_dbl(
      m,
      ~ estimate_asymptotic_capacity(., patch_num) * estimate_ratio_area(.)
    )) %>%
    select(-sampled_patch, -m)
}

#' Calculate the patch importance of sub-metapopulations
#' @param M landscape/connectivity matrix
#' @param subpatch_num_min minimum number of samples
#' @param subpatch_num_max maximum number of samples
#' @param Nrand_num number of replicates with a given number of sampled patches
#' @return patch importance for sub-metapopulations
#' @export
sample_sub_metapopulation_mc <- function(M,
                                         subpatch_num_min = 2,
                                         subpatch_num_max = patch_num / 2,
                                         Nrand_num = 20) {
  samplable_patches <- sample(1:patch_num, patch_num)

  expand_grid(
    subpatch_num = c(10:subpatch_num_max),
    Nrand = 1:Nrand_num
  ) %>%
    mutate(sampled_patch = map(subpatch_num, ~ sample(samplable_patches, .))) %>%
    mutate(m = map(sampled_patch, ~ M[., .])) %>%
    mutate(mc_empirical = map2(.data$m, .data$sampled_patch, function(m, sampled_patch) {
      calculate_metapopulation_capacity(m, patch_mc = T)$patch_mc %>%
        enframe(
          name = "patch_label",
          value = "mc_estimated"
        ) %>%
        mutate(patch_label = sampled_patch)
    })) %>%
    select(-m, -sampled_patch) %>%
    unnest(mc_empirical)
}
