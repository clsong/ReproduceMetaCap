#' Reproduce Figure 2 in the main text
#'
#' @export
plot_figure_2 <- function() {
  # generate a random metapopulation ----------------------------------------
  patch_num <- 100

  area_type <- "uniform"
  area_para <- 1

  location_distribution <- "cluster_inhomogeneous"

  kernel_type <- "exponential"
  charactersitic_distance <- 1

  self_regulation <- F
  extinction_rate <- 1
  dispersal_rate <- 1

  set.seed(1010)
  metapopulation <- generate_metapopulation(
    patch_num,
    kernel_type,
    charactersitic_distance,
    location_distribution,
    self_regulation,
    dispersal_rate,
    extinction_rate,
    area_type,
    area_para,
  )

  M <- metapopulation$landscape_matrix
  area <- metapopulation$patch_area

  subpatch_num_max <- 50
  set.seed(1234)
  sampled_subpatches <- sample(1:patch_num, subpatch_num_max)

  # plot location -----------------------------------------------------------
  ggthemr("fresh", layout = "minimal")
  p_location <-
    tibble(
      x = metapopulation$patch_location[, 1],
      y = metapopulation$patch_location[, 2],
      area = metapopulation$patch_area,
      sampled = ifelse(1:patch_num %in% sampled_subpatches, "sampled", "unsampled")
    ) %>%
    ggplot(aes(
      x = x, y = y,
      group = sampled,
      color = sampled,
      fill = sampled
    )) +
    geom_point(aes(size = area), shape = 21) +
    scale_size(range = c(0, 3)) +
    labs(
      x = "",
      y = ""
    ) +
    scale_fill_manual(
      values = c("#64ADC2", "white")
    ) +
    scale_color_manual(
      values = c("#64ADC2", "black")
    ) +
    guides(size = "none") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      aspect.ratio = 1,
      panel.border = element_rect(fill = NA),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      line = element_blank(),
      legend.margin = margin(-40, 0, 0, 0)
    )

  # plot capa ---------------------------------------------------------------

  ggthemr("fresh")

  Nrand_num <- 20
  set.seed(1010)
  df_capa <- expand_grid(
    subpatch_num = c(2:subpatch_num_max),
    Nrand = 1:Nrand_num
  ) %>%
    mutate(sampled_patch = map(subpatch_num, ~ sample(sampled_subpatches, .))) %>%
    rowwise() %>%
    mutate(capacity = calculate_metapopulation_capacity(M[sampled_patch, sampled_patch])) %>%
    mutate(area_var = mean(area[subpatch_num]^2)) %>%
    select(subpatch_num, capacity, area_var)

  naive_estimate <- lm(capacity ~ subpatch_num, df_capa) %>%
    tidy()
  p_capa <-
    df_capa %>%
    ggplot(aes(subpatch_num, capacity)) +
    geom_point() +
    labs(
      x = "Number of sampled subpatches",
      y = "Metapopulation capacity \u03BB"
    ) +
    geom_point(
      data = tibble(
        subpatch_num = patch_num,
        capacity = calculate_metapopulation_capacity(M)
      ),
      color = "#F75C2F",
      size = 3
    ) +
    geom_segment(
      aes(
        x = 2,
        y = 2 * naive_estimate$estimate[2] + naive_estimate$estimate[1],
        xend = subpatch_num_max,
        yend = subpatch_num_max * naive_estimate$estimate[2] + naive_estimate$estimate[1]
      ),
      color = "gray"
    ) +
    geom_segment(
      aes(
        x = subpatch_num_max,
        y = subpatch_num_max * naive_estimate$estimate[2] + naive_estimate$estimate[1],
        xend = patch_num,
        yend = patch_num * naive_estimate$estimate[2] + naive_estimate$estimate[1]
      ),
      color = "gray",
      linetype = "dashed"
    )

  # plot mc ---------------------------------------------------------------
  set.seed(1010)
  df_mc <- expand_grid(
    subpatch_num = c(2:subpatch_num_max),
    Nrand = 1:Nrand_num
  ) %>%
    mutate(sampled_patch = map(subpatch_num, ~ sample(sampled_subpatches, .))) %>%
    mutate(mc = map(sampled_patch, function(sampled_patch) {
      M[sampled_patch, sampled_patch] %>%
        calculate_metapopulation_capacity(patch_mc = T) %>%
        {
          .$patch_mc
        } %>%
        enframe(
          name = "patch_label",
          value = "mc"
        ) %>%
        mutate(patch_label = sampled_patch)
    })) %>%
    select(-sampled_patch) %>%
    unnest(mc)

  mc_true_normalized <-
    calculate_metapopulation_capacity(M, patch_mc = T)$patch_mc %>%
    enframe(
      name = "patch_label",
      value = "true_mc"
    ) %>%
    mutate(true_mc = true_mc * patch_num) %>%
    filter(patch_label %in% sampled_subpatches) %>%
    left_join(
      df_mc %>%
        filter(subpatch_num == subpatch_num_max) %>%
        group_by(patch_label) %>%
        summarise(mc = mean(mc)) %>%
        ungroup() %>%
        select(patch_label, mc) %>%
        mutate(mc = mc * subpatch_num_max) %>%
        rename(pred_mc = mc)
    )

  p_mc <-
    df_mc %>%
    group_by(subpatch_num, patch_label) %>%
    summarise(
      mc = mean(mc),
      .groups = "drop"
    ) %>%
    ggplot(aes(subpatch_num, mc * subpatch_num)) +
    geom_line(aes(group = patch_label)) +
    labs(
      x = "Number of sampled subpatches",
      y = "Normalized patch importance \u03C9"
    ) +
    geom_segment(
      aes(
        x = subpatch_num_max,
        y = pred_mc,
        xend = patch_num,
        yend = true_mc
      ),
      color = "gray",
      linetype = "dashed",
      data = mc_true_normalized
    ) +
    geom_point(
      aes(x = patch_num, y = true_mc),
      data = mc_true_normalized,
      color = "#F75C2F"
    )

  # export plot -------------------------------------------------------------
  p_location + (p_capa + p_mc) + theme(element_text(size = 20)) +
    plot_annotation(tag_levels = "A")
}
