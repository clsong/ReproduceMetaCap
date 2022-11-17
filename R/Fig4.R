#' Reproduce Figure 4 in the main text
#'
#' @export
plot_figure_4 <- function() {
  ggthemr("fresh")

  p_capacity <-
    empirical_capacity_predicted %>%
    select(-capacity_true, -sub_metapopulation, -capacity_predicted) %>%
    mutate(sample_prob = as.factor(sample_prob)) %>%
    group_by(prediction_method, sample_prob) %>%
    filter(error < quantile(error, .975) &
      error > quantile(error, .025)) %>%
    ungroup() %>%
    filter(prediction_method != "capacity_nonlinear") %>%
    mutate(
      prediction_method =
        ifelse(prediction_method == "capacity_linear",
          "Regression based",
          "Analytic based"
        )
    ) %>%
    filter(prediction_method == "Analytic based") %>%
    ggplot(aes(x = error, y = sample_prob)) +
    geom_density_ridges2(aes(fill = sample_prob),
      stat = "binline", bins = 40,
      scale = .9,
      color = "white"
    ) +
    labs(
      x = "Prediction error of\nmetapopoulation capacity",
      y = "Percentage of\nsampled patches"
    ) +
    scale_fill_manual(
      values = paletteer_dynamic("cartography::grey.pal", 5)
    ) +
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 12)
    ) +
    geom_hline(yintercept = "0.1", color = "black", size = .1) +
    geom_hline(yintercept = "0.2", color = "black", size = .1) +
    geom_hline(yintercept = "0.3", color = "black", size = .1) +
    geom_hline(yintercept = "0.4", color = "black", size = .1) +
    geom_hline(yintercept = "0.5", color = "black", size = .1)


  p_mc <-
    empirical_mc_predicted %>%
    mutate(sample_prob = as.factor(sample_prob)) %>%
    group_by(sample_prob) %>%
    filter(cor > quantile(cor, .05)) %>%
    ungroup() %>%
    ggplot(aes(x = cor, y = sample_prob)) +
    geom_density_ridges2(aes(fill = sample_prob),
      stat = "binline", bins = 40,
      scale = .9,
      color = "white"
    ) +
    labs(
      x = "Correlation between sampled and \nempirical patch importance",
      y = "Percentage of\nsampled patches"
    ) +
    # xlim(c(NA, 1))+
    scale_fill_manual(
      values = paletteer_dynamic("cartography::grey.pal", 5)
    ) +
    geom_hline(yintercept = "0.1", color = "black", linewidth = .1) +
    geom_hline(yintercept = "0.2", color = "black", linewidth = .1) +
    geom_hline(yintercept = "0.3", color = "black", linewidth = .1) +
    geom_hline(yintercept = "0.4", color = "black", linewidth = .1) +
    geom_hline(yintercept = "0.5", color = "black", linewidth = .1)

  p_capacity + p_mc +
    plot_annotation(tag_levels = list(LETTERS[2:3])) &
    theme_ridges(
      grid = F,
      center_axis_labels = TRUE,
      font_size = 12
    ) +
      theme(
        panel.grid.major.y = element_line(
          colour = "black",
          size = .1
        ),
        legend.position = "none"
      )
}
