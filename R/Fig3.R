#' Reproduce Figure 3 in the main text
#'
#' @export
plot_figure_3 <- function() {
  ggthemr("fresh")
  p_capacity <-
    capacity_predicted %>%
    group_by(prediction_method) %>%
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
    ggplot(aes(error, group = prediction_method)) +
    geom_vline(xintercept = 0, color = "black", linewidth = .3) +
    geom_histogram(aes(y = ..density.., fill = prediction_method), bins = 30) +
    facet_grid(~prediction_method) +
    scale_fill_manual(
      values = c("#EFBB24", "#5DAC81")
    ) +
    labs(
      x = "Prediction error of metapopoulation capacity",
      y = "Density"
    ) +
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 12)
    )
  p_mc <- mc_predicted %>%
    filter(cor > quantile(cor, .025)) %>%
    ggplot(aes(cor)) +
    geom_histogram(aes(y = ..density.., ), bins = 30, fill = "#5DAC81") +
    labs(
      x = "Correlation between sampled and empirical patch importance",
      y = "Density"
    ) +
    theme(
      legend.position = "none"
    )

  p_capacity / p_mc
}
