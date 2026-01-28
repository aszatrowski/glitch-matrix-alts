library(ggplot2)

cov_df <- readr::read_csv(snakemake@input$covariances_csv) |>
  dplyr::mutate(
    gcov_bin = cut(
      Genetic_Covariance,
      # breaks = seq(min(Genetic_Covariance), max(Genetic_Covariance), by = 0.05),
      breaks = c(-Inf, seq(-0.25, 0.25, by = 0.05), 0.5, 0.75, Inf),
      include.lowest = TRUE,
      ordered_result = TRUE
    )
  ) |>
  dplyr::group_by(gcov_bin) |>
  dplyr::summarize(
    mean_gcov = mean(Genetic_Covariance),
    mean_y = mean(Y),
    .groups = "drop"
  )

p <- ggplot(cov_df, aes(x = gcov_bin, y = mean_y)) +
  geom_point() +
  geom_line(group = 1) +
  labs(
    title = "Distribution of Y by Genetic Covariance Bins",
    subtitle = paste("Arch:", snakemake@wildcards$arch),
    x = "Genetic Covariance",
    y = "Phenotypic Covariance"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 

ggsave(
  filename = snakemake@output$covariance_plot,
  plot = p,
  width = 10,
  height = 6
)