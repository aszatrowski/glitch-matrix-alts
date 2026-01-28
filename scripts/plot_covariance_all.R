library(ggplot2)

cov_df <- readr::read_csv(snakemake@input$covariances_csv, show_col_types = FALSE) |>
  dplyr::filter(Genetic_Covariance <= 1.125)
he_regression <- lm(Y ~ Genetic_Covariance, data = cov_df)

p <- ggplot(cov_df, aes(x = Genetic_Covariance, y = Y)) +
  geom_point() +
  geom_abline(
    slope = coef(he_regression)[2],
    intercept = coef(he_regression)[1],
    color = "blue"
  ) +
  geom_label(
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.5,
    label = paste0("h^2 = ", round(coef(he_regression)[2], 3))
  ) +
  labs(
    title = "Phenotypic Covariance vs Genetic Covariance",
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