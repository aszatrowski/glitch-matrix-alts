library(ggplot2)

cov_df <- readr::read_csv(snakemake@input$covariances_csv, show_col_types = FALSE) |>
  dplyr::filter(Genetic_Covariance < 1.125 & Genetic_Covariance > -1.5)

he_regression <- lm(Y ~ Genetic_Covariance, data = cov_df)
he_est <- round(coef(he_regression)[2], 3)
he_label <- bquote(h^2[HE] == .(round(he_est, 3)))
p <- ggplot(cov_df, aes(x = Genetic_Covariance, y = Y)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Phenotypic Covariance vs Genetic Covariance",
    subtitle = bquote(h[HE]^2 == .(round(he_est, 3))),
    x = "Genetic Covariance",
    y = "Phenotypic Covariance"
  )

ggsave(
  filename = snakemake@output$covariance_plot,
  plot = p,
  width = 10,
  height = 6
)