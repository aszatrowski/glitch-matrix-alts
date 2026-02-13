library(ggplot2)

cov_df <- readr::read_csv(snakemake@input$covariances_csv, show_col_types = FALSE) |>
  dplyr::slice_sample(n = 50000)

he_regression <- lm(phenotype_covariance ~ genotype_covariance, data = cov_df)
he_est <- round(coef(he_regression)[2], 3)
he_label <- bquote(h^2[HE] == .(round(he_est, 3)))
title_text <- bquote(h^2 == .(snakemake@wildcards$h2) ~ b^2 == .(snakemake@wildcards$b2))

p <- ggplot(cov_df, aes(x = genotype_covariance, y = phenotype_covariance)) +
  geom_density_2d_filled(alpha = 0.5) +
  geom_abline(slope = he_est, intercept = 0, color = "blue", linetype = "solid") +
  labs(
    title = title_text,
    subtitle = bquote(
      .(paste("Arch:", snakemake@wildcards$arch)) ~ h[HE]^2 == .(round(he_est, 3))
    ),
    x = "Genetic Covariance",
    y = "Phenotypic Covariance"
  )

ggsave(
  filename = snakemake@output$covariance_plot,
  plot = p,
  width = 10,
  height = 6
)