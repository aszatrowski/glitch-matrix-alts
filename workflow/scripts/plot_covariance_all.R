library(ggplot2)
library(dplyr)
cov_df <- arrow::read_parquet(snakemake@input$covariances_csv) |>
  dplyr::as_tibble() |>
  dplyr::filter(genotype_covariance > -0.2 & genotype_covariance < 0.65)

arch <- snakemake@wildcards$arch
he_regression <- lm(phenotype_covariance ~ genotype_covariance, data = cov_df)
he_est <- round(coef(he_regression)[2], 3)
parental_coef <- snakemake@wildcards$parental_coef
title_text <- bquote(h^2 == .(snakemake@wildcards$h2) ~ b^2 == .(snakemake@wildcards$b2) ~ c[m]^2 == .(parental_coef))

p <- ggplot(cov_df, aes(x = genotype_covariance, y = phenotype_covariance)) +
  stat_bin2d(bins = 40, aes(fill = after_stat(log10(count)))) +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  # geom_abline(slope = he_est, intercept = coef(he_regression[1]), color = "cyan", linetype = "solid", linewidth = 0.8) +
  labs(
    title = title_text,
    subtitle = bquote(
      .(paste("Arch:", snakemake@wildcards$arch)) ~ h[HE]^2 == .(round(he_est, 3))
    ),
    x = "Genetic Covariance",
    y = "Phenotypic Covariance",
    fill = expression(log[10] * "(count)")
  ) +
  theme_bw()

ggsave(
  filename = snakemake@output$covariance_plot,
  plot = p,
  width = 10,
  height = 6
)