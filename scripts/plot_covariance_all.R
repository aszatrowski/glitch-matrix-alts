library(ggplot2)

cov_df <- readr::read_csv(snakemake@input$covariances_csv, show_col_types = FALSE) |>
  dplyr::filter(genotype_covariance < 1.5 & genotype_covariance > -1.5)

he_regression <- lm(phenotype_covariance ~ genotype_covariance, data = cov_df)
he_est <- round(coef(he_regression)[2], 3)
he_label <- bquote(h^2[HE] == .(round(he_est, 3)))
title_text <- bquote(h^2 == .(snakemake@wildcards$h2) ~ b^2 == .(snakemake@wildcards$b2))

p <- ggplot(cov_df, aes(x = genotype_covariance, y = phenotype_covariance)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
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