library(ggplot2)
suppressPackageStartupMessages(library("dplyr"))

binwidth <- snakemake@params$binwidth
min_obs_in_bin <- snakemake@params$min_obs_in_bin 
arch <- snakemake@wildcards$arch

x_axis_cuts_relhat <- seq(-1, 1, binwidth)
cov_df <- arrow::read_parquet(snakemake@input$covariances_csv) |> dplyr::as_tibble() |>
  dplyr::group_by(gc_bins = cut(genotype_covariance, x_axis_cuts_relhat)) |>
  dplyr::mutate(n=n()) |>
  dplyr::filter(n > min_obs_in_bin) |>
  dplyr::summarize(across(-contains("id"), function(x) mean(x, na.rm=TRUE))) |>
  dplyr::filter(genotype_covariance > -0.2 & genotype_covariance < 0.65)

he_regression <- lm(phenotype_covariance ~ genotype_covariance, data = cov_df)
he_est <- round(coef(he_regression)[2], 3)
he_intercept <- coef(he_regression)[1]
he_label <- bquote(h[HE]^2 == .(he_est))

parental_coef <- snakemake@wildcards$parental_coef
title_text <- bquote(h^2 == .(snakemake@wildcards$h2) ~ b^2 == .(snakemake@wildcards$b2) ~ c[m] == .(parental_coef))

p <- ggplot(cov_df, aes(x = genotype_covariance, y = phenotype_covariance)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept = he_intercept, slope = he_est, color = "blue") +
  labs(
    title = title_text,
    subtitle = he_label,
    x = "Genetic Covariance",
    y = "Phenotypic Covariance"
  ) +
  theme_classic()

ggsave(
  filename = snakemake@output$covariance_plot,
  plot = p,
  width = 10,
  height = 6
)