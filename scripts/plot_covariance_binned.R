library(ggplot2)
library(dplyr)

binwidth <- snakemake@params$binwidth
min_obs_in_bin <- snakemake@params$min_obs_in_bin 
arch <- "M-P equal VCT"

x_axis_cuts_relhat <- seq(-1, 1, binwidth)
cov_df <- readr::read_csv(snakemake@input$covariances_csv, show_col_types = FALSE) |>
  dplyr::group_by(gc_bins = cut(genotype_covariance, x_axis_cuts_relhat)) |>
  dplyr::mutate(n=n()) |>
  dplyr::filter(n > min_obs_in_bin) |>
  dplyr::summarize(across(-contains("id"), function(x) mean(x, na.rm=TRUE)))

he_regression <- lm(phenotype_covariance ~ genotype_covariance, data = cov_df)
he_est <- round(coef(he_regression)[2], 3)
he_intercept <- coef(he_regression)[1]
he_label <- bquote(h[HE]^2 == .(he_est))
title_text <- bquote(paste(.(arch), ", ", h^2 == .(snakemake@wildcards$h2), " ", b^2 == .(snakemake@wildcards$b2)))

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