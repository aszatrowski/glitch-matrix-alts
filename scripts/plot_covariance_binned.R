library(ggplot2)
library(dplyr)
x_axis_cuts_relhat <- seq(-1, 1, 0.01)
cov_df <- readr::read_csv(snakemake@input$covariances_csv, show_col_types = FALSE) |>
  dplyr::group_by(GCC = cut(Genetic_Covariance, x_axis_cuts_relhat)) |>
  dplyr::mutate(n=n()) |>
  dplyr::summarize(across(-contains("id"), function(x) mean(x, na.rm=TRUE)))

title_text <- bquote(h^2 == .(snakemake@wildcards$h2) ~ b^2 == .(snakemake@wildcards$b2))
p <- ggplot(cov_df, aes(x = Genetic_Covariance, y = Y)) +
  geom_point() +
  geom_line() +
  labs(
    title = title_text,
    subtitle = paste("Arch:", snakemake@wildcards$arch),
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