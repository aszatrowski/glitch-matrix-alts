library(ggplot2)
covariances <- readr::read_csv("data/gcta_covariances.csv") |>
  dplyr::filter(Genetic_Covariance >= 0.2)

covariances_grouped <- covariances |>
  dplyr::mutate(
    geno_bin = cut(
      Genetic_Covariance,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE,
      labels = paste0("Geno ", seq(0, 0.9, by = 0.1), "-", seq(0.1, 1, by = 0.1))
    )
  ) |>
  dplyr::group_by(geno_bin) |>
  dplyr::summarise(
    mean_geno = mean(Genetic_Covariance),
    mean_pheno = mean(Y),
    .groups = "drop"
  )
he_regression <- lm(Y ~ Genetic_Covariance, data = covariances)
ggplot(covariances, aes(x = Genetic_Covariance, y = Y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    title = "Mean Phenotypic Covariance vs Mean Genetic Covariance",
    subtitle = bquote(hat(h^2) == .(round(coef(he_regression)[2], 3))),
    x = bquote("Genetic Covariance," ~ pi),
    y = "Phenotypic Covariance"
  ) +
  theme_minimal()
ggsave("results/mean_covariances.png", width = 8, height = 6)