library(data.table)
library(ggplot2)

binwidth <- snakemake@params$binwidth
min_obs_in_bin <- snakemake@params$min_obs_in_bin 
h2_val <- snakemake@wildcards$h2
arch <- "M-P equal VCT"

# Read all parquet files and extract b2 from filename, binding into single table
cov_paths <- snakemake@input$covariances_csv_list
cov_list <- lapply(seq_along(cov_paths), function(i) {
  path <- cov_paths[i]
  # Extract b2 from filename: h2_{h2}_b2_{b2}_covmatrix_merged.parquet
  b2_match <- regmatches(path, regexpr("b2_[^_]+", path))
  b2_val <- sub("b2_", "", b2_match)
  
  # Read and add b2 column
  dt <- arrow::read_parquet(path) |> as.data.table()
  dt[, b2 := b2_val]
  return(dt)
})

# Bind all tables
cov_dt <- rbindlist(cov_list)

# Bin by genotype covariance and compute means within each bin, preserving b2
x_axis_cuts_relhat <- seq(-1, 1, binwidth)
binned_dt <- cov_dt[,
  .(
    genotype_covariance = mean(genotype_covariance),
    phenotype_covariance = mean(phenotype_covariance),
    n = .N
  ),
  by = .(
    b2 = b2,
    gc_bins = cut(genotype_covariance, x_axis_cuts_relhat)
  )
][n > min_obs_in_bin & genotype_covariance > -0.2 & genotype_covariance < 0.65]

# Convert b2 to factor and order for consistent coloring
binned_dt[, b2 := factor(b2, levels = sort(unique(as.numeric(b2))))]

# Create plot with color by b2, point and line for each b2
p <- ggplot(binned_dt, aes(x = genotype_covariance, y = phenotype_covariance, color = b2)) +
  geom_point(size = 2) +
  geom_line() +
  scale_color_viridis_d() +
  labs(
    title = bquote("Phenotype covariance by" ~ b^2),
    x = "Genetic Covariance",
    y = "Phenotypic Covariance",
    color = bquote(b^2)
  ) +
  theme_classic() +
  theme(legend.position = "right")

ggsave(
  filename = snakemake@output$covariance_plot,
  plot = p,
  width = 12,
  height = 6
)
