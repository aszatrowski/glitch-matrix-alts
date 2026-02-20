library(data.table)
library(ggplot2)

binwidth <- snakemake@params$binwidth
min_obs_in_bin <- snakemake@params$min_obs_in_bin 
h2_val <- snakemake@wildcards$h2
arch <- snakemake@wildcards$arch 
parental_coef <- snakemake@wildcards$parental_coef

# Read all parquet files and extract b2 from filename, binding into single table
cov_paths <- snakemake@input$covariances_csv_list
cov_list <- lapply(seq_along(cov_paths), function(i) {
  path <- cov_paths[i]
  # Extract b2 from filename: h2_{h2}_b2_{b2}_covmatrix_merged.parquet
  # Use the basename and a capture to reliably extract numeric b2 (handles dots, signs, exponents)
  fname <- basename(path)
  m <- regexec("_b2_([-0-9.eE]+)(?:_|\\.|$)", fname)
  reg <- regmatches(fname, m)
  if (length(reg) && length(reg[[1]]) >= 2) {
    b2_val <- reg[[1]][2]
  } else {
    # Fallback: simple match (older pattern)
    b2_match <- regmatches(fname, regexpr("b2_[^_.]+", fname))
    b2_val <- sub("b2_", "", b2_match)
  }

  # Coerce to numeric and warn if parsing failed (NA)
  b2_num <- as.numeric(b2_val)
  if (is.na(b2_num)) {
    warning(sprintf("Could not parse b2 from filename '%s' (extracted: '%s')", fname, b2_val))
  }

  # Read parquet and add numeric b2 column
  dt <- arrow::read_parquet(path) |> as.data.table()
  dt[, b2 := b2_num]
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
  scale_x_continuous(limits = c(-0.05, 0.65)) +
  labs(
    title = bquote("Phenotype covariance by" ~ b^2),
    subtitle = bquote(paste("arch:" ~ .(arch), ", ", h^2 == .(h2_val) ~ c[m] == .(parental_coef))),
    caption = paste(
      "Replicates:", snakemake@params$replicates,
      "| Variants:", snakemake@params$m_variants,
      "| Causal:", snakemake@params$n_causal
    ),
    x = "Genetic Covariance",
    y = "Phenotypic Covariance",
    color = bquote(b^2)
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(
  filename = snakemake@output$covariance_plot,
  plot = p,
  width = 12,
  height = 6
)
