library(ggplot2)
suppressPackageStartupMessages(library("dplyr"))
library(data.table)

binwidth <- snakemake@params$binwidth
min_obs_in_bin <- snakemake@params$min_obs_in_bin 
arch <- "M-P equal VCT"

x_axis_cuts_relhat <- seq(-1, 1, binwidth)
import_and_bin <- function(parquet, b2_value) {
    cov_df <- arrow::read_parquet(snakemake@input$covariances_csv) |>
        dplyr::as_tibble() |>
        dplyr::group_by(gc_bins = cut(genotype_covariance, x_axis_cuts_relhat)) |>
        dplyr::mutate(n=n()) |>
        dplyr::filter(n > min_obs_in_bin) |>
        dplyr::summarize(across(-contains("id"), function(x) mean(x, na.rm=TRUE))) |>
        dplyr::mutate(b2 = b2_value)
}

cov_df_all_b2 <- 
