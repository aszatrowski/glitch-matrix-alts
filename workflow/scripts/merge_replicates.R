library(data.table)

import_melt_cov <- function(cov_path, iids_path, replicate_num, value_name) {
  # Read in covariance matrix and IID files, convert to long format with replicate number
  cov_dt <- fread(cov_path, nThread = snakemake@threads)
  iids_dt <- fread(iids_path, nThread = 1) |>
    dplyr::mutate(
      xftsim_id = paste0(
        # generations are 0-indexed
        snakemake@params$generations - 1,
        "..",
        IID,
        ".",
        `#FID`
      )
    )
  
  setnames(cov_dt, old = names(cov_dt), new = iids_dt$xftsim_id)
  
  cov_dt[, xftsim_id1 := iids_dt$xftsim_id]
  
  cov_long <- melt(cov_dt,
                   id.vars = "xftsim_id1",
                   variable.name = "xftsim_id2",
                   value.name = value_name)
  
  # Canonicalize ID ordering for consistent joins
  cov_long[, c("xftsim_id1", "xftsim_id2") := .(
    pmin(as.character(xftsim_id1), as.character(xftsim_id2)),
    pmax(as.character(xftsim_id1), as.character(xftsim_id2))
  )]
  
  # remove upper-triangle 0 entries that were blank in plink matrix
  cov_long[get(value_name) != 0, replicate := replicate_num]
  return(cov_long)
}

import_melt_multiple <- function(paths, iids_paths, value_name) {
  rbindlist(lapply(seq_along(paths), function(i) {
    import_melt_cov(paths[i], iids_paths[i], replicate_num = i, value_name)
  }))
}

grm_file_paths <- snakemake@input[['grm_replicates']]
iid_paths <- snakemake@input[['grm_replicate_iids']]
pcov_paths <- snakemake@input[['phenotype_covariance_replicates']]

grm_all <- import_melt_multiple(grm_file_paths, iid_paths, "genotype_covariance")
pcov_all <- rbindlist(lapply(seq_along(pcov_paths), function(i) {
  pcov_dt <- arrow::read_parquet(pcov_paths[i]) |> as.data.table()
  pcov_dt[, replicate := i]
}))
grm_pcov_averaged <- grm_all[pcov_all, on = .(xftsim_id1, xftsim_id2, replicate)][,
  .(genotype_covariance = mean(genotype_covariance),
    phenotype_covariance = mean(phenotype_covariance)),
  by = .(xftsim_id1, xftsim_id2)
]
# remove self covariances
grm_pcov_averaged <- grm_pcov_averaged[xftsim_id1 != xftsim_id2, ]
# grm_pcov_averaged <- grm_pcov_averaged[genotype_covariance < 0.75, ]

arrow::write_parquet(grm_pcov_averaged, snakemake@output[["merged_replicates"]])
