library(dplyr)

# returns list of data frames, one per replicate
replicates <- purrr::map(snakemake@input[["replicates"]], readr::read_csv, show_col_types = FALSE)
for (i in seq_along(replicates)) {
  # append replicate identifier to all columns except pedigree_id1 and pedigree_id2
  replicates[[i]] <- replicates[[i]] |>
    dplyr::mutate(rep = i - 1)
}
merged_replicates <- dplyr::bind_rows(replicates) |>
  dplyr::select(pedigree_id1, pedigree_id2, Y, Genetic_Covariance) |>
  dplyr::group_by(pedigree_id1, pedigree_id2) |>
  dplyr::summarize(
    dplyr::across(c(Y, Genetic_Covariance), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )
  
readr::write_csv(merged_replicates, snakemake@output[["merged_replicates"]])