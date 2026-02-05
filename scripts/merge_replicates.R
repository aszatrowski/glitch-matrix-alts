library(dplyr)

# returns list of data frames, one per replicate
replicates <- purrr::map(snakemake@input[["replicates"]], readr::read_csv, show_col_types = FALSE)
for (i in seq_along(replicates)) {
  # append replicate identifier to all columns except pedigree_id1 and pedigree_id2
  # Posted by akrun, modified by community. See post 'Timeline' for change history
  # Retrieved 2026-02-05, License - CC BY-SA 4.0
  # Source - https://stackoverflow.com/a/64188680
  replicates[[i]] <- dplyr::rename_with(
    replicates[[i]],
    ~ paste(., "rep", i - 1, sep = "_"),  # 0-indexed to match your rep wildcard
    !c("pedigree_id1", "pedigree_id2"), # columns to exclude from renaming
  )
}

merged_replicates <- purrr::reduce(replicates, dplyr::full_join, by = c("pedigree_id1", "pedigree_id2"))
print(snakemake@output[["merged_replicates"]])
readr::write_csv(merged_replicates, snakemake@output[["merged_replicates"]])