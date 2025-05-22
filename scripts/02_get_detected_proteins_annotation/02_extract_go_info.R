################################################################################
#Loading Inputs from Snakemake
################################################################################
uniprot_annotated_protein_info = snakemake@input[["uniprot_annotated_protein_info"]]
go_annotations_out = snakemake@output[["go_annotations"]]

# Loading annotated protein info
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info)

# Defining columns in output
shared_columns = c("protein_id","organism_type","superkingdom","kingdom",
                   "phylum","class","order","family","genus","species")

################################################################################
# Saving GO data from the annotations
################################################################################
go_annotations = uniprot_annotated_protein_info |>
  dplyr::select("go",dplyr::all_of(shared_columns)) |>
  # Split by semicolon and unnest into separate rows
  dplyr::mutate(go = stringr::str_split(go, "; ")) |>
  tidyr::unnest(go) |>
  # Extract description and GO term
  dplyr::mutate(
    description = stringr::str_extract(go, "^[^\\[]+"), # Extract the part before the bracket
    go_term = stringr::str_extract(go, "\\[GO:[^\\]]+\\]") # Extract the GO term inside the brackets
  ) |>
  # Clean up the columns (remove extra spaces and brackets)
  dplyr::mutate(
    description = stringr::str_trim(description),
    go_term = stringr::str_replace_all(go_term, "\\[|\\]", "")# Remove square brackets from GO term
  ) |>
  # Select the relevant columns
  dplyr::select(description, go_term,protein_id,dplyr::all_of(shared_columns))

# Write to file.
readr::write_delim(go_annotations,go_annotations_out,delim = "\t")




