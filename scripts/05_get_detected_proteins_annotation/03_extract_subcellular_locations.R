################################################################################
#Loading Inputs from Snakemake
################################################################################
uniprot_annotated_protein_info = snakemake@input[["uniprot_annotated_protein_info"]]
subcellular_locations_out = snakemake@output[["subcellular_locations"]]

# #test
# uniprot_annotated_protein_info="user_input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt"
# subcellular_locations_out="user_input/00_database_resources/detected_protein_resources/04_subcellular_locations.txt"
################################################################################
# Getting Cellular Location Information
################################################################################
# Defining columns in output
shared_columns = c("protein_id","organism_type","superkingdom","kingdom",
                   "phylum","class","order","family","genus","species")

uniprot_annotated_protein_info =readr::read_delim(uniprot_annotated_protein_info)


subcellular_locations = uniprot_annotated_protein_info |>
  dplyr::select("cc_subcellular_location",dplyr::all_of(shared_columns)) |>
  dplyr::mutate(cc_subcellular_location = stringr::str_split(cc_subcellular_location, "; ")) |>  # Split into list
  tidyr::unnest(cc_subcellular_location) |>  # Expand to long format
  dplyr::mutate(
    subcellular_location = stringr::str_extract(cc_subcellular_location, "^[^{]+"),  # Extract location (before {)
    evidence = stringr::str_extract(cc_subcellular_location, "\\{([^}]+)\\}")  # Extract text inside curly braces
  ) |>
  dplyr::mutate(
    evidence = stringr::str_remove_all(evidence, "[{}]"),  # Clean up curly braces
    cc_subcellular_location = gsub("SUBCELLULAR LOCATION: ", "",
                                   cc_subcellular_location),
    cc_subcellular_location = gsub("\\{.*","",cc_subcellular_location),
    cc_subcellular_location = gsub("\\.","",cc_subcellular_location),
    note = stringr::str_extract(cc_subcellular_location, "Note.*"),
    cc_subcellular_location = gsub("Note.*","",cc_subcellular_location)) |>
  dplyr::select(cc_subcellular_location, evidence,note, protein_id,
                dplyr::all_of(shared_columns))  # Keep relevant columns

# Writing to file.
readr::write_delim(subcellular_locations,
                   subcellular_locations_out,
                   delim = "\t")
