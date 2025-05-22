# Loading files
uniprot_annotated_protein_info = snakemake@input[["uniprot_annotated_protein_info"]]
detected_taxonomy =snakemake@output[["detected_taxonomy"]]

# Testing
#uniprot_annotated_protein_info="user_input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt"
#detected_taxonomy="user_input/00_database_resources/detected_protein_resources/06_detected_taxonomy.txt"

# Select only relevant columns
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info) |>
  dplyr::select(protein_id,
                organism_type,
                superkingdom,
                kingdom,
                phylum,
                class,
                order,
                family,
                genus,
                species)

readr::write_delim(uniprot_annotated_protein_info,detected_taxonomy)
