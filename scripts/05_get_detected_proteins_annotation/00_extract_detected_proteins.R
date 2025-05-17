# Get input and output files from Snakemake workflow
protein_info_df <- snakemake@input[["protein_info_df"]]
protein_info_fasta <- snakemake@input[["protein_info_fasta"]]
report_pg_matrix <- snakemake@input[["report_pg_matrix"]]
detected_protein_info_df <- snakemake@output[["detected_protein_info_df"]]
detected_protein_info_fasta<- snakemake@output[["detected_protein_info_fasta"]]

# Make the directory if it doesn't exist
if(!dir.exists(dirname(detected_protein_info_df))){
dir.create(dirname(detected_protein_info_df))
}

# Reading in protein information
protein_info = readr::read_tsv(protein_info_df)

# Determining what the detected proteins from DIA-NN were detected
report_pg_matrix = readr::read_tsv(report_pg_matrix)
detected_uniprot_ids = report_pg_matrix |>
  dplyr::pull("Protein.Group") |>
  strsplit(";") |>
  unlist() |>
  unique()

# Filtering detected protein information
detected_protein_info = protein_info |>
  dplyr::filter(protein_id %in% detected_uniprot_ids)

# Writing to file - data frame
readr::write_tsv(detected_protein_info,detected_protein_info_df)

### Filtering Fasta file to only include detected proteins ###
# Reading in fasta file
sequences <- Biostrings::readAAStringSet(protein_info_fasta)
# Extract UniProt IDs from FASTA headers
uniprot_ids <- names(sequences) |> stringr::str_extract("(?<=\\|)[A-Z0-9]+(?=\\|)")
# Create a named vector mapping UniProt IDs to full headers
header_map <- setNames(names(sequences), uniprot_ids)

# Filtering the header map to include only the ids detected in study.
detected_headers = header_map[names(header_map) %in% detected_protein_info$protein_id]

# Producing a filtered_fasta
filtered_fasta <- sequences[names(sequences) %in% detected_headers]

# Writing Filtered Fasta file
Biostrings::writeXStringSet(filtered_fasta,detected_protein_info_fasta)
