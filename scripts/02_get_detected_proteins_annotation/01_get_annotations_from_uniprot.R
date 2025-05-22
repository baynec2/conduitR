################################################################################
# Extract Protein Info From Uniprot For Detected Proteins
################################################################################
# Get input and output files from Snakemake workflow
detected_protein_info <- snakemake@input[["detected_protein_info"]]
uniprot_annotated_protein_info <- snakemake@output[["uniprot_annotated_protein_info"]]

detected_protein_ids = readr::read_delim(detected_protein_info) |>
  dplyr::pull("protein_id")

detected_protein_info = readr::read_delim(detected_protein_info)

# Annotating with the relevant information.
uniprot_annotations <- conduitR::get_annotations_from_uniprot(detected_protein_ids)

# Appending uniprot annotations to protein info
protein_info <- dplyr::inner_join(detected_protein_info, uniprot_annotations,
                                 by = c("protein_id" = "accession"))
# Writing to file
readr::write_delim(protein_info,uniprot_annotated_protein_info)
