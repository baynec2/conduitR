# Loading files
uniprot_annotated_protein_info = snakemake@input[["uniprot_annotated_protein_info"]]
kegg_annotations_out =snakemake@output[[1]]


# Testing
#uniprot_annotated_protein_info="user_input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt"
#kegg_annotations="user_input/00_database_resources/detected_protein_resources/05_kegg_annotations.txt"

# Reading in detected protein info
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info)

# Defining columns that all outputs will share.
shared_columns = c("protein_id","organism_type","superkingdom","kingdom",
                   "phylum","class","order","family","genus","species")

################################################################################
# Getting KEGG Information for each KEGG ID
################################################################################
# Pulling out the KEGG ids
kegg_ids = uniprot_annotated_protein_info |>
  dplyr::filter(!is.na(xref_kegg)) |>
  dplyr::pull(xref_kegg)


# Getting the terms for each annotation
kegg_annotations = conduitR::get_kegg_in_batches(kegg_ids)

# Joining back to original annotation file
kegg_annotations_final = uniprot_annotated_protein_info |>
  tidyr::separate_rows(sep =";",xref_kegg) |>
  dplyr::left_join(kegg_annotations,by = c("xref_kegg" = "kegg_id")) |>
  dplyr::select(kegg_pathway,kegg_id = xref_kegg,ko,code,
                dplyr::all_of(shared_columns)) |>
  dplyr::filter(kegg_id != "")

# Writing to File
readr::write_delim(kegg_annotations_final,kegg_annotations_out,delim = "\t")
