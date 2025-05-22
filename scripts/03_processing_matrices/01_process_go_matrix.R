# Getting inputs and outputs from snakemake.
report_pg_matrix=snakemake@input[["report_pg_matrix"]]
go_annotations = snakemake@input[["go_annotations"]]
go_annotations_matrix_fp = snakemake@output[["go_annotations_matrix"]]
go_annotations_matrix_taxa_fp = snakemake@output[["go_annotations_taxa_matrix"]]

#
#  # Testing
#  report_pg_matrix="output/03_diann_output/report.pg_matrix.tsv"
#  go_annotations = "user_input/00_database_resources/detected_protein_resources/03_go_annotations.txt"
#  go_annotations_matrix_fp = "output/05_output_files/go_matrix.tsv"
#  go_annotations_matrix_taxa_fp = "output/05_output_files/go_taxa_matrix.tsv"

report_pg_matrix = readr::read_tsv(report_pg_matrix)
go_annotations = readr::read_delim(go_annotations)

# Finding first and last column containing data
first_col = which(names(report_pg_matrix) == "First.Protein.Description") + 1
last_col = length(report_pg_matrix)

selected_columns = names(report_pg_matrix)[first_col:last_col]

# Dividing each intensity by the number of proteins it maps to
report_pg_matrix_mod <- report_pg_matrix |>
  dplyr::mutate(
    # Median centering to adjust for sample loading.
    dplyr::across(first_col:last_col,
                  ~ .x - median(.x, na.rm = TRUE)),
    # Adjusting negative values to be positive. Offsetting by smallest value.
    dplyr::across(first_col:last_col,
                  ~ .x + abs(min(.x, na.rm = TRUE)) + 1)  # Adjusting to make all values positive
  ) |>
  dplyr::mutate(master_protein_id = gsub(";.*","",Protein.Group),.before = "Protein.Group")

combined = dplyr::inner_join(go_annotations,
                             report_pg_matrix_mod,
                             by = c("protein_id" = "master_protein_id"))

go_annotations_matrix = combined |>
  dplyr::mutate(go_id = paste0(organism_type,"_",description)) |>
  dplyr::group_by(go_id,organism_type,description) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(selected_columns),
                   ~ sum(.x, na.rm = TRUE))
  )

# Writing to file
readr::write_tsv(go_annotations_matrix,go_annotations_matrix_fp)


go_annotations_matrix_taxa = combined |>
  dplyr::mutate(go_id = paste0(organism_type,"_", superkingdom,"_",kingdom,"_",
                               phylum,"_",class,"_",order,"_",family,"_",genus,
                               "_",species,"_",description),.before = 1

  ) |>
  dplyr::group_by(go_id,organism_type,description,superkingdom,kingdom,phylum,
                  class,order,family,genus,species) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(selected_columns),
                                 ~ sum(.x, na.rm = TRUE)))

# Writing to file
readr::write_tsv(go_annotations_matrix_taxa,go_annotations_matrix_taxa_fp)

