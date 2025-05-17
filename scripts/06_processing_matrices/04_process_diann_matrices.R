# Defining files
# Inputs
report_pg_matrix = snakemake@input[[report_pg_matrix]]
report_pr_matrix = snakemake@input[[report_pr_matrix]]
# Outputs
protein_group_matrix = snakemake@output[[protein_group_matrix]]
precursor_matrix = snakemake@output[[precursor_matrix]]
peptide_matrix = snakemake@output[[peptide_matrix]]

# # # # Testing
# report_pg_matrix="output/03_diann_output/report.pg_matrix.tsv"
# report_pr_matrix ="output/03_diann_output/report.pr_matrix.tsv"
# protein_group_matrix = "output/05_output_files/protein_group_matrix.tsv"
# precursor_matrix = "output/05_output_files/precursor_matrix.tsv"
# peptide_matrix = "output/05_output_files/peptide_matrix.tsv"

# Converting precursor format - ID needs to be first.
precursor_matrix_out = readr::read_tsv(report_pr_matrix) |>
  dplyr::select(Precursor.Id,dplyr::everything())

readr::write_tsv(precursor_matrix_out,precursor_matrix)

# Creating peptide format - summing peptides with same sequence (including mods)
first_col = names(precursor_matrix_out)[which(names(precursor_matrix_out) == "Precursor.Charge") + 1]
last_col = names(precursor_matrix_out[ncol(precursor_matrix_out)])

peptide_matrix_out = precursor_matrix_out |>
  dplyr::group_by(Modified.Sequence,Protein.Group,Protein.Ids,Protein.Names,
                  First.Protein.Description,Stripped.Sequence) |>
  dplyr::summarise(dplyr::across(first_col:last_col,
                   ~ sum(.x, na.rm = TRUE))
  )

readr::write_tsv(peptide_matrix_out,peptide_matrix)

# Protein Group format stays the same
file.copy(report_pg_matrix,protein_group_matrix)



