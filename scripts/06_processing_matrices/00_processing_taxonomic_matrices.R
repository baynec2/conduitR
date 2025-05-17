# Getting inputs and outputs from snakemake.
report_pr_matrix= snakemake@input[["report_pr_matrix"]]
protein_info=snakemake@input[["protein_info"]]

superkingdom_matrix=snakemake@output[["superkingdom_matrix"]]
kingdom_matrix=snakemake@output[["kingdom_matrix"]]
phylum_matrix=snakemake@output[["phylum_matrix"]]
class_matrix=snakemake@output[["class_matrix"]]
order_matrix=snakemake@output[["order_matrix"]]
family_matrix=snakemake@output[["family_matrix"]]
genus_matrix=snakemake@output[["genus_matrix"]]
species_matrix=snakemake@output[["species_matrix"]]


# # # Testing
# report_pr_matrix="output/03_diann_output/report.pr_matrix.tsv"
# protein_info="user_input/00_database_resources/02_protein_info.txt"
# #
# superkingdom_matrix="output/05_output_files/superkingdom_matrix.tsv"
# kingdom_matrix="output/05_output_files/kingdom_matrix.tsv"
# phylum_matrix="output/05_output_files/phylum_matrix.tsv"
# class_matrix="output/05_output_files/class_matrix.tsv"
# order_matrix="output/05_output_files/order_matrix.tsv"
# family_matrix="output/05_output_files/family_matrix.tsv"
# genus_matrix="output/05_output_files/genus_matrix.tsv"
# species_matrix="output/05_output_files/species_matrix.tsv"


# Reading in files
report_pr_matrix = readr::read_tsv(report_pr_matrix) |>
  dplyr::mutate(n_proteins_mapped_to = stringr::str_count(Protein.Ids, ";") + 1,.after = Protein.Ids)


# Dividing each intensity by the number of proteins it maps to
report_pr_matrix_mod <- report_pr_matrix |>
  dplyr::mutate(
    # Median centering to adjust for sample loading.
    dplyr::across(12:ncol(report_pr_matrix),
                  ~ .x - median(.x, na.rm = TRUE)),
    # Adjusting negative values to be positive. Offsetting by smallest value.
    dplyr::across(12:ncol(report_pr_matrix),
                  ~ .x + abs(min(.x, na.rm = TRUE)) + 1)  # Adjusting to make all values positive
  ) |>
  tidyr::separate_rows("Protein.Ids", sep = ";")

# Reading in protein info file
protein_info = readr::read_delim(protein_info)

# Combining protein info file with report_pr_matrix
combined_data = dplyr::inner_join(protein_info,report_pr_matrix_mod,
                                  by = c("protein_id" = "Protein.Ids"))

# Finding the first and last column holding the data
first_col = which(names(combined_data) == "Precursor.Id") + 1
last_col = length(combined_data)

# Selecting only the necessary columns
combined_data = combined_data |>
  dplyr::select("organism_type","superkingdom","kingdom","phylum","class","order","family",
                "genus","species",first_col:last_col)

# Defining the taxonomic ranks of interest
taxonomic_ranks = c("superkingdom","kingdom","phylum","class","order","family",
                    "genus","species")

# Loop through taxonomic ranks, save resulting matrices
for (i in taxonomic_ranks) {
  # Summarize by organism_type and superkingdom, applying sum to the relevant columns
  summed_data <- combined_data |>
    dplyr::group_by(organism_type, !!dplyr::sym(i)) |>
    dplyr::summarise(
      dplyr::across(10:ncol(combined_data)-2, .fns = \(x) sum(x, na.rm = TRUE))
    ) |>
    dplyr::mutate(taxa_id = paste0(organism_type,"_",!!dplyr::sym(i)),.before = 1)
readr::write_tsv(summed_data,paste0("output/05_output_files/",i,"_matrix.tsv"))
}





