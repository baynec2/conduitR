################################################################################
# Processing Taxonomic Matrices
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running 00_processing_taxonomic_matrices.R script")

conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[2]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[2]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[3]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[4]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[5]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[6]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[7]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[8]]))


# Getting inputs and outputs from snakemake.
# Input files
report_pr_matrix= snakemake@input[["report_pr_matrix"]]
protein_info=snakemake@input[["protein_info"]]
# Output files
domain_matrix=snakemake@output[["domain_matrix"]]
kingdom_matrix=snakemake@output[["kingdom_matrix"]]
phylum_matrix=snakemake@output[["phylum_matrix"]]
class_matrix=snakemake@output[["class_matrix"]]
order_matrix=snakemake@output[["order_matrix"]]
family_matrix=snakemake@output[["family_matrix"]]
genus_matrix=snakemake@output[["genus_matrix"]]
species_matrix=snakemake@output[["species_matrix"]]

# Reading in files
conduitR::log_with_timestamp(paste0("Reading in report_pr_matrix from ", report_pr_matrix))
report_pr_matrix = readr::read_tsv(report_pr_matrix) |>
  dplyr::mutate(n_proteins_mapped_to = stringr::str_count(Protein.Ids, ";") + 1,.after = Protein.Ids)

conduitR::log_with_timestamp("Dividing each intensity by the number of proteins it maps to")
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
  dplyr::select("organism_type","domain","kingdom","phylum","class","order","family",
                "genus","species",first_col:last_col)

# Defining the taxonomic ranks of interest
taxonomic_ranks = c("domain","kingdom","phylum","class","order","family",
                    "genus","species")

# Loop through taxonomic ranks, save resulting matrices
for (i in taxonomic_ranks) {
  conduitR::log_with_timestamp(paste0("Processing ", i, " matrix"))
  # Summarize by organism_type and superkingdom, applying sum to the relevant columns
  summed_data <- combined_data |>
    dplyr::group_by(organism_type, !!dplyr::sym(i)) |>
    dplyr::summarise(
      dplyr::across(10:ncol(combined_data)-2, .fns = \(x) sum(x, na.rm = TRUE))
    ) |>
    dplyr::mutate(taxa_id = paste0(organism_type,"_",!!dplyr::sym(i)),.before = 1)
  conduitR::log_with_timestamp(paste0("Writing ", i, " matrix to ", paste0("output/05_output_files/",i,"_matrix.tsv")))
  readr::write_tsv(summed_data,paste0("output/05_output_files/",i,"_matrix.tsv"))
}

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed 00_processing_taxonomic_matrices.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)
