# Inputs
protein_info = snakemake@input[["protein_info"]]
detected_protein_info = snakemake@input[["detected_protein_info"]]

# Defining outputs
database_taxonomy=snakemake@output[["database_taxonomy"]]
database_metrics=snakemake@output[["database_metrics"]]
detected_protein_taxonomy = snakemake@output[["detected_protein_taxonomy"]]
detected_protein_metrics = snakemake@output[["detected_protein_taxonomy"]]
combined_metrics = snakemake@output[["combined_metrics"]]


# Testing
protein_info = "output/00_database_resources/02_protein_info.txt"
detected_protein_info = "output/00_database_resources/detected_protein_resources/00_detected_protein_info.txt"
database_taxonomy="output/05_output_files/database_taxonomy.tsv"
database_metrics="output/05_output_files/database_metrics.tsv"
detected_protein_taxonomy = "output/05_output_files/detected_protein_taxonomy.tsv"
detected_protein_metrics = "output/05_output_files/detected_protein_metrics.tsv"
combined_metrics = "output/05_output_files/combined_metrics.tsv"

################################################################################
# Database information
################################################################################

protein_info = readr::read_tsv(protein_info)

# Generating database_taxonomy
database_taxonomy_out = protein_info |>
  dplyr::select(protein_id,organism_type,superkingdom,kingdom,phylum,class,
                order,family,genus,species)

# Writing to file
readr::write_tsv(database_taxonomy_out,database_taxonomy)

# Generating metrics
database_metrics_out= protein_info |>
  dplyr::select(-organism_name,-organism_id,-sequence) |>
  tidyr::pivot_longer(1:10,names_to = "metric") |>
  dplyr::group_by(metric) |>
  dplyr::summarise(n_in_db = dplyr::n_distinct(value)) |>
  dplyr::ungroup() |>
  dplyr::arrange(factor(metric, levels = c("organism_type","superkingdom","kingdom","phylum","class","order","family","genus","species","protein_id")))


# Writing to file
readr::write_tsv(database_metrics_out,database_metrics)

################################################################################
# Detected protein Information
################################################################################

detected_protein_info = readr::read_tsv(detected_protein_info)

# Generating database_taxonomy
detected_protein_taxonomy_out = detected_protein_info |>
  dplyr::select(protein_id,organism_type,superkingdom,kingdom,phylum,class,
                order,family,genus,species)

# Writing to file
readr::write_tsv(detected_protein_taxonomy_out,detected_protein_taxonomy)

# Generating metrics
detected_protein_metrics_out= detected_protein_info |>
  dplyr::select(-organism_name,-organism_id,-sequence) |>
  tidyr::pivot_longer(1:10,names_to = "metric") |>
  dplyr::group_by(metric) |>
  dplyr::summarise(n_detected = dplyr::n_distinct(value)) |>
  dplyr::ungroup() |>
  dplyr::arrange(factor(metric, levels = c("organism_type","superkingdom","kingdom","phylum","class","order","family","genus","species","protein_id")))


# Writing to file
readr::write_tsv(detected_protein_metrics_out,detected_protein_metrics)

# Combining metrics into one file.
combined_metrics_out = dplyr::inner_join(database_metrics_out,detected_protein_metrics_out,
                            by = "metric") |>
  dplyr::mutate(per_detected = round((n_detected/n_in_db) * 100)) |>
  dplyr::ungroup() |>
  dplyr::arrange(factor(metric, levels = c("organism_type","superkingdom","kingdom","phylum","class","order","family","genus","species","protein_id")))

readr::write_tsv(combined_metrics_out,combined_metrics)



