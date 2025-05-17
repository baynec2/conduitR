# Get input and output files from Snakemake workflow
taxonomy_data_fp <- snakemake@input[[1]]
output_file_fp <- snakemake@output[[1]]

# #Testing
# taxonomy_data_fp = "user_input/00_database_resources/01_taxonomy.txt"
# output_file_fp = "user_input/00_database_resources/03_database_taxonomic_tree.pdf"
# Creating the plot
taxonomy_data_fp = readr::read_delim(taxonomy_data_fp) |>
  conduitR::plot_taxa_tree()

# Saving the plot
ggplot2::ggsave(output_file_fp,
                taxonomy_data_fp,
                height = 6,
                width = 6,
                units = "in")
