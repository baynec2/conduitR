################################################################################
# Getting Taxonomy From Organism IDs
################################################################################
# Get input and output files from Snakemake workflow
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Read organism IDs from the input file
organism_txt <- readr::read_delim(input_file,
                                  col_types = "cc")

organism_ids = organism_txt |>
  dplyr::pull(organism_id)

# Pull all taxonomy information
taxonomy = conduitR::get_ncbi_taxonomy(organism_ids) |>
  dplyr::inner_join(organism_txt,by = "organism_id") |>
  dplyr::select("organism_type","organism_id","domain","kingdom","phylum",
                "class","order","family","genus","species"
    )

# Writing to file.
readr::write_delim(taxonomy,output_file)
