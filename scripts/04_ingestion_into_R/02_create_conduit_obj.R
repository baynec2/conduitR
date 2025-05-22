# Defining inputs
qf = snakemake@input[["qf"]]
database_taxonomy = snakemake@input[["database_taxonomy"]]
database_metrics = snakemake@input[["database_metrics"]]
detected_protein_taxonomy = snakemake@input[["detected_protein_taxonomy"]]
detected_protein_metrics = snakemake@input[["detected_protein_metrics"]]

# Defining outputs
conduit_obj = "output/05_output_files/conduit_output.rds"

# # Testing
# qf = "output/05_output_files/qf.rds"
# database_taxonomy = "output/05_output_files/database_taxonomy.tsv"
# database_metrics = "output/05_output_files/database_metrics.tsv"
# detected_protein_taxonomy = "output/05_output_files/detected_protein_taxonomy.tsv"
# detected_protein_metrics = "output/05_output_files/detected_protein_metrics.tsv"
# conduit_obj = "output/05_output_files/conduit_output.rds"


conduit = conduitR::create_conduit_obj(qf,
                                              database_taxonomy,
                                              database_metrics,
                                              detected_protein_taxonomy,
                                              detected_protein_metrics)

# Saving to File
saveRDS(conduit,conduit_obj)
