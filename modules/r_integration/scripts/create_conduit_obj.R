################################################################################
# Creating Conduit Object
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages
start_time <- Sys.time()

conduitR::log_with_timestamp("Running create_conduit_obj.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[2]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[3]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[4]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[5]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

# Defining inputs
qf = snakemake@input[["qf"]]
database_taxonomy = snakemake@input[["database_taxonomy"]]
database_metrics = snakemake@input[["database_metrics"]]
detected_protein_taxonomy = snakemake@input[["detected_protein_taxonomy"]]
detected_protein_metrics = snakemake@input[["detected_protein_metrics"]]

# Defining outputs
conduit_obj = snakemake@output[["conduit_obj"]]

conduitR::log_with_timestamp("Creating Conduit object")
conduit = conduitR::create_conduit_obj(qf,
                                       database_taxonomy,
                                       database_metrics,
                                       detected_protein_taxonomy,
                                       detected_protein_metrics)
conduitR::log_with_timestamp(paste0("Saving Conduit object to ", conduit_obj))
# Saving to File
saveRDS(conduit,conduit_obj)

end_time <- Sys.time()
conduitR::log_with_timestamp("create_conduit_obj.R script completed. Time taken: %.2f minutes", 
                              as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)