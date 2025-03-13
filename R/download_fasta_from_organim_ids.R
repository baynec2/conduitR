#' download_fasta_from_organism_ids
#'
#' Given a vector of organism ids, pull down a list of proteomes from uniprot
#' for each organism. The first one from the list is then downloaded as a fasta
#' file. This is then repeated for each organism id and concatenated into a
#' combined fasta file.
#'
#' This is particularly useful for generating a metaproteomics .fasta database
#'
#' @param organism_ids vector of organism ids
#' @param destination_fp path to the destination fasta file. Must be .fasta
#' @param additonal_organism_ids vector of any additional organism ids to include
#'
#' @returns a downloaded .fasta file to user specified file path.
#' @export
#'
#' @examples
#' Making a database with Human and E.coli
#'get_fasta_from_organism_ids(organism_ids = c(9606, 818))
download_fasta_from_organism_ids <- function(organism_ids,
                                        parallel = FALSE,
                                        destination_fp = paste0(getwd(),"/",
                                                                Sys.Date(),
                                                                ".fasta")) {
  # Check if file exists - if so delete it.
  if (file.exists(destination_fp)) {
    unlink(destination_fp)
    cat("File deleted at", destination_fp, "\n")
  }
  # Making sure organism_ids are unique
  organism_ids <- unique(organism_ids)
  ##############################################################################
  # Generating List of Proteome IDs
  ##############################################################################
  cat("Getting Proteome IDs cooresponding to organism IDs")
  proteome_ids = get_proteome_ids_from_organism_ids(organism_ids,
                                                    parallel = parallel )$`Proteome Id`
  cat("Proteome IDs Sucessfully Retrieved")
  ##############################################################################
  # Downloading all FASTA files into the temp directory
  ##############################################################################
  fasta_dir = paste0(dirname(destination_fp),"/temp")
  if(!dir.exists(fasta_dir)){
  dir.create(fasta_dir)
  }else{
    message("temp/ dir already existed, old version was deleted")
    unlink(fasta_dir)
    dir.create(fasta_dir)
  }
   if (parallel) {
    future::plan(future::multisession, workers = future::availableCores() - 1)
    furrr::future_map_dfr(proteome_ids,
                          get_fasta_file,
                          fasta_dir = fasta_dir ,
                          .progress = TRUE)
    future::plan(future::sequential) # Reset to sequential processing
  } else {
    purrr::map_dfr(proteome_ids,
               get_fasta_file,
               fasta_dir = fasta_dir,
               .progress = TRUE)
  }
  ##############################################################################
  # Concatenating all files in temp dir
  ##############################################################################
  cat("Concatinating fasta files...")
  concatenate_fasta_files(fasta_dir, destination_fp)
  # Delete the temp directory and its contents
    unlink(fasta_dir, recursive = TRUE)
    cat("Temporary directory deleted.\n")
  }
