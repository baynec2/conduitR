#' Download and Combine FASTA Files from UniProt Organism IDs
#'
#' Downloads proteome FASTA files from UniProt for a list of organism IDs and combines
#' them into a single FASTA database. This function is particularly useful for creating
#' custom protein sequence databases for metaproteomics analysis. It automatically
#' selects reference proteomes when available and handles parallel processing for
#' efficient downloads.
#'
#' @param organism_ids A numeric vector of NCBI taxonomy IDs (e.g., 9606 for human,
#'   562 for E. coli). These IDs are used to identify organisms in UniProt.
#' @param parallel Logical indicating whether to use parallel processing (default: FALSE).
#'   When TRUE, uses all available CPU cores minus one for downloading FASTA files.
#' @param destination_fp Character string specifying the path where the combined
#'   FASTA file should be saved (default: creates a file named with the current date
#'   in the working directory).
#'
#' @return A combined FASTA file containing protein sequences from all specified
#'   organisms. The file includes:
#'   \itemize{
#'     \item Protein sequences from reference proteomes when available
#'     \item Alternative proteomes for organisms without reference proteomes
#'     \item Standard UniProt FASTA headers with protein and organism information
#'   }
#'   The function also provides progress messages about:
#'   \itemize{
#'     \item File deletion if a file exists at the destination
#'     \item Proteome ID retrieval status
#'     \item FASTA file download progress
#'     \item File concatenation status
#'   }
#'
#' @export
#'
#' @examples
#' # Download human and E. coli proteomes:
#' # download_fasta_from_organism_ids(c(9606, 562))
#' 
#' # Download with custom file path and parallel processing:
#' # download_fasta_from_organism_ids(
#' #   organism_ids = c(9606, 562, 10090),  # Human, E. coli, Mouse
#' #   parallel = TRUE,
#' #   destination_fp = "custom_database.fasta"
#' # )
#' 
#' # Use the resulting FASTA file with other functions:
#' # extract_fasta_info("database.fasta")
#' # prepare_taxonomy_matricies("taxonomy.tsv", "peptides.tsv", "database.fasta")
#'
#' @note
#' This function:
#' \itemize{
#'   \item Automatically removes duplicate organism IDs
#'   \item Prioritizes reference proteomes when available
#'   \item Creates a temporary directory for intermediate files
#'   \item Handles API rate limits and connection errors
#'   \item Cleans up temporary files after completion
#' }
#' 
#' For large numbers of organisms or when downloading frequently, consider:
#' \itemize{
#'   \item Using parallel processing for faster downloads
#'   \item Specifying a custom destination path
#'   \item Monitoring available disk space
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[get_proteome_ids_from_organism_ids]{get_proteome_ids_from_organism_ids}}
#'     for retrieving proteome IDs
#'   \item \code{\link[get_fasta_file]{get_fasta_file}} for downloading individual
#'     FASTA files
#'   \item \code{\link[extract_fasta_info]{extract_fasta_info}} for parsing the
#'     resulting FASTA file
#' }
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
