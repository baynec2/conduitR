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
#' @param fasta_destination_fp Character string specifying the path where the combined
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
#' #   fasta_destination_fp = "custom_database.fasta"
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
                                             proteome_id_destination_fp = paste0(
                                               getwd(), "/",
                                               Sys.Date(),
                                               ".txt"
                                             ),
                                             fasta_destination_fp = paste0(
                                               getwd(), "/",
                                               Sys.Date(),
                                               ".fasta"
                                             )) {
  # Check if proteome id file exists - if so delete it.
  if (file.exists(proteome_id_destination_fp)) {
    unlink(proteome_id_destination_fp)
    log_with_timestamp(paste0("File deleted at", proteome_id_destination_fp))
  }
  # Check if fasta file exists - if so delete it.
  if (file.exists(fasta_destination_fp)) {
    unlink(fasta_destination_fp)
    log_with_timestamp(paste0("File deleted at", fasta_destination_fp))
  }
  ##############################################################################
  # Generating List of Proteome IDs
  ##############################################################################
  log_with_timestamp("Getting Proteome IDs cooresponding to organism IDs")

  # Only looking at unique organism ids
  organism_ids <- unique(organism_ids)

  # Getting proteome_ids
  proteome_id_df <- get_proteome_ids_from_organism_ids(organism_ids,
    parallel = parallel
  ) |>
    dplyr::mutate(selected_proteome_id = dplyr::case_when(is.na(redundant_to) ~
      proteome_id, TRUE ~ redundant_to))

  # Figuring out what level of taxonomy the proteome ids are at
  taxon_df <- organism_ids |>
    purrr::map_dfr(.f = get_parent_taxonomy_id, .progress = TRUE)

  proteome_taxa_level <- dplyr::left_join(proteome_id_df,
    taxon_df,
    by = c("organism_id" = "child_id")
  )
  ################################################################################
  # Dealing with Stain level taxonomy
  ################################################################################
  strain <- proteome_taxa_level |>
    dplyr::filter(child_rank == "strain")

  if(nrow(strain > 0)){
  # Defining reference proteome classifications
  good_pt_filter <- c(
    "Representative proteome",
    "Reference and representative proteome"
  )

  # User supplied strains that are representative are good proteome ids
  log_with_timestamp("Determining what strains have good proteomes")

  good_strain_proteome_df <- strain |>
    dplyr::filter(proteome_type %in% good_pt_filter)

  # Strain proteome ids that are not reference or representative or NA
  bad_strain_proteome_ids_df <- strain |>
    dplyr::filter(proteome_type %!in% good_pt_filter)

  # Checking parent taxonomy
  parents_of_bad_ids <- bad_strain_proteome_ids_df |>
    dplyr::pull(parent_id) |>
    unique()

  log_with_timestamp("Getting species level proteomes of strains with bad proteomes")
  # Getting proteomes at the parent id
  if(length(parents_of_bad_ids > 1)){
  parent_proteome_df <- get_proteome_ids_from_organism_ids(parents_of_bad_ids,
    parallel = parallel
  ) |>
    dplyr::mutate(
      selected_proteome_id = dplyr::case_when(is.na(redundant_to) ~
        proteome_id, TRUE ~ redundant_to),
      child_rank = "species"
    )
  } else {
    parent_proteome_df <- dplyr::tibble(
      organism_id = integer(),
      proteome_id = character(),
      proteome_type = character(),
      redundant_to = character(),
      parent_id = integer(),
      child_rank = character(),
      selected_proteome_id = character()
    )
}
  # Figuring out which have representative proteomes at the species level
  good_parent_df <- parent_proteome_df |>
    dplyr::filter(proteome_type %in% good_pt_filter)

  log_with_timestamp("Determining what species level proteomes are better than those at the strain level")

  # Figuring out what proteomes to just use the strain level info from.
  strain_level_ids_okay_df <- bad_strain_proteome_ids_df |>
    dplyr::filter(parent_id %!in% good_parent_df$organism_id) |>
    dplyr::mutate(selected_proteome_id = dplyr::case_when(is.na(redundant_to) ~
      proteome_id, TRUE ~ redundant_to))

  # Defining the ids
  strain_level_ids_okay <- strain_level_ids_okay_df |>
    dplyr::filter(!is.na(proteome_id)) |>
    dplyr::pull(selected_proteome_id)

  # Selecting the ids to search
  strain_resolved_ids <- c(
    good_strain_proteome_df$proteome_id,
    good_parent_df$proteome_id,
    strain_level_ids_okay
  )
  } else {
    strain_resolved_ids <- c()
    good_strain_proteome_df <- dplyr::tibble()
    good_parent_df <- dplyr::tibble()
    strain_level_ids_okay_df <- dplyr::tibble()
  }

  ################################################################################
  # Dealing with Other taxonomic levels (Species or subspecies)
  ################################################################################
  log_with_timestamp("Dealing with non-strain level proteomes")

  species_df <- proteome_taxa_level |>
    dplyr::filter(child_rank != "strain")

  if(nrow(species_df) >0){

  species_id <- species_df |>
    dplyr::filter(!is.na(selected_proteome_id)) |>
    dplyr::pull(selected_proteome_id)

  }else{
    species_id <- c()
  }

  # Proteome ids to search
  proteome_ids_to_search <- c(strain_resolved_ids, species_id)

  n_ids_provided <- length(organism_ids)

  log_with_timestamp("Proteome IDs Sucessfully Retrieved")
  log_with_timestamp(glue::glue(
    "{n_ids_provided} unique taxonomy IDs were provided. ",
    "UniProtKB proteomes were successfully retrieved for {length(proteome_ids_to_search)} taxa, ",
    "while {n_ids_provided - length(proteome_ids_to_search)} had no corresponding proteome available."
  ))

  ##############################################################################
  # Downloading all FASTA files into the temp directory
  ##############################################################################
  fasta_dir <- paste0(dirname(fasta_destination_fp), "/temp")
  if (!dir.exists(fasta_dir)) {
    dir.create(fasta_dir)
  } else {
    log_with_timestamp("temp/ dir already existed, old version was deleted")
    unlink(fasta_dir)
    dir.create(fasta_dir)
  }
  if (parallel) {
    future::plan(future::multisession, workers = future::availableCores() - 1)
    furrr::future_map_dfr(proteome_ids_to_search,
      get_fasta_file,
      fasta_dir = fasta_dir,
      .progress = TRUE
    )
    future::plan(future::sequential) # Reset to sequential processing
  } else {
    purrr::map_dfr(
      proteome_ids_to_search,
      function(id, ...) {
        get_fasta_file(id, ...)
      },
      fasta_dir = fasta_dir,
      .progress = TRUE
    )
  }
  ##############################################################################
  # Retry for any files that were not downloaded the first time
  ##############################################################################

  fasta_files <- list.files(fasta_dir)
  n_downloaded_proteomes <- length(fasta_files)

  if (n_downloaded_proteomes != length(proteome_ids_to_search)) {
    warning(glue::glue(
      "{length(proteome_ids_to_search) - n_downloaded_proteomes} proteomes were not successfully downloaded."
    ))
  }


  ##############################################################################
  # Concatenating all files in temp dir
  ##############################################################################

  log_with_timestamp(glue::glue(
    "Concatenating FASTA files for {n_downloaded_proteomes} downloaded proteomes..."
  ))

  concatenate_fasta_files(fasta_dir, fasta_destination_fp)

  # Assembling proteome id dataframe
  final_proteome_df <- dplyr::bind_rows(
    good_strain_proteome_df,
    good_parent_df,
    strain_level_ids_okay_df,
    species_df
  )

  success <- gsub(".fasta", "", fasta_files)

  # annotating with whether or not the proteome was downloaded
  annotated_downloads <- final_proteome_df |>
    dplyr::mutate(download_info = dplyr::case_when(
      selected_proteome_id %!in% success ~ "not_downloaded",
      TRUE ~ proteome_type
    ))

  # Writing proteome ids that coorespond with taxa ids to file path.
  log_with_timestamp(paste0(
    "Writing downloaded proteome information to ",
    proteome_id_destination_fp
  ))

  readr::write_delim(annotated_downloads, file = proteome_id_destination_fp)

  # Delete the temp directory and its contents
  unlink(fasta_dir, recursive = TRUE)
  log_with_timestamp("Temporary directory deleted.\n")
}
