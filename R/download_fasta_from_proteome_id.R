download_fasta_from_proteome_id <- function(proteome_ids,
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
  # Replicating any duplicates if they exist.
  proteome_ids_to_search <- unique(proteome_ids)
  if(length(proteome_ids_to_search) > length(proteome_ids)){
    warning("The protein ids supplied contain duplicates")
  }
  ################################################################################
  # Downloading all FASTA files into the temp directory
  ################################################################################
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

    # Checking to make sure the correct number of proteome ids were downloaded.
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

    # Assembling proteome id dataframe (organism_id included)
    final_proteome_df <- get_proteome_taxids_and_types(proteome_ids_to_search)

    success <- gsub(".fasta", "", fasta_files)

    # annotating with whether or not the proteome was downloaded
    annotated_downloads <- final_proteome_df |>
      dplyr::mutate(download_info = dplyr::case_when(
        proteome_id %!in% success ~ "not_downloaded",
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
}

# Useage

proteome_ids <- readr::read_delim(("../Desktop/proteome_ids.txt"))$proteome_id
