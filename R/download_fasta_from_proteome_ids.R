#' Download and Concatenate FASTA Files for Multiple Proteomes
#'
#' For each UniProt proteome ID, downloads the proteome FASTA from UniProt
#' (UniProtKB or UniParc), saves per-proteome FASTA files in a temporary
#' directory, then concatenates them into a single FASTA file. Also writes a
#' delimited file with proteome IDs and associated taxonomy/type information.
#'
#' @param proteome_ids Character vector of UniProt proteome IDs (e.g.
#'   `"UP000005640"`). Duplicates are removed with a warning.
#' @param parallel Logical. If `TRUE`, download proteomes in parallel using
#'   `future::multisession` with `availableCores() - 1` workers. Default
#'   `FALSE`.
#' @param proteome_id_destination_fp Character. Path for the output file
#'   containing proteome IDs and metadata (default: `getwd()` plus current
#'   date and `.txt`).
#' @param fasta_destination_fp Character. Path for the concatenated FASTA
#'   output (default: `getwd()` plus current date and `.fasta`).
#'
#' @return When `parallel` is `FALSE`, returns the result of
#'   `concatenate_fasta_files` after writing the proteome info file; the temp
#'   directory is removed. When `parallel` is `TRUE`, returns the combined
#'   result of `get_fasta_file` across proteomes (no concatenation or
#'   proteome info file in that path). In both cases, FASTA and metadata files
#'   are written to the specified paths.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Download two proteomes and write one FASTA + one metadata file
#' download_fasta_from_proteome_ids(
#'   c("UP000005640", "UP000000625"),
#'   fasta_destination_fp = "my_proteomes.fasta",
#'   proteome_id_destination_fp = "my_proteomes.txt"
#' )
#'
#' # Parallel downloads (no concatenation; check temp dir for FASTA files)
#' download_fasta_from_proteome_ids(
#'   c("UP000005640", "UP000000625"),
#'   parallel = TRUE
#' )
#' }
download_fasta_from_proteome_ids <- function(proteome_ids,
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
    download_results <- purrr::map_dfr(
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

    # Annotate using per-proteome source (uniprotkb vs uniparc) from get_fasta_file
    annotated_downloads <- final_proteome_df |>
      dplyr::left_join(
        download_results |> dplyr::select("proteome_id", "source"),
        by = "proteome_id"
      ) |>
      dplyr::mutate(download_info = dplyr::case_when(
        is.na(.data$source) | .data$source == "not_downloaded" ~ "not_downloaded",
        .data$source == "uniparc" ~ "uniparc",
        TRUE ~ .data$proteome_type
      )) |>
      dplyr::select(-"source")

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
