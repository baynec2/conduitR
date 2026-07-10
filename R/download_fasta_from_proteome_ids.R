#' Download and Concatenate FASTA Files for Multiple Proteomes
#'
#' For each UniProt proteome ID, downloads the proteome FASTA from UniProt
#' (UniProtKB or UniParc), saves per-proteome FASTA files in a temporary
#' directory, then concatenates them into a single FASTA file. Also writes a
#' delimited file with proteome IDs and associated taxonomy/type information.
#'
#' Before concatenating (non-parallel path), the per-proteome download manifest
#' is checked for completeness: the function errors if the total delivered
#' sequence count is zero or grossly short of the summed expected `proteinCount`,
#' and warns when some proteomes failed to download completely. This prevents an
#' empty or truncated search database from silently flowing downstream into
#' DIA-NN (see issue #20). No output FASTA is written when the gate errors.
#'
#' @param proteome_ids Character vector of UniProt proteome IDs (e.g.
#'   `"UP000005640"`). Duplicates are removed with a warning.
#' @param parallel Logical. If `TRUE`, download proteomes in parallel across
#'   proteomes using `future::multisession` with `availableCores() - 1` workers.
#'   This is purely a speed knob: the completeness gate, concatenation, and
#'   metadata output are identical to the sequential path. Default `FALSE`.
#' @param proteome_id_destination_fp Character. Path for the output file
#'   containing proteome IDs and metadata (default: `getwd()` plus current
#'   date and `.txt`).
#' @param fasta_destination_fp Character. Path for the concatenated FASTA
#'   output (default: `getwd()` plus current date and `.fasta`).
#'
#' @return Invisibly, the annotated proteome metadata tibble (proteome IDs with
#'   taxonomy/type and per-proteome download source). As side effects, writes
#'   the concatenated FASTA to `fasta_destination_fp` and the proteome metadata
#'   to `proteome_id_destination_fp`, and removes the temporary directory.
#'   Behaviour is identical whether `parallel` is `TRUE` or `FALSE`. Errors
#'   (without writing the FASTA) if the completeness gate fails.
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
#' # Same result, downloaded in parallel across proteomes (faster)
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
  # Always clean up the temp directory when the function exits, including when
  # the completeness gate below aborts the build.
  on.exit(unlink(fasta_dir, recursive = TRUE), add = TRUE)

  ##############################################################################
  # Download every proteome into the temp directory. `parallel` is only a speed
  # knob: it selects the download strategy (across-proteome parallelism via
  # furrr vs sequential purrr) and produces the SAME per-proteome manifest
  # either way, which then flows through the same gate -> concatenate -> write
  # tail below. (Pages within a single proteome are fetched sequentially by
  # get_fasta_file regardless, since UniProt uses cursor pagination.)
  ##############################################################################
  if (parallel) {
    future::plan(future::multisession, workers = max(1L, future::availableCores() - 1L))
    on.exit(future::plan(future::sequential), add = TRUE) # reset even on error
    download_results <- furrr::future_map(
      proteome_ids_to_search,
      get_fasta_file,
      fasta_dir = fasta_dir,
      .progress = TRUE
    ) |> dplyr::bind_rows()
  } else {
    download_results <- purrr::map(
      proteome_ids_to_search,
      get_fasta_file,
      fasta_dir = fasta_dir,
      .progress = TRUE
    ) |> dplyr::bind_rows()
  }

  # Hard-fail gate (issue #20): refuse to build an empty or grossly truncated
  # search database. Runs BEFORE concatenation so a bad database is never
  # written to fasta_destination_fp; warns on partial (some-but-not-all)
  # failures and errors on an empty or grossly short total.
  check_search_db_completeness(download_results)

  fasta_files <- list.files(fasta_dir)
  n_downloaded_proteomes <- length(fasta_files)

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
      .default = .data$proteome_type
    )) |>
    dplyr::select(-"source")

  # Writing proteome ids that coorespond with taxa ids to file path.
  log_with_timestamp(paste0(
    "Writing downloaded proteome information to ",
    proteome_id_destination_fp
  ))

  readr::write_delim(annotated_downloads, file = proteome_id_destination_fp)

  invisible(annotated_downloads)
}
