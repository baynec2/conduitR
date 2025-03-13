#' annotate_uniprot_ids
#'
#' @param uniprot_ids
#' @param columns
#' @param batch_size
#' @param parallel
#'
#' @returns
#' @export
#'
#' @examples
annotate_uniprot_ids <- function(uniprot_ids,
                                 columns = NULL,
                                 batch_size = 150,
                                 parallel = TRUE) {
  # Validate uniprot IDs
  uniprot_ids_filtered <- validate_uniprot_accession_ids(uniprot_ids)

  # Split into batches of batch size.
  batches <- split(
    uniprot_ids_filtered,
    ceiling(seq_along(uniprot_ids_filtered) / batch_size)
  )

  # Enable parallelization if requested
  if (parallel) {
    future::plan(future::multisession, workers = parallel::detectCores() - 1) # Use max cores - 1
    results <- furrr::future_map_dfr(batches, get_uniprot_data,
      columns = columns,
      batch_size = batch_size,
      .progress = TRUE
    )
    future::plan(future::sequential) # Reset to sequential processing
  } else {
    results <- purrr::map_dfr(batches, get_uniprot_data,
      columns = columns,
      batch_size = batch_size,
      .progress = TRUE
    )
  }

  # Convert to tibble
  output <- tibble::as_tibble(results)

  # Fix column names if provided
  if (!is.null(columns)) {
    col_names <- unlist(strsplit(columns, ","))
    names(output) <- col_names
  }

  return(output)
}
