#' Annotate UniProt Accession IDs with Protein Information
#'
#' Retrieves detailed protein information from UniProt for a list of accession IDs.
#' This function efficiently fetches data in batches and supports parallel processing
#' for improved performance with large datasets. It automatically validates the input
#' IDs and handles API rate limits.
#'
#' @param uniprot_ids Character vector containing UniProt accession IDs to annotate.
#'   IDs are automatically validated using \code{\link[validate_uniprot_accession_ids]{validate_uniprot_accession_ids}}.
#' @param columns Character string specifying which UniProt fields to retrieve,
#'   separated by commas. If NULL, returns all available fields. Common fields include:
#'   \itemize{
#'     \item accession: UniProt accession number
#'     \item id: Entry name
#'     \item protein_name: Full protein name
#'     \item gene_primary: Primary gene name
#'     \item organism_name: Full organism name
#'     \item organism_id: NCBI taxonomy ID
#'     \item go: Gene Ontology terms
#'     \item xref_kegg: KEGG database cross-references
#'     \item xref_eggnog: EggNOG database cross-references
#'     \item cc_subcellular_location: Subcellular location
#'     \item cc_tissue_specificity: Tissue specificity
#'     \item xref_cazy: CAZy database cross-references
#'     \item xref_pfam: Pfam database cross-references
#'     \item xref_interpro: InterPro database cross-references
#'   }
#'   For a complete list of available fields, see the UniProt REST API documentation.
#' @param batch_size Integer specifying the number of IDs to process in each batch
#'   (default: 150). This helps manage API rate limits and memory usage.
#' @param parallel Logical indicating whether to use parallel processing (default: TRUE).
#'   When TRUE, uses all available CPU cores minus one for processing.
#'
#' @return A tibble containing the requested UniProt information for each valid
#'   accession ID. The columns correspond to the requested fields, and rows
#'   represent individual proteins. Invalid IDs are automatically filtered out,
#'   and a message indicates how many were removed.
#'
#' @export
#'
#' @examples
#' # Basic usage with default fields:
#' # protein_info <- annotate_uniprot_ids(c("P01308", "P01325"))
#' 
#' # Request specific fields:
#' # detailed_info <- annotate_uniprot_ids(
#' #   c("P01308", "P01325"),
#' #   columns = "accession,protein_name,gene_primary,go,xref_kegg"
#' # )
#' 
#' # Process a large list of IDs in parallel:
#' # all_proteins <- annotate_uniprot_ids(
#' #   protein_list,
#' #   batch_size = 200,
#' #   parallel = TRUE
#' # )
#' 
#' # Use the results with other functions:
#' # kegg_info <- get_kegg_in_batches(detailed_info$xref_kegg)
#' # eggnog_info <- get_eggnog_function(detailed_info$xref_eggnog)
#'
#' @note
#' This function:
#' \itemize{
#'   \item Automatically validates UniProt accession IDs
#'   \item Handles API rate limits through batching
#'   \item Supports parallel processing for improved performance
#'   \item Preserves the order of valid IDs in the output
#' }
#' For large datasets, consider:
#' \itemize{
#'   \item Adjusting batch_size based on your network connection
#'   \item Using parallel = TRUE for faster processing
#'   \item Requesting only necessary columns to reduce data transfer
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[validate_uniprot_accession_ids]{validate_uniprot_accession_ids}} for ID validation
#'   \item \code{\link[get_kegg_in_batches]{get_kegg_in_batches}} for KEGG annotation
#'   \item \code{\link[get_eggnog_function]{get_eggnog_function}} for EggNOG annotation
#' }
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
