#' Get Protein Data from UniProt API
#'
#' Retrieves protein information from UniProt's REST API for a batch of accession IDs.
#' This function is typically used internally by \code{\link[annotate_uniprot_ids]{annotate_uniprot_ids}}
#' to fetch data in batches, but can also be used directly for custom queries.
#'
#' @param ids Character vector containing UniProt accession IDs to query. These should
#'   be valid UniProt accession IDs (e.g., "P01308", "A0A023GPI8"). The function
#'   will process these IDs in a single API request.
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
#' @param batch_size Integer specifying the maximum number of IDs to include in a
#'   single API request (default: 150). This helps manage API rate limits and
#'   response size.
#'
#' @return A tibble containing the requested UniProt information for each accession ID.
#'   The columns correspond to the requested fields, and rows represent individual
#'   proteins. If the API request fails, returns NULL and prints an error message
#'   with the status code.
#'
#' @export
#'
#' @examples
#' # Get basic information for a single protein:
#' # insulin_data <- get_uniprot_data("P01308")
#' 
#' # Request specific fields for multiple proteins:
#' # protein_data <- get_uniprot_data(
#' #   c("P01308", "P01325"),
#' #   columns = "accession,protein_name,gene_primary,go"
#' # )
#' 
#' # Use with other functions:
#' # annotated_data <- annotate_uniprot_ids(
#' #   c("P01308", "P01325"),
#' #   columns = "accession,protein_name,go"
#' # )
#'
#' @note
#' This function:
#' \itemize{
#'   \item Makes a single API request to UniProt
#'   \item Handles API errors gracefully
#'   \item Returns data in a tabular format
#'   \item Is designed for batch processing
#' }
#' 
#' Important considerations:
#' \itemize{
#'   \item The function requires an internet connection
#'   \item API rate limits may affect request success
#'   \item Large batches may timeout or return incomplete data
#'   \item Invalid IDs are silently excluded from results
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[annotate_uniprot_ids]{annotate_uniprot_ids}} for
#'     processing multiple batches with parallel support
#'   \item \code{\link[validate_uniprot_accession_ids]{validate_uniprot_accession_ids}}
#'     for validating UniProt accession IDs
#'   \item \code{\link[get_kegg_in_batches]{get_kegg_in_batches}} for
#'     retrieving KEGG annotations
#'   \item \code{\link[get_eggnog_function]{get_eggnog_function}} for
#'     retrieving EggNOG annotations
#' }
get_uniprot_data <- function(ids, columns = NULL, batch_size = 150) {
  # Base URL for request
  base_url <- "https://rest.uniprot.org/uniprotkb/search?"
  query <- paste0(paste0("accession:", ids, collapse = " OR "))
  req <- httr2::request(base_url) |>
    httr2::req_url_query(
      query = query,
      format = "tsv",
      fields = columns,
      size = batch_size # Max limit per request
    ) |>
    httr2::req_perform()

  # Parse response
  if (httr2::resp_status(req) == 200) {
    httr2::resp_body_string(req) |>
      readr::read_tsv(show_col_types = FALSE)
  } else {
    message("Request failed, status code: ", httr2::resp_status(req))
    return(NULL)
  }
}
