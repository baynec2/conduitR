#' Get Proteome ID for a Single Organism
#'
#' Retrieves the proteome ID for a single organism from UniProt using its NCBI
#' taxonomy ID. This function is typically used internally by
#' \code{\link[get_proteome_ids_from_organism_ids]{get_proteome_ids_from_organism_ids}}
#' for batch processing of multiple organisms.
#'
#' @param id A single NCBI taxonomy ID (e.g., 9606 for human, 562 for E. coli).
#'   This ID is used to query UniProt's proteome database.
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item Proteome Id: UniProt proteome identifier
#'     \item Organism: Scientific name of the organism
#'     \item Organism Id: NCBI taxonomy ID
#'     \item Protein count: Number of proteins in the proteome
#'   }
#'   If the request fails, returns NULL and prints an error message.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Human (NCBI taxid 9606)
#' human_proteome <- get_proteome_id_from_organism_id(9606)
#'
#' # E. coli (562)
#' ecoli_proteome <- get_proteome_id_from_organism_id(562)
#'
#' # Batch: use get_proteome_ids_from_organism_ids instead
#' proteomes <- get_proteome_ids_from_organism_ids(c(9606, 562))
#' }
#'
#' @note
#' This function:
#' \itemize{
#'   \item Makes a single API request to UniProt
#'   \item Returns only the first proteome if multiple are found
#'   \item Handles API errors gracefully
#'   \item Is designed for use in batch processing
#' }
#'
#' For processing multiple organisms, use
#' \code{\link[get_proteome_ids_from_organism_ids]{get_proteome_ids_from_organism_ids}}
#' which includes additional features like:
#' \itemize{
#'   \item Parallel processing
#'   \item Reference proteome prioritization
#'   \item Duplicate ID handling
#'   \item Missing ID reporting
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[get_proteome_ids_from_organism_ids]{get_proteome_ids_from_organism_ids}}
#'     for processing multiple organisms
#'   \item \code{\link[get_all_reference_proteomes]{get_all_reference_proteomes}}
#'     for retrieving reference proteome information
#'   \item \code{\link[get_fasta_file]{get_fasta_file}} for downloading
#'     proteome FASTA files
#' }
get_proteome_id_from_organism_id <- function(id) {
  base_url <- "https://rest.uniprot.org/proteomes/search"
  proteome_types <- c(1, 2, 3, 4)  # 1 = reference, 2 = other, 3 = redundant, 4 = excluded

  for (ptype in proteome_types) {
    query <- paste0("organism_id:", id, " AND proteome_type:", ptype)

    resp <- httr2::request(base_url) |>
      httr2::req_url_query(
        query = query,
        format = "json",
        sort = "cpd asc, busco desc, protein_count desc",
        size = 1
      ) |>
      httr2::req_headers(accept = "application/json") |>
      httr2::req_perform()

    if (httr2::resp_status(resp) != 200) {
      warning("Request failed for ", id, " (status: ", httr2::resp_status(resp), ")")
      next
    }

    out <- httr2::resp_body_json(resp)

    # Skip empty results
    if (length(out$results) == 0) next

    # Take first result only
    result <- out$results[[1]]

    null_to_na <- function(x) if (is.null(x)) NA_character_ else x

    # Return structured output
    return(tibble::tibble(
      proteome_id = result$id,
      organism_id = id,
      organism = result$taxonomy$scientificName,
      protein_count = result$proteinCount,
      proteome_type = result$proteomeType,
      redundant_to = null_to_na(result$redundantTo),
      genome_assembly_id = result$genomeAssembly$assemblyId,
      genome_assembly_level = result$genomeAssembly$level,
      annotation_score = result$annotationScore

    ))
  }

  # If all attempts fail
  warning("No proteomes found for taxonomy id: ", id)
  return(tibble::tibble(
    proteome_id = NA_character_,
    organism_id = id,
    organism = NA_character_,
    protein_count = NA_integer_,
    proteome_type = NA_character_,
    redundant_to = NA_character_,
    genome_assembly_id = NA_character_,
    genome_assembly_level = NA_character_,
    annotation_score = NA_integer_
  ))
}
