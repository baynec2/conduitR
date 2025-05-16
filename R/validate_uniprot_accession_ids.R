#' Validate UniProt Accession IDs
#'
#' Validates a vector of UniProt accession IDs using regular expressions to ensure
#' they match the official UniProt format. Invalid IDs are filtered out, and a
#' message reports how many were removed.
#'
#' @param uniprot_ids Character vector containing UniProt accession IDs to validate.
#'   The function checks for two valid formats:
#'   \itemize{
#'     \item Format 1: [OPQ][0-9][A-Z0-9]{3}[0-9] (e.g., "P12345", "O12345")
#'     \item Format 2: [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2} (e.g., "A0A023GPI8")
#'   }
#'
#' @return A character vector containing only the valid UniProt accession IDs.
#'   If any invalid IDs were found, a message is printed indicating how many
#'   were removed.
#'
#' @export
#'
#' @examples
#' # Validate a list of UniProt IDs:
#' # valid_ids <- validate_uniprot_accession_ids(c("P12345", "invalid", "A0A023GPI8"))
#' 
#' # The function will print a message if invalid IDs are found:
#' # "1 invalid UniProt IDs removed."
#' 
#' # The validated IDs can be used with other functions:
#' # protein_info <- extract_fasta_info("uniprot.fasta")
#' # valid_proteins <- protein_info[protein_info$protein_id %in% valid_ids,]
validate_uniprot_accession_ids <- function(uniprot_ids) {
  # Define the regex pattern for valid UniProt accessions
  uniprot_regex <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"

  # Validate UniProt accessions
  uniprot_ids_filtered <- uniprot_ids[stringr::str_detect(uniprot_ids, uniprot_regex)]
  n_invalid <- length(uniprot_ids) - length(uniprot_ids_filtered)
  if (n_invalid > 0) message(paste0(n_invalid, " invalid UniProt IDs removed."))

  return(uniprot_ids_filtered)
}
