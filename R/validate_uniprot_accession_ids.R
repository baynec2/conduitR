#' Validate UniProt Accession IDs
#'
#' Validates a vector of UniProt accession IDs using regular expressions to ensure
#' they match the official UniProt format. Invalid IDs are filtered out, and a
#' message reports how many were removed.
#'
#' @param uniprot_ids Character vector containing UniProt accession IDs to validate.
#'   The function checks for two valid formats via regex: (1) Swiss-Prot/TrEMBL
#'   style starting with O, P, or Q followed by digits and alphanumerics (e.g.
#'   "P12345", "O12345"); (2) TrEMBL style starting with A-N or R-Z, then digits
#'   and blocks of 4 characters (e.g. "A0A023GPI8").
#'
#' @return A character vector containing only the valid UniProt accession IDs.
#'   If any invalid IDs were found, a message is printed indicating how many
#'   were removed.
#'
#' @export
#'
#' @examples
#' # Valid IDs are kept; invalid ones removed (message printed)
#' validate_uniprot_accession_ids(c("P12345", "invalid", "A0A023GPI8"))
#' # Returns: c("P12345", "A0A023GPI8")
#'
#' # All valid: no message
#' validate_uniprot_accession_ids(c("P12345", "O12345", "A0A023GPI8"))
validate_uniprot_accession_ids <- function(uniprot_ids) {
  # Define the regex pattern for valid UniProt accessions
  uniprot_regex <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"

  # Validate UniProt accessions
  uniprot_ids_filtered <- uniprot_ids[stringr::str_detect(uniprot_ids, uniprot_regex)]
  n_invalid <- length(uniprot_ids) - length(uniprot_ids_filtered)
  if (n_invalid > 0) message(paste0(n_invalid, " invalid UniProt IDs removed."))

  return(uniprot_ids_filtered)
}
