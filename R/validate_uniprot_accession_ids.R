#' validate_uniprot_accession_ids
#'
#' validates uniprot accession ids using regular expression. Filters out ones
#' that don't match.
#'
#' @param uniprot_ids
#'
#' @returns
#' @export
#'
#' @examples
validate_uniprot_accession_ids <- function(uniprot_ids) {
  # Define the regex pattern for valid UniProt accessions
  uniprot_regex <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"

  # Validate UniProt accessions
  uniprot_ids_filtered <- uniprot_ids[stringr::str_detect(uniprot_ids, uniprot_regex)]
  n_invalid <- length(uniprot_ids) - length(uniprot_ids_filtered)
  if (n_invalid > 0) message(paste0(n_invalid, " invalid UniProt IDs removed."))

  return(uniprot_ids_filtered)
}
