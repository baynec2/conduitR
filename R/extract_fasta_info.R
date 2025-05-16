#' Extract Information from FASTA Files
#'
#' Parses a FASTA file and extracts key information including protein IDs,
#' organism names, organism IDs, and protein sequences. This function is
#' particularly useful for processing UniProt FASTA files and extracting
#' metadata from the headers.
#'
#' @param fasta_file Character string specifying the path to a FASTA file.
#'   The file should contain protein sequences with headers in UniProt format,
#'   which typically include:
#'   \itemize{
#'     \item Protein ID (e.g., "sp|P12345|PROT_NAME")
#'     \item Organism name (OS=...)
#'     \item Organism ID (OX=...)
#'   }
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item protein_id: UniProt protein accession ID
#'     \item organism_name: Scientific name of the organism
#'     \item organism_id: NCBI taxonomy ID of the organism
#'     \item sequence: Amino acid sequence of the protein
#'   }
#'   Missing values (NA) are used for any information that cannot be extracted
#'   from the headers.
#'
#' @export
#'
#' @examples
#' # Extract information from a UniProt FASTA file:
#' # protein_info <- extract_fasta_info("uniprot_proteins.fasta")
#' 
#' # The extracted information can be used for:
#' # - Creating taxonomy databases
#' # - Mapping proteins to organisms
#' # - Analyzing sequence properties
#' 
#' # Note: The function expects UniProt-style headers. For other formats,
#' # the extraction of organism information may not work as expected.
extract_fasta_info <- function(fasta_file) {
  # Read the FASTA file
  fasta_data <- Biostrings::readAAStringSet(fasta_file)
  sequences <- as.character(fasta_data) # Convert the sequences to character strings

  # Extract headers and sequences
  headers <- names(fasta_data)

  # Extract organism name (OS=...)
  organism_name <- ifelse(
    grepl("OS=", headers),
    sub(".*OS=([^=]+?)\\s*OX=.*", "\\1", headers),
    NA
  )
  organism_name <- trimws(organism_name) # Remove extra spaces

  # Extract organism ID (OX=...)
  organism_id <- ifelse(
    grepl("OX=", headers),
    sub(".*OX=([0-9]+).*", "\\1", headers),
    NA
  )
  organism_id <- as.integer(organism_id) # Convert to integer

  # Extract protein ID (first field in UniProt-style headers)
  protein_id <- ifelse(
    grepl("\\|", headers),
    sub("^.*\\|([^|]+)\\|.*$", "\\1", headers),
    NA
  )


  # Create a DataFrame
  df <- data.frame(
    protein_id = protein_id,
    organism_name = organism_name,
    organism_id = organism_id,
    sequence = sequences,
    stringsAsFactors = FALSE
  )

  return(df)
}
