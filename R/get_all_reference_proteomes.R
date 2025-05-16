#' Get Current Reference Proteomes from UniProt
#'
#' Retrieves a list of all current reference proteomes from UniProt's FTP server.
#' This list is updated every 8 weeks according to UniProt's release schedule.
#' The function downloads and parses the README file from UniProt's reference
#' proteomes directory.
#'
#' @return A tibble containing the following columns:
#'   \itemize{
#'     \item Proteome_ID: UniProt proteome identifier
#'     \item Tax_ID: NCBI taxonomy ID
#'     \item OSCODE: Organism code
#'     \item SUPERREGNUM: Super-regnum classification
#'     \item #(1): Number of entries in main FASTA (canonical sequences)
#'     \item #(2): Number of entries in additional FASTA (isoforms)
#'     \item #(3): Number of entries in gene2acc mapping file
#'     \item Species Name: Scientific name of the organism
#'   }
#'
#' @export
#'
#' @examples
#' # Get the current list of reference proteomes:
#' # proteomes <- get_all_reference_proteomes()
#' 
#' # Filter for specific organisms:
#' # human_proteome <- proteomes[proteomes$`Species Name` == "Homo sapiens",]
#' 
#' # Use with other functions:
#' # concatenate_fasta_files("path/to/proteomes", "combined.fasta")
#' 
#' # Note: This function requires an internet connection to access
#' # UniProt's FTP server. The data is updated every 8 weeks.
#' 
#' @references
#' UniProt Reference Proteomes Documentation:
#' https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README
get_all_reference_proteomes <- function() {
  current_reference_proteomes <- readr::read_tsv("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README",
    skip = 153, skip_empty_rows = TRUE
  )
  if (sum(names(current_reference_proteomes) == c(
    "Proteome_ID", "Tax_ID", "OSCODE", "SUPERREGNUM", "#(1)", "#(2)", "#(3)",
    "Species Name"
  )) == 8) {
    return(current_reference_proteomes)
  } else {
    stop("Unexpected column names are returned from the file")
  }
}
