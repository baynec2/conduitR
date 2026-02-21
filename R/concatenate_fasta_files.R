#' Concatenate Multiple FASTA Files
#'
#' Combines multiple FASTA files from a directory into a single output file.
#' This is useful for merging protein sequence databases or combining
#' multiple reference proteomes into a single file.
#'
#' @param fasta_dir Character string specifying the path to a directory
#'   containing FASTA files to be concatenated
#' @param destination_fp Character string specifying the path where the
#'   concatenated FASTA file should be saved
#'
#' @return Invisibly returns NULL. The function creates a new file at
#'   destination_fp containing all sequences from the input files.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Combine all FASTA files in a folder into one file
#' concatenate_fasta_files("path/to/fasta_files", "combined.fasta")
#'
#' # Then extract metadata from the combined file
#' protein_info <- extract_fasta_info("combined.fasta")
#' }
concatenate_fasta_files <- function(fasta_dir, destination_fp) {
  files <- list.files(fasta_dir, full.names = TRUE) # Get all file paths.
  # Read all lines from each file
  all_fasta <- purrr::map(files, readLines) |>
    unlist() |> # Flatten the list of character vectors
    paste(collapse = "\n") # Concatenate all lines into a single string
  # Write the concatenated content to the output file
  writeLines(all_fasta, destination_fp)
}
