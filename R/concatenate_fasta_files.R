#' concatenate_fasta_files
#' concatenate fasta files.
#'
#' @param directory
#' @param output_file
#'
#' @returns
#' @export
#'
#' @examples
concatenate_fasta_files <- function(fasta_dir, destination_fp) {
  files <- list.files(fasta_dir, full.names = TRUE) # Get all file paths.
  # Read all lines from each file
  all_fasta <- purrr::map(files, readLines) |>
    unlist() |> # Flatten the list of character vectors
    paste(collapse = "\n") # Concatenate all lines into a single string
  # Write the concatenated content to the output file
  writeLines(all_fasta, destination_fp)
}
