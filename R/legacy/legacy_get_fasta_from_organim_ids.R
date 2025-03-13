#' get_fasta_from_organism_ids
#'
#' Given a vector of organism ids, pull down a list of proteomes from uniprot
#' for each organism. The first one from the list is then downloaded as a fasta
#' file. This is then repeated for each organism id and concatintated into a
#' combined fasta file.
#'
#' This is particularly useful for generating a metaproteomics .fasta database
#'
#' @param organism_ids vector of organism ids
#' @param destination_fp path to the destination fasta file. Must be .fasta
#' @param additonal_organism_ids vector of any additional organism ids to include
#'
#' @returns a downloaded .fasta file to user specified file path.
#' @export
#'
#' @examples
#' Making a database with Human and E.coli
#'get_fasta_from_organism_ids(organism_ids = c(9606, 562))
legacy_get_fasta_from_organism_ids <- function(organism_ids,
                                        destination_fp = paste0(getwd(),"/",
                                                                Sys.Date(),
                                                                ".fasta")) {
# Check if file exists - if so delete it.
  if (file.exists(destination_fp)) {
    unlink(destination_fp)
    cat("File deleted at", destination_fp, "\n")
  }
  # Making sure organism_ids are unique
  organism_ids <- unique(organism_ids)

  # Downloading list of reference proteomes
  rp = get_all_reference_proteomes() |>
    dplyr::pull(.data$Proteome_ID)

  # Create temp directory at the destination of database
  temp_dir <- paste0(dirname(destination_fp), "/temp/")
  dir.create(temp_dir)

  ##############################################################################
  # Download a FASTA file for each organism id
  ##############################################################################
  for (i in seq_along(organism_ids)) {
    organism_id <- organism_ids[i]
    message("Started downloading fasta for organism id: ",organism_id)
    tryCatch({
      # UniProt API request to get the Proteome ID
      request_url <- paste0(
        "https://rest.uniprot.org/proteomes/search?query=organism_id:",
        organism_id, "&format=tsv"
      )

      # Requesting data and parsing it
      response <- httr2::request(request_url) |> httr2::req_perform()
      org_to_rp <- read.csv(text = httr2::resp_body_string(response), header = TRUE, sep = "\t")

      # Taking first proteome id
      if(org_to_rp$Proteome.Id[1] %in% rp){
        first_rp <- org_to_rp |>
          dplyr::filter(.data$Proteome.Id %in% rp) |>
          dplyr::pull(.data$Proteome.Id)
        message("Reference proteome found for this organism!")
      } else {
        # If we can't find a reference, take the first one returned.
        first_rp <- org_to_rp$Proteome.Id[1]
        message("Reference proteome not found for this organism")
      }

      # Proceed with downloading the FASTA file if the proteome ID is found
      if (!is.na(first_rp)) {
        fasta_url <- paste0(
          "https://rest.uniprot.org/uniprotkb/stream?query=proteome:",
          first_rp, "&format=fasta"
        )

        # Downloading the reference proteome
        fasta_fp <- file.path(temp_dir, paste0(first_rp, ".fasta"))
        download_response <- httr2::request(fasta_url) |> httr2::req_perform()

        # Save FASTA content to file
        if (httr2::resp_status(download_response) < 400) {
          writeLines(httr2::resp_body_string(download_response), fasta_fp)
          message("Fasta File for organism:", organism_id," + Proteome ID: ",first_rp, " sucessfully downloaded")


        } else {
          message("Failed to download FASTA for Proteome: ", first_rp)
          message("Error code: ", httr2::resp_status(download_response))
        }
      }

      Sys.sleep(0.5)  # Throttle requests to avoid hitting rate limits
    }, error = function(e) {
      message(paste("Error processing organism ID", organism_id, ":",
                    e$message))
    })
  }

  ##############################################################################
  # Concatenate all FASTA files in the temp directory into one file
  ##############################################################################
  fasta_files <- list.files(temp_dir, pattern = "\\.fasta$", full.names = TRUE)
  if (length(fasta_files) > 0) {
    # Loop through each FASTA file and append its content
    for (i in fasta_files) {
      fasta_content <- readLines(i)
      if (file.exists(destination_fp)) {
        readr::write_lines(fasta_content, destination_fp, append = TRUE)
      } else {
        readr::write_lines(fasta_content, destination_fp)
      }
    }
    cat("FASTA files concatenated into:", destination_fp, "\n")

    # Delete the temp directory and its contents
    unlink(temp_dir, recursive = TRUE)
    cat("Temporary directory deleted.\n")
  } else {
    cat("No FASTA files found in the temporary directory.\n")
  }
}
