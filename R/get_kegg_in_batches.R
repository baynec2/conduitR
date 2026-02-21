#' Fetch KEGG Gene/Pathway Information in Batches
#'
#' Queries the KEGG API for functional information (gene, orthology, pathway,
#' etc.) for a set of KEGG IDs. IDs can be semicolon-separated; they are split
#' and queried in batches to respect API limits.
#'
#' @param kegg_ids Character vector of KEGG IDs in format like `"hsa:1234"` or
#'   `"hsa:123;hsa:456"`. Semicolons are split and whitespace is trimmed.
#' @param batch_size Integer. Number of IDs per API batch (default: 10).
#'
#' @return A tibble with KEGG-derived columns (e.g. gene ID, organism, name,
#'   orthology, pathway). Failed or missing entries are omitted.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' insulin <- get_kegg_in_batches("hsa:3630")
#' get_kegg_in_batches(c("hsa:3630", "hsa:3640"), batch_size = 5)
#' }
get_kegg_in_batches <- function(kegg_ids, batch_size = 10) {
  # original vector (e.g. kegg_ids <- c("hsa:123", "hsa:456;hsa:789", ...))
  split_ids <- stringr::str_split(kegg_ids, ";", simplify = FALSE) |> unlist()
  # remove leading/trailing whitespace
  split_ids <- stringr::str_trim(split_ids)
  split_ids <- split_ids[split_ids != ""]

  # Initiating list
  all_results <- list()
  # Split IDs into batches
  batches <- split(split_ids, ceiling(seq_along(split_ids) / batch_size))
  for (i in seq_along(batches)) {
    message(sprintf("Processing batch %d of %d of KEGG Ids", i, length(batches)))
    batch <- batches[[i]]
    # Fetch KEGG entries
    entries <- lapply(batch, function(id) {
      tryCatch(KEGGREST::keggGet(id), error = function(e) NULL)
    })
    # Flatten list and remove NULLs
    entries <- unlist(entries, recursive = FALSE)
    entries <- entries[!sapply(entries, is.null)]
    # Process each entry and store results
    if (length(entries) > 0) {
      results <- purrr::map_dfr(entries, function(entry) {
        if (is.null(entry)) {
          return(NULL)
        } # Handle missing entries gracefully

        gene_num <- entry$ENTRY # Gene ID
        org_id <- names(entry$ORGANISM) # Organism ID
        gene_id <- paste0(org_id, ":", gene_num) # Full gene ID
        gene_description = entry$NAME

        # Extract KO ID safely, if more than one ko they will be comma separated
        ko <- if (!is.null(entry$ORTHOLOGY)) paste(names(entry$ORTHOLOGY),collapse = ";") else NA_character_

        ko_description <- if (!is.null(entry$ORTHOLOGY)) paste(entry$ORTHOLOGY,collapse = ";") else NA_character_

        kegg_pathway <- if (!is.null(entry$PATHWAY)) {
          entry$PATHWAY
        } else {
          NA_character_
        }
        kegg_pathway_id <- if (!is.null(entry$PATHWAY)) {
         names(entry$PATHWAY)
        } else {
          NA_character_
        }

        # Return extracted info as a tibble
        tibble::tibble(
          kegg_id = gene_id,
          gene_description = gene_description,
          kegg_pathway_id = kegg_pathway_id,
          kegg_pathway = kegg_pathway,
          ko = ko,
          ko_description = ko_description,
          org_id = org_id
          # brite_info = brite_info
        )
      })

      all_results[[i]] <- results
    }

    # Avoid hitting KEGG rate limits (3x per second)
    # 0.34 led to IP getting blocked.
    Sys.sleep(1)
  }
  # Combine all batch results into a single dataframe
  out <- dplyr::bind_rows(all_results) |>
    dplyr::select(kegg_id,gene_description,kegg_pathway_id,kegg_pathway, ko, ko_description,org_id)

  return(out)
}
