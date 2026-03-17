#' Add KEGG Gene IDs to QFeatures rowData
#'
#' Stores UniProt KEGG cross-reference IDs (e.g. \code{"hsa:1234"}) as a
#' list-column in the \code{protein_groups} rowData. Unlike
#' \code{\link{add_annotation_to_qf}}, this function does \emph{not} create an
#' aggregated assay — it only adds the gene-ID list-column so that
#' \code{\link{plot_kegg_pathway}} can join expression data onto pathway nodes.
#'
#' Call this function \strong{before} \code{\link{build_conduit_obj}} so that
#' the IDs are available when \code{\link{perform_limma_analysis}} joins rowData
#' into the top-table.
#'
#' @param qf A QFeatures object containing a \code{protein_groups} assay.
#' @param uniprot_annotation A data frame or tibble returned by
#'   \code{\link{get_annotations_from_uniprot}} (or any frame with columns
#'   \code{accession} and \code{xref_kegg}).
#' @param id_column Unquoted name of the protein-ID column in
#'   \code{uniprot_annotation} (default: \code{accession}).
#' @param kegg_col Unquoted name of the KEGG cross-reference column in
#'   \code{uniprot_annotation} (default: \code{xref_kegg}).
#' @param out_col Name of the list-column to create in rowData
#'   (default: \code{"xref_kegg"}).
#'
#' @return The QFeatures object with an \code{xref_kegg} list-column added to
#'   the \code{protein_groups} rowData.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' uniprot_ann <- get_annotations_from_uniprot(protein_ids)
#' qf <- add_kegg_ids_to_qf(qf, uniprot_ann)
#' # xref_kegg is now in rowData(qf[["protein_groups"]])
#' }
add_kegg_ids_to_qf <- function(qf,
                                uniprot_annotation,
                                id_column  = accession,
                                kegg_col   = xref_kegg,
                                out_col    = "xref_kegg") {

  # Build a per-protein-group list of KEGG gene IDs
  pg_df <- SummarizedExperiment::rowData(qf[["protein_groups"]])[, "Protein.Group", drop = FALSE] |>
    tibble::as_tibble()

  # Expand semicolon-delimited Protein.Groups to individual protein IDs
  pg_expanded <- pg_df |>
    dplyr::mutate(protein_id = strsplit(Protein.Group, ";")) |>
    tidyr::unnest(protein_id) |>
    dplyr::mutate(protein_id = trimws(protein_id))

  # Extract the KEGG ID column from the annotation table
  kegg_df <- uniprot_annotation |>
    dplyr::select(
      protein_id = {{ id_column }},
      kegg_raw   = {{ kegg_col }}
    ) |>
    dplyr::filter(!is.na(kegg_raw), kegg_raw != "") |>
    dplyr::mutate(
      kegg_ids = stringr::str_extract_all(kegg_raw, "[^;\\s]+")
    ) |>
    dplyr::select(protein_id, kegg_ids) |>
    tidyr::unnest(kegg_ids)

  # Join back to protein groups and collect into lists
  pg_kegg <- pg_expanded |>
    dplyr::left_join(kegg_df, by = "protein_id") |>
    dplyr::group_by(Protein.Group) |>
    dplyr::summarise(
      !!out_col := list(unique(stats::na.omit(kegg_ids))),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      !!out_col := purrr::map(.data[[out_col]], ~ if (length(.x) == 0) "Noterm" else .x)
    )

  # Re-order to match rowData row order
  ordered <- pg_kegg[match(rownames(qf[["protein_groups"]]), pg_kegg$Protein.Group), ]

  # Assign into rowData
  SummarizedExperiment::rowData(qf[["protein_groups"]])[[out_col]] <- ordered[[out_col]]

  return(qf)
}
