#' Attach an Annotation Column to a QFeatures Object
#'
#' Extracts annotation terms (e.g. GO, Pfam, KEGG) from a long-format
#' annotation table, attaches them as a list-column on the
#' `protein_groups` assay's rowData, and registers the column as a
#' functional aggregation target. Aggregation
#' itself is deferred — call [aggregate_assay_by_annotation()] later when
#' you want an aggregated assay.
#'
#' @param qf A QFeatures object containing an assay named `"protein_groups"`.
#' @param id_column Unquoted name of the column in `conduit_annotations`
#'   that holds protein-group identifiers (default: `Protein.Group`).
#' @param conduit_annotations A data frame with one row per protein-group /
#'   annotation pair (long format), or a wide table with one row per
#'   protein-group and the annotation in a delimited cell.
#' @param column_name Unquoted name of the column in `conduit_annotations`
#'   that holds the annotation terms to extract.
#' @param regex Regular expression used to extract individual terms from
#'   each cell of `column_name`. Defaults to splitting on `;`.
#'
#' @return The input QFeatures object with one extra list-column on
#'   `rowData(qf[["protein_groups"]])` and the column registered as a
#'   functional aggregation target. No new assays are created.
#'
#' @seealso [aggregate_assay_by_annotation()], [aggregation_targets()]
#' @export
add_annotation_to_qf <- function(qf,
                                 id_column = Protein.Group,
                                 conduit_annotations,
                                 column_name = go,
                                 regex = "[^;]+(?=;|$)") {

  annotation_list_df <- conduit_annotations |>
    dplyr::select("Protein.Group" = {{ id_column }}, {{ column_name }}) |>
    dplyr::mutate({{ column_name }} := stringr::str_extract_all(
      {{ column_name }}, regex))

  pg_df <- SummarizedExperiment::rowData(qf[["protein_groups"]])[, "Protein.Group", drop = FALSE] |>
    tibble::as_tibble()

  combined <- dplyr::left_join(pg_df, annotation_list_df, by = "Protein.Group")

  column_name_string <- rlang::as_string(rlang::ensym(column_name))

  protein_identifiers <- rownames(qf[["protein_groups"]])
  annotation_list_df_ordered <- combined[match(protein_identifiers,
                                               combined$Protein.Group), ]

  SummarizedExperiment::rowData(qf[["protein_groups"]])[[column_name_string]] <-
    annotation_list_df_ordered[[column_name_string]]

  qf <- register_aggregation_target(qf,
                                    target = column_name_string,
                                    kind = "functional",
                                    from = "protein_groups")
  qf
}
