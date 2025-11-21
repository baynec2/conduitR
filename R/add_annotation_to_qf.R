#'Add Annotations to a QFeatures Object

#' This function adds annotation terms (e.g., GO, Pfam, InterPro) to a QFeatures
#' object.It extracts terms from a UniProt annotation file, stores them as a
#' list-column in the protein assay, and creates an adjacency matrix for
#' aggregation while keeping assay links intact.
#'
#' @param qf A QFeatures object containing a protein assay named "protein_groups".
#' @param uniprot_annotation A data frame or tibble containing protein annotations,
#'   with a column for protein IDs and a column for the terms of interest.
#' @param column_name Unquoted name of the column in \code{uniprot_annotation} that
#'   contains the annotation terms to extract.
#' @param regex A regular expression pattern used to extract annotation terms into a list.
#'   For example: \code{xref_pfam = "PF\\d{5}"}, \code{xref_interpro = "IPR\\d{6}"}.
#'
#' @return A QFeatures object with:
#'   \itemize{
#'     \item A new list-column in the protein assay containing the extracted annotation terms.
#'     \item An adjacency matrix stored in \code{adjacencyMatrix()} for aggregation.
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add GO annotations
#' qf <- add_annotation_to_qf(
#'   qf,
#'   uniprot_annotation = uniprot_df,
#'   column_name = go,
#'   regex = "GO:\\d{7}"
#' )
#'
#' # Add Pfam annotations
#' qf <- add_annotation_to_qf(
#'   qf,
#'   uniprot_annotation = uniprot_df,
#'   column_name = xref_pfam,
#'   regex = "PF\\d{5}"
#' )
#' }
add_annotation_to_qf = function(qf,
                                id_column = Protein.Group,
                                conduit_annotations,
                                column_name = go,
                                regex = "[^;]+(?=;|$)"
                                ){
  # Prepare the annotation data for selected column
  annotation_list_df <- conduit_annotations |>
    dplyr::select("Protein.Group" = {{ id_column }}, {{ column_name }}) |>
    dplyr::mutate({{ column_name }} := stringr::str_extract_all({{ column_name }}, regex))

  # Extracting protein groups from QF.
  pg_df <- SummarizedExperiment::rowData(qf[["protein_groups"]])[,"Protein.Group", drop = FALSE] |>
    tibble::as_tibble()

  # Combining the data frames for the ones that have different matches
  combined = dplyr::left_join(pg_df,
                              annotation_list_df,
                              by = "Protein.Group")

  # Replacing NA with Noterm
  combined <- combined |>
    dplyr::mutate({{ column_name }} := purrr::map({{ column_name }}, ~ {
      if (is.null(.x) || all(is.na(.x))) {
        "Noterm"
      } else {
        .x
      }
    }))

  # Add the list-column to the QFeatures object's rowData
  #    - First, reorder the annotation data to match the protein assay
  protein_identifiers <- rownames(qf[["protein_groups"]])
  annotation_list_df_ordered <- combined[match(protein_identifiers,
                                               combined$Protein.Group), ]
 # Extracting the column name as the string
 column_name_string = rlang::as_string(rlang::ensym(column_name))

 # Adding the specified annotation data to the rowData.
 SummarizedExperiment::rowData(qf[["protein_groups"]])[[column_name_string]] <-
   annotation_list_df_ordered[[column_name_string]]

  # Creating an adjacency matrix
  protein_annotation_list <- SummarizedExperiment::rowData(qf[["protein_groups"]])[[column_name_string]]
  all_proteins <- rownames(qf[["protein_groups"]])
  all_annotation_terms <- unique(unlist(protein_annotation_list))

  # Create a data frame in long format from the list-column
  protein_annotation_df <- data.frame(
    protein = rep(all_proteins, lengths(protein_annotation_list)),
    annotation_term = unlist(protein_annotation_list)
  )

  # Create a sparse adjacency matrix
  adjacency_matrix <- with(protein_annotation_df, {
    Matrix::sparseMatrix(
      i = match(protein, all_proteins),
      j = match(annotation_term, all_annotation_terms),
      x = 1,
      dims = c(length(all_proteins), length(all_annotation_terms)),
      dimnames = list(all_proteins, all_annotation_terms)
    )
  })
  # There is a problem with multiple adjacency matrices being stored in rowData.
  # Identify adjacency matrix columns in rowData (list-columns of matrices)
  is_adj <- sapply(SummarizedExperiment::rowData(qf[["protein_groups"]]), function(x) is.matrix(x) || inherits(x, "dgCMatrix"))

  # Keep only non-adjacency columns
  SummarizedExperiment::rowData(qf[["protein_groups"]]) <- SummarizedExperiment::rowData(qf[["protein_groups"]])[ , !is_adj, drop = FALSE]

  QFeatures::adjacencyMatrix(qf[["protein_groups"]]) <- adjacency_matrix
  # Define the name of the adjacency matrix column
  adjacency_matrix_name = paste0(column_name_string,"_adjacency_matrix")
  # Rename the adjacency matrix column
  colnames(SummarizedExperiment::rowData(qf[["protein_groups"]]))[colnames(SummarizedExperiment::rowData(qf[["protein_groups"]])) == "adjacencyMatrix"] <- adjacency_matrix_name

  # Aggregate the QFeatures object using the adjacency matrix
  qf <- QFeatures::aggregateFeatures(
    qf,
    i = "protein_groups",
    name = rlang::as_string(rlang::ensym(column_name)),
    fcol = adjacency_matrix_name,
    fun = MsCoreUtils::colSumsMat,
    na.rm = TRUE)

  # Define adjacency columns after aggregation / renaming
  adj_cols <- grepl("_adjacency_matrix$",
                    colnames(SummarizedExperiment::rowData(
                      qf[["protein_groups"]])
                      )
                    )

  # Keep only non-adjacency columns
  SummarizedExperiment::rowData(qf[["protein_groups"]]) <-
    SummarizedExperiment::rowData(qf[["protein_groups"]])[ ,
                                                           !adj_cols
                                                           , drop = FALSE]

  return(qf)
}

