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
                                uniprot_annotation,
                                column_name = go,
                                regex = "GO:\\d{7}"
                                ){
  # 1. Prepare the GO annotation data
  #    - Filter for relevant columns
  #    - Handle multiple GO terms by splitting the string into a list
  annotation_list_df <- uniprot_annotation |>
    dplyr::select(protein_id, {{ column_name }}) |>
    dplyr::mutate({{ column_name }} := stringr::str_extract_all({{ column_name }}, regex))

  # There is an annotation problem, uniprot annotates protein ids, we have
  # protein groups in our data. Need to map the two. Here we are making a
  # mapping between Protein Group and Protein ID, and counting how many Protein
  # IDs are in a Protein Group.
  protein_count <- SummarizedExperiment::rowData(qf[["protein_groups"]])[,"Protein.Group", drop = FALSE] |>
    tibble::as_tibble() |>
    dplyr::mutate(count = stringr:::str_count(Protein.Group,pattern =";") + 1) |>
    dplyr::group_by(Protein.Group,Protein.Group.Full = Protein.Group) |>
    tidyr::separate_rows(Protein.Group,sep = ";") |>
    dplyr::select(Protein.Group = Protein.Group.Full,protein_id = Protein.Group,count)

  # Dealing with proteins belonging to protein groups that have different GO terms
  # For the ones that have different matches
  combined = annotation_list_df |>
    dplyr::left_join(protein_count,annotation_list_df, by = "protein_id") |>
    dplyr::group_by(Protein.Group,count, {{ column_name }}) |>
    # Counting how many different annotation terms there are for a protein group
    dplyr::summarise(group_count = dplyr::n(),.groups = "drop") |>
    dplyr::mutate(difference = count - group_count) |>
    dplyr::group_by(Protein.Group) |>
    dplyr::filter(difference == min(difference)) |>
    # In event of tie, arbitrarily taking the first one
    dplyr::slice(1) |>
    dplyr::select(Protein.Group,{{ column_name }})

  # Replacing NA with No
  combined <- combined |>
    dplyr::mutate({{ column_name }} := purrr::map({{ column_name }}, ~ {
      if (is.null(.x) || all(is.na(.x))) {
        "Noterm"
      } else {
        .x
      }
    }))

  # # Checking rows, they should all match
  # `%!in%` = Negate(`%in%`)
  # length(combined$Protein.Group)
  # go_t = table(combined$Protein.Group) |> as.data.frame()
  # length(rowData(qf[["protein_groups"]])$Protein.Group)
  # go_t = table(combined$Protein.Group) |> as.data.frame()

  # 2. Add the GO list-column to the QFeatures object's rowData
  #    - First, reorder the annotation data to match the protein assay
  protein_identifiers <- rownames(qf[["protein_groups"]]) # Assumes "proteins" is the protein assay name
  annotation_list_df_ordered <- combined[match(protein_identifiers, combined$Protein.Group), ]


 column_name_string = rlang::as_string(rlang::ensym(column_name))

  #    - Add the list-column to the QFeatures object.
  #      Note: If you have protein groups, your protein_id might not be unique.
  #      This code assumes unique IDs per row.
  # Replaced this code with tidy eval
  #rowData(qf[["protein_groups"]])$go_terms <- annotation_list_df_ordered$go_id
 SummarizedExperiment::rowData(qf[["protein_groups"]])[[column_name_string]] <- annotation_list_df_ordered[[column_name_string]]


  # 3. Create the many-to-many adjacency matrix
  #    - This handles the case where one protein maps to multiple GO terms.
  protein_annotation_list <- rowData(qf[["protein_groups"]])[[column_name_string]]
  all_proteins <- rownames(qf[["protein_groups"]])
  all_annotation_terms <- unique(unlist(protein_annotation_list))

  #    - Create a data frame in long format from the list-column
  protein_annotation_df <- data.frame(
    protein = rep(all_proteins, lengths(protein_annotation_list)),
    annotation_term = unlist(protein_annotation_list)
  )

  #    - Create the sparse adjacency matrix
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
  is_adj <- sapply(rowData(qf[["protein_groups"]]), function(x) is.matrix(x) || inherits(x, "dgCMatrix"))

  # Keep only non-adjacency columns
  rowData(qf[["protein_groups"]]) <- rowData(qf[["protein_groups"]])[ , !is_adj, drop = FALSE]

  adjacencyMatrix(qf[["protein_groups"]]) <- adjacency_matrix
  # Define the name of the adjacency matrix column
  adjacency_matrix_name = paste0(column_name_string,"_adjacency_matrix")
  # Rename the adjacency matrix column
  colnames(rowData(qf[["protein_groups"]]))[colnames(rowData(qf[["protein_groups"]])) == "adjacencyMatrix"] <- adjacency_matrix_name

  # 4. Aggregate the QFeatures object using the adjacency matrix
  qf <- aggregateFeatures(
    qf,
    i = "protein_groups",
    name = rlang::as_string(rlang::ensym(column_name)),
    fcol = adjacency_matrix_name,
    fun = MsCoreUtils::colSumsMat,
    na.rm = TRUE)

  return(qf)

}

