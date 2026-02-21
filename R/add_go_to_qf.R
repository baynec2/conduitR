#' Add Gene Ontology Terms to a QFeatures Object
#'
#' Maps Gene Ontology (GO) terms from a UniProt annotation table onto the
#' protein group assay, adds a GO list-column to rowData, builds a
#' protein–GO adjacency matrix, and aggregates the assay by GO terms to create
#' a new assay `go_terms`. Assumes the protein assay is named `"protein_groups"`.
#'
#' @param qf A `QFeatures` object with an assay named `"protein_groups"` and
#'   rowData containing `Protein.Group`.
#' @param uniprot_annotation A data frame or tibble with at least columns
#'   `protein_id` and `go` (GO IDs, e.g. from `get_annotations_from_uniprot`).
#'
#' @return The same `QFeatures` object with GO terms in rowData of
#'   `protein_groups` and a new assay `go_terms` (aggregated by GO).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ann <- get_annotations_from_uniprot(rownames(qf[["protein_groups"]]))
#' qf <- add_go_to_qf(qf, ann)
#' }
add_go_to_qf = function(qf,
                        uniprot_annotation){
  # 1. Prepare the GO annotation data
  #    - Filter for relevant columns
  #    - Handle multiple GO terms by splitting the string into a list
  go_list_df <- uniprot_annotation |>
    dplyr::select(protein_id, go) |>
    dplyr::mutate(go_id = stringr::str_extract_all(go, "GO:\\d{7}")) |>
    dplyr::select(-go)

  # There is an annotation problem, uniprot annotates protein ids, we have
  # protein groups. Need to map the two.
  protein_count <- rowData(qf[["protein_groups"]]) |>
    as.data.frame() |>
    dplyr::mutate(count = stringr:::str_count(Protein.Group,pattern =";") + 1) |>
    dplyr::group_by(Protein.Group,Protein.Group.Full = Protein.Group) |>
    tidyr::separate_rows(Protein.Group,sep = ";") |>
    dplyr::select(Protein.Group = Protein.Group.Full,protein_id = Protein.Group,count)


  # Filtering to only include the terms that match between all proteins in group
  go_combined = go_list_df |>
    dplyr::left_join(protein_count,go, by = "protein_id") |>
    dplyr::group_by(Protein.Group,count,go_id) |>
    dplyr::summarise(group_count = dplyr::n(),.groups = "drop") |>
    dplyr::mutate(difference = count - group_count) |>
    #This would only contain the ones that have matches for all proteins in PG
    dplyr::group_by(Protein.Group) |>
    dplyr::filter(difference == min(difference)) |>
    dplyr::select(Protein.Group,go_id)

  # Won't work with NAs. Getting rid of here.
  go_combined$go_id[is.na(go_combined$go_id)] <- "NoGO"

  # # Checking rows, they should all match
  # `%!in%` = Negate(`%in%`)
  # length(go_combined$Protein.Group)
  # go_t = table(go_combined$Protein.Group) |> as.data.frame()
  # length(rowData(qf[["protein_groups"]])$Protein.Group)
  # go_t = table(go_combined$Protein.Group) |> as.data.frame()

  # 2. Add the GO list-column to the QFeatures object's rowData
  #    - First, reorder the annotation data to match the protein assay
  protein_identifiers <- rownames(qf[["protein_groups"]]) # Assumes "proteins" is the protein assay name
  go_list_df_ordered <- go_combined[match(protein_identifiers, go_combined$Protein.Group), ]

  #    - Add the list-column to the QFeatures object.
  #      Note: If you have protein groups, your protein_id might not be unique.
  #      This code assumes unique IDs per row.
  rowData(qf[["protein_groups"]])$go_terms <- go_list_df_ordered$go_id

  # 3. Create the many-to-many adjacency matrix
  #    - This handles the case where one protein maps to multiple GO terms.
  protein_go_list <- rowData(qf[["protein_groups"]])$go_terms
  all_proteins <- rownames(qf[["protein_groups"]])
  all_go_terms <- unique(unlist(protein_go_list))

  #    - Create a data frame in long format from the list-column
  protein_go_df <- data.frame(
    protein = rep(all_proteins, lengths(protein_go_list)),
    go_term = unlist(protein_go_list)
  )

  #    - Create the sparse adjacency matrix
  go_adjacency_matrix <- with(protein_go_df, {
    Matrix::sparseMatrix(
      i = match(protein, all_proteins),
      j = match(go_term, all_go_terms),
      x = 1,
      dims = c(length(all_proteins), length(all_go_terms)),
      dimnames = list(all_proteins, all_go_terms)
    )
  })


  adjacencyMatrix(qf[["protein_groups"]]) <- go_adjacency_matrix

  # 4. Aggregate the QFeatures object using the adjacency matrix
  qf <- aggregateFeatures(
    qf,
    i = "protein_groups",
    name = "go_terms",
    fcol = "adjacencyMatrix",
    fun = MsCoreUtils::colSumsMat,
    na.rm = TRUE)


return(qf)

}
