#' get_annotations_from_uniprot
#'
#' given the pg_matrix output of diann, assign go terms to each protein group
#' and summarise.
#'
#' @param pg_matrix_fp pg_matrix output of diann
#' @param outdir out directory
#'
#' @returns tsv at output directory location
#' @export
#'
#' @examples
#' insulin <- get_annotations_from_uniprot("P01308")
#' t <- get_annotations_from_uniprot("P47340")
get_annotations_from_uniprot <- function(uniprot_ids,
                                         columns = "accession,go,xref_kegg,xref_eggnog,cc_subcellular_location,cc_tissue_specificity,xref_cazy,xref_pfam,xref_interpro"
 ) {

  annotated_data <- annotate_uniprot_ids(uniprot_ids,
    columns = columns,
    batch_size = 150
  ) |>
    # Cleaning up go and Kegg and eggnog data.
    dplyr::mutate(
      go = dplyr::na_if(go, " "),
      go = dplyr::na_if(go, ""),
      xref_kegg = dplyr::na_if(xref_kegg, " "),
      xref_kegg = dplyr::na_if(xref_kegg, ""),
      xref_eggnog = dplyr::na_if(xref_eggnog, " "),
      xref_eggnog = dplyr::na_if(xref_eggnog, "")
    )

  return(annotated_data)
}
