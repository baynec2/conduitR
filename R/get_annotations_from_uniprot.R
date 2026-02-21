#' Get Annotations from UniProt for a Set of Accession IDs
#'
#' Fetches UniProt annotations (GO, KEGG, EggNOG, subcellular location,
#' tissue specificity, CAZy, Pfam, InterPro) for the given IDs and returns a
#' tibble with empty/missing values cleaned (NA for blank or single space).
#'
#' @param uniprot_ids Character vector of UniProt accession IDs.
#' @param columns Character. Comma-separated list of UniProt fields to request
#'   (default includes accession, protein_name, go, xref_kegg, xref_eggnog,
#'   cc_subcellular_location, cc_tissue_specificity, xref_cazy, xref_pfam,
#'   xref_interpro).
#'
#' @return A tibble with one row per protein and columns corresponding to the
#'   requested UniProt fields.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ann <- get_annotations_from_uniprot(c("P01308", "P47340"))
#' get_annotations_from_uniprot("P01308", columns = "accession,protein_name,go")
#' }
get_annotations_from_uniprot <- function(uniprot_ids,
                                         columns = "accession,protein_name,go,xref_kegg,xref_eggnog,cc_subcellular_location,cc_tissue_specificity,xref_cazy,xref_pfam,xref_interpro"
 ) {

  annotated_data <- annotate_uniprot_ids(uniprot_ids,
    columns = columns,
    batch_size = 100
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
