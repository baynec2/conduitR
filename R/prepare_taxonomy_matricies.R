#' Build Taxonomy-Level Matrices from DIA-NN and FASTA Data
#'
#' Joins full taxonomy, DIA-NN protein group matrix, and FASTA-derived protein
#' IDs to produce one TSV matrix per taxonomic level (domain, kingdom, phylum,
#' class, order, family, genus, species). Each file contains log2-transformed
#' intensities by sample and taxon.
#'
#' @param full_taxonomy_fp Character. Path to the full taxonomy file (e.g.
#'   organism_id and taxonomic columns).
#' @param pr_group_matrix_fp Character. Path to the DIA-NN protein group matrix
#'   (Protein.Ids and sample intensity columns).
#' @param fasta_fp Character. Path to the FASTA used in the experiment (for
#'   protein_id to organism_id mapping via `extract_fasta_info`).
#' @param output_dir Character. Directory to write `*_matrix.tsv` files
#'   (default: `getwd()`).
#'
#' @return No return value; writes one TSV per taxonomic level named
#'   `{level}_matrix.tsv` in `output_dir`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' prepare_taxonomy_matricies(
#'   full_taxonomy_fp = "taxonomy.tsv",
#'   pr_group_matrix_fp = "pg_matrix.tsv",
#'   fasta_fp = "combined.fasta",
#'   output_dir = "taxonomy_matrices"
#' )
#' }
prepare_taxonomy_matricies <- function(full_taxonomy_fp,
                                       pr_group_matrix_fp,
                                       fasta_fp,
                                       output_dir = getwd()) {
  # Reading in the files
  full_taxonomy <- readr::read_delim(full_taxonomy_fp)
  peptide_groups <- readr::read_tsv(pr_group_matrix_fp)
  fasta_info <- extract_fasta_info(fasta_fp)

  # Joining data
  peptide_and_org <- dplyr::inner_join(peptide_groups, fasta_info,
    by = c("Protein.Ids" = "protein_id")
  )

  peptide_and_taxa <- dplyr::inner_join(peptide_and_org, full_taxonomy,
    by = c("organism_id" = "tax_id")
  )


  taxa_units <- c(
    "domain", "kingdom", "phylum", "class", "order", "family",
    "genus", "species"
  )

  # Initialize an empty list to store results

  for (i in taxa_units) {
    df <- peptide_and_taxa %>%
      tidyr::pivot_longer(
        cols = contains("raw"), names_to = "colData",
        values_to = "intensity"
      ) %>%
      dplyr::select(i, Protein.Ids, colData, intensity) %>%
      dplyr::group_by(colData, !!dplyr::sym(i)) %>%
      dplyr::summarise(
        intensity = sum(intensity, na.rm = TRUE) + 1,
        .groups = "drop"
      ) %>%
      dplyr::mutate(log2 = log2(intensity)) %>%
      dplyr::select(-intensity) %>%
      tidyr::pivot_wider(names_from = colData, values_from = log2)

    readr::write_tsv(df, paste0(output_dir, "/", i, "_matrix.tsv"))
  }
}
