#' Add Taxonomy / LCA Annotations to a QFeatures Object
#'
#' Assigns a Lowest Common Ancestor (LCA) taxonomy to each peptide using the
#' peptide's `Protein.Ids` (the full set of proteins matching the peptide as
#' reported by DIA-NN), then derives a per-protein-group label from the
#' peptide LCAs. Per-rank taxonomy columns are written to both the
#' `peptides` and `protein_groups` rowData, and the eight ranks are
#' registered as aggregation targets — call
#' [aggregate_assay_by_annotation()] to materialise a per-rank assay on
#' demand.
#'
#' The protein-group label is a *strict* LCA of its constituent peptides: a
#' group is resolved at a rank only if every peptide agrees there and none is
#' ambiguous (`NA`) at that rank. An ambiguous peptide therefore pulls the
#' group label up to the most specific rank the whole group shares, rather
#' than abstaining — so a small resolved minority cannot promote a
#' mostly-ambiguous group to a confident (e.g. species) call. The
#' `protein_groups` rowData additionally records `n_peptides`,
#' `n_species_resolved`, and `species_resolved_fraction` so downstream code
#' can flag or filter low-confidence group calls.
#'
#' Peptide-centric summation is the canonical aggregation path for taxonomy
#' (intensity attributes to a rank only at the resolution the peptide
#' evidence supports), so each rank is registered with `from = "peptides"`.
#'
#' @param qf A QFeatures object containing assays named `peptides` (with a
#'   `Protein.Ids` column in rowData) and `protein_groups` (with rownames
#'   matching `Protein.Group`).
#' @param uniprot_annotation A data frame with columns `protein_id`,
#'   `domain`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`,
#'   `species`.
#'
#' @return The input QFeatures object with:
#'   \itemize{
#'     \item Per-rank taxonomy and `lca` columns added to `peptides` rowData.
#'     \item Per-rank taxonomy and `lca` columns (strict LCA) added to
#'           `protein_groups` rowData, plus `n_peptides`,
#'           `n_species_resolved`, and `species_resolved_fraction`.
#'     \item Each rank registered as a taxonomic aggregation target sourced
#'           from `peptides`.
#'   }
#'
#' @seealso [aggregate_assay_by_annotation()], [aggregation_targets()]
#' @export
add_taxonomy_to_qf <- function(qf,
                               uniprot_annotation) {

  ranks <- c("domain", "kingdom", "phylum", "class",
             "order", "family", "genus", "species")

  if (!"peptides" %in% names(qf)) {
    stop("qf must contain a 'peptides' assay")
  }
  if (!"protein_groups" %in% names(qf)) {
    stop("qf must contain a 'protein_groups' assay")
  }
  pep_rd <- SummarizedExperiment::rowData(qf[["peptides"]])
  if (!"Protein.Ids" %in% colnames(pep_rd)) {
    stop("peptides rowData must contain a 'Protein.Ids' column; see diann_to_qfeatures()")
  }
  if (!"Protein.Group" %in% colnames(pep_rd)) {
    stop("peptides rowData must contain a 'Protein.Group' column; see diann_to_qfeatures()")
  }

  taxonomy <- uniprot_annotation |>
    dplyr::select(protein_id, dplyr::all_of(ranks))

  peptide_tax <- tibble::tibble(
      Stripped.Sequence = rownames(pep_rd),
      Protein.Ids = pep_rd$Protein.Ids
    ) |>
    tidyr::separate_rows(Protein.Ids, sep = ";") |>
    dplyr::rename(protein_id = Protein.Ids) |>
    dplyr::left_join(taxonomy, by = "protein_id")

  peptide_lca <- peptide_tax |>
    dplyr::group_by(Stripped.Sequence) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(ranks),
                    ~ if (length(unique(.)) == 1) unique(.) else NA_character_),
      .groups = "drop"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(lca = dplyr::coalesce(species, genus, family, order, class, phylum, kingdom, domain)) |>
    dplyr::ungroup()

  peptide_lca_ordered <- peptide_lca[match(rownames(pep_rd), peptide_lca$Stripped.Sequence), ]
  peptide_lca_ordered <- peptide_lca_ordered |> dplyr::select(-Stripped.Sequence)
  SummarizedExperiment::rowData(qf[["peptides"]]) <- cbind(pep_rd, peptide_lca_ordered)

  group_lca <- tibble::tibble(
      Protein.Group = pep_rd$Protein.Group,
      domain = peptide_lca_ordered$domain,
      kingdom = peptide_lca_ordered$kingdom,
      phylum = peptide_lca_ordered$phylum,
      class = peptide_lca_ordered$class,
      order = peptide_lca_ordered$order,
      family = peptide_lca_ordered$family,
      genus = peptide_lca_ordered$genus,
      species = peptide_lca_ordered$species
    ) |>
    dplyr::group_by(Protein.Group) |>
    dplyr::summarise(
      # Diagnostics (computed before the across() so they bind to the raw
      # per-peptide columns): group size and how many peptides carry a
      # resolved species-level LCA. species_resolved_fraction lets downstream
      # code flag low-confidence group calls.
      n_peptides = dplyr::n(),
      n_species_resolved = sum(!is.na(species)),
      # Strict LCA: a group is resolved at a rank only if *every* constituent
      # peptide agrees there and none abstain. Ambiguous (NA) peptides count as
      # dissent (unique() keeps NA -> length > 1 -> NA), so a small resolved
      # minority can no longer promote a mostly-ambiguous group; the label
      # falls back up to the most specific rank the whole group shares. This
      # matches the peptide-level rule above. The previous na.omit() dropped
      # abstainers and over-promoted the group (see issue #16).
      dplyr::across(dplyr::all_of(ranks),
                    ~ if (length(unique(.)) == 1) unique(.) else NA_character_),
      .groups = "drop"
    ) |>
    dplyr::mutate(species_resolved_fraction = n_species_resolved / n_peptides) |>
    dplyr::rowwise() |>
    dplyr::mutate(lca = dplyr::coalesce(species, genus, family, order, class, phylum, kingdom, domain)) |>
    dplyr::ungroup()

  pg_rd <- SummarizedExperiment::rowData(qf[["protein_groups"]])
  group_lca_ordered <- group_lca[match(rownames(pg_rd), group_lca$Protein.Group), ]
  group_lca_ordered <- group_lca_ordered |> dplyr::select(-Protein.Group)
  SummarizedExperiment::rowData(qf[["protein_groups"]]) <- cbind(pg_rd, group_lca_ordered)

  for (r in ranks) {
    qf <- register_aggregation_target(qf, target = r,
                                      kind = "taxonomic",
                                      from = "peptides")
  }
  qf
}
