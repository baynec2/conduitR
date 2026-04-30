#' Add Taxonomy / LCA Annotations to a QFeatures Object
#'
#' This function assigns a Lowest Common Ancestor (LCA) taxonomy to each
#' peptide in a QFeatures object using the peptide's `Protein.Ids` (the full
#' set of proteins matching the peptide, as reported by DIA-NN), and then
#' derives a per-protein-group label from the peptide LCAs. Per-rank
#' abundance assays (`domain` ... `species`) are aggregated from the
#' peptides assay so that intensity attributes to each rank only at the
#' resolution the peptide evidence supports.
#'
#' @param qf A QFeatures object containing assays named `peptides` (with a
#'   `Protein.Ids` column in rowData) and `protein_groups` (with rownames
#'   matching `Protein.Group`).
#' @param uniprot_annotation A data frame or tibble containing taxonomy
#'   information. Must include columns: `protein_id`, `domain`, `kingdom`,
#'   `phylum`, `class`, `order`, `family`, `genus`, `species`.
#'
#' @return A QFeatures object with:
#'   \itemize{
#'     \item Per-rank taxonomy and `lca` columns added to `peptides` rowData.
#'     \item Per-rank taxonomy and `lca` columns added to `protein_groups`
#'           rowData (derived from the constituent peptides' LCAs).
#'     \item New assays for each rank (`domain` ... `species`) aggregated
#'           from the peptides assay.
#'   }
#'
#' @export
add_taxonomy_to_qf = function(qf,
                              uniprot_annotation){

  ranks <- c("domain","kingdom","phylum","class","order","family","genus","species")

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

  # Per-peptide LCA: split Protein.Ids, join taxonomy, collapse to consensus
  # at each rank (NA where peptide's matching proteins disagree at that rank).
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

  # Order to match peptide assay rownames and write
  peptide_lca_ordered <- peptide_lca[match(rownames(pep_rd), peptide_lca$Stripped.Sequence), ]
  peptide_lca_ordered <- peptide_lca_ordered |> dplyr::select(-Stripped.Sequence)
  SummarizedExperiment::rowData(qf[["peptides"]]) <- cbind(pep_rd, peptide_lca_ordered)

  # Derive per-protein-group LCA from the constituent peptides' rank values.
  # Same "all-equal-or-NA" collapse, grouped by Protein.Group.
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
      dplyr::across(dplyr::all_of(ranks),
                    ~ if (length(unique(stats::na.omit(.))) == 1) unique(stats::na.omit(.)) else NA_character_),
      .groups = "drop"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(lca = dplyr::coalesce(species, genus, family, order, class, phylum, kingdom, domain)) |>
    dplyr::ungroup()

  pg_rd <- SummarizedExperiment::rowData(qf[["protein_groups"]])
  group_lca_ordered <- group_lca[match(rownames(pg_rd), group_lca$Protein.Group), ]
  group_lca_ordered <- group_lca_ordered |> dplyr::select(-Protein.Group)
  SummarizedExperiment::rowData(qf[["protein_groups"]]) <- cbind(pg_rd, group_lca_ordered)

  # Build per-rank abundance assays by aggregating peptides by per-peptide LCA.
  # Peptides with NA at a given rank drop out for that rank's assay; column
  # sums shrink at finer ranks, which is the honest accounting.
  for (i in ranks){
    rank_vals <- SummarizedExperiment::rowData(qf[["peptides"]])[[i]]
    if (all(is.na(rank_vals))) {
      log_with_timestamp("Skipping aggregation for rank '%s' — all values are NA", i)
      next
    }
    qf <- QFeatures::aggregateFeatures(
      qf,
      i = "peptides",
      fcol = i,
      name = i,
      fun = colSums,
      na.rm = TRUE
    )
  }
  return(qf)
}
