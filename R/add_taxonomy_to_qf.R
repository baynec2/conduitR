#' Add Taxonomy / LCA Annotations to a QFeatures Object
#'
#' This function assigns a Lowest Common Ancestor (LCA) taxonomy to each
#' protein group in a QFeatures object. It uses the taxonomy of each protein
#' within a group and collapses them to the most specific shared taxonomic level
#' (domain → species). The full taxonomy columns and the computed LCA are added
#' to the protein assay's rowData.
#'
#' @param qf A QFeatures object containing a protein assay named "protein_groups".
#' @param uniprot_annotation A data frame or tibble containing taxonomy information.
#'   Must include columns: \code{domain}, \code{kingdom}, \code{phylum}, \code{class},
#'   \code{order}, \code{family}, \code{genus}, and \code{species}.
#'
#' @return A QFeatures object with:
#'   \itemize{
#'     \item Taxonomy columns added to the protein assay rowData.
#'     \item A new \code{lca} column representing the lowest common ancestor for each protein group.
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' qf <- add_lca_to_qf(
#'   qf,
#'   uniprot_annotation = uniprot_taxonomy_df
#' )
#' # Access LCA column
#' rowData(qf[["protein_groups"]])$lca
#' }
add_taxonomy_to_qf = function(qf,
                              uniprot_annotation){

  # Load the Annotation data
  taxonomy <- uniprot_annotation |>
    dplyr::select(protein_id,organism_id,organism_type,domain,kingdom,phylum,class,
           order,family,genus,species)

  protein_count <- SummarizedExperiment::rowData(qf[["protein_groups"]])[,"Protein.Group", drop = FALSE] |>
    tibble::as_tibble() |>
    dplyr::mutate(count = stringr:::str_count(Protein.Group,pattern =";") + 1) |>
    dplyr::group_by(Protein.Group,Protein.Group.Full = Protein.Group) |>
    tidyr::separate_rows(Protein.Group,sep = ";") |>
    dplyr::select(Protein.Group = Protein.Group.Full,protein_id = Protein.Group,count)

  # Merge taxonomy info for each protein in each group
  protein_tax <- protein_count |>
    dplyr::left_join(taxonomy, by = c("protein_id" = "protein_id"))

  # Defining a function to get lca
  get_lca <- function(tax_df) {
    # tax_df: one protein group with multiple rows
    lca <- sapply(tax_df, function(col) {
      unique_vals <- unique(col)
      if (length(unique_vals) == 1) return(unique_vals)
      else return(NA)  # stop at the level where there is disagreement
    })
    return(lca)
  }

  lca_per_group <- protein_tax |>
    dplyr::group_by(Protein.Group) |>
    dplyr::summarise(dplyr::across(domain:species, ~ if(length(unique(.)) == 1) unique(.) else NA),
              .groups = "drop") |>
    dplyr::rowwise() |>
    dplyr::mutate(lca = dplyr::coalesce(species, genus, family, order, class, phylum, kingdom, domain)) |>
    dplyr::ungroup()

  # Make sure the rows match the assay rownames
  lca_ordered <- lca_per_group[match(rownames(qf[["protein_groups"]]), lca_per_group$Protein.Group), ]

  # Drop the Protein.Group column if you like (optional)
  lca_ordered <- lca_ordered |>  dplyr::select(-Protein.Group)

  # Assign all columns to rowData
  rowData(qf[["protein_groups"]]) <- cbind(rowData(qf[["protein_groups"]]), lca_ordered)

  # Define the taxonomic ranks
  ranks <- c("domain","kingdom","phylum","class","order","family","genus","species")

  for(i in ranks){
    qf <- QFeatures::aggregateFeatures(
      qf,
      i = "protein_groups",
      fcol = i,   # grouping column in rowData
      name = i,
      fun = colSums,
      na.rm = TRUE
    )
  }
  return(qf)
}
