#' Create an UpSet Plot Showing Feature Overlap Between Groups
#'
#' Generates an UpSet plot to visualize the overlap of detected features
#' (e.g., proteins, peptides) across experimental groups. UpSet plots scale
#' cleanly beyond three groups where Venn diagrams become unreadable.
#'
#' A feature is considered detected in a group if it has a non-\code{NA} value
#' in at least one sample belonging to that group.
#'
#' @param qf A QFeatures object containing the data.
#' @param assay_name Character string specifying which assay to use.
#' @param group_by Character string specifying a column in \code{colData(qf)}
#'   to use for grouping samples.
#'
#' @return A \code{ComplexHeatmap} UpSet plot object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_upset(qf, "protein_groups_log2_MinDet_none", group_by = "treatment")
#' }
plot_upset <- function(qf, assay_name, group_by) {

  se <- qf[[assay_name]]
  mat <- SummarizedExperiment::assay(se)
  col_data <- SummarizedExperiment::colData(se)

  groups <- levels(factor(col_data[[group_by]]))

  sets <- lapply(groups, function(grp) {
    samples_in_grp <- rownames(col_data)[col_data[[group_by]] == grp]
    sub_mat <- mat[, samples_in_grp, drop = FALSE]
    rownames(sub_mat)[rowSums(!is.na(sub_mat)) > 0]
  })
  names(sets) <- groups

  comb_mat <- ComplexHeatmap::make_comb_mat(sets)
  ComplexHeatmap::UpSet(comb_mat)
}
