#' Replace Zero Values with NA in QFeatures Object
#'
#' Replaces all zero values with NA in all assays of a QFeatures object.
#' This is typically done as a preprocessing step before log transformation
#' and imputation, as zeros in proteomics data often represent missing values
#' rather than true zero abundances.
#'
#' @param qf A QFeatures object containing the data to process
#'
#' @return A QFeatures object with all zero values replaced by NA in all assays.
#'   The structure and metadata of the object remain unchanged.
#'
#' @export
#'
#' @examples
#' # Replace zeros with NA in all assays:
#' # qf_no_zeros <- replace_zero_with_na(qfeatures_obj)
#' 
#' # This is often used as a preprocessing step:
#' # qf_processed <- qfeatures_obj |>
#' #   replace_zero_with_na() |>
#' #   add_log_imputed_norm_assays()
#' 
#' # The resulting object can be used with plotting functions:
#' # plot_missing_val_heatmap(qf_no_zeros, "protein")
replace_zero_with_na <- function(qf) {
  for (assay_name in names(qf)) {
    expr_mat <- SummarizedExperiment::assay(qf[[assay_name]])  # Extract assay data
    expr_mat[expr_mat == 0] <- NA  # Replace 0s with NA
    SummarizedExperiment::assay(qf[[assay_name]]) <- expr_mat  # Update the assay
  }
  return(qf)
}
