#' replace 0 with NA.
#'
#' @param qf qfeatures object
#'
#' @returns
#' @export
#'
#' @examples
replace_zero_with_na <- function(qf) {
  for (assay_name in names(qf)) {
    expr_mat <- SummarizedExperiment::assay(qf[[assay_name]])  # Extract assay data
    expr_mat[expr_mat == 0] <- NA  # Replace 0s with NA
    SummarizedExperiment::assay(qf[[assay_name]]) <- expr_mat  # Update the assay
  }
  return(qf)
}
