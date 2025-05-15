#' find_possible_contrast_terms
#' Helper function to assist the user in finding contrast terms that are
#' possible given the specified formula.
#'
#' @param qf
#' @param assay_name
#' @param formula
#'
#' @returns
#' @export
#'
#' @examples
find_possible_contrast_terms <- function(qf,
                                         assay_name,
                                         formula) {
  # Extract assay and metadata
  se <- qf[[assay_name]]
  exprs <- SummarizedExperiment::assay(se)

  # Ensure all character/logical columns are factors
  SummarizedExperiment::colData(se)[] <- lapply(
    SummarizedExperiment::colData(se),
    function(x) if (is.character(x) || is.logical(x)) as.factor(x) else x
  )

  colData <- SummarizedExperiment::colData(se)

  # Build design matrix from the provided formula
  design <- model.matrix(formula, data = colData)

  contrast_terms = colnames(design)

  # Return the terms that can be used to make a formula.
  return(contrast_terms)
}
