#' Find Possible Contrast Terms for Limma Analysis
#'
#' Helper function that identifies valid contrast terms that can be used in limma analysis
#' based on the provided experimental design formula. This is useful for determining
#' which comparisons are possible in differential expression analysis.
#'
#' @param qf A QFeatures object containing the proteomics data
#' @param assay_name Character string specifying which assay to use for the analysis
#' @param formula A formula object specifying the experimental design (e.g., ~group + batch)
#'
#' @return A character vector containing the names of all possible contrast terms that
#'   can be used in limma analysis. These terms correspond to the columns of the
#'   design matrix that would be created from the formula.
#'
#' @export
#'
#' @examples
#' # Find possible contrasts for a simple group comparison:
#' # terms <- find_possible_contrast_terms(qfeatures_obj, "protein", ~group)
#'
#' # Find possible contrasts for a more complex design:
#' # terms <- find_possible_contrast_terms(qfeatures_obj, "protein", ~group + batch)
#'
#' # The returned terms can be used to construct contrasts for limma analysis:
#' # contrast <- paste(terms[2], "-", terms[1])  # e.g., "groupB - groupA"
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

  # SANITIZE the design matrix column names
  colnames(design) <- make.names(colnames(design))

  contrast_terms = colnames(design)

  # Return the terms that can be used to make a formula.
  return(contrast_terms)
}
