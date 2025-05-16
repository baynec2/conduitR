#' Perform Differential Expression Analysis Using Limma
#'
#' Conducts differential expression analysis on proteomics data using the limma package.
#' This function handles the entire limma workflow including design matrix creation,
#' model fitting, contrast testing, and empirical Bayes moderation.
#'
#' @param qf A QFeatures object containing the proteomics data to analyze
#' @param assay_name Character string specifying which assay to use for the analysis
#' @param formula A formula object specifying the experimental design (e.g., ~group + batch)
#' @param contrast Character string specifying the contrast to test (e.g., "groupB - groupA")
#'
#' @return A list containing:
#'   \itemize{
#'     \item top_table: A data frame with differential expression results including:
#'       \itemize{
#'         \item logFC: log fold change
#'         \item AveExpr: average expression
#'         \item t: moderated t-statistic
#'         \item P.Value: raw p-value
#'         \item adj.P.Val: Benjamini-Hochberg adjusted p-value
#'         \item neg_log10.adj.P.Val: -log10 of adjusted p-value
#'       }
#'     \item design: The design matrix used in the analysis
#'     \item coefficients: Estimated model coefficients
#'     \item model_terms: Names of the model terms (design matrix columns)
#'     \item contrast_matrix: The contrast matrix used for testing
#'   }
#'
#' @export
#'
#' @examples
#' # Basic usage with a simple group comparison:
#' # results <- perform_limma_analysis(qfeatures_obj, "protein", ~group, "groupB - groupA")
#' 
#' # More complex design with batch effect:
#' # results <- perform_limma_analysis(qfeatures_obj, "protein", ~group + batch, "groupB - groupA")
#' 
#' # Multiple group comparison:
#' # results <- perform_limma_analysis(qfeatures_obj, "protein", ~treatment, "treatmentB - treatmentA")
perform_limma_analysis <- function(qf,
                                   assay_name,
                                   formula,
                                   contrast) {
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

  # Fit the linear model
  fit <- limma::lmFit(exprs, design)

  # Automatically infer levels from the design matrix
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast,
    levels = colnames(design)  # This is key
  )

  # Apply contrasts and empirical Bayes
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  # Output full table of results
  top_table <- limma::topTable(fit2, adjust.method = "BH", number = Inf) |>
    # Converting for shiny app.
    dplyr::mutate(neg_log10.adj.P.Val = -log10(adj.P.Val),.after = adj.P.Val)


  # Need to have some ways that plot_volcano can figure out what the logFC is
  # Relative to.
  return(list(
    top_table = top_table,
    design = design,
    coefficients = fit2$coefficients,
    model_terms = colnames(design),
    contrast_matrix = contrast_matrix
  ))
}


