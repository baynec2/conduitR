#' perform_limma_analysis
#'
#' Performs limma differential expression analysis on a specified assay within a QFeatures object.
#'
#' @param qf_object A QFeatures object
#' @param assay_name Name of the assay to analyze
#' @param formula A formula used to build the design matrix from colData
#' @param contrast A character string specifying the contrast (e.g., "groupB - groupA")
#'
#' @return A list containing:
#'   \item{top_table}{The result of limma::topTable()}
#'   \item{design}{The design matrix}
#'   \item{coefficients}{Estimated model coefficients}
#'   \item{model_terms}{Names of the model terms (design matrix columns)}
#'
#' @export
#'
#' @examples
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


