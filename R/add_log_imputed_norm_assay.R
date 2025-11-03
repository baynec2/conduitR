#' Add Log-Transformed, Imputed, and Normalized Assays
#'
#' Performs a series of data transformations on all assays in a QFeatures object:
#' 1. Replaces zeros with NA values
#' 2. Applies log transformation with a pseudocount
#' 3. Imputes missing values
#' 4. Optionally normalizes the data
#' Each step creates a new assay in the QFeatures object.
#'
#' @param qf A QFeatures object containing the data to transform
#' @param base Numeric value specifying the base for log transformation (default: 2)
#' @param impute_method Character string specifying the imputation method to use.
#'   Must be one of the methods supported by QFeatures::impute() (default: "min")
#' @param norm_method Character string specifying the normalization method to use.
#'   Must be one of the methods supported by QFeatures::normalize(), or "none"
#'   to skip normalization (default: "none")
#'
#' @return A QFeatures object with additional assays:
#'   \itemize{
#'     \item {assay_name}_log{base}: Log-transformed data
#'     \item {assay_name}_log{base}_imputed: Log-transformed and imputed data
#'     \item {assay_name}_log{base}_imputed_norm_{norm_method}: Normalized data
#'       (only if norm_method is not "none")
#'   }
#'
#' @export
#'
#' @examples
#' # Basic usage with default settings:
#' # qf_transformed <- add_log_imputed_norm_assays(qfeatures_obj)
#'
#' # Using log10 transformation and min imputation:
#' # qf_transformed <- add_log_imputed_norm_assays(qfeatures_obj, base = 10)
#'
#' # With normalization:
#' # qf_transformed <- add_log_imputed_norm_assays(qfeatures_obj,
#' #                                             norm_method = "center.scale")
#'
#' # The transformed assays can be used with plotting functions:
#' # plot_density(qf_transformed, "protein", color = "group")
add_log_imputed_norm_assay <- function(qf,
                                       assay = "protein_groups",
                                       base = 2,
                                       impute_method = "min",
                                       norm_method = "none",
                                       min_n = 1) {

  # Allowed imputation methods
  allowed_methods <- c(MsCoreUtils::imputeMethods(), "no_missing")
  if (!(impute_method %in% allowed_methods)) {
    stop("Supplied impute method not allowed.")
  }

  # Names of resulting assays
  log_name <- paste0(assay, "_log", base)
  impute_name <- paste0(log_name, "_imputed")
  norm_name <- paste0(impute_name, "_norm")

  # Step 1: Replace zeros with NA
  qf <- replace_zero_with_na(qf)

  # Step 2: Log-transform
  qf <- QFeatures::logTransform(qf,
                                base = base,
                                pc = 1,
                                i = assay,
                                name = log_name)

  # Step 3: Impute
  if (impute_method == "no_missing") {
    se <- qf[[log_name]]
    se_no_na <- se[rowSums(is.na(SummarizedExperiment::assay(se))) == 0, ]
    qf[[impute_name]] <- se_no_na
  } else {
    qf <- QFeatures::impute(qf,
                            method = impute_method,
                            i = log_name,
                            name = impute_name)
  }

  # Step 4: Normalize
  if (norm_method == "none") {
    # If no normalization has been applied, we will just propogate that to the
    # norm value. This will make the logic easier for the shiny app
    qf[[norm_name]] <- qf[[impute_name]]

  } else {
    qf <- QFeatures::normalize(qf,
                               method = norm_method,
                               i = impute_name,
                               name = norm_name)
  }


  # Step 5: Filter by .n
  if (!is.null(min_n)) {
    se <- qf[[norm_name]]
    se <- se[SummarizedExperiment::rowData(se)$.n >= min_n, ]
    qf[[norm_name]] <- se
  }

  return(qf)
}
