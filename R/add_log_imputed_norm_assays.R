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
add_log_imputed_norm_assays <- function(qf,
                                        base = 2,
                                        impute_method = "min",
                                        norm_method = "none") {

  # Step 1: Replace any zeros with NA
  qf <- replace_zero_with_na(qf)

  # Step 2: Log2 Transform All Assays. Add a pseudocount of 1 to prevent log(0)
  qf <- QFeatures::logTransform(qf,
                                base = base,
                                pc = 1,
                                i = names(qf),
                                name = paste0(names(qf), "_log", base))

  # Identify log2-transformed assays
  log2_assays = names(qf)[grepl(".*log2", names(qf))]

  # Step 3: Impute Missing Values
  qf <- QFeatures::impute(qf,
                          method = impute_method,
                          i = log2_assays,
                          name = paste0(log2_assays, "_imputed"))

  # Step 4: Normalize Log-transformed Imputed Values
  if (norm_method != "none") {
    # Identify imputed assays
    imputed_assays = names(qf)[grepl(".*imputed", names(qf))]

    # Perform normalization only if imputed assays exist
    if (length(imputed_assays) > 0) {
      qf <- QFeatures::normalize(qf,
                                 method = norm_method,
                                 i = imputed_assays,
                                 name = paste0(imputed_assays, "_norm_", norm_method))
    }
  }

  return(qf)
}
