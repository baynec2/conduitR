#' Add Log-Transformed, Imputed, and Normalized Assays
#'
#' Applies the log-transform, imputation, and normalization pipeline to every
#' assay currently in a QFeatures object. Each step creates a new assay:
#' \enumerate{
#'   \item Zeros replaced with NA
#'   \item Log transformation
#'   \item Imputation
#'   \item Normalization (or copy if norm_method = "none")
#' }
#'
#' @param qf A QFeatures object.
#' @param base Numeric log base (default: 2).
#' @param impute_method Imputation method supported by \code{QFeatures::impute()} (default: "MinDet").
#' @param norm_method Normalization method supported by \code{QFeatures::normalize()}, or
#'   "none" to skip (default: "none").
#' @param min_n Minimum number of non-NA values required to keep a row (default: 1).
#'
#' @return A QFeatures object with additional derived assays for each original assay.
#' @export
#'
#' @examples
#' # qf_transformed <- add_log_imputed_norm_assays(qfeatures_obj)
#' # qf_transformed <- add_log_imputed_norm_assays(qfeatures_obj, base = 10)
add_log_imputed_norm_assays <- function(qf,
                                        base = 2,
                                        impute_method = "MinDet",
                                        norm_method = "none",
                                        min_n = 1) {
  for (assay in names(qf)) {
    qf <- add_log_imputed_norm_assay(qf,
                                     assay = assay,
                                     base = base,
                                     impute_method = impute_method,
                                     norm_method = norm_method,
                                     min_n = min_n)
  }
  return(qf)
}

#' Add Log-Transformed, Imputed, and Normalized Assays for a Single Assay
#'
#' Transforms one assay in a QFeatures object through log transformation,
#' imputation, and optional normalization.
#'
#' @param qf A QFeatures object.
#' @param assay Name of the assay to transform (default: "protein_groups").
#' @param base Numeric log base (default: 2).
#' @param impute_method Imputation method supported by \code{QFeatures::impute()} (default: "MinDet").
#' @param norm_method Normalization method or "none" (default: "none").
#' @param min_n Minimum number of non-NA values required to keep a row (default: 1).
#'
#' @return A QFeatures object with additional \code{{assay}_log{base}},
#'   \code{{assay}_log{base}_imputed}, and \code{{assay}_log{base}_imputed_norm} assays.
#' @export
add_log_imputed_norm_assay <- function(qf,
                                       assay = "protein_groups",
                                       base = 2,
                                       impute_method = "MinDet",
                                       norm_method = "none",
                                       min_n = 1) {

  # Allowed imputation methods
  allowed_methods <- c(MsCoreUtils::imputeMethods(), "no_missing")
  if (!(impute_method %in% allowed_methods)) {
    stop("Supplied impute method not allowed.")
  }

  # Names of resulting assays — encode method in name so provenance is self-documenting
  log_name    <- paste0(assay, "_log", base)
  impute_name <- paste0(log_name, "_", impute_method)
  norm_name   <- paste0(impute_name, "_", norm_method)

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


  # Step 5: Filter by .n (only when the column exists in rowData)
  if (!is.null(min_n) && ".n" %in% colnames(SummarizedExperiment::rowData(qf[[norm_name]]))) {
    se <- qf[[norm_name]]
    se <- se[SummarizedExperiment::rowData(se)$.n >= min_n, ]
    qf[[norm_name]] <- se
  }

  return(qf)
}
