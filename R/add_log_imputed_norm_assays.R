#' add_log_imputed_norm_assays
#'
#'
#' @param qf
#' @param impute_method
#' @param norm_method
#'
#' @returns
#' @export
#'
#' @examples
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
