#' Convert QFeatures Object to Tidy Data Format
#'
#' Transforms a QFeatures object into a long-format tibble by combining assay data,
#' sample metadata (colData), and feature metadata (rowData). This function is
#' particularly useful for preparing data for visualization and statistical analysis
#' using tidyverse functions.
#'
#' @param qf A QFeatures object containing the proteomics data and metadata
#' @param assay_name Character string specifying which assay to extract from the
#'   QFeatures object. This should be the name of an existing assay in the object.
#'
#' @return A tibble in long format containing:
#'   \itemize{
#'     \item rowid: Feature identifiers (e.g., protein IDs)
#'     \item file: Sample identifiers
#'     \item value: Quantitative measurements (e.g., intensities, abundances)
#'     \item All columns from colData (sample metadata)
#'     \item All columns from rowData (feature metadata)
#'   }
#'   The data is organized such that each row represents a single measurement
#'   (feature-sample combination), making it ideal for:
#'   \itemize{
#'     \item Creating visualizations with ggplot2
#'     \item Performing statistical analyses
#'     \item Data filtering and transformation
#'     \item Integration with other tidyverse functions
#'   }
#'
#' @export
#'
#' @examples
#' # Basic usage to prepare data for visualization:
#' # tidy_data <- tidy_conduit(qfeatures_obj, "protein")
#' # 
#' # # Create a boxplot of protein abundances by group:
#' # tidy_data |>
#' #   ggplot2::ggplot(ggplot2::aes(x = group, y = value)) +
#' #   ggplot2::geom_boxplot()
#' 
#' # Use with other analysis functions:
#' # - Calculate summary statistics
#' # tidy_data |>
#' #   dplyr::group_by(group) |>
#' #   dplyr::summarise(mean_abundance = mean(value))
#' 
#' # - Filter and transform data
#' # tidy_data |>
#' #   dplyr::filter(value > 0) |>
#' #   dplyr::mutate(log_value = log2(value))
#'
#' @note
#' This function is designed to work seamlessly with other functions in the package:
#' \itemize{
#'   \item Use with \code{\link[plot_relative_abundance]{plot_relative_abundance}}
#'     for abundance visualization
#'   \item Use with \code{\link[plot_density]{plot_density}} for distribution plots
#'   \item Use with \code{\link[plot_selected_features]{plot_selected_features}}
#'     for feature-specific plots
#' }
#' 
#' Important considerations:
#' \itemize{
#'   \item The function maintains all metadata from both samples and features
#'   \item Memory usage increases with the size of the QFeatures object
#'   \item For very large datasets, consider filtering before tidying
#'   \item The output format is optimized for tidyverse operations
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[QFeatures]{QFeatures}} for the input data structure
#'   \item \code{\link[plot_relative_abundance]{plot_relative_abundance}} for
#'     visualizing relative abundances
#'   \item \code{\link[plot_density]{plot_density}} for distribution plots
#'   \item \code{\link[plot_selected_features]{plot_selected_features}} for
#'     feature-specific visualizations
#' }
tidy_conduit = function(qf,
                        assay_name){

  # Extracting the assay data
  assay_data <- SummarizedExperiment::assay(qf, assay_name)  |>
    tibble::rownames_to_column("rowid")

  # Finding the ncol of assay data. This is the last col in pivot_longer()
  ncol = ncol(assay_data)

  # Pivoting to wide format
  assay_data = assay_data |>
    tidyr::pivot_longer(2:ncol,
                        names_to = "file",
                        values_to = "value")

  #Extracting the rowdata
  row_data = SummarizedExperiment::rowData(qf[[assay_name]]) |>
    as.data.frame() |>
    tibble::rownames_to_column("rowid")

  # Extracting the colData
  colData = SummarizedExperiment::colData(qf) |>
    as.data.frame() |>
    tibble::rownames_to_column("file")

  # Combining the data sources
  out = colData |>
    dplyr::inner_join(assay_data, by = "file") |>
    dplyr::inner_join(row_data, by = "rowid")

  return(out)

}
