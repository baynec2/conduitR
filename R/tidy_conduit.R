#' tidy_conduit
#'
#' tidy data from a QFeatures object.This function combines,
#' assay data, colData, and rowData into a tbl_df in long format.
#'
#' @param qf
#' @param assay_name
#'
#' @returns
#' @export
#'
#' @examples
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
