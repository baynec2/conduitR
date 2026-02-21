#' Convert DIA-NN Parquet Output into a QFeatures Object
#'
#' Reads a DIA-NN parquet file, applies standard quality filters (Q.Value,
#' Lib.Q.Value, Lib.PG.Q.Value, Global.PG.Q.Value, PG.Q.Value), and builds a
#' QFeatures object with three assays: precursors, peptides (aggregated by
#' stripped sequence), and protein groups (using DIA-NN's PG.MaxLFQ), with
#' assay links between them.
#'
#' @param diann_parquet_fp Character. Path to the DIA-NN output parquet file.
#'
#' @return A `QFeatures` object with assays `precursors`, `peptides`, and
#'   `protein_groups`, and assay links from precursors to peptides and from
#'   peptides to protein groups.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' qf <- diann_to_qfeatures("path/to/report.parquet")
#' QFeatures::assayNames(qf)
#' # "precursors" "peptides" "protein_groups"
#' dim(qf)[["protein_groups"]]
#' }
diann_to_qfeatures = function(diann_parquet_fp){
  # Read DIANN parquet file
  diann_parquet <- arrow::read_parquet(diann_parquet_fp) |>
    # Filtering per DIA-NN User Manual
    dplyr::filter(Q.Value <= 0.01,
                  Lib.Q.Value <= 0.01,
                  Lib.PG.Q.Value <= 0.01,
                  Global.PG.Q.Value <= 0.01,
                  PG.Q.Value <= 0.05)

  # Define columns to keep: only shared values across runs
  cols_to_keep <- c("Run","Precursor.Id","Modified.Sequence",
                    "Stripped.Sequence","Precursor.Charge",
                    "Precursor.Lib.Index","Decoy","Proteotypic","Precursor.Mz",
                    "Protein.Ids","Protein.Group","Protein.Names","Genes",
                    "Global.Q.Value","Lib.Q.Value","Precursor.Normalised"
                    )

  # Keep only precursor-relevant columns and pivot to wide format
  precursors_df <- diann_parquet |>
    dplyr::select(dplyr::any_of(cols_to_keep)) |>
    tidyr::pivot_wider(
      names_from = "Run",
      values_from = "Precursor.Normalised"
    )
  # Define the end of the rowData
  rowdata_end_index = which(colnames(precursors_df) == "Lib.Q.Value")
  # Create SummarizedExperiment from Precursor data
  rowdat_cols <- colnames(precursors_df)[1:rowdata_end_index]

  rowdat <- precursors_df |> dplyr::select(dplyr::any_of(rowdat_cols))
  assay_mat <- precursors_df |>
    dplyr::select(-dplyr::any_of(rowdat_cols)) |>
    as.matrix()
  rownames(assay_mat) <- rowdat$Precursor.Id  # or another unique ID column

  se_prec <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = assay_mat),
    rowData = rowdat
  )

  # Create a QFeatures object containing the precursor
  qf <- QFeatures::QFeatures(list(precursors = se_prec))

  # Aggregate precursors to peptides
  qf <- QFeatures::aggregateFeatures(
    qf,
    i = "precursors",
    fcol = "Stripped.Sequence",   # grouping column in rowData
    name = "peptides",
    fun = colSums,
    na.rm = TRUE
  )

  # Add protein groups to the QFeatures assay
  # Instead of aggregating peptides to protein groups, we will use the values
  # already calculated by DIA-NN
  protein_groups = diann_parquet |>
    dplyr::select(Run,Protein.Group,Protein.Names,Genes,
                  PG.MaxLFQ,Global.PG.Q.Value) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(names_from = Run,
                       values_from = PG.MaxLFQ)

  # Counting peptide sequences per protein_group
  n_peptides_detected <- diann_parquet |>
    dplyr::select(Protein.Group,Stripped.Sequence) |>
    dplyr::distinct() |>
    dplyr::group_by(Protein.Group) |>
    dplyr::summarise(.n = dplyr::n()) |>
    dplyr::ungroup()

  protein_groups <- dplyr::left_join(protein_groups,
                                     n_peptides_detected,
                                     by = "Protein.Group")


  # Defining protein group row data
  protein_rowdata <- protein_groups |>
    dplyr::select(Protein.Group,Protein.Names,
                  Genes, Global.PG.Q.Value,.n) |>
    dplyr::mutate(Protein.Group2 = Protein.Group) |>
    tibble::column_to_rownames("Protein.Group2")



  # Creating SummarizedExperiment for protein groups
  protein_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = as.matrix(protein_groups |>
                                        dplyr::select(-Protein.Names,
                                                     -Global.PG.Q.Value,
                                                     -Genes,
                                                     -Protein.Group,
                                                     -.n))),
    rowData = protein_rowdata
  )

  # Adding protein groups assay to QFeatures object
  qf <- QFeatures::addAssay(qf,protein_se,name = "protein_groups")

  # Adding assay links from peptides to protein groups
  qf <- QFeatures::addAssayLink(qf,
                                from = "peptides",
                                to = "protein_groups",
                                varFrom = "Protein.Group",
                                varTo = "Protein.Group")

  return(qf)

}


