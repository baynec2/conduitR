#' Convert DIANN Parquet file into a QFeatures object containing precursors,
#' peptides, and proteins, all with assay links.
#'
#' @param diann_parquet_fp
#'
#' @returns a QFeatures object
#' @export
#'
#' @examples
diann_to_qfeatures = function(diann_parquet_fp){

  # Read DIANN parquet file
  diann_parquet <- arrow::read_parquet(diann_parquet_fp) |>
    # Filtering per DIA-NN User Manual
    dplyr::filter(Q.Value <= 0.01,
                  Lib.Q.Value <= 0.01,
                  Lib.PG.Q.Value <= 0.01,
                  Global.PG.Q.Value <= 0.01,
                  PG.Q.Value <= 0.05)

  # Define columns to remove- these are irrelevant for precursor level
  cols_to_remove <- c(
    "PG.TopN", "PG.MaxLFQ", "Genes.TopN",
    "Genes.MaxLFQ", "Genes.MaxLFQ.Unique", "PG.MaxLFQ.Quality",
    "Genes.MaxLFQ.Quality", "Genes.MaxLFQ.Unique.Quality",
    "Peptidoform.Q.Value", "Global.Peptidoform.Q.Value",
    "Lib.Peptidoform.Q.Value","PTM.Site.Confidence",
    "Site.Occupancy.Probabilities", "Protein.Sites",
    "Lib.PTM.Site.Confidence", "Translated.Q.Value", "Channel.Q.Value",
    "PG.Q.Value", "PG.PEP", "GG.Q.Value",
    "Protein.Q.Value", "Global.PG.Q.Value", "Lib.PG.Q.Value"
  )

  # Keep only precursor-relevant columns and pivot to wide format
  precursors_df <- diann_parquet |>
    dplyr::select(-dplyr::any_of(cols_to_remove)) |>
    tidyr::pivot_wider(
      names_from = "Run",
      values_from = "Precursor.Normalised"
    )

  # Create SummarizedExperiment from Precursor data
  rowdat_cols <- setdiff(colnames(precursors_df),
                         colnames(diann_parquet)[grepl("Run",
                                                       colnames(diann_parquet))])
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
    fun = colSums
  )

  # Add protein groups to the QFeatures assay
  # Instead of aggregating peptides to protein groups, we will use the values
  # already calculated by DIA-NN
  protein_groups = diann_parquet |>
    dplyr::select(Run,Protein.Group,Protein.Names,Genes,
                  PG.MaxLFQ.Quality,PG.MaxLFQ,Global.PG.Q.Value) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(names_from = Run,
                       values_from = PG.MaxLFQ)

  # Defining protein group row data
  protein_rowdata <- protein_groups |>
    dplyr::select(Protein.Group,Protein.Names,
                  Genes,PG.MaxLFQ.Quality, Global.PG.Q.Value) |>
    dplyr::mutate(Protein.Group2 = Protein.Group) |>
    tibble::column_to_rownames("Protein.Group2")

  # Creating SummarizedExperiment for protein groups
  protein_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = as.matrix(protein_groups |>
                                        dplyr::select(-Protein.Names,
                                                     -Genes, -PG.MaxLFQ.Quality,
                                                     -Global.PG.Q.Value,
                                                     -Protein.Group))),
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


