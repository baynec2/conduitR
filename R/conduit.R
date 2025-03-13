#' conduit Class
#'
#' A class to store QFeatures and additional data about the database used for
#' detection, as well as the detected proteins.
#'
#' @slot qfeatures A QFeatures object containing proteomics data.
#' @slot experiment_metadata A tibble containing the experimental metadata.
#' @slot database_taxonomy A tibble containing taxonomy information for the database.
#' @slot database_metrics A tibble containing metrics related to the database.
#' @slot detected_protein_taxonomy A tibble containing taxonomy information for detected proteins.
#' @slot detected_protein_metrics A tibble containing metrics for detected proteins.
#' @export
setClass(Class = "conduit",
         slots = list(
           QFeatures = "QFeatures",
           database_taxonomy = "tbl_df",
           database_metrics = "tbl_df",
           detected_protein_taxonomy = "tbl_df",
           detected_protein_metrics = "tbl_df"
         ))

#' Initialize a Conduit object
#' @param qf A QFeatures object.
#' @param experiment_metadata A tibble containing the experiment metadata.
#' @param database_taxonomy A tibble containing the database taxonomy.
#' @param database_metrics A tibble containing the database metrics.
#' @param detected_protein_taxonomy A tibble containing the detected protein taxonomy.
#' @param detected_protein_metrics A tibble containing the detected protein metrics.
#' @return A Conduit object.
#' @export
setMethod("initialize", "conduit",
          function(.Object, QFeatures,
                   database_taxonomy = NULL,
                   database_metrics = NULL,
                   detected_protein_taxonomy = NULL,
                   detected_protein_metrics = NULL) {

            .Object@QFeatures <- QFeatures
            .Object@database_taxonomy <- database_taxonomy
            .Object@database_metrics <- database_metrics
            .Object@detected_protein_taxonomy <- detected_protein_taxonomy
            .Object@detected_protein_metrics <- detected_protein_metrics

            return(.Object)
          })

#' Show method for Conduit object
#'
#' @param object A Conduit object.
#' @export
setMethod("show", "conduit",
          function(object) {
            cat("conduit Summary:\n")
            cat("Number of Assays in QFeatures:", length(SummarizedExperiment::assays(object@QFeatures)), "\n")
            cat("Detected Proteins Taxonomy:", if (!is.null(object@detected_protein_taxonomy)) "Available" else "Not Available", "\n")
            cat("Detected Proteins Metrics:", if (!is.null(object@detected_protein_metrics)) "Available" else "Not Available", "\n")
          })
