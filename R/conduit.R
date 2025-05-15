#' conduit Class
#'
#' A class to store QFeatures and additional data about the database used for
#' detection, as well as the detected proteins.
#'
#' @slot QFeatures A QFeatures object containing proteomics data.
#' @slot combined_metrics A tibble containing metrics about how many proteins were detected at various levels.
#' @slot database_protein_taxonomy A tibble containing taxonomy information for the database.
#' @slot detected_protein_taxonomy A tibble containing taxonomy information for detected proteins.
#' @slot detected_protein_metrics A tibble containing metrics for detected proteins.
#' @export
setClass(Class = "conduit",
         slots = list(
           QFeatures = "QFeatures",
           combined_metrics = "tbl_df",
           database_protein_taxonomy = "tbl_df",
           detected_protein_taxonomy = "tbl_df"
           )
         )

#' Initialize a Conduit object
#' @param QFeatures A QFeatures object.
#' @param combined_metrics A tibble containing proteins and their assigned taxonomy, as used in database
#' @param database_protein_taxonomy A tibble containing the database metrics.
#' @param detected_protein_taxonomy A tibble containing the detected protein taxonomy.
#' @return A Conduit object.
#' @export
setMethod("initialize", "conduit",
          function(.Object,
                   QFeatures,
                   combined_metrics = NULL,
                   database_protein_taxonomy = NULL,
                   detected_protein_taxonomy = NULL) {
            .Object@QFeatures <- QFeatures
            .Object@combined_metrics <- combined_metrics
            .Object@database_protein_taxonomy <- database_protein_taxonomy
            .Object@detected_protein_taxonomy <- detected_protein_taxonomy
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
            cat("Combined (in Database + Detected) Metrics:", if (!is.null(object@combined_metrics)) "Available" else "Not Available", "\n")
            cat("Database Proteins Taxonomy:", if (!is.null(object@database_protein_taxonomy)) "Available" else "Not Available", "\n")
            cat("Detected Proteins Taxonomy:", if (!is.null(object@detected_protein_taxonomy)) "Available" else "Not Available", "\n")
          })
