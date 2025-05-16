#' Conduit Class for Proteomics Data Analysis
#'
#' A comprehensive class for storing and analyzing proteomics data, combining QFeatures
#' objects with taxonomic and detection metrics. This class serves as the central data
#' structure for the conduit package, enabling integrated analysis of protein abundance,
#' taxonomy, and detection statistics.
#'
#' @slot QFeatures A QFeatures object containing proteomics data, including:
#'   \itemize{
#'     \item Protein abundance measurements across samples
#'     \item Sample metadata (colData)
#'     \item Feature metadata (rowData)
#'     \item Multiple assays (e.g., raw, normalized, log-transformed data)
#'   }
#' @slot combined_metrics A tibble containing comprehensive metrics about protein detection,
#'   including:
#'   \itemize{
#'     \item Number of proteins detected at each taxonomic level
#'     \item Detection rates and coverage statistics
#'     \item Quality metrics for protein identification
#'   }
#' @slot database_protein_taxonomy A tibble containing taxonomy information for all proteins
#'   in the database, including:
#'   \itemize{
#'     \item Complete taxonomic classification (domain to species)
#'     \item Organism type and source information
#'     \item Protein identifiers and annotations
#'   }
#' @slot detected_protein_taxonomy A tibble containing taxonomy information specifically for
#'   detected proteins, including:
#'   \itemize{
#'     \item Taxonomic classification of detected proteins
#'     \item Detection statistics for each taxon
#'     \item Quality metrics for identified proteins
#'   }
#'
#' @export
#'
#' @examples
#' # Create a conduit object from workflow outputs:
#' # conduit_obj <- create_conduit_obj(
#' #   QFeatures.rds = "results/qfeatures.rds",
#' #   combined_metrics.tsv = "results/metrics.tsv",
#' #   database_protein_taxonomy.tsv = "results/database_taxonomy.tsv",
#' #   detected_protein_taxonomy.tsv = "results/detected_taxonomy.tsv"
#' # )
#' 
#' # Access and analyze the data:
#' # - View summary of the object
#' # show(conduit_obj)
#' 
#' # - Access QFeatures data
#' # assays <- assays(conduit_obj@QFeatures)
#' # metadata <- colData(conduit_obj@QFeatures)
#' 
#' # - Analyze protein detection rates
#' # rates <- calc_percent_proteins_detected(conduit_obj)
#' 
#' # - Create taxonomic visualizations
#' # plot_taxa_tree(conduit_obj@database_protein_taxonomy)
#' # plot_sunburst(conduit_obj@detected_protein_taxonomy)
#'
#' @note
#' The conduit class is designed to work seamlessly with other functions in the package:
#' \itemize{
#'   \item Use \code{\link[create_conduit_obj]{create_conduit_obj}} to create new instances
#'   \item Use \code{\link[calc_percent_proteins_detected]{calc_percent_proteins_detected}}
#'     for detection rate analysis
#'   \item Use \code{\link[plot_taxa_tree]{plot_taxa_tree}} and
#'     \code{\link[plot_sunburst]{plot_sunburst}} for taxonomic visualization
#'   \item Use \code{\link[perform_limma_analysis]{perform_limma_analysis}} for
#'     differential expression analysis
#' }
#' For large datasets, consider:
#' \itemize{
#'   \item Using appropriate filtering before creating the object
#'   \item Managing memory usage when working with multiple assays
#'   \item Utilizing the package's visualization functions for efficient data exploration
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[create_conduit_obj]{create_conduit_obj}} for creating new instances
#'   \item \code{\link[QFeatures]{QFeatures}} for details about the QFeatures class
#'   \item \code{\link[plot_taxa_tree]{plot_taxa_tree}} and
#'     \code{\link[plot_sunburst]{plot_sunburst}} for visualization
#'   \item \code{\link[calc_percent_proteins_detected]{calc_percent_proteins_detected}}
#'     for detection rate analysis
#' }
setClass(Class = "conduit",
         slots = list(
           QFeatures = "QFeatures",
           combined_metrics = "tbl_df",
           database_protein_taxonomy = "tbl_df",
           detected_protein_taxonomy = "tbl_df"
           )
         )

#' Initialize a Conduit Object
#'
#' Creates a new instance of the conduit class, combining proteomics data with
#' taxonomic and detection metrics into a single integrated object.
#'
#' @param QFeatures A QFeatures object containing proteomics data and metadata
#' @param combined_metrics A tibble containing metrics about protein detection and
#'   taxonomy (default: NULL)
#' @param database_protein_taxonomy A tibble containing taxonomy information for
#'   all proteins in the database (default: NULL)
#' @param detected_protein_taxonomy A tibble containing taxonomy information for
#'   detected proteins (default: NULL)
#'
#' @return A new conduit object containing the provided data
#'
#' @export
#'
#' @examples
#' # Create a new conduit object:
#' # new_conduit <- new("conduit",
#' #   QFeatures = qfeatures_obj,
#' #   combined_metrics = metrics_tibble,
#' #   database_protein_taxonomy = db_taxonomy,
#' #   detected_protein_taxonomy = detected_taxonomy
#' # )
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

#' Show Method for Conduit Object
#'
#' Displays a summary of the conduit object's contents, including the number of
#' assays and availability of various data components.
#'
#' @param object A conduit object to summarize
#'
#' @export
#'
#' @examples
#' # Display summary of a conduit object:
#' # show(conduit_obj)
setMethod("show", "conduit",
          function(object) {
            cat("conduit Summary:\n")
            cat("Number of Assays in QFeatures:", length(SummarizedExperiment::assays(object@QFeatures)), "\n")
            cat("Combined (in Database + Detected) Metrics:", if (!is.null(object@combined_metrics)) "Available" else "Not Available", "\n")
            cat("Database Proteins Taxonomy:", if (!is.null(object@database_protein_taxonomy)) "Available" else "Not Available", "\n")
            cat("Detected Proteins Taxonomy:", if (!is.null(object@detected_protein_taxonomy)) "Available" else "Not Available", "\n")
          })
