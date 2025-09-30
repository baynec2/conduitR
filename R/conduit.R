#' Conduit Class for Integrated Proteomics Analysis
#'
#' The `conduit` class provides a unified container for proteomics data,
#' integrating quantitative measurements, taxonomic information, and detection metrics.
#' It serves as the core data structure of the **conduit** package, enabling
#' downstream analysis and visualization of protein abundance, taxonomy, and
#' quality statistics.
#'
#' @slot QFeatures A [QFeatures::QFeatures-class] object containing proteomics data,
#'   including:
#'   \itemize{
#'     \item Protein abundance measurements across samples
#'     \item Sample metadata (`colData`)
#'     \item Feature metadata (`rowData`)
#'     \item Multiple assays (e.g., raw, normalized, or transformed data)
#'   }
#' @slot metrics A named list of tibbles, each containing a distinct set of metrics
#'   from the experiment. Examples include file-level statistics reported by DIA-NN,
#'   summaries of detected taxonomy, or taxonomic coverage. This structure allows for
#'   flexible storage of multiple, heterogeneous metric types and is designed to be
#'   extensible as new metrics are added.
#'   \itemize{
#'     \item diann_stats: File-level statistics from DIA-NN
#'     \item detected_taxonomy_summary: Summaries of detected taxa
#'     \item protein_coverage_by_taxa: Number of proteins detected per taxon
#'     relative to the total number of proteins in the database for that taxon.
#'     The percentage represents coverage of the taxon's proteome.
#'   }
#' @slot database A tibble containing taxonomic information for all proteins in
#'   the reference database, including:
#'   \itemize{
#'     \item Complete taxonomic classification (domain through species)
#'     \item Organism type and source information
#'     \item Protein identifiers and annotations
#'   }
#' @slot annotations A tibble (long format) containing taxonomy and detection
#'   information for identified protein groups, enriched with functional annotations
#'   retrieved from external sources (e.g., GO, KEGG, Pfam). This slot complements
#'   the `database` slot by providing additional metadata specific to detected
#'   proteins, beyond what is present in the original FASTA database.
#'   Includes:
#'   \itemize{
#'     \item Taxonomic classification of detected proteins
#'     \item Detection statistics for each taxon or protein group
#'     \item Functional annotations retrieved via API calls (e.g., GO, KEGG, Pfam)
#'     \item Quality metrics for protein identification
#'   }
#'   \itemize{
#'     \item Taxonomic classification of detected proteins
#'     \item Detection statistics by taxon
#'     \item Quality metrics for protein identification
#'   }
#'
#' @details
#' The `conduit` class is designed to combine raw quantitative proteomics data
#' with functional and taxonomic context, streamlining analyses that link protein
#' abundance with microbial community structure.
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
#' # Summarize the object
#' # show(conduit_obj)
#'
#' # Access QFeatures data
#' # assays(conduit_obj@QFeatures)
#' # colData(conduit_obj@QFeatures)
#'
#' # Analyze detection metrics
#' # calc_percent_proteins_detected(conduit_obj)
#'
#' # Visualize taxonomy
#' # plot_taxa_tree(conduit_obj@database)
#' # plot_sunburst(conduit_obj@annotations)
setClass(
  Class = "conduit",
  slots = list(
    QFeatures   = "QFeatures",
    metrics     = "list",
    database    = "tbl_df",
    annotations = "tbl_df"
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
                   metrics = NULL,
                   database = NULL,
                   annotations = NULL) {
            .Object@QFeatures   <- QFeatures
            .Object@metrics     <- metrics
            .Object@database    <- database
            .Object@annotations <- annotations
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
            cat("Number of Assays in QFeatures:",
                length(SummarizedExperiment::assays(object@QFeatures)), "\n")

            cat("Metrics:",
                if (!is.null(object@metrics)) "Available" else "Not Available", "\n")

            cat("Database Taxonomy:",
                if (!is.null(object@database)) "Available" else "Not Available", "\n")

            cat("Annotations:",
                if (!is.null(object@annotations)) "Available" else "Not Available", "\n")
          })
