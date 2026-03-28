#' Create a provenance list for a conduit object
#'
#' Constructs and validates a provenance list capturing reproducibility
#' metadata for a [conduit-class] object. The returned list is suitable for
#' passing to the \code{provenance} argument of [build_conduit_obj()] or
#' assigning to the \code{@@provenance} slot directly.
#'
#' @param workflow_version character(1). Version string of conduit-ascent
#'   (the Snakemake workflow that generated the data), e.g. \code{"1.2.3"}.
#' @param uniprotkb_release character(1). UniProtKB release identifier used
#'   as the reference database, e.g. \code{"2024_05"}. Defaults to the result
#'   of [get_uniprotkb_release()], which queries the UniProt REST API.
#' @param config Optional named \code{list} of \code{tbl_df} objects, where
#'   each element represents a distinct configuration source (e.g.
#'   \code{snakemake_yaml}, \code{diann_run_cfg}, \code{runtime}). Every
#'   tibble must have columns \code{parameter} (character) and \code{value}
#'   (character). Pass \code{NULL} (the default) to omit configuration details.
#' @param generated_date Date. The date the conduit object was generated.
#'   Defaults to today via [Sys.Date()].
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{workflow_version}{character(1)}
#'     \item{generated_date}{Date}
#'     \item{uniprotkb_release}{character(1)}
#'     \item{config}{named list of tbl_df, or NULL}
#'   }
#' @export
#'
#' @examples
#' prov <- create_provenance(
#'   workflow_version  = "1.2.3",
#'   uniprotkb_release = "2024_05",
#'   config = list(
#'     snakemake_yaml = tibble::tibble(parameter = "fdr", value = "0.01"),
#'     runtime        = tibble::tibble(parameter = "snakemake_version", value = "8.4.6")
#'   )
#' )
create_provenance <- function(workflow_version,
                              uniprotkb_release = get_uniprotkb_release(),
                              config = NULL,
                              generated_date = Sys.Date()) {
  stopifnot(
    is.character(workflow_version),
    length(workflow_version) == 1
  )
  stopifnot(
    is.character(uniprotkb_release),
    length(uniprotkb_release) == 1
  )
  if (!is.null(config)) {
    stopifnot(is.list(config), !is.data.frame(config))
    for (nm in names(config)) {
      stopifnot(
        is.data.frame(config[[nm]]),
        all(c("parameter", "value") %in% names(config[[nm]]))
      )
    }
  }
  list(
    workflow_version  = workflow_version,
    generated_date    = as.Date(generated_date),
    uniprotkb_release = uniprotkb_release,
    config            = config
  )
}
