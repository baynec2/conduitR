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
#' @param config Optional \code{tbl_df} (or \code{data.frame}) with columns
#'   \code{parameter} (character) and \code{value} (character) containing
#'   key-value pairs from the Snakemake configuration. Pass \code{NULL} (the
#'   default) to omit configuration details.
#' @param generated_date Date. The date the conduit object was generated.
#'   Defaults to today via [Sys.Date()].
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{workflow_version}{character(1)}
#'     \item{generated_date}{Date}
#'     \item{uniprotkb_release}{character(1)}
#'     \item{config}{tbl_df or NULL}
#'   }
#' @export
#'
#' @examples
#' prov <- create_provenance(
#'   workflow_version  = "1.2.3",
#'   uniprotkb_release = "2024_05",
#'   config = tibble::tibble(parameter = "fdr", value = "0.01")
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
    stopifnot(
      is.data.frame(config),
      all(c("parameter", "value") %in% names(config))
    )
  }
  list(
    workflow_version  = workflow_version,
    generated_date    = as.Date(generated_date),
    uniprotkb_release = uniprotkb_release,
    config            = config
  )
}
