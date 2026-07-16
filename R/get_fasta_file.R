#' Download Proteome FASTA File from UniProt
#'
#' Downloads a proteome FASTA file from UniProt using their REST API. This function
#' retrieves protein sequences for a specified proteome ID and saves them to a local
#' file. The downloaded file contains protein sequences with UniProt-style headers
#' including protein IDs, organism information, and other metadata.
#'
#' Sequences are fetched from the paginated `/search` endpoint (following the
#' `Link: rel="next"` cursor) rather than `/stream`. `/search` reports an
#' `X-Total-Results` count per page and lets a transient failure retry a single
#' page instead of the whole proteome, so completeness can be validated inline:
#' a download is accepted only if the delivered sequence count matches the
#' endpoint's `X-Total-Results` and is not grossly short of the proteome's
#' `proteinCount` (from `/proteomes/{id}`). UniProtKB is tried first; if it is
#' empty or short, the UniParc `/search` endpoint is tried as a fallback. This
#' replaces the previous `/stream` implementation, whose empty/error 200
#' responses could silently yield an empty or truncated database (see issue #20).
#'
#' UniParc `/search` deflines are bare (`>UPI... status=active`) with no organism
#' annotation. Because that would leave every UniParc protein without a taxon
#' downstream, each defline is stamped with `OS=`/`OX=` taken from the proteome's
#' record, inserted after the leading `UPI` accession so the first-token
#' `protein_id` (used by DIA-NN and `extract_fasta_info()`) is preserved.
#'
#' @param proteome_id Character string specifying the UniProt proteome identifier
#'   (e.g., "UP000005640" for human proteome). These IDs can be found in the
#'   UniProt proteomes database.
#' @param fasta_dir Character string specifying the directory where the FASTA file
#'   should be saved. Defaults to the current working directory.
#' @param max_tries Integer. Maximum attempts per HTTP request (each page is
#'   retried independently) before giving up. Default 5.
#' @param wait_seconds Numeric. Base backoff between retries, scaled linearly by
#'   attempt number (attempt `i` waits `wait_seconds * i`). Default 30.
#' @param completeness_tol Numeric in (0, 1]. A download is accepted only if the
#'   delivered sequence count is at least `completeness_tol` times the proteome's
#'   `proteinCount`. Guards against silent truncation while tolerating minor
#'   count drift (isoforms/redundancy). Default 0.9. Ignored when `proteinCount`
#'   is unavailable, in which case the endpoint's own `X-Total-Results` is the
#'   sole completeness signal.
#'
#' @return A tibble containing download information with the following columns:
#'   \itemize{
#'     \item proteome_id: The input proteome identifier
#'     \item resp_status: HTTP response status code (200 for success)
#'     \item source: Either "uniprotkb", "uniparc", or "not_downloaded"
#'     \item n_sequences: Number of sequences delivered (0 if none)
#'     \item expected: Expected protein count from the UniProt proteomes record
#'       (`NA` if it could not be fetched)
#'   }
#'   The function also saves a FASTA file named \code{<proteome_id>.fasta} in the
#'   specified directory when a complete proteome is found in UniProtKB or
#'   UniParc. A partial or empty download writes no file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Download human proteome to current directory
#' result <- get_fasta_file("UP000005640")
#'
#' # Download to a specific directory
#' result <- get_fasta_file("UP000005640", fasta_dir = "data/proteomes")
#'
#' # Check download status
#' if (result$source != "not_downloaded") {
#'   message("Download successful: ", result$n_sequences, " sequences")
#' } else {
#'   message("Download failed")
#' }
#' }
#'
#' @note
#' This function:
#' \itemize{
#'   \item Requires an internet connection to access UniProt
#'   \item Uses the UniProt REST API
#'   \item Creates a new file or overwrites existing files
#'   \item Handles API errors with informative messages
#' }
#' The UniProt API is free to use but has rate limits. For bulk downloads,
#' consider implementing appropriate delays between requests.
get_fasta_file <- function(proteome_id,
                           fasta_dir = getwd(),
                           max_tries = 5,
                           wait_seconds = 30,
                           completeness_tol = 0.9) {
  # Log the proteome id and the time to console.
  log_with_timestamp(paste0("Processing proteome ID: ", proteome_id))

  fasta_fp <- file.path(fasta_dir, paste0(proteome_id, ".fasta"))

  # Count FASTA headers (records) in a body of concatenated FASTA text.
  count_fasta_headers <- function(txt) {
    if (!nzchar(txt)) return(0L)
    m <- gregexpr("(?m)^>", txt, perl = TRUE)[[1]]
    if (length(m) == 1L && m[1] == -1L) 0L else length(m)
  }

  # Shared retry/timeout policy: retry transient HTTP statuses AND connection-
  # level failures (retry_on_failure), with linear backoff. Applied per page,
  # so a single flaky page is retried rather than the whole proteome.
  add_resilience <- function(req) {
    req |>
      httr2::req_timeout(seconds = 60) |>
      httr2::req_retry(
        max_tries = max_tries,
        backoff = function(i) wait_seconds * i,
        retry_on_failure = TRUE,
        is_transient = function(resp) {
          httr2::resp_status(resp) %in% c(429, 500, 502, 503, 504)
        }
      )
  }

  # Proteome record: protein count (to catch silent truncation, e.g. 5 of 6064)
  # plus taxon id and organism name (to backfill OS=/OX= onto bare UniParc
  # deflines below). All fields optional; missing values degrade to NA. Expected
  # count NA -> the endpoint's own X-Total-Results is the sole completeness signal.
  fetch_proteome_meta <- function() {
    na <- list(expected = NA_integer_, taxon_id = NA_integer_,
               organism_name = NA_character_)
    tryCatch({
      resp <- httr2::request(paste0("https://rest.uniprot.org/proteomes/", proteome_id)) |>
        add_resilience() |>
        httr2::req_error(is_error = function(r) FALSE) |>
        httr2::req_perform()
      if (httr2::resp_status(resp) != 200) return(na)
      j <- httr2::resp_body_json(resp)
      pc <- j$proteinCount
      tx <- j$taxonomy$taxonId
      on <- j$taxonomy$scientificName
      list(
        expected      = if (is.null(pc)) NA_integer_ else as.integer(pc),
        taxon_id      = if (is.null(tx)) NA_integer_ else as.integer(tx),
        organism_name = if (is.null(on)) NA_character_ else as.character(on)
      )
    }, error = function(e) na)
  }

  # UniParc /search deflines are bare (">UPI... status=active") with no organism
  # annotation, so extract_fasta_info() cannot map the protein to a taxon and all
  # UniParc taxonomy is lost downstream. Restore OS=/OX= from the proteome record,
  # inserting them AFTER the leading UPI accession so the first-token protein_id
  # (what DIA-NN and extract_fasta_info key on) is unchanged. A proteome-scoped
  # download shares one taxon, so a single taxon_id/organism_name applies to every
  # defline. Organism scientific names contain no regex-replacement metacharacters.
  stamp_uniparc_headers <- function(txt, taxon_id, organism_name) {
    if (is.na(taxon_id)) return(txt)  # no taxon to stamp; leave body unchanged
    ins <- if (!is.na(organism_name) && nzchar(organism_name)) {
      sprintf(" OS=%s OX=%s", organism_name, taxon_id)
    } else {
      sprintf(" OX=%s", taxon_id)
    }
    gsub("(?m)^(>\\S+)", paste0("\\1", ins), txt, perl = TRUE)
  }

  # Paginated FASTA fetch from a /search endpoint, following the rel="next"
  # cursor to completion. max_reqs = Inf is REQUIRED: httr2 defaults to 20,
  # which would silently cap large proteomes at 20 * 500 = 10,000 sequences.
  fetch_search_fasta <- function(base_url) {
    req <- httr2::request(base_url) |>
      httr2::req_url_query(
        query = paste0("proteome:", proteome_id),
        format = "fasta",
        size = 500
      ) |>
      add_resilience()

    resps <- tryCatch(
      httr2::req_perform_iterative(
        req,
        next_req = httr2::iterate_with_link_url(rel = "next"),
        max_reqs = Inf,
        on_error = "return"
      ),
      error = function(e) {
        log_with_timestamp("Search request failed (%s): %s", proteome_id, conditionMessage(e))
        list()
      }
    )

    ok <- httr2::resps_successes(resps)
    if (length(ok) == 0) {
      return(list(status = NA_integer_, n = 0L, total = NA_integer_, body = ""))
    }
    total <- suppressWarnings(as.integer(httr2::resp_header(ok[[1]], "x-total-results")))
    body <- paste0(
      vapply(ok, function(r) {
        if (httr2::resp_has_body(r)) httr2::resp_body_string(r) else ""
      }, character(1)),
      collapse = ""
    )
    list(
      status = httr2::resp_status(ok[[length(ok)]]),
      n = count_fasta_headers(body),
      total = total,
      body = body
    )
  }

  proteome_meta <- fetch_proteome_meta()
  expected <- proteome_meta$expected

  result_row <- function(status, source, n) {
    tibble::tibble(
      proteome_id = proteome_id,
      resp_status = as.integer(status),
      source = source,
      n_sequences = as.integer(n),
      expected = as.integer(expected)
    )
  }

  # Accept a result only if it delivered everything the endpoint reports
  # (n >= X-Total-Results) AND is not grossly short of the proteomes record's
  # proteinCount. Either signal being absent (NA) falls back to the other.
  is_complete <- function(res) {
    res$n > 0 &&
      (is.na(res$total) || res$n >= res$total) &&
      (is.na(expected) || res$n >= completeness_tol * expected)
  }

  # 1) UniProtKB paginated search
  kb <- fetch_search_fasta("https://rest.uniprot.org/uniprotkb/search")
  if (is_complete(kb)) {
    writeLines(kb$body, fasta_fp)
    log_with_timestamp("Proteome %s: %d sequences downloaded (UniProtKB)", proteome_id, kb$n)
    return(result_row(kb$status, "uniprotkb", kb$n))
  }
  if (kb$n > 0) {
    log_with_timestamp(
      "Proteome %s: UniProtKB delivered %d of %s expected (incomplete); trying UniParc",
      proteome_id, kb$n, ifelse(is.na(expected), "?", as.character(expected))
    )
  }

  # 2) UniParc paginated search (fallback). Bare deflines are stamped with the
  #    proteome's OS=/OX= so extract_fasta_info() can recover the taxon.
  up <- fetch_search_fasta("https://rest.uniprot.org/uniparc/search")
  if (is_complete(up)) {
    up_body <- stamp_uniparc_headers(up$body, proteome_meta$taxon_id,
                                     proteome_meta$organism_name)
    writeLines(up_body, fasta_fp)
    log_with_timestamp("Proteome %s: %d sequences downloaded (UniParc)", proteome_id, up$n)
    return(result_row(up$status, "uniparc", up$n))
  }

  # 3) Neither source produced a complete proteome. Report the shortfall
  #    loudly and write nothing -- the caller's completeness gate decides
  #    whether the overall database is usable.
  best_n <- max(kb$n, up$n)
  log_with_timestamp(
    "Proteome %s: FAILED -- best delivered %d of %s expected (UniProtKB=%d, UniParc=%d)",
    proteome_id, best_n, ifelse(is.na(expected), "?", as.character(expected)), kb$n, up$n
  )
  result_row(
    if (!is.na(kb$status)) kb$status else up$status,
    "not_downloaded",
    best_n
  )
}
