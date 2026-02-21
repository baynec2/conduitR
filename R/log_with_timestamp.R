#' Log a Message with a Timestamp
#'
#' Prints a message to the console prefixed with the current time. Uses
#' `sprintf`-style formatting for the message when additional arguments
#' are passed.
#'
#' @param message Character string to log; may include `sprintf` placeholders
#'   (e.g. `"%s"`, `"%d"`).
#' @param ... Values passed to `sprintf(message, ...)` when formatting.
#'
#' @return Invisibly returns nothing; output is printed to the console.
#'
#' @export
#'
#' @examples
#' log_with_timestamp("Processing started")
#' log_with_timestamp("Proteome %s: %d proteins", "UP000005640", 20000)
log_with_timestamp <- function(message, ...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, sprintf(message, ...)))
}