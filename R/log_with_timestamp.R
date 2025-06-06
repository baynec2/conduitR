#' Log a message with a timestamp
#'
#' This function logs a message with a timestamp to the console.
#'
#' @param message A character string containing the message to log.
#' @param ... Additional arguments to pass to the message formatting.
#'
#' @return None. The function logs the message to the console.
#'
#' @export
log_with_timestamp <- function(message, ...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, sprintf(message, ...)))
}