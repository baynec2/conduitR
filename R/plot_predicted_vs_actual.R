#' Plot Predicted vs Actual Values for Regression Results
#'
#' Creates a scatter plot of predicted versus actual values from a regression model,
#' with a diagonal reference line representing perfect prediction and model
#' performance metrics in the subtitle.
#'
#' @param predict_regression_list A list containing prediction results, including:
#'   \itemize{
#'     \item test_predictions: A data frame with `.pred` and the actual outcome column for test data
#'     \item training_predictions: A data frame with `.pred` and the actual outcome column for training data
#'     \item outcome: The name of the outcome variable
#'   }
#' @param data_set Character string specifying which dataset to plot. Must be either
#'   "test" or "training". Defaults to "test".
#'
#' @return A ggplot object containing a scatter plot with:
#'   \itemize{
#'     \item X-axis showing actual values
#'     \item Y-axis showing predicted values
#'     \item Points for each observation
#'     \item A dashed diagonal reference line (perfect prediction)
#'     \item RMSE and R-squared in the subtitle
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- predict_regression(qf, "protein_groups", "concentration")
#' plot_predicted_vs_actual(res, data_set = "test")
#' plot_predicted_vs_actual(res, data_set = "training")
#' }
plot_predicted_vs_actual <- function(predict_regression_list,
                                     data_set = "test") {

  if (data_set == "test") {
    preds <- predict_regression_list$test_predictions
  } else if (data_set == "training") {
    preds <- predict_regression_list$training_predictions
  } else {
    stop("data_set must be 'test' or 'training'")
  }

  outcome <- predict_regression_list$outcome

  rmse_val <- yardstick::rmse(preds, truth = !!rlang::sym(outcome), estimate = .pred)$.estimate
  rsq_val  <- yardstick::rsq(preds,  truth = !!rlang::sym(outcome), estimate = .pred)$.estimate

  p1 <- ggplot2::ggplot(preds, ggplot2::aes(x = !!rlang::sym(outcome), y = .pred)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::labs(
      title    = paste("Predicted vs Actual -", data_set, "Data"),
      subtitle = paste0("RMSE = ", round(rmse_val, 3), "   R\u00b2 = ", round(rsq_val, 3)),
      x        = paste("Actual", outcome),
      y        = paste("Predicted", outcome)
    )

  return(p1)
}
