#' Plot ROC Curve for Classification Results
#'
#' Creates a Receiver Operating Characteristic (ROC) curve plot for classification model predictions,
#' showing the trade-off between true positive rate and false positive rate.
#'
#' @param predict_classification_list A list containing prediction results, including:
#'   \itemize{
#'     \item test_predictions: A data frame with prediction probabilities and true outcomes for test data
#'     \item training_predictions: A data frame with prediction probabilities and true outcomes for training data
#'     \item outcome: The name of the outcome variable
#'   }
#' @param data_set Character string specifying which dataset to plot. Must be either "test" or "training".
#'   Defaults to "test".
#'
#' @return A ggplot object containing the ROC curve plot with AUC score in the subtitle.
#'
#' @export
#'
#' @examples
#' # Assuming you have a classification model output:
#' # plot_roc(model_predictions, data_set = "test")
#' # plot_roc(model_predictions, data_set = "training")
plot_roc = function(predict_classification_list,
                    data_set = "test"){

  if(data_set == "test"){
  # Generate ROC data
  roc_data <- yardstick::roc_curve(
    data = predict_classification_list$test_predictions,
    truth = predict_classification_list$outcome,
    names(predict_classification_list$test_predictions)[2]  # assumes binary classification with levels 0 and 1, and `.pred_1` is the prob for class 1
  )

  roc_auc <- yardstick::roc_auc(
    data = predict_classification_list$test_predictions,
    truth = predict_classification_list$outcome,
    names(predict_classification_list$test_predictions)[2]  # assumes binary classification with levels 0 and 1, and `.pred_1` is the prob for class 1
  )$.estimate

  } else if(data_set == "training"){
    # Generate ROC data
    roc_data <- yardstick::roc_curve(
      data = predict_classification_list$training_predictions,
      truth = predict_classification_list$outcome,
      names(predict_classification_list$training_predictions)[2]  # assumes binary classification with levels 0 and 1, and `.pred_1` is the prob for class 1
    )

    roc_auc <- yardstick::roc_auc(
      data = predict_classification_list$training_predictions,
      truth = predict_classification_list$outcome,
      names(predict_classification_list$training_predictions)[2]  # assumes binary classification with levels 0 and 1, and `.pred_1` is the prob for class 1
    )$.estimate


  } else{
    stop("data_set must be test or training")

  }

  # Plot via autoplot
  p1 = ggplot2::autoplot(roc_data)+
    ggplot2::labs(
    title = paste("ROC Curve -", data_set, "Data"),
    subtitle = paste("AUC = ", round(roc_auc, 2)
                     )
    )

  return(p1)
}
