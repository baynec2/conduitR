#' Plot Feature Importance from Classification or Regression Model
#'
#' Creates a horizontal bar plot showing the importance of features (e.g., proteins)
#' from a trained machine learning model. Features are ordered by their importance score.
#'
#' @param predict_classification_list A list containing model results, including:
#'   \itemize{
#'     \item importance: A data frame with columns 'feature' and 'importance' containing
#'       feature importance scores
#'   }
#' @param start Integer specifying the starting rank of features to plot (default: 1)
#' @param end Integer specifying the ending rank of features to plot (default: 10)
#'
#' @return A ggplot object containing a horizontal bar plot with:
#'   \itemize{
#'     \item Y-axis showing feature names, ordered by importance
#'     \item X-axis showing importance scores
#'     \item Horizontal bars representing feature importance
#'   }
#'
#' @export
#'
#' @examples
#' # Plot top 10 most important features:
#' # plot_feature_importance(model_results)
#' 
#' # Plot features ranked 5-15:
#' # plot_feature_importance(model_results, start = 5, end = 15)
#' 
#' # Plot only top 5 features:
#' # plot_feature_importance(model_results, end = 5)
plot_feature_importance <- function(predict_classification_list,
                                    start = 1,
                                    end = 10){


  data = predict_classification_list$importance |>
    dplyr::slice(start:end)


  p1 = data |>
    ggplot(aes(reorder(feature,importance),importance))+
    geom_col()+
    coord_flip()+
    xlab("Feature")

  p1


}
