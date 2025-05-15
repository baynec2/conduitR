#' plot feature importance
#'
#' @param predict_classification_list
#' @param start
#' @param end
#'
#' @returns
#' @export
#'
#' @examples
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
