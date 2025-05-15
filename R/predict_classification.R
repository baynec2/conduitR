#' Classification Model Prediction with Optional Tuning
#'
#' Fits a classification model using lasso regression, random forest, or XGBoost
#' on a `QFeatures` object to a specified outcome variable in colData.
#' Optionally performs hyperparameter tuning using cross-validation.
#'
#' @param qf A `QFeatures` object.
#' @param assay_name Name of the assay within `qf` to use.
#' @param outcome Column name of the classification outcome variable.
#' @param train_percent Percentage of data to use for training. Default is 70.
#' @param model_type One of `"lasso_regression"`, `"random_forest"`, or `"xgboost"`.
#' @param v Number of folds for cross-validation. Default is 5. If set to 1, tuning is skipped.
#' @param grid_size Number of tuning grid points. Default is 20.
#'
#' @return A list with fitted model, test and training predictions, confusion matrix, and variable importance.
#' @export
#'
#' @examples
#' # predict_classification(qf, "myassay", "Group", model_type = "random_forest")

  predict_classification <- function(qf,
                                     assay_name,
                                     outcome,
                                     train_percent = 70,
                                     model_type = "lasso_regression",
                                     v = 5,
                                     grid_size = 20) {

    supported_models <- c("lasso_regression","random_forest","xgboost")

    if(model_type %in% supported_models){
      print(paste0( model_type," model_type is supported"))
    }else{stop("model_type must be one of ", supported_models)}

    if(v >1 ){tune = TRUE}else if(v == 1){tune = FALSE}else{stop("v must be >= 1")}
    if (!inherits(qf, "QFeatures")) stop("qf must be a QFeatures object.")
    if (!assay_name %in% names(qf)) stop("Invalid assay name.")

    tidy <- conduitR::tidy_conduit(qf, assay_name)

    selected_data <- tidy |>
      dplyr::select(sample, !!rlang::sym(outcome), rowid, value) |>
      tidyr::pivot_wider(names_from = rowid, values_from = value) |>
      dplyr::mutate(!!rlang::sym(outcome) := as.factor(!!rlang::sym(outcome)))


    split <- rsample::initial_split(selected_data,
                                    prop = train_percent / 100,
                                    strata = !!rlang::sym(outcome))
    train <- rsample::training(split)
    test <- rsample::testing(split)

    if (is.numeric(train[[outcome]])) {
      stop("This function only supports classification.")
    }

    rec <- recipes::recipe(as.formula(paste(outcome, "~ .")), data = train) |>
      recipes::update_role(sample, new_role = "id") |>
      recipes::step_zv()

    if (model_type %in% c("lasso_regression", "xgboost")) {
      rec <- rec |> recipes::step_normalize(recipes::all_numeric_predictors())
    }

    # Model spec with or without tuning
    if (tune) {
      model_spec <- switch(model_type,
                           "lasso_regression" = parsnip::logistic_reg(penalty = tune(), mixture = 1) |>
                             parsnip::set_engine("glmnet") |>
                             parsnip::set_mode("classification"),

                           "random_forest" = parsnip::rand_forest(mtry = tune(), min_n = tune(), trees = 500) |>
                             parsnip::set_engine("ranger", importance = "impurity") |>
                             parsnip::set_mode("classification"),

                           "xgboost" = parsnip::boost_tree(
                             trees = 1000,
                             tree_depth = tune(),
                             learn_rate = tune(),
                             mtry = tune(),
                             loss_reduction = tune(),
                             sample_size = tune()
                           ) |>
                             parsnip::set_engine("xgboost") |>
                             parsnip::set_mode("classification")
      )
    } else {
      model_spec <- switch(model_type,
                           "lasso_regression" = parsnip::logistic_reg(penalty = 0.1, mixture = 1) |>
                             parsnip::set_engine("glmnet") |>
                             parsnip::set_mode("classification"),

                           "random_forest" = parsnip::rand_forest(trees = 500) |>
                             parsnip::set_engine("ranger", importance = "impurity") |>
                             parsnip::set_mode("classification"),

                           "xgboost" = parsnip::boost_tree(trees = 100, mtry = 3) |>
                             parsnip::set_engine("xgboost") |>
                             parsnip::set_mode("classification")
      )
    }

    wf <- workflows::workflow() |>
      workflows::add_recipe(rec) |>
      workflows::add_model(model_spec)

    # If tuning is enabled, set up cross-validation and grid search
    if (tune) {
      cv_folds <- rsample::vfold_cv(train, v = v, strata = !!rlang::sym(outcome))

      # Finalize parameters (e.g., mtry) after preparing recipe
      finalized_params <- dials::finalize(
        dials::parameters(wf),
        train
      )

      grid <- dials::grid_space_filling(
        finalized_params,
        size = grid_size
      )

      tuned <- tune::tune_grid(
        wf,
        resamples = cv_folds,
        grid = grid,
        metrics = yardstick::metric_set(yardstick::roc_auc)
      )

      best_params <- tune::select_best(tuned, metric = "roc_auc")
      wf <- tune::finalize_workflow(wf, best_params)
    }

    # Generating the fit
    fit <- parsnip::fit(wf, data = train)

    # Getting predictions on test set
    preds <- dplyr::bind_cols(
      predict(fit, test, type = "prob"),
      predict(fit, test),
      test[, outcome, drop = FALSE]
    )

    # Getting predictions on training set
    training_predictions <- dplyr::bind_cols(
      predict(fit, train, type = "prob"),
      predict(fit, train),
      train[, outcome, drop = FALSE]
    )

    # Return importance
      if (model_type == "random_forest") {
        rf_fit <- workflows::extract_fit_parsnip(fit)
        importance <- tibble::enframe(rf_fit$fit$variable.importance,
                                      name = "feature", value = "importance") |>
          dplyr::arrange(desc(importance))
      } else if (model_type == "xgboost") {
        xgb_fit <- workflows::extract_fit_parsnip(fit)
        importance <- tibble::as_tibble(xgboost::xgb.importance(model = xgb_fit$fit)) |>
          dplyr::rename(feature = Feature, importance = Gain) |>
          dplyr::arrange(desc(importance))
      } else if (model_type == "lasso_regression") {
        lasso_fit <- workflows::extract_fit_parsnip(fit)
        coefs <- as.matrix(coef(lasso_fit$fit, s = lasso_fit$spec$args$penalty))[, 1]
        importance <- tibble::tibble(feature = names(coefs), importance = coefs) |>
          dplyr::filter(feature != "(Intercept)") |>
          dplyr::arrange(desc(abs(importance)))
      }

    # Create confusion matrix
    confusion_matrix <- yardstick::conf_mat(
      data = preds,
      truth = !!rlang::sym(outcome),
      estimate = .pred_class
    )

   output <-  list(
      fit = fit,
      confusion_matrix = confusion_matrix,
      training_predictions = training_predictions,
      test_predictions = preds,
      outcome = outcome,
      test = test,
      importance = importance
    )

   return(output)
  }
