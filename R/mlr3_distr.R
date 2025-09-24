# mlr3 classes related to distributional regression

library(mlr3)
# library(R6)

# example learner, ngr from crch package
LearnerRegrNGR = R6::R6Class(
  "LearnerRegrNGR",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      param_set = paradox::ps(
        cp             = paradox::p_dbl(0, 1, default = 0.01, tags = "train"),
        maxcompete     = paradox::p_int(0L, default = 4L, tags = "train"),
        maxdepth       = paradox::p_int(1L, 30L, default = 30L, tags = "train"),
        maxsurrogate   = paradox::p_int(0L, default = 5L, tags = "train"),
        minbucket      = paradox::p_int(1L, tags = "train"),
        minsplit       = paradox::p_int(1L, default = 20L, tags = "train"),
        surrogatestyle = paradox::p_int(0L, 1L, default = 0L, tags = "train"),
        usesurrogate   = paradox::p_int(0L, 2L, default = 2L, tags = "train"),
        xval           = paradox::p_int(0L, default = 10L, tags = "train")
      )
      param_set$set_values(xval = 10L)
      super$initialize(
        id = "regr.ngr",
        feature_types = c("logical", "integer", "numeric", "factor", "ordered"),
        predict_types = "distr",
        packages = "crch",
        param_set = param_set,
        properties = c("weights", "missings"),
        label = "Non-homogeneous Gaussian Regression"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pv = self$param_set$get_values(tags = "train")
      if ("weights" %in% task$properties) {
        pv$weights = task$weights$weight
      }
      # crch::crch(y ~ . | ., data = cbind(x, y = y), method = 'boosting', type = 'crps',
      #      maxit = 1000, nu = 0.05, mstop = 'aic')
      # invoke(
      #   crch::crch,
      #   formula = task$formula(),
      #   data = task$data(),
      #   .args = pv
      # )
      if ('formula' %in% names(task$extra_args)) {
        formula = task$extra_args$formula
      } else {
        formula = task$formula()
      }
      crch::crch(formula, data = task$data(), method = 'boosting',
                 type = 'crps', maxit = 1000, nu = 0.05, mstop = 'aic')
    },
    .predict = function(task) {
      pv = self$param_set$get_values(tags = "predict")

      # ensure same column order in train and predict
      # newdata = mlr3extralearners:::ordered_features(task, self)
      # response = mlr3misc::invoke(predict, self$model, newdata = task$data(), .args = pv)
      # need to wrap in a `VectorDistribution` (see `LearnerRegr`)
      params = mlr3misc::invoke(predict, self$model, newdata = task$data(),
                                .args = c(type = 'parameter', pv))
      names(params) = c('mean', 'sd')
      distrs = distr6::VectorDistribution$new(distribution = "Normal",
                                              params = params)
      # list(response = unname(response), distr = distrs)
      list(distr = distrs)
    }
  )
)

# # quick test
# # task = tsk("mtcars", formula = mpg ~ . | .)
# # task = as_task_regr(x = mpg ~ . | ., mtcars)
# # task = as_task_regr(x = mpg ~ ., mtcars)
# # task = as_task_regr(mtcars, target = 'mpg', formula = mpg ~ . | .)
# # task = as_task_regr(mtcars, target = 'mpg', extra_args = list(formula = mpg ~ . | .))
# # task = Task$new('mtcars', 'regr', mtcars, extra_args = list(formula = mpg ~ . | .))
# task = TaskRegr$new('mtcars', mtcars, 'mpg', extra_args = list(formula = mpg ~ . | .))
# learner = LearnerRegrNGR$new(predict_type = 'distr')
# learner$train(task)
# p = learner$predict(task)

# # get distribution type from prediction
# p$distr$modelTable$Distribution[1]

# # get parameters
# unlist(p$distr$getParameterValue('mean'))

# p$score(msr("regr.mse"))

LearnerRegrDrf = R6::R6Class(
  "LearnerRegrDrf",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      param_set = paradox::ps(
          num.trees                    = paradox::p_int(1L, default = 500L, tags = c("train", "predict", "hotstart")),
          # num.features ?
          sample.fraction              = paradox::p_dbl(0L, 1L, tags = "train"),
          mtry                         = paradox::p_int(lower = 1L, special_vals = list(NULL), tags = "train"),
          mtry.ratio                   = paradox::p_dbl(lower = 0, upper = 1, tags = "train"),
          min.node.size                = paradox::p_int(1L, default = 15L, special_vals = list(NULL), tags = "train"),
          honesty                      = paradox::p_lgl(default = TRUE, tags = "train"),
          honesty.fraction             = paradox::p_dbl(lower = 0, upper = 1, default = .5, tags = "train"),
          honesty.prune.leaves         = paradox::p_lgl(default = TRUE, tags = "train"),
          alpha                        = paradox::p_dbl(lower = 0, upper = .25, default = .05, tags = "train"),
          ci.group.size                = paradox::p_int(default = 2L, tags = "train"),
          num.threads                  = paradox::p_int(lower = 1L, special_vals = list(NULL), tags = "train")
      )
      super$initialize(
        id = "regr.drf",
        feature_types = c("logical", "integer", "numeric", "factor", "ordered"),
        predict_types = "distr",
        packages = "drf",
        param_set = param_set,
        properties = c("weights", "missings"),
        label = "Non-homogeneous Gaussian Regression"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pv = self$param_set$get_values(tags = "train")
      pv = mlr3learners:::convert_ratio(pv, "mtry", "mtry.ratio", length(task$feature_names))
      if ("weights" %in% task$properties) {
        pv$sample.weights = task$weights$weight
      }
      # mlr3misc::invoke(
      #     drf::drf,
      #     X = task$data(cols = task$feature_names),
      #     Y = task$data(cols = task$target_names),
      #     .args = pv
      #     )
      # drf::drf(X = task$data(cols = task$feature_names),
      #          Y = task$data(cols = task$target_names))
      args = list(X = task$data(cols = task$feature_names),
                  Y = task$data(cols = task$target_names))
      args = c(args, pv)
      tryCatch({
        do.call(drf::drf, args)
      }, error = function(e) {
        print('Error message:')
        print(as.character(e))
        warning(as.character(e))
      })
    },
    .predict = function(task) {
      pv = self$param_set$get_values(tags = "predict")
      newdata = task$data(cols = task$feature_names)
      n_pred = nrow(newdata)
      if (inherits(self$model, 'character')) {
        # print('Model is an error message:')
        # print(self$model)
        train_params = self$param_set$get_values(tags = "train")
        param_txt = paste(capture.output(print(train_params)), collapse=' ')
        warning_txt = paste0('Model is an error message. Parameters: ', param_txt)
        warning(warning_txt)
        # for now just return bad predictions
        params = data.frame(mean = rep(0, nrow(newdata)),
                            sd = sqrt(.Machine$double.xmax))
        distrs = distr6::VectorDistribution$new(distribution = "Normal",
                                                params = params)
        return(list(distr = distrs))
      }
      means = predict(self$model, newdata = newdata, functional = 'mean')
      sds = predict(self$model, newdata = newdata, functional = 'sd')
      # need to wrap in a `VectorDistribution` (see `LearnerRegr`)
      params = data.frame(mean = means, sd = sds)
      distrs = distr6::VectorDistribution$new(distribution = "Normal",
                                              params = params)
      list(distr = distrs)
    }
  )
)

# average performance of multiple quantiles
MeasureRegrCRPS = R6::R6Class(
    "MeasureRegrCRPS",
    inherit = MeasureRegr,
    public = list(
        initialize = function() {
          super$initialize(
                    id = "regr.crps",
                    packages = 'scoringRules',
                    range = c(0, Inf),
                    minimize = TRUE,
                    predict_type = 'distr',
                    label = "Mean CRPS"
                )
        }
    ),
    private = list(
        .score = function(prediction, ...) {
          # Note that the predictions will be out of order if used after
          # `Prediction$filter()`. See
          # https://github.com/mlr-org/mlr3/issues/1400
          dist_types = unique(prediction$distr$modelTable$Distribution)
          if (length(dist_types) > 1) stop('All predictions must use the same distribution type')
          
          # first, get the type of distribution or sampling
          dist_type = prediction$distr$modelTable$Distribution[1]

          # then, use the corresponding CRPS method
          if (dist_type == 'Normal') {
            # check for wrong order
            resp1 = prediction$response
            resp2 = unlist(prediction$distr$getParameterValue('mean'))
            in_order = all(resp1 == resp2)
            if (!in_order) stop('Inconsistent responses, likely due to `filter` bug')
            mean = unlist(prediction$distr$getParameterValue('mean'))
            sd = unlist(prediction$distr$getParameterValue('sd'))
            pred_crps = scoringRules::crps_norm(prediction$truth, mean = mean, sd = sd)
          } else if (dist_type == 'Empirical') {
            sample_mat = prediction$distr$wrappedModels() %>%
              lapply(function(x) x$parameters()$values$data$samples) %>%
              do.call(rbind, .)
            pred_crps = scoringRules::crps_sample(prediction$truth, sample_mat)
          } else {
            stop(paste0('Distribution type ', dist_type, ' not implemented'))
          }
          mean(pred_crps)
        }
    )
)
mlr_measures$add("regr.crps", function() MeasureRegrCRPS$new())

# # quick test
# p$score(msr("regr.crps"))

# needed to scale the Y values
PipeOpTargetTrafoScaleRangeDistr = R6::R6Class("PipeOpTargetTrafoScaleRangeDistr",
  inherit = PipeOpTargetTrafoScaleRange,
  public = list(),
  private = list(
    .invert = function(prediction, predict_phase_state) {
      orig_means = unlist(prediction$distr$getParameterValue('mean'))
      orig_sds = unlist(prediction$distr$getParameterValue('sd'))
      new_means = (orig_means - self$state$offset) / self$state$scale
      new_sds = orig_sds / self$state$scale
      distrs = data.frame(mean = new_means, sd = new_sds) %>%
        distr6::VectorDistribution$new(distribution = "Normal",
                                       params = .)
      PredictionRegr$new(row_ids = prediction$row_ids,
                         truth = predict_phase_state$truth, distr = distrs)
    }
  )
)
