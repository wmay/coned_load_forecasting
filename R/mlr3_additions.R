# extra classes for mlr3

library(mlr3)
library(mlr3pipelines)
# library(R6)

# CV using contiguous blocks rather than random folds
ResamplingBlockCV = R6::R6Class(
    "ResamplingBlockCV",
    inherit = ResamplingCV,
    private = list(
        .sample = function(ids, ...) {
          n = length(ids)
          data.table(
              row_id = ids,
              fold = floor(self$param_set$values$folds * (seq_len(n) - 1) / n) + 1L,
              key = "fold"
          )
        }
    )
)
mlr_resamplings$add("block_cv", function() ResamplingBlockCV$new())

# average performance of multiple quantiles
MeasureRegrQuantiles = R6::R6Class(
    "MeasureRegrQuantiles",
    inherit = MeasureRegr,
    public = list(
        initialize = function() {
          super$initialize(
                    id = "regr.mean_pinball",
                    range = c(0, Inf),
                    minimize = TRUE,
                    predict_type = "quantiles",
                    man = "mlr3measures::pinball",
                    label = "Mean Pinball Loss"
                )
        }
    ),
    private = list(
        .score = function(prediction, ...) {
          truth = prediction$truth
          quantiles_mat = prediction$quantiles
          levels = colnames(quantiles_mat) %>%
            sub('^q', '', .) %>%
            as.numeric
          pinball_losses = sapply(seq_along(levels), function(i) {
            mlr3measures::pinball(truth, quantiles_mat[, i], alpha = levels[i])
          })
          mean(pinball_losses)
        }
    )
)
mlr_measures$add("regr.mean_pinball", function() MeasureRegrQuantiles$new())

# A model that uses the predictors as the prediction. Useful for benchmarking
# NWP models
LearnerRegrIdent = R6::R6Class(
  "LearnerRegrIdentity",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      super$initialize(
        id = "regr.ident",
        feature_types = c("logical", "integer", "numeric", "factor", "ordered"),
        predict_types = c('response', "distr"),
        packages = character(),
        properties = c("weights", "missings"),
        label = "Identity"
      )
    }
  ),
  private = list(
      .train = function(task) {
        TRUE # no model needed
      },
      .predict = function(task) {
        newdata = task$data(cols = task$feature_names)
        if (self$predict_type == 'response') {
          # it's just the first predictor column
          list(response = as.data.frame(newdata)[, 1])
        } else if (self$predict_type == 'distr') {
          samples = newdata
          distlist = split(samples, 1:nrow(samples)) %>%
            lapply(function(x) {
              if (all(is.na(x))) {
                print(x)
                stop('all samples are missing')
              }
              if (any(is.na(x))) warning('task contains missing samples')
              distr6::Empirical$new(samples = na.omit(x))
            })
          distrs = distr6::VectorDistribution$new(distlist = distlist)
          list(distr = distrs)
        }
      }
  )
)

# distributional regression with ranger
LearnerRegrRangerDist = R6::R6Class(
  "LearnerRegrRangerDist",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      ps = ps(
          always.split.variables       = paradox::p_uty(tags = "train"),
          holdout                      = paradox::p_lgl(default = FALSE, tags = "train"),
          importance                   = paradox::p_fct(c("none", "impurity", "impurity_corrected", "permutation"), tags = "train"),
          keep.inbag                   = paradox::p_lgl(default = FALSE, tags = "train"),
          max.depth                    = paradox::p_int(default = NULL, lower = 1L, special_vals = list(NULL), tags = "train"),
          min.bucket                   = paradox::p_int(1L, default = 1L, tags = "train"),
          min.node.size                = paradox::p_int(1L, default = 5L, special_vals = list(NULL), tags = "train"),
          mtry                         = paradox::p_int(lower = 1L, special_vals = list(NULL), tags = "train"),
          mtry.ratio                   = paradox::p_dbl(lower = 0, upper = 1, tags = "train"),
          na.action                    = paradox::p_fct(c("na.learn", "na.omit", "na.fail"), default = "na.learn", tags = "train"),
          node.stats                   = paradox::p_lgl(default = FALSE, tags = "train"),
          num.random.splits            = paradox::p_int(1L, default = 1L, tags = "train", depends = quote(splitrule == "extratrees")),
          num.threads                  = paradox::p_int(1L, default = 1L, tags = c("train", "predict", "threads")),
          num.trees                    = paradox::p_int(1L, default = 500L, tags = c("train", "predict", "hotstart")),
          oob.error                    = paradox::p_lgl(default = TRUE, tags = "train"),
          poisson.tau                  = paradox::p_dbl(default = 1, tags = "train", depends = quote(splitrule == "poisson")),
          regularization.factor        = paradox::p_uty(default = 1, tags = "train"),
          regularization.usedepth      = paradox::p_lgl(default = FALSE, tags = "train"),
          replace                      = paradox::p_lgl(default = TRUE, tags = "train"),
          respect.unordered.factors    = paradox::p_fct(c("ignore", "order", "partition"), tags = "train"),
          sample.fraction              = paradox::p_dbl(0L, 1L, tags = "train"),
          save.memory                  = paradox::p_lgl(default = FALSE, tags = "train"),
          scale.permutation.importance = paradox::p_lgl(default = FALSE, tags = "train", depends = quote(importance == "permutation")),
          se.method                    = paradox::p_fct(c("jack", "infjack"), default = "infjack", tags = "predict"), # FIXME: only works if predict_type == "se". How to set dependency?
          seed                         = paradox::p_int(default = NULL, special_vals = list(NULL), tags = c("train", "predict")),
          split.select.weights         = paradox::p_uty(default = NULL, tags = "train"),
          splitrule                    = paradox::p_fct(c("variance", "extratrees", "maxstat", "beta", "poisson"), default = "variance", tags = "train"),
          verbose                      = paradox::p_lgl(default = TRUE, tags = c("train", "predict")),
          write.forest                 = paradox::p_lgl(default = TRUE, tags = "train")
      )
      ps$set_values(num.threads = 1L)
      super$initialize(
        id = "regr.rangerdist",
        feature_types = c("logical", "integer", "numeric", "factor", "ordered"),
        predict_types = "distr",
        packages = "ranger",
        param_set = ps,
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
        pv$case.weights = task$weights$weight
      }
      mlr3misc::invoke(
          ranger::ranger,
          formula = task$formula(),
          data = task$data(),
          quantreg = TRUE,
          .args = pv
      )
    },
    .predict = function(task) {
      pv = self$param_set$get_values(tags = "predict")
      newdata = task$data(cols = task$feature_names)
      pred = predict(self$model, newdata, type = "quantiles",
                     what = function(x) c(mean(x), sd(x)))
      # need to wrap in a `VectorDistribution` (see `LearnerRegr`)
      params = as.data.frame(pred$predictions)
      names(params) = c('mean', 'sd')
      distrs = distr6::VectorDistribution$new(distribution = "Normal",
                                              params = params)
      list(distr = distrs)
    }
  )
)

PipeOpSubtractColumn = R6::R6Class('PipeOpSubtractColumn',
  inherit = mlr3pipelines::PipeOpTargetTrafo,
  public = list(
      initialize = function(id = 'subtract.column', param_vals = list()) {
      ps = ps(
          column = paradox::p_uty(tags = c("train", "predict"))
      )
      # param_vals = list(column = column)
      # ps$values = list(trafo = identity, inverter = identity)
      super$initialize(id = id, param_set = ps, param_vals = param_vals)
    }
  ),
  private = list(
    .transform = function(task, phase) {
      new_target = task$data(cols = task$target_names) -
        task$data(cols = self$param_set$values$column)
      # if (!is.data.frame(new_target) && !is.matrix(new_target)) {
      #   stopf("Hyperparameter 'trafo' must be a function returning a 'data.frame', 'data.table', or 'matrix', not '%s'.", class(new_target)[[1L]])
      # }
      task$cbind(new_target)
      convert_task(task, target = colnames(new_target), drop_original_target = TRUE)
    },
    .train_invert = function(task) {
      # return a predict_phase_state object (can be anything)
      list(truth = task$truth(),
           subtracted = task$data(cols = self$param_set$values$column))
    },
    .invert = function(prediction, predict_phase_state) {
      type = prediction$task_type
      # need to adjust the means of the predicted distributions
      means = unlist(prediction$distr$getParameterValue('mean'))
      means = means + predict_phase_state$subtracted[[1]]
      sds = unlist(prediction$distr$getParameterValue('sd'))
      # need to wrap in a `VectorDistribution` (see `LearnerRegr`)
      params = data.frame(mean = means, sd = sds)
      distrs = distr6::VectorDistribution$new(distribution = "Normal",
                                              params = params)
      new_pred = list(distr = distrs)
      mlr3misc::invoke(get(mlr_reflections$task_types[type, mult = "first"]$prediction)$new,
                       row_ids = prediction$row_ids,
                       truth = predict_phase_state$truth, .args = new_pred)
    }
  )
)
