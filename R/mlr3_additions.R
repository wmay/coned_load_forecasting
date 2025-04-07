# extra classes for mlr3

library(mlr3)
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
