# Evaluate daily peak load forecasting methods

# pkg_deps = c('ncmeta', 'scoringRules')
# install.packages(pkg_deps)
# remotes::install_github("mlr-org/mlr3temporal")
# distr6_repos =  c(CRAN = 'https://cloud.r-project.org',
#                   raphaels1 = 'https://raphaels1.r-universe.dev')
# install.packages('distr6', repos = distr6_repos)
# install.packages("mlr3proba", repos = "https://mlr-org.r-universe.dev")
# library(timeDate) # holidays
library(magrittr)
library(stars) # also requires ncmeta
library(future)
library(future.apply)
library(mlr3)
library(torch)
library(mlr3torch)
# library(mlr3learners)
library(mlr3tuning)
library(mlr3pipelines)
library(mlr3temporal)
# library(mlr3proba)
# library(ggplot2)
source('R/mlr3_additions.R') # block CV, mean pinball loss, etc.
source('R/mlr3_distr.R')

# see notes in `?benchmark` and `?mlr_tuners_random_search`
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

all_eff_tmps = read.csv('results/process_station_data/eff_tmp.csv')
all_loads = read.csv('results/process_load_data/loads.csv')

# if coordinates are not in the same order for every variable, `read_mdim` works
# fine while `read_ncdf` messes up the coordinates
gefs_3day = read_ncdf('results/process_nwp_data/gefs_3day_wmean.nc')

# get values from stars object at the centroid of a network
networks = readRDS('results/maps/coned_networks_cleaned.rds')
# add centroid coordinates
network_coords = networks %>%
  st_centroid %>%
  st_coordinates %>%
  as.data.frame
row.names(network_coords) = networks$id
names(network_coords) = c('nad83_x', 'nad83_y')
stations = readRDS('results/station_data/stations.rds')

# note that the day returned here is the valid day
extract_values = function(network, days_ahead) {
  network_idx = which(attr(gefs_3day, 'dimensions')$network$values == network)
  day_idx = days_ahead + 1
  gefs_3day %>%
    dplyr::slice('edt9pm_day', day_idx) %>%
    dplyr::slice('network', network_idx) %>%
    as.data.frame %>%
    # getElement('time') %>%
    # class
    transform(valid_day = as.Date(time) + days_ahead) %>%
    subset(select = -c(time, network))
}

# add observed effective temps to short-term TV forecasts
fill_incomplete_tv = function(tv, day, network, days_ahead) {
  # return(tv)
  eff_tmps_col =
    if (network == 'system') 'system.new' else paste0('network.', network)
  system_eff_tmp = all_eff_tmps[, c('day', eff_tmps_col)]
  names(system_eff_tmp)[2] = 'eff_temp'
  if (days_ahead == 1) {
    tv + .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  } else if (days_ahead == 0) {
    tv + .2 * system_eff_tmp$eff_temp[match(day - 1, system_eff_tmp$day)] +
      .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  }
}

prepare_load_task = function(network, days_ahead) {
  # Note: netcdf files use forecast day, while `system_tv` uses valid day. Must
  # be consistent!
  is_system = startsWith(network, 'system')
  extract_name = if (is_system) 'system' else network
  out = extract_values(extract_name, days_ahead) %>%
    subset(select = -eff_temp) %>% # 3-day eff_temp is same as TV
    # the time measure, to allow for a nonstationary model
    transform(fct_run = as.integer(valid_day - days_ahead))
  # sometimes soilw is all NA? for now just skip it if so
  if (all(is.na(out$soilw))) {
    # warning('All soilw values are missing')
    out$soilw = NULL
  }
  # For days_ahead < 2, need to add observed (past) portion of TV to the
  # forecast portion. TV_sd is correct though, since observations have 0 sd
  if (days_ahead < 2) {
    out$TV = with(out, fill_incomplete_tv(TV, valid_day, network, days_ahead))
  }
  out$doy = as.POSIXlt(out$valid_day)$yday
  # get load data
  load_col = paste0('reading.', network)
  id = paste('load', network, days_ahead, sep = '_')
  out %>%
    transform(load = all_loads[match(valid_day, all_loads$day), load_col]) %>%
    subset(select = -valid_day) %>%
    na.omit %>%
    as_task_regr('load', id = id)
}

# create a task that combines lead times, as in https://doi.org/10.1002/qj.4701
prepare_load_task_combined = function(network) {
  task = prepare_load_task(network, 0)
  task$cbind(data = data.frame(days_ahead = rep(0, task$nrow)))
  # append the other lead times
  for (i in 1:7) {
    task_i = prepare_load_task(network, i)
    task_i$cbind(data = data.frame(days_ahead = rep(i, task_i$nrow)))
    task$rbind(data = task_i$data())
  }
  task$id = paste('load', network, sep = '_')
  # make sure resampling is grouped by day_ahead
  task$set_col_roles('days_ahead', add_to = 'stratum')
  task
}

prepare_load_task_all_networks = function(network_ids = networks$id) {
  network_tasks = network_ids %>%
    lapply(prepare_load_task_combined) %>%
    lapply(function(x) {
      # add the network as a predictor
      network_id = sub('load_', '', x$id)
      # remove soilw since some networks don't have it
      if ('soilw' %in% x$feature_names) {
        x$select(setdiff(x$feature_names, 'soilw'))
      }
      x$cbind(data.frame(network = rep(network_id, x$nrow)))
    })
  out_task = network_tasks[[1]]
  for (i in 2:length(network_tasks)) {
    out_task$rbind(network_tasks[[i]]$data())
  }
  # add NAD83 coordinates as predictors
  out_task$cbind(network_coords[out_task$data()$network, ])
  # convert network column to factor
  network_fct = factor(out_task$data()$network, levels = network_ids)
  out_task$select(out_task$feature_names[out_task$feature_names != 'network'])
  out_task$cbind(data.frame(network = network_fct))
  out_task$id = paste('network_loads')
  out_task
}

# Adding distributional regression, and fixing a bug-- see
# https://github.com/mlr-org/mlr3pipelines/issues/989
PipeOpTargetTrafoScaleRangeDistr = R6::R6Class("PipeOpTargetTrafoScaleRangeDistr",
  inherit = PipeOpTargetTrafoScaleRange,
  public = list(
    initialize = function(id = "targettrafoscalerangedistr", param_vals = list()) {
      super$initialize(id = id, param_vals = list())
    }
  ),
  private = list(
    .transform = function(task, phase) {
      x = task$data(cols = task$target_names)
      new_target = self$state$offset + x * self$state$scale
      data.table::setnames(new_target, paste0(colnames(new_target), ".scaled"))
      task$cbind(new_target)
      # address the internal valid task
      if (!is.null(task$internal_valid_task)) {
        x = task$internal_valid_task$data(cols = task$internal_valid_task$target_names)
        new_target = self$state$offset + x * self$state$scale
        data.table::setnames(new_target, paste0(colnames(new_target), ".scaled"))
        task$internal_valid_task$cbind(new_target)
      }
      convert_task(task, target = colnames(new_target), drop_original_target = TRUE, drop_levels = FALSE)
    },
    .invert = function(prediction, predict_phase_state) {
      orig_means = unlist(prediction$distr$getParameterValue('mean'))
      orig_sds = unlist(prediction$distr$getParameterValue('sd'))
      offsets = self$state$offset
      scales = self$state$scale
      # new_means = (orig_means / scales) - offsets
      new_means = (orig_means - offsets) / scales
      new_sds = orig_sds / scales
      # new_truth = (prediction$truth * scales) + offsets
      # new_truth = (prediction$truth - offsets) / scales
      distrs = data.frame(mean = new_means, sd = new_sds) %>%
        distr6::VectorDistribution$new(distribution = "Normal",
                                       params = .)
      PredictionRegr$new(row_ids = prediction$row_ids,
                         truth = predict_phase_state$truth, distr = distrs)
    }
  )
)

categ_scale_target = function(task, scale_params) {
  old_target = task$data(cols = task$target_names)
  network = task$data()$network
  # scale each network's loads
  new_target = (old_target - scale_params[network, 'mean']) /
    scale_params[network, 'sd']
  # I'm honestly confused about how the duplicate columns are handled but it
  # seems to work
  task$cbind(new_target)
  convert_task(task, target = colnames(new_target), drop_original_target = TRUE,
               drop_levels = FALSE)
}

LearnerTorchModuleDistr = R6::R6Class("LearnerTorchModuleDistr",
  inherit = LearnerTorchModule,
  public = list(
      #packages = 'torch',
      initialize = function(module_generator = NULL, param_set = NULL,
                            ingress_tokens = NULL, task_type, properties = NULL,
                            optimizer = NULL, loss = NULL, callbacks = list(),
                            packages = character(0), feature_types = NULL,
                            predict_types = NULL) {
        super$initialize(module_generator = module_generator,
                         param_set = param_set, ingress_tokens = ingress_tokens,
                         task_type = task_type, properties = properties,
                         optimizer = optimizer, loss = loss,
                         callbacks = callbacks, packages = packages,
                         feature_types = feature_types,
                         predict_types = predict_types)
      }
  ),
  private = list(
      .encode_prediction = function(predict_tensor, task) {
        # see `encode_prediction_default` for example
        if (self$predict_type == "response") {
          list(response = as.numeric(predict_tensor[, 1]))
        } else if (self$predict_type == "distr") {
          predict_mat = predict_tensor %>%
            as.numeric %>%
            matrix(ncol = 2, dimnames = list(NULL, c('mean', 'sd')))
          # check the numbers for problems
          if (any(is.na(predict_mat))) {
            # print(tensor_mat)
            # stop('NA tensor values')
            warning('NA tensor values will be replaced with zero (mean) and max double (sd)')
          }
          # for some reason `VectorDistribution` won't accept double.xmax?
          max_double = sqrt(.Machine$double.xmax)
          params = predict_mat %>%
            as.data.frame %>%
            # deal with bad numbers
            transform(mean = replace(mean, is.na(mean), 0),
                      sd = replace(sd, is.na(sd), .Machine$double.xmax)) %>%
            # original number is log sd
            transform(sd = exp(sd)) %>%
            # deal with too large numbers
            transform(sd = replace(sd, is.infinite(sd), max_double)) %>%
            transform(sd = replace(sd, sd > max_double, max_double))
          tryCatch({
            distrs = distr6::VectorDistribution$new(distribution = "Normal",
                                                    params = params)
          }, error = function(e) {
            print(summary(params))
            print(apply(params, 2, range, na.rm = TRUE))
            print(apply(params, 2, function(x) table(is.na(x))))
            stop(e)
          })
          list(distr = distrs, response = params$mean)
        } else {
          stopf("Invalid predict_type for task_type 'regr'.")
        }
      }
  )
)

# can use this to flexibly transform target variable. Also see
# https://github.com/mlr-org/mlr3pipelines/issues/956
PipeOpTargetCategScale = R6::R6Class("PipeOpTargetCategScale",
  inherit = PipeOpTargetTrafo,
  public = list(
    packages = c('magrittr'),
    initialize = function(id = "targetscale") super$initialize(id = id)
  ),
  private = list(
    .transform = function(task, phase) {
      require('magrittr', quietly = TRUE)
      old_target = task$data(cols = task$target_names)
      network = task$data()$network
      if (phase == 'train') {
        # get the scale factors
        self$state$scale_params = old_target %>%
          by(network, function(x) c(mean = mean(x$load), sd = sd(x$load))) %>%
          sapply(identity) %>%
          t
      }
      internal_valid = NULL
      if (!is.null(task$internal_valid_task)) {
        # preserve internal validation task
        internal_valid = task$internal_valid_task
      }
      # including the function definition here because the future package
      # doesn't automatically export it
      categ_scale_target = function(task, scale_params) {
        old_target = task$data(cols = task$target_names)
        network = task$data()$network
        # scale each network's loads
        new_target = (old_target - scale_params[network, 'mean']) /
          scale_params[network, 'sd']
        # I'm honestly confused about how the duplicate columns are handled but it
        # seems to work
        task$cbind(new_target)
        convert_task(task, target = colnames(new_target), drop_original_target = TRUE,
                     drop_levels = FALSE)
      }
      out = categ_scale_target(task, self$state$scale_params)
      # out = convert_task(task, target = colnames(new_target), new_type = private$.new_task_type, drop_original_target = TRUE, drop_levels = FALSE)
      if (!is.null(internal_valid)) {
        out$internal_valid_task =
          categ_scale_target(internal_valid, self$state$scale_params)
      }
      out
    },
    .train_invert = function(task) {
      # return a predict_phase_state object (can be anything)
      list(truth = task$truth(), network = task$data()$network)
    },
    .invert = function(prediction, predict_phase_state) {
      orig_means = unlist(prediction$distr$getParameterValue('mean'))
      orig_sds = unlist(prediction$distr$getParameterValue('sd'))
      offsets = self$state$scale_params[predict_phase_state$network, 'mean']
      scales = self$state$scale_params[predict_phase_state$network, 'sd']
      new_means = (orig_means * scales) + offsets
      new_sds = orig_sds * scales
      new_truth = predict_phase_state$truth
      distrs = data.frame(mean = new_means, sd = new_sds) %>%
        distr6::VectorDistribution$new(distribution = "Normal",
                                       params = .)
      PredictionRegr$new(row_ids = prediction$row_ids,
                         truth = new_truth, distr = distrs)
    }
  )
)
# mlr_pipeops$add("targetmutate", PipeOpTargetMutate)

nn_normal_crps_loss <- nn_module(
    "nn_normal_crps_loss",
    inherit = torch:::nn_loss,
    initialize = function() {},
    forward = function(output, target) {
      # following
      # https://github.com/slerch/ppnn/blob/7af9af3cfb754b8b6da54d2fe5d917cd54e32b9d/nn_postprocessing/nn_src/losses.py#L14C1-L47C24
      mu = output[, 1, drop = FALSE]
      log_sigma = output[, 2, drop = FALSE]
      sigma = log_sigma$exp()
      loc = (target - mu) / sigma
      phi = (1 / sqrt(2 * pi)) * torch_exp(-loc$square() / 2)
      Phi = .5 * (1 + torch_erf(loc / sqrt(2)))
      crps =  sigma * (loc * (2 * Phi - 1) + 2 * phi - 1 / sqrt(pi))
      crps$mean()
    }
)

# ########################
# NEXT STEP: Multiple networks
# ########################

drn_emb = nn_module("drn_emb",
  initialize = function(task, size_hidden, size_emb) {
    require(torch, quietly = TRUE)
    n_networks = length(levels(task$data()$network))
    self$network_emb = nn_embedding(n_networks, size_emb)
    n_features = task$n_features - 1 + size_emb
    self$first = nn_linear(n_features, size_hidden)
    self$second = nn_linear(size_hidden, 2)
  },
  # argument x corresponds to the ingress token x
  forward = function(x_num, x_cat) {
    self$network_emb(x_cat[, 1]) %>%
      list(x_num) %>%
      torch_cat(dim = 2) %>%
      self$first() %>%
      nnf_relu %>%
      self$second()
  }
)

CallbackSetBatchHistory = R6::R6Class("CallbackSetBatchHistory",
  inherit = mlr3torch::CallbackSetHistory,
  lock_objects = FALSE,
  public = list(
    packages = 'mlr3torch',
    #' @description
    #' Initializes lists where the train and validation metrics are stored.
    on_begin = function() {
      self$train = list(list(epoch = numeric(0)))
      self$valid = list(list(epoch = numeric(0)))
      self$batch_idx = 0
    },
    #' @description
    #' Converts the lists to data.tables.
    state_dict = function() {
      train = data.table::rbindlist(self$train, fill = TRUE)
      colnames(train)[-(1:2)] = paste0("train.", colnames(train)[-(1:2)])
      valid = data.table::rbindlist(self$valid, fill = TRUE)
      colnames(valid)[-(1:2)] = paste0("valid.", colnames(valid)[-(1:2)])
      state = if (nrow(valid) == 0 && nrow(train) == 0) {
        data.table::data.table(epoch = numeric(0))
      } else if (nrow(valid) == 0) {
        train
      } else if (nrow(train) == 0) {
        valid
      } else {
        merge(train, valid, by = c("epoch", 'step'))
      }
      if (is.null(self$prev_state)) {
        state
      } else {
        rbind(state, self$prev_state)
      }
    },
    #' @description
    #' Add the latest training scores to the history.
    on_batch_end = function() {
      self$batch_idx = self$batch_idx + 1
      if (self$batch_idx %% 5 != 0) return()
      if (length(self$ctx$last_loss)) {
        self$train[[length(self$train) + 1]] = c(
          list(epoch = self$ctx$epoch, step = self$ctx$step),
          self$ctx$last_loss
        )
        # we need to calculate valid loss
        ctx = self$ctx
        task_valid = self$ctx$task_valid
        ctx$network$eval()
        pred_tensor = mlr3torch:::torch_network_predict_valid(ctx)
        truth_tensor = task_valid$data(cols = task_valid$target_names) %>%
          as.matrix %>%
          torch_tensor
        last_scores_valid = ctx$loss_fn(pred_tensor, truth_tensor)$item()
        ctx$network$train()
        self$valid[[length(self$valid) + 1]] = c(
          list(epoch = self$ctx$epoch, step = self$ctx$step),
          last_scores_valid
        )
      }
    },
    #' @description
    #' Add the latest training scores to the history.
    on_before_valid = function() {},
    #' @description
    #' Add the latest validation scores to the history.
    on_epoch_end = function() {}
  )
)
mlr3torch_callbacks$add("batchhistory", function() {
  TorchCallback$new(
    callback_generator = CallbackSetBatchHistory,
    param_set = ps(),
    id = "batchhistory",
    label = "Batch History",
    man = "mlr3torch::mlr_callback_set.history"
  )
})

get_network_msrs = function(network, task) {
  print(network)
  categ_scale_target # help future know to export this
  network_rows = which(task$data(cols = 'network')$network == network)
  network_row_ids = task$row_ids[network_rows]
  n_task = task$data(network_row_ids) %>%
    as_task_regr(target = 'load')
  print(n_task$nrow)
  n_pred = gl$predict(n_task)
  0:7 %>%
    sapply(function(i) {
      # print(i)
      which_i = which(n_task$data(cols = 'days_ahead')$days_ahead == i)
      rows_i = n_pred$row_ids[which_i]
      # print(length(rows_i))
      pred_i = n_pred$clone()$filter(rows_i)
      # fix the filter bug! https://github.com/mlr-org/mlr3/issues/1400
      pred_i$data$distr = n_pred$distr[which_i]
      pred_i$score(msrs(c('regr.crps', 'regr.mae', 'regr.mape')))
    }) %>%
    t %>%
    as.data.frame %>%
    transform(network = network, days_ahead = 0:7) %>%
    subset(select = c(network, days_ahead, regr.crps, regr.mae, regr.mape))
}

drn_emb_params = ps(
    size_hidden = p_int(lower = 1L, tags = "train"),
    size_emb = p_int(lower = 1L, tags = "train")
)
l2 = LearnerTorchModuleDistr$new(
    module_generator = drn_emb,
    param_set = drn_emb_params,
    task_type = 'regr',
    ingress_tokens = list(x_num = ingress_num(), x_cat = ingress_categ()),
    optimizer = NULL,
    loss = nn_normal_crps_loss,
    # callbacks = t_clbks(c("history", "batchhistory")),
    callbacks = t_clbks(c('batchhistory')),
    predict_types = 'distr'
)
l2$param_set$set_values(
  epochs = 10, batch_size = 1024, device = "cpu",
  size_hidden = to_tune(8:64),
  size_emb = to_tune(1:32),
  opt.lr = to_tune(.001, .05, TRUE),
  opt.weight_decay = to_tune(.005, .05, TRUE)
)
tt = l2 %>%
  PipeOpLearner$new() %>%
  ppl("targettrafo", trafo_pipeop = PipeOpTargetCategScale$new(), graph = .)
gl = po('fixfactors') %>>% po('scale') %>>% tt %>%
  as_learner
set_validate(gl, validate = 'test')

# time series resampling, testing on summer months
ResamplingLoadForecastCV = R6::R6Class("ResamplingLoadForecastCV",
  inherit = ResamplingForecastCV,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      super$initialize()
      # maybe should set `folds` to 2, and change id
    }
  ),
  active = list(
    #' @template field_iters
    iters = function(rhs) {
      # as.integer(self$param_set$values$folds)
      4
    }
  ),
  private = list(
    .sample = function(ids, task) {
      row_dates = as.Date(task$data(rows = ids)$fct_run)
      # Rolling window cross validation, in 1-month segments (where we have
      # adequate data). Remember to edit `iters` after changing this.
      train_ends = c(seq(as.Date('2024-07-01'), as.Date('2024-09-01'), by = 'month'),
                     as.Date('2025-06-01'))
      test_ends = c(seq(as.Date('2024-08-01'), as.Date('2024-10-01'), by = 'month'),
                    as.Date('2025-07-01'))
      train_ids = lapply(as.list(train_ends), function(x) ids[row_dates < x])
      test_ids = lapply(seq_along(train_ends), function(i) {
        ids[row_dates >= train_ends[i] & row_dates < test_ends[i]]
      })
      # drop periods without sufficient data with a warning
      # ...
      stopifnot(all(lengths(train_ids) > 0))
      stopifnot(all(lengths(test_ids) > 0))
      list(train = train_ids, test = test_ids)
    }
  )
)
# resamplings[["load_forecast_cv"]] = ResamplingForecastCV

task_tune = prepare_load_task_all_networks()

# multicore will raise a cuda error when using GPU. Or it might freeze for no
# apparent reason
# future::plan('multicore', workers = 4)
n_workers = 3
future::plan('multisession', workers = n_workers)
# everything has to be limited to 1 thread, otherwise it breaks multicore
data.table::setDTthreads(1)
torch_set_num_threads(1)
torch_set_num_interop_threads(1)

# Get the best hyperparameters. For now, I'm skipping the formal benchmarking
system.time({
  tune_res <- tune(
      tuner = tnr("random_search", batch_size = n_workers),
      task = task_tune,
      learner = gl,
      resampling = ResamplingLoadForecastCV$new(),
      measure = msr('regr.crps'),
      # measure = msr("internal_valid_score",
      #               select = "regr.mlplss.regr.logloss", minimize = TRUE),
      term_evals = 12,
      store_models = T
  )
})


# # m = model, s = sample I think?
# plot_tune_hist = function(tune_res, m, s) {
#   t1 = tune_res$archive$learners(m)[[s]]$base_learner()
#   with(t1$model$callbacks$batchhistory, {
#     # multiply by 25, because history is collected every 25 batches
#     x = 25 * seq_len(length(train.V1))
#     # ylim = c(0, max(c(train.V1, valid.regr.crps)))
#     ylim = range(c(train.V1, valid.V1))
#     plot(x, train.V1, type = 'l', ylim = ylim)
#     lines(x, valid.V1, col = 'blue')
#     grid(nx = NA, ny = NULL)
#   })
# }

# par(mfrow = c(3, 4))
# for (m in 1:3) for (s in 1:4) plot_tune_hist(tune_res, m, s)

#    regr.module.size_hidden regr.module.size_emb regr.module.opt.lr regr.module.opt.weight_decay regr.crps
#                     <char>               <char>              <num>                        <num>     <num>
# 1:                      20                    8          -3.617886                     -5.04026  3.441419

# will need to get accuracy by lead time
# res = 0:7 %>%
#   sapply(function(i) {
#     days_scale = attr(scale_res, 'scaled:scale')['days_ahead']
#     days_offset = attr(scale_res, 'scaled:center')['days_ahead']
#     days_i = (i - days_offset) / days_scale
#     task_i = task_sc$clone()$filter(which(task_sc$data()$days_ahead == days_i))
#     pred_i = l1$predict(task_i)
#     pred_i$score(msrs(c('regr.crps', 'regr.mae'))) * load_scale
#   }) %>%
#   t %>%
#   as.data.frame %>%
#   transform(days_ahead = 0:7)

# train model with the best hyperparameters

l2 = LearnerTorchModuleDistr$new(
    module_generator = drn_emb,
    param_set = drn_emb_params,
    task_type = 'regr',
    ingress_tokens = list(x_num = ingress_num(), x_cat = ingress_categ()),
    optimizer = NULL,
    loss = nn_normal_crps_loss,
    # callbacks = t_clbks(c("history", "batchhistory")),
    callbacks = t_clbks(c('batchhistory')),
    predict_types = 'distr'
)
l2$param_set$set_values(
  epochs = 10, batch_size = 1024, device = "cpu",
  size_hidden = tune_res$result_learner_param_vals$regr.module.size_hidden,
  size_emb = tune_res$result_learner_param_vals$regr.module.size_emb,
  opt.lr = tune_res$result_learner_param_vals$regr.module.opt.lr,
  opt.weight_decay = tune_res$result_learner_param_vals$regr.module.opt.weight_decay
)
tt = l2 %>%
  PipeOpLearner$new() %>%
  ppl("targettrafo", trafo_pipeop = PipeOpTargetCategScale$new(), graph = .)
gl = po('fixfactors') %>>% po('scale') %>>% tt %>%
  as_learner

# gl$param_set$values = tune_res$result_learner_param_vals
# set_validate(gl, validate = NULL)

set_validate(gl, validate = .2)
gl$train(task_tune)

# save model -- see
# https://mlr3torch.mlr-org.com/reference/mlr_learners_torch.html#saving-a-learner

gl$marshal()

saveRDS(gl, 'results/forecast_load/network_model.rds')


# OK now training the system model

drn_vanilla = nn_module("drn_vanilla",
  initialize = function(task, size_hidden) {
    require(torch, quietly = TRUE)
    n_features = task$n_features
    self$first = nn_linear(n_features, size_hidden)
    self$second = nn_linear(size_hidden, 2)
  },
  # argument x corresponds to the ingress token x
  forward = function(x_num) {
    x_num %>%
      self$first() %>%
      nnf_relu %>%
      self$second()
  }
)

drn_vanilla_params = ps(
    size_hidden = p_int(lower = 1L, tags = "train")
)
learner_sys = LearnerTorchModuleDistr$new(
    module_generator = drn_vanilla,
    param_set = drn_vanilla_params,
    task_type = 'regr',
    ingress_tokens = list(x_num = ingress_num()),
    optimizer = NULL,
    loss = nn_normal_crps_loss,
    # callbacks = t_clbks(c("history", "batchhistory")),
    callbacks = t_clbks(c('batchhistory')),
    predict_types = 'distr'
)
# data set is around 70 times smaller than the network dataset, therefore we
# want something like 70x as many epochs
learner_sys$param_set$set_values(
  epochs = 200, batch_size = 1024, device = "cpu",
  size_hidden = to_tune(8:64),
  opt.lr = to_tune(.001, .05, TRUE),
  opt.weight_decay = to_tune(.01, .5, TRUE)
)
tt_sys = learner_sys %>%
  PipeOpLearner$new() %>%
  ppl("targettrafo", trafo_pipeop = PipeOpTargetTrafoScaleRangeDistr$new(),
      graph = .)
gl_sys = po('scale') %>>% tt_sys %>%
  as_learner
set_validate(gl_sys, validate = 'test')

task_tune_sys = prepare_load_task_combined('system')
# task_tune_sys2 = task_tune_sys$clone()$select(c('TV', 'fct_run'))
# plot(task_tune_sys$data()$TV, task_tune_sys$data()$load)

# Get the best hyperparameters. For now, I'm skipping the formal benchmarking
system.time({
  tune_res_sys <- tune(
      tuner = tnr("random_search", batch_size = n_workers),
      task = task_tune_sys,
      learner = gl_sys,
      resampling = ResamplingLoadForecastCV$new(),
      measure = msrs('regr.crps'),
      # measure = msr("internal_valid_score",
      #               select = "regr.mlplss.regr.logloss", minimize = TRUE),
      term_evals = 12,
      store_models = T
  )
})

# train system model with the best hyperparameters

learner_sys = LearnerTorchModuleDistr$new(
    module_generator = drn_vanilla,
    param_set = drn_vanilla_params,
    task_type = 'regr',
    ingress_tokens = list(x_num = ingress_num()),
    optimizer = NULL,
    loss = nn_normal_crps_loss,
    # callbacks = t_clbks(c("history", "batchhistory")),
    callbacks = t_clbks(c('batchhistory')),
    predict_types = 'distr'
)
learner_sys$param_set$set_values(
  epochs = 10, batch_size = 1024, device = "cpu",
  size_hidden = tune_res_sys$result_learner_param_vals[[1]]$regr.module.size_hidden,
  opt.lr = tune_res_sys$result_learner_param_vals[[1]]$regr.module.opt.lr,
  opt.weight_decay = tune_res_sys$result_learner_param_vals[[1]]$regr.module.opt.weight_decay
)
tt_sys = learner_sys %>%
  PipeOpLearner$new() %>%
  ppl("targettrafo", trafo_pipeop = PipeOpTargetTrafoScaleRangeDistr$new(),
      graph = .)
gl_sys = po('scale') %>>% tt_sys %>%
  as_learner

# gl$param_set$values = tune_res$result_learner_param_vals
# set_validate(gl, validate = NULL)

set_validate(gl_sys, validate = .2)
gl_sys$train(task_tune_sys)

# save model -- see
# https://mlr3torch.mlr-org.com/reference/mlr_learners_torch.html#saving-a-learner
gl_sys$marshal()

saveRDS(gl_sys, 'results/forecast_load/system_model.rds')
