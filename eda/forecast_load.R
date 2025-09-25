# Evaluate daily peak load forecasting methods

# Issues to be addressed:

# - Feature engineering:
#   - how to summarize/aggregate hourly forecasts for daily predictions
#   - including lagged values
#   - summarizing/aggregating gridded values
# - Statistical issues:
#   - regularization due to small data size
#   - forecasting for multiple nearby areas (spatially autocorrelated errors)
#   - serially autocorrelated errors for longer forecast times

# A helpful trick: can subtract the lagged portion of TV from the output so it
# doesn't need to be predicted (and add it back afterward)

# NEXT STEPS:
# - get confidence intervals/p-values (randomize data etc.) May want to shuffle
#   minibatches. See https://arxiv.org/pdf/1206.5533, p. 6

# pkg_deps = c('ncmeta', 'crch', 'ranger', 'drf', 'RandomForestsGLS',
#              'scoringRules')
# install.packages(pkg_deps)
# remotes::install_github("mlr-org/mlr3temporal")
# distr6_repos =  c(CRAN = 'https://cloud.r-project.org',
#                   raphaels1 = 'https://raphaels1.r-universe.dev')
# install.packages('distr6', repos = distr6_repos)
# install.packages("mlr3proba", repos = "https://mlr-org.r-universe.dev")
setwd('..')
# library(timeDate) # holidays
library(magrittr)
library(stars) # also requires ncmeta
library(future)
# library(future.batchtools)
library(mlr3)
library(torch)
library(mlr3torch)
library(mlr3learners)
library(mlr3tuning)
# library(mlr3tuningspaces)
library(mlr3pipelines)
library(mlr3temporal)
library(mlr3proba)
source('R/mlr3_additions.R') # block CV, mean pinball loss, etc.
source('R/mlr3_distr.R')

# see notes in `?benchmark` and `?mlr_tuners_random_search`
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

all_eff_tmps = read.csv('results/process_station_data/eff_tmp.csv')
# all_tvs = read.csv('results/process_station_data/tv.csv')
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
# gam_fct = readRDS('results/load_curve_forecast/gam_forecasts.rds')

# get_valid_day = function(nc, days_ahead) {
#   # time is generally read by R as a date, due to the 1 day differences
#   attr(nc, 'dimensions')$time$values + days_ahead
# }

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

plot_loss = function(learner, loss, skip = 0, scale = NULL, train = TRUE) {
  with(learner$model$callbacks, {
    h = as.data.frame(history)
    if (skip) {
      h = tail(h, -skip)
    }
    tloss = h[, paste0('train.', loss)]
    vloss = h[, paste0('valid.', loss)]
    if (!is.null(scale)) {
      tloss = tloss * scale
      vloss = vloss * scale
    }
    if (train) {
      ylim = range(c(tloss, vloss))
      plot(h$epoch, tloss, type = 'l', ylim = ylim, xlab = 'Epoch', ylab = loss)
      points(h$epoch, vloss, type = 'l', col = 'blue')
    } else {
      plot(h$epoch, vloss, type = 'l', xlab = 'Epoch', ylab = loss)
    }
  })
  grid(nx = NA, ny = NULL)
}

nn_normal_nll_loss <- nn_module(
    "nn_normal_nll_loss",
    inherit = torch:::nn_loss,
    initialize = function() {},
  forward = function(output, target) {
    # output: [batch, 2] -> mean, log_var
    # NOTE: `target` is a 2D tensor, and operations with a 1D array will lead to
    # an outer product-- that's not what we want!! So we must preserve the 2nd
    # dimensions here
    mean <- output[, 1, drop = FALSE]
    log_sd <- output[, 2, drop = FALSE]
    sd <- torch_exp(log_sd)
    -distr_normal(mean, sd)$log_prob(target)$mean()
  }
)

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


LearnerTorchMLPDistr <- R6::R6Class(
  "LearnerTorchMLPDistr",
  inherit = LearnerTorch,
  public = list(
    # initialize = function(task_type, optimizer = NULL, loss = NULL, callbacks = list()) {
    #   out = super$initialize(task_type, optimizer, loss, callbacks)
    #   # ps <- self$param_set$clone(deep = TRUE)
    #   # ps$values$out_features <- 2  # mean + logvar
    #   # self$param_set$values$out_features <- 2  # mean + logvar
    #   self$predict_types = c('response', "distr")
    #   self$predict_type = "distr"
    #   out
    # },
    initialize = function(task_type, optimizer = NULL, loss = NULL, callbacks = list()) {
      # check_activation = mlr3misc::crate(function(x) check_class(x, "nn_module"))
      # check_activation_args = mlr3misc::crate(function(x) check_list(x, names = "unique"))
      # check_neurons = mlr3misc::crate(function(x) checkmate::check_integerish(x, any.missing = FALSE, lower = 1))
      # check_shape = mlr3misc::crate(function(x) check_shape(x, null_ok = TRUE, len = 2L))
      param_set = ps(
        neurons         = p_uty(tags = c("train", "predict")),
        p               = p_dbl(lower = 0, upper = 1, tags = "train"),
        n_layers        = p_int(lower = 1L, tags = "train"),
        activation      = p_uty(tags = c("required", "train")),
        activation_args = p_uty(tags = c("required", "train")),
        shape           = p_uty(tags = "train")
      )
      param_set$set_values(
        activation = nn_relu,
        activation_args = list(),
        neurons = integer(0),
        p = 0.5
      )
      super$initialize(
        task_type = task_type,
        id = paste0(task_type, ".mlplss"),
        label = "Multi Layer Perceptron",
        param_set = param_set,
        optimizer = optimizer,
        callbacks = callbacks,
        loss = loss,
        predict_types = c('distr', 'response'), # NEW!
        man = "mlr3torch::mlr_learners.mlp",
        feature_types = c("numeric", "integer", "lazy_tensor"),
        jittable = TRUE
      )
    }
  ),
  private = list(
      .encode_prediction = function(predict_tensor, task) {
        # see `encode_prediction_default` for example
        if (self$predict_type == "response") {
          list(response = as.numeric(predict_tensor[, 1]))
        } else if (self$predict_type == "distr") {
          params = predict_tensor %>%
            as.numeric %>%
            matrix(ncol = 2, dimnames = list(NULL, c('mean', 'sd'))) %>%
            as.data.frame %>%
            # original number is log sd
            transform(sd = exp(sd)) %>%
            # deal with bad numbers
            # transform(mean = replace(mean, is.na(mean), 0)) %>%
            # transform(sd = replace(sd, sd < 0, 0)) %>%
            transform(sd = replace(sd, is.infinite(sd), sqrt(.Machine$double.xmax)))
          distrs = distr6::VectorDistribution$new(distribution = "Normal",
                                                  params = params)
          list(distr = distrs, response = params$mean)
        } else {
          stopf("Invalid predict_type for task_type 'regr'.")
        }
      },
    .network = function(task, param_vals) {
      # d_out = output_dim_for(task)
      d_out = 2
      d_in = private$.ingress_tokens(task, param_vals)[[1L]]$shape[2L]
      network = mlr3misc::invoke(mlr3torch:::make_mlp, .args = param_vals, d_in = d_in, d_out = d_out)
      network
    }
  )
)

LearnerTorchModuleDistr = R6::R6Class("LearnerTorchModuleDistr",
  inherit = LearnerTorchModule,
  public = list(
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

drn_one_layer = nn_module("drn_one_layer",
  initialize = function(task, size_hidden) {
    self$first = nn_linear(task$n_features, size_hidden)
    self$second = nn_linear(size_hidden, 2)
  },
  # argument x corresponds to the ingress token x
  forward = function(x) {
    x = self$first(x)
    x = nnf_relu(x)
    self$second(x)
  }
)

# try out the custom DRN
# drn_params = ps(
#     size_hidden = p_int(lower = 1L, tags = "train"),
#     opt.lr = to_tune(1e-3, .1, TRUE),
#     opt.weight_decay = to_tune(1e-4, 1, TRUE)
# )
scale_res = prepare_load_task('system', 0)$data() %>%
                                         scale
load_scale = attr(scale_res, 'scaled:scale')['load']
task_sc = prepare_load_task('system', 0)$data() %>%
                                       scale %>%
                                       as_task_regr(target = 'load')
drn_params = ps(
    size_hidden = p_int(lower = 1L, tags = "train")
)
l1 = LearnerTorchModuleDistr$new(
    module_generator = drn_one_layer,
    param_set = drn_params,
    task_type = 'regr',
    ingress_tokens = list(input = ingress_num(shape = 1)),
    optimizer = NULL,
    loss = nn_normal_crps_loss,
    callbacks = t_clbk("history"),
    predict_types = 'distr'
)
l1$param_set$set_values(
  epochs = 200, batch_size = task_sc$nrow, device = "cpu",
  size_hidden = 15,
  opt.lr = 0.005, opt.weight_decay = .05,
  # validate = .2,
  # Measures to track
  measures_valid = msrs(c('regr.crps', 'regr.mae')),
  measures_train = msrs(c('regr.crps', 'regr.mae'))
)
set_validate(l1, .2)

l1$train(task_sc)

plot_loss(l1, 'regr.crps', 150, load_scale)
plot_loss(l1, 'regr.mae', 150, load_scale)


# ok good let's try tuning. can't use `tune` here because it doesn't work well
# with mlr3torch validation tracking

l2 = LearnerTorchModuleDistr$new(
    module_generator = drn_one_layer,
    param_set = drn_params,
    task_type = 'regr',
    ingress_tokens = list(input = ingress_num(shape = 1)),
    optimizer = NULL,
    loss = nn_normal_crps_loss,
    callbacks = t_clbk("history"),
    predict_types = 'distr'
)
l2$param_set$set_values(
  epochs = 10, batch_size = 128, device = "cpu",
  size_hidden = to_tune(5, 25),
  # opt.lr = .005,
  opt.lr = to_tune(.005, .1, TRUE),
  opt.weight_decay = to_tune(1e-2, .1, TRUE)
)

# let's do it manually

# the sampling is done by paradox::SamplerUnif
s1 = SamplerUnif$new(l2$param_set$search_space())
# sam_bak = sam
sam = s1$sample(10)

# future::plan('multicore', workers = 2)
# not sure why, but multicore isn't working
future::plan('multisession', workers = 4)
library(future.apply)

results = future_apply(sam$data, 1, function(x) {
  # make sure `future` knows about regr.crps
  mlr_measures$add("regr.crps", function() MeasureRegrCRPS$new())
  params_i = x %>%
    as.list %>%
    transform(opt.lr = exp(opt.lr),
              opt.weight_decay = exp(opt.weight_decay))
  l_i = LearnerTorchModuleDistr$new(
                                    module_generator = drn_one_layer,
                                    param_set = drn_params,
                                    task_type = 'regr',
                                    ingress_tokens = list(input = ingress_num(shape = 1)),
                                    optimizer = NULL,
                                    loss = nn_normal_crps_loss,
                                    callbacks = t_clbk("history"),
                                    predict_types = 'distr'
                                )
  l_i$param_set$set_values(
                    epochs = 1000, batch_size = task_sc$nrow, device = "cpu",
                    measures_valid = msrs(c('regr.crps', 'regr.mae')),
                    measures_train = msrs(c('regr.crps', 'regr.mae'))
                )
  do.call(l_i$param_set$set_values, params_i)
  set_validate(l_i, .2)
  l_i$train(task_sc)
  l_i
}, simplify = FALSE, future.seed = TRUE)

lapply(results, function(x) {
  tail(x$model$callbacks$history, 1)
}) %>%
  do.call(rbind, .)

par(mfrow = c(3, 4))
for (r in results) plot_loss(r, 'regr.crps', 300, scale = load_scale)

par(mfrow = c(3, 4))
for (r in results) plot_loss(r, 'regr.crps', 300, scale = load_scale, train = F)

par(mfrow = c(3, 4))
for (r in results) plot_loss(r, 'regr.mae', 100, scale = load_scale)

par(mfrow = c(3, 4))
for (r in results) plot_loss(r, 'regr.mae', 700, scale = load_scale, train = F)

tune_res = results %>%
  sapply(function(x) {
    x$model$callbacks$history[, c('valid.regr.crps', 'valid.regr.mae')] %>%
      tail(1) %>%
      unlist
  }) %>%
  t %>%
  as.data.frame %>%
  cbind(sam$data) %>%
  transform(opt.lr = exp(opt.lr),
            opt.weight_decay = exp(opt.weight_decay),
            valid.regr.crps = valid.regr.crps * load_scale,
            valid.regr.mae = valid.regr.mae * load_scale)
tune_res = tune_res[, c(3:5, 1:2)]
tune_res = tune_res[order(tune_res$valid.regr.crps), ]
# looks like we want weight decay .02, lr around .01-.02, hidden size matters
# little but maybe 16?

# how repeatable is this?
results2 = future_apply(sam$data, 1, function(x) {
  # make sure `future` knows about regr.crps
  mlr_measures$add("regr.crps", function() MeasureRegrCRPS$new())
  # params_i = list(size_hidden = 16, opt.lr = .01, opt.weight_decay = .02)
  l_i = LearnerTorchModuleDistr$new(
                                    module_generator = drn_one_layer,
                                    param_set = drn_params,
                                    task_type = 'regr',
                                    ingress_tokens = list(input = ingress_num(shape = 1)),
                                    optimizer = NULL,
                                    loss = nn_normal_crps_loss,
                                    callbacks = t_clbk("history"),
                                    predict_types = 'distr'
                                )
  l_i$param_set$set_values(
                    epochs = 1000, batch_size = task_sc$nrow, device = "cpu",
                    size_hidden = 16, opt.lr = .01, opt.weight_decay = .025,
                    measures_valid = msrs(c('regr.crps', 'regr.mae')),
                    measures_train = msrs(c('regr.crps', 'regr.mae'))
                )
  # do.call(l_i$param_set$set_values, params_i)
  set_validate(l_i, .2)
  l_i$train(task_sc)
  l_i
}, simplify = FALSE, future.seed = TRUE)

tune_res2 = results2 %>%
  sapply(function(x) {
    x$model$callbacks$history[, c('valid.regr.crps', 'valid.regr.mae')] %>%
      tail(1) %>%
      unlist
  }) %>%
  t %>%
  as.data.frame %>%
  cbind(sam$data) %>%
  transform(opt.lr = exp(opt.lr),
            opt.weight_decay = exp(opt.weight_decay),
            valid.regr.crps = valid.regr.crps * load_scale,
            valid.regr.mae = valid.regr.mae * load_scale)
tune_res2 = tune_res2[, c(3:5, 1:2)]
tune_res2 = tune_res2[order(tune_res2$valid.regr.crps), ]

par(mfrow = c(3, 4))
for (r in results2) plot_loss(r, 'regr.crps', 100, scale = load_scale)

par(mfrow = c(3, 4))
for (r in results2) plot_loss(r, 'regr.crps', 100, scale = load_scale, train = F)

par(mfrow = c(3, 4))
for (r in results2) plot_loss(r, 'regr.mae', 100, scale = load_scale)

par(mfrow = c(3, 4))
for (r in results2) plot_loss(r, 'regr.mae', 700, scale = load_scale, train = F)



# res_hist = results %>%
#   lapply(function(x) tail(x$model$callbacks$history, 1)) %>%
#   do.call(rbind, .) %>%
#   cbind(sam$data, iter = 1:nrow(sam$data)) %>%
#   transform(opt.lr = exp(opt.lr),
#             opt.weight_decay = exp(opt.weight_decay)) %>%
#   as.data.frame %>%
#   subset(select = -epoch)
# res_hist = res_hist[order(res_hist$valid.regr.logloss), ]
# res_cols = c('iter', 'opt.lr', 'opt.weight_decay', 'size_hidden',
#              'valid.regr.logloss')
# res_hist[, res_cols]

cor(res_hist[, res_cols[-1]])
res_hist[, res_cols[-1]] %>%
  subset(valid.regr.logloss <= 100) %>%
  pairs
# looks like we want weight decay .05, lr around .02-.03, hidden size matters
# little but maybe 16?

# ########################
# NEXT STEP: Multiple lead times
# ########################

scale_res = prepare_load_task_combined('system')$data() %>%
                                         scale
load_scale = attr(scale_res, 'scaled:scale')['load']
task_sc = prepare_load_task_combined('system')$data() %>%
                                       scale %>%
                                       as_task_regr(target = 'load')
drn_params = ps(
    size_hidden = p_int(lower = 1L, tags = "train")
)
l1 = LearnerTorchModuleDistr$new(
    module_generator = drn_one_layer,
    param_set = drn_params,
    task_type = 'regr',
    ingress_tokens = list(input = ingress_num(shape = 1)),
    optimizer = NULL,
    loss = nn_normal_crps_loss,
    callbacks = t_clbk("history"),
    predict_types = 'distr'
)
l1$param_set$set_values(
  epochs = 200, batch_size = 512, device = "cpu",
  size_hidden = 16,
  opt.lr = 0.01, opt.weight_decay = .025,
  # validate = .2,
  # Measures to track
  measures_valid = msrs(c('regr.crps', 'regr.mae')),
  measures_train = msrs(c('regr.crps', 'regr.mae'))
)
set_validate(l1, .2)

system.time(l1$train(task_sc))

plot_loss(l1, 'regr.crps', 50, load_scale)
plot_loss(l1, 'regr.mae', 0, load_scale)



# get the loss by lead time
res = 0:7 %>%
  sapply(function(i) {
    days_scale = attr(scale_res, 'scaled:scale')['days_ahead']
    days_offset = attr(scale_res, 'scaled:center')['days_ahead']
    days_i = (i - days_offset) / days_scale
    task_i = task_sc$clone()$filter(which(task_sc$data()$days_ahead == days_i))
    pred_i = l1$predict(task_i)
    pred_i$score(msrs(c('regr.crps', 'regr.mae'))) * load_scale
  }) %>%
  t %>%
  as.data.frame %>%
  transform(days_ahead = 0:7)
# looks good!



# ########################
# NEXT STEP: Multiple networks
# ########################

drn_emb = nn_module("drn_emb",
  initialize = function(task, size_hidden, size_emb) {
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

# can use this to flexibly transform target variable
PipeOpTargetCategScale = R6::R6Class("PipeOpTargetCategScale",
  inherit = PipeOpTargetTrafo,
  public = list(
    initialize = function(id = "targetscale") super$initialize(id = id)
  ),
  private = list(
    .transform = function(task, phase) {
      # print(self$state)
      old_target = task$data(cols = task$target_names)
      network = task$data()$network
      if (phase == 'train') {
        # get the scale factors
        self$state$scale_params = old_target %>%
          by(network, function(x) c(mean = mean(x$load), sd = sd(x$load))) %>%
          sapply(identity) %>%
          t
      }
      # scale each network's loads
      new_target = (old_target - self$state$scale_params[network, 'mean']) /
        self$state$scale_params[network, 'sd']
      # new_target = self$param_set$values$trafo(task$data(cols = task$target_names))
      if (!is.data.frame(new_target) && !is.matrix(new_target)) {
        stopf("Hyperparameter 'trafo' must be a function returning a 'data.frame', 'data.table', or 'matrix', not '%s'.", class(new_target)[[1L]])
      }
      task$cbind(new_target)
      convert_task(task, target = colnames(new_target), new_type = private$.new_task_type, drop_original_target = TRUE, drop_levels = FALSE)
    },
    .train_invert = function(task) {
      # return a predict_phase_state object (can be anything)
      list(network = task$data()$network)
    },
    .invert = function(prediction, predict_phase_state) {
      orig_means = unlist(prediction$distr$getParameterValue('mean'))
      orig_sds = unlist(prediction$distr$getParameterValue('sd'))
      # print(head(self$state$scale_params))
      offsets = self$state$scale_params[predict_phase_state$network, 'mean']
      scales = self$state$scale_params[predict_phase_state$network, 'sd']
      new_means = (orig_means * scales) + offsets
      new_sds = orig_sds * scales
      new_truth = (prediction$truth * scales) + offsets
      # print(head(data.frame(mean = new_means, sd = new_sds)))
      distrs = data.frame(mean = new_means, sd = new_sds) %>%
        distr6::VectorDistribution$new(distribution = "Normal",
                                       params = .)
      PredictionRegr$new(row_ids = prediction$row_ids,
                         truth = new_truth, distr = distrs)
    }
  )
)
# mlr_pipeops$add("targetmutate", PipeOpTargetMutate)

# networks_task = prepare_load_task_all_networks()
# because that's a lot
networks_task = prepare_load_task_all_networks(sample(networks$id, 5))

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
    callbacks = t_clbk("history"),
    predict_types = 'distr'
)
l2$param_set$set_values(
  epochs = 10, batch_size = 1024, device = "cpu",
  size_hidden = 16,
  size_emb = 2,
  opt.lr = 0.005, opt.weight_decay = .025,
  # Measures to track
  measures_valid = msrs(c('regr.crps', 'regr.mae')),
  measures_train = msrs(c('regr.crps', 'regr.mae'))
)
# set_validate(l2, .2)
tt = l2 %>%
  PipeOpLearner$new() %>%
  ppl("targettrafo", trafo_pipeop = PipeOpTargetCategScale$new(), graph = .)
gl = po("scale") %>>% tt %>%
  as_learner
set_validate(gl$base_learner(), validate = .2)
# internal validation with pipeline notes here:
# https://mlr3book.mlr-org.com/chapters/chapter15/predsets_valid_inttune.html

system.time(gl$train(networks_task))
#    user  system elapsed 
# 117.173   0.016 116.058 # CPU work
#    user  system elapsed 
# 115.089   0.150 114.114 # CPU work, intel mkl
 #   user  system elapsed 
 # 64.964   0.863  62.511 # GPU home
 #   user  system elapsed 
 # 62.633   0.083  58.601 # CPU home (lol)
# I think the home laptop is using turbo and switching CPU cores, which won't be
# sustainable if I run many in parallel

#    user  system elapsed 
# 106.096   0.047 105.326 # 5 batches of 1024
#     user   system  elapsed 
# 5593.180    5.111 5612.521 # 50 batches of 1024
# runtime seems to scale non-linearly with dataset size. I don't get why

# rt50 ~= rt5 * 50, but why?

plot_loss(gl$base_learner(), 'regr.crps', 0)
plot_loss(gl$base_learner(), 'regr.mae', 0)

pred = gl$predict(networks_task2)

pred %>%
  as.data.table %>%
  head %>%
  subset(select = -distr)

pred$distr$getParameterValue('mean') %>%
  head %>%
  unlist

pred$distr$getParameterValue('sd') %>%
  head %>%
  unlist

get_network_msrs = function(network, task) {
  print(network)
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

res2 = networks_task$data()$network %>%
  unique %>%
  as.character %>%
  # head %>%
  lapply(get_network_msrs, task = networks_task) %>%
  do.call(rbind, .)

res %>%
  aggregate(regr.crps ~ days_ahead, FUN = mean)

res %>%
  aggregate(regr.mape ~ days_ahead, FUN = mean)

res %>%
  aggregate(regr.mae ~ days_ahead, FUN = mean)

library(ggplot2)

ggplot(res2, aes(days_ahead, regr.mae, group = network)) + geom_line()

ggplot(res2, aes(days_ahead, regr.mape * 100, group = network)) +
  geom_line(color = '#00000066') +
  stat_summary(aes(group = 1), geom = 'line', colour = 'blue', linewidth = 1.5) +
  ylim(c(0, NA))


# get the loss by lead time
res = networks_task_sc$data()$days_ahead %>%
  unique %>%
  sort %>%
  sapply(function(i) {
    task_i = networks_task_sc$clone()$filter(which(networks_task_sc$data()$days_ahead == i))
    pred_i = lt$predict(task_i)
    pred_i$score(msrs(c('regr.crps', 'regr.mae')))
  }) %>%
  t %>%
  as.data.frame %>%
  transform(days_ahead = 0:7)
# looks good!





# ok let's try the pipe and do a benchmark, to get results in MW
l3 = LearnerTorchModuleDistr$new(
    module_generator = drn_one_layer,
    param_set = drn_params,
    task_type = 'regr',
    ingress_tokens = list(input = ingress_num(shape = 1)),
    optimizer = NULL,
    loss = nn_normal_nll_loss,
    predict_types = 'distr'
)
l3$param_set$set_values(
    epochs = 1000, batch_size = 128, device = "cpu",
    size_hidden = 16,
    opt.lr = .025,
    opt.weight_decay = .05
)
pos = po("scale")
task_sc = pos$train(list(task1))$output
tlearner = ppl("targettrafo",
               trafo_pipeop = PipeOpTargetTrafoScaleRangeDistr$new(param_vals = list(lower = -1)),
               graph = PipeOpLearner$new(l3)) %>%
  as_learner

bgrid = benchmark_grid(
    tasks = task_sc,
    #learners = ngr,
    learners = c(tlearner),
    resamplings = rsmp('forecast_holdout', ratio = 5 / 6)
)
system.time(bres <- benchmark(bgrid))
bres$aggregate(msrs(c('regr.crps', 'regr.mae', 'regr.mape', 'regr.logloss')))
#    regr.crps regr.mae  regr.mape regr.logloss
#        <num>    <num>      <num>        <num>
# 1:  195.5224 254.9854 0.03550228     12.15634
# not bad!

# This is a small dataset-- can we just use a traditional optimizer like
# `optim_lbfgs`? Will need to change .train method-- see
# https://torch.mlverse.org/docs/reference/optim_lbfgs#details



# I wonder if we can do the thing where we incrementally train? -- No, shouldn't
# do that for benchmarking. Better to have independent results




## Part 1: fit a model, make sure results are reasonable and validation loss has
## stabilized

task1 = prepare_load_task('6B', 0, use_gam = FALSE)

useful_msrs = c('regr.mae', 'regr.crps', 'regr.mape', 'regr.logloss')

# Example usage
# task <- tsk("mtcars")
learner0 <- LearnerTorchMLPDistr$new(
  'regr',
  optimizer = t_opt("adam", lr = 0.005, weight_decay = .05),
  loss = nn_normal_nll_loss,
  callbacks = t_clbk("history")
)
# learner$encapsulate("none")
learner0$param_set$set_values(
  epochs = 1000, batch_size = 128, device = "cpu",
  neurons = 15, p = 0,
  #validate = 'predefined',
  # Measures to track
  measures_valid = msrs(c("regr.rmse", "regr.crps")),
  measures_train = msrs(c("regr.rmse", "regr.crps"))
)
# learner0$validate = .2
# set_validate(learner0, validate = 'predefined')
pos = po("scale")
task_sc = pos$train(list(task1))$output
task_sc$col_roles[['original_ids']] = 'fct_run'
# task_sc$internal_valid_task = sample(task_sc$nrow, floor(.2 * task_sc$nrow))
# learner$train(task_sc)
tlearner = ppl("targettrafo",
               trafo_pipeop = PipeOpTargetTrafoScaleRangeDistr$new(param_vals = list(lower = -1)),
               graph = PipeOpLearner$new(learner0)) %>%
  as_learner
# set_validate(tlearner, validate = .2, ids = "regr.mlplss")
# set_validate(tlearner, validate = 'predefined')
tlearner$graph$pipeops$regr.mlplss$learner$validate = .2
tlearner$train(task_sc)
# scaling seems to be super important here!

pred <- tlearner$predict(task_sc)
pred$score(msrs(c('regr.mae', 'regr.crps', 'regr.mape')))
# pred$distr$distributions[[1]]
learner = tlearner$pipeops$regr.mlplss$learner_model

head(learner$model$callbacks$history)
tail(learner$model$callbacks$history)
# tail(tlearner$pipeops$regr.mlp$learner_model$model$callbacks$history)

with(learner$model$callbacks$history, {
  par(mfrow = c(1, 2))
  plot(epoch, train.regr.crps, type = 'l', ylim = c(0, max(train.regr.crps)))
  points(epoch, valid.regr.crps, type = 'l', col = 'blue')
  plot(epoch, train.regr.rmse, type = 'l', ylim = c(0, max(train.regr.rmse)))
  points(epoch, valid.regr.rmse, type = 'l', col = 'blue')
})

with(learner$model$callbacks$history, {
  par(mfrow = c(1, 2))
  plot(epoch, train.regr.crps, type = 'l', ylim = c(0, max(train.regr.crps)))
  # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
  plot(epoch, log(train.regr.crps), type = 'l', ylim = range(log(train.regr.crps)))
  points(epoch, log(valid.regr.crps), type = 'l', col = 'blue')
})

learner$model$callbacks$history %>%
  tail(60) %>%
  with({
    plot(epoch, log(train.regr.crps), type = 'l')
    points(epoch, log(valid.regr.crps), type = 'l', col = 'blue')
  })

learner$model$callbacks$history %>%
  tail(500) %>%
  with({
    plot(epoch, train.regr.crps, type = 'l')
    # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
  })

# quick benchmarking
# t_measures = msrs(c("regr.crps", 'regr.mae', 'regr.mape'))

# ginstance = ti(
#     task = tsk(task_sc),
#     learner = tlearner,
#     resampling = rsmp('forecast_holdout', ratio = 5 / 6),
#     measures = msrs("regr.crps")
# )

learner0 <- LearnerTorchMLPDistr$new(
                                    'regr',
                                    optimizer = t_opt("adam"),
                                    loss = nn_normal_nll_loss,
                                    callbacks = t_clbk("history")
                                )
learner0$param_set$set_values(
                      epochs = 1000, batch_size = 128, device = "cpu",
                      neurons = to_tune(5:50),
                      p = 0,
                      opt.lr = to_tune(1e-3, .1, TRUE),
                      opt.weight_decay = to_tune(1e-4, 1, TRUE),
                      measures_train = msrs(c("regr.mae", "regr.crps", "regr.logloss"))
                   )
task_sc = prepare_load_task('6B', 0, use_gam = FALSE)$data() %>%
  scale %>%
  as_task_regr(target = 'load')
# set_validate(learner0, 'test')
# set_validate(learner0, .2)

future::plan('multicore', workers = 3)
future::plan('multisession', workers = 3)
future::plan('sequential')

tinstance = tune(
  tuner = tnr("random_search"),
  task = task_sc,
  learner = learner0,
  resampling = rsmp("holdout", ratio = 5 / 6),
  measure = msr("regr.logloss"),
  # measure = msr("internal_valid_score",
  #               select = "regr.mlplss.regr.logloss", minimize = TRUE),
  term_evals = 10
)

l1 = tinstance$archive$benchmark_result$learners$learner[[1]]

# this isn't doing what I want, so let's just use a for loop

SamplerUnif$new(learner0$param_set)


head(l1$callbacks$history)
tail(learner$model$callbacks$history)

# based on that
learner0 <- LearnerTorchMLPDistr$new(
                                    'regr',
                                    optimizer = t_opt("adam"),
                                    loss = nn_normal_nll_loss,
                                    callbacks = t_clbk("history")
                                )
learner0$param_set$set_values(
                      epochs = 1000, batch_size = 128, device = "cpu",
                      neurons = 2,
                      # p = 0,
                      opt.lr = .1,
                      opt.weight_decay = 0,
                      measures_train = msrs(c("regr.mae", "regr.crps"))# ,
                      # measures_valid = msrs(c("regr.mae", "regr.crps"))
                   )
pos = po("scale")
# what if we change the features?
# task1$feature_names = c("TV", "TV_sd")
task1$col_roles$feature = c("TV", "TV_sd")
task_sc = pos$train(list(task1))$output
tlearner = ppl("targettrafo",
               trafo_pipeop = PipeOpTargetTrafoScaleRangeDistr$new(param_vals = list(lower = -1)),
               graph = PipeOpLearner$new(learner0)) %>%
  as_learner

# tlearner$graph$pipeops$regr.mlplss$learner$validate = .2
tlearner$train(task_sc)
# scaling seems to be super important here!
pred <- tlearner$predict(task_sc)
pred$score(msrs(c('regr.crps', 'regr.mae', 'regr.mape')))

plot(pred$truth, pred$response)
abline(0, 1)
abline(v = 8300)
# it's just predicting the mean load!

learner = tlearner$pipeops$regr.mlplss$learner_model

head(learner$model$callbacks$history)
tail(learner$model$callbacks$history)

with(learner$model$callbacks$history, {
  par(mfrow = c(1, 2))
  plot(epoch, train.regr.crps, type = 'l', ylim = c(0, max(train.regr.crps)))
  points(epoch, valid.regr.crps, type = 'l', col = 'blue')
  plot(epoch, train.regr.mae, type = 'l', ylim = c(0, max(train.regr.mae)))
  points(epoch, valid.regr.mae, type = 'l', col = 'blue')
})

with(learner$model$callbacks$history, {
  par(mfrow = c(1, 2))
  plot(epoch, train.regr.crps, type = 'l', ylim = c(0, max(train.regr.crps)))
  # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
  plot(epoch, log(train.regr.crps), type = 'l', ylim = range(log(train.regr.crps)))
  points(epoch, log(valid.regr.crps), type = 'l', col = 'blue')
})

learner$model$callbacks$history %>%
  tail(400) %>%
  with({
    plot(epoch, log(train.regr.crps), type = 'l')
    points(epoch, log(valid.regr.crps), type = 'l', col = 'blue')
  })

learner$model$callbacks$history %>%
  tail(500) %>%
  with({
    plot(epoch, train.regr.crps, type = 'l')
    # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
  })




#### What the absolute fuck???
task2 = prepare_load_task('6B', 0, use_gam = FALSE)$data() %>%
  scale %>%
  subset(select = c(TV, load)) %>%
  as_task_regr(target = 'load')

learner2 <- LearnerTorchMLPDistr$new(
  'regr',
  optimizer = t_opt("adam", lr = 0.01, weight_decay = 0),
  # optimizer = optim_sgd,
  loss = nn_normal_nll_loss,
  callbacks = t_clbk("history")
)
learner2$param_set$set_values(
                      epochs = 300, batch_size = 128, device = "cpu",
                      neurons = integer(),
                      # p = 0,
                      # opt.lr = .05,
                      # opt.weight_decay = 0,
                      measures_train = msrs(c("regr.mae", "regr.crps", 'regr.logloss')),
                      measures_valid = msrs(c("regr.mae", "regr.crps", 'regr.logloss'))
                   )

# tlearner$graph$pipeops$regr.mlplss$learner$validate = .2
set_validate(learner2, .2)
learner2$train(task2)
learner = learner2

pred <- learner2$predict(task2)
pred$score(msrs(c('regr.crps', 'regr.mae', 'regr.rmse', 'regr.mape')))

plot(pred$truth, pred$response)
abline(0, 1)
# it's just predicting the mean load!

plot(task2$data()$TV, pred$response)
abline(0, 1)
# it's just predicting the mean load!

plot(task2$data()$TV, pred$se)

head(learner$model$callbacks$history)
tail(learner$model$callbacks$history)

with(learner$model$callbacks$history, {
  y_range = range(c(train.regr.logloss, valid.regr.logloss))
  plot(epoch, train.regr.logloss, type = 'l', ylim = y_range)
  points(epoch, valid.regr.logloss, type = 'l', col = 'blue')
})

learner$model$callbacks$history %>%
  tail(250) %>%
  with({
    y_range = range(c(train.regr.logloss, valid.regr.logloss))
    plot(epoch, train.regr.logloss, type = 'l', ylim = y_range)
    points(epoch, valid.regr.logloss, type = 'l', col = 'blue')
  })

with(learner$model$callbacks$history, {
  y_range = range(c(train.regr.crps, valid.regr.crps))
  plot(epoch, train.regr.crps, type = 'l', ylim = y_range)
  points(epoch, valid.regr.crps, type = 'l', col = 'blue')
})

learner$model$callbacks$history %>%
  tail(1500) %>%
  with({
    y_range = range(c(train.regr.crps, valid.regr.crps))
    plot(epoch, train.regr.crps, type = 'l', ylim = y_range)
    points(epoch, valid.regr.crps, type = 'l', col = 'blue')
  })

learner$model$callbacks$history %>%
  tail(1500) %>%
  with({
    plot(epoch, train.regr.crps, type = 'l')
  })

# with(learner$model$callbacks$history, {
#   par(mfrow = c(1, 2))
#   plot(epoch, train.regr.crps, type = 'l', ylim = c(0, max(train.regr.crps)))
#   points(epoch, valid.regr.crps, type = 'l', col = 'blue')
#   plot(epoch, train.regr.mae, type = 'l', ylim = c(0, max(train.regr.mae)))
#   points(epoch, valid.regr.mae, type = 'l', col = 'blue')
# })






bgrid = benchmark_grid(
    tasks = task_sc,
    #learners = ngr,
    learners = c(tlearner),
    resamplings = rsmp('forecast_holdout', ratio = 5 / 6)
)
system.time(bres <- progressr::with_progress(benchmark(bgrid, store_models = TRUE)))
bres$aggregate(c(msr('regr.mape'), msr('regr.mae')))


g1 = po("scale") %>>%
  learner %>%
  as_learner
g1$train(task1)

pred <- g1$predict(task1)

head(g1$graph_model$pipeops$regr.mlp$learner_model$model$callbacks$history)
tail(g1$graph_model$pipeops$regr.mlp$learner_model$model$callbacks$history)

with(g1$graph_model$pipeops$regr.mlp$learner_model$model$callbacks$history, {
  par(mfrow = c(1, 2))
  plot(epoch, train.regr.crps, type = 'l', ylim = c(0, max(train.regr.crps)))
  # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
  plot(epoch, train.regr.rmse, type = 'l', ylim = c(0, max(train.regr.rmse)))
})

with(g1$graph_model$pipeops$regr.mlp$learner_model$model$callbacks$history, {
  plot(epoch, train.regr.crps, type = 'l', ylim = c(0, max(train.regr.crps)))
  # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
})



# Ok forget all the GAM stuff for now. Just running a neural network. Including
# time predictor to allow for nonstationarity. Will have to run incremental
# updates over time. For this version, neither NGR nor random forest make much
# sense (nonlinear main effect, and requires extrapolation)

# just doing an MLP now, but will want to change to DRN. Also why does MLP only
# support dropout? Probably need to extend base torch class to support
# distributional regression
lrn("regr.mlp", ...)

# needed for evaluating nonstationary model
rsmp("forecast_cv")

learner_mlp = lrn("regr.mlp",
  # defining network parameters
  activation     = nn_relu,
  neurons        = integer(),
  p = 0,
  # training parameters
  batch_size     = 128,
  epochs         = 300,
  device         = "cpu",
  # Proportion of data to use for validation
  validate = NULL,
  # Defining the optimizer, loss, and callbacks
  optimizer      = t_opt("adam", lr = 0.05, weight_decay = 0),
  #optimizer      = t_opt("sgd", lr = 0.05, weight_decay = 0),
  loss           = t_loss("mse"),
  callbacks      = t_clbk("history"), # this saves the history in the learner
  # Measures to track
  measures_valid = msrs(c("regr.rmse")),
  measures_train = msrs(c("regr.rmse")),
  # predict type (required by logloss)
  predict_type = "response"
)
set_validate(learner_mlp, .2)

# task0 = prepare_load_task_combined('6B', use_gam = FALSE)
# task0 = prepare_load_task('6B', 0, use_gam = FALSE)
task0 = prepare_load_task('6B', 0, use_gam = FALSE)$data() %>%
  scale %>%
  subset(select = c(TV, load)) %>%
  as_task_regr(target = 'load')

learner_mlp$train(task0)

pred <- learner_mlp$predict(task0)
pred$score(msrs(c('regr.rmse', 'regr.mape', 'regr.mae')))

plot(pred$truth, pred$response)
abline(0, 1)

head(learner_mlp$model$callbacks$history)
tail(learner_mlp$model$callbacks$history)

with(learner_mlp$model$callbacks$history, {
  y_range = range(c(train.regr.rmse, valid.regr.rmse))
  plot(epoch, train.regr.rmse, type = 'l', ylim = y_range)
  points(epoch, valid.regr.rmse, type = 'l', col = 'blue')
})

with(learner_mlp$model$callbacks$history, plot(epoch, train.regr.rmse, type = 'b'))

g0 = po("scale") %>>% learner_mlp %>%
  as_learner

g0$train(task0)

head(learner_mlp$model$callbacks$history)

with(g0$graph_model$pipeops$regr.mlp$learner_model$model$callbacks$history, {
  plot(epoch, train.regr.rmse, type = 'l', ylim = c(0, max(train.regr.rmse)))
  # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
})

g0$graph_model$pipeops$regr.mlp$learner_model$model$callbacks$history %>%
  tail(-20) %>%
  with({
    plot(epoch, train.regr.rmse, type = 'l')
    # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
  })

# get accuracy per lead time
pg0 = g0$predict(task0)

# get accuracy by lead time
data.frame(row_ids = pg0$row_ids, response = pg0$response,
           truth = pg0$truth) %>%
  transform(days_ahead = task0$data(rows = row_ids)$days_ahead) %>%
  transform(abs_err = abs(response - truth)) %>%
  aggregate(abs_err ~ days_ahead, data = ., FUN = mean)
  # head

bgrid = benchmark_grid(
    tasks = task0,
    #learners = ngr,
    learners = c(g0),
    resamplings = rsmp('forecast_holdout', ratio = 5 / 6)
)
system.time(bres <- progressr::with_progress(benchmark(bgrid, store_models = TRUE)))
bres$aggregate(c(msr('regr.mape'), msr('regr.mae')))
# bres$aggregate(c(msr('regr.crps'), msr('regr.mape'), msr('regr.mae')))

# get the loss by lead time
bres$obs_loss() %>%
  transform(days_ahead = task0$data(rows = row_ids)$days_ahead,
            abs_err = abs(response - truth)) %>%
  aggregate(abs_err ~ days_ahead, FUN = mean, data = .)

bench_mlp = bres$score()$learner[[1]]$graph_model$pipeops$regr.mlp$learner_model

with(bench_mlp$model$callbacks$history, {
  plot(epoch, train.regr.rmse, type = 'l', ylim = c(0, max(train.regr.rmse)))
  # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
})

bench_mlp$model$callbacks$history %>%
  tail(-20) %>%
  with({
    plot(epoch, train.regr.rmse, type = 'l')
    # points(epoch, train.regr.rmse, type = 'b', col = 'blue')
  })

# https://mlr3torch.mlr-org.com/reference/mlr_learners_torch.html#saving-a-learner:
# In order to save a LearnerTorch for later usage, it is necessary to call the
# $marshal() method on the Learner before writing it to disk, as the object will
# otherwise not be saved correctly. After loading a marshaled LearnerTorch into
# R again, you then need to call $unmarshal() to transform it into a useable
# state.
