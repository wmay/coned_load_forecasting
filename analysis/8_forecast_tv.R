# use ensemble MOS with neural networks to predict ConEd's TV

# Mostly following Rasp 2018 (10.1175/MWR-D-18-0187.1) and Lerch 2022
# (https://arxiv.org/abs/2204.05102). May also try some sort of RNN for the load
# forecasting part

# Results of this Masters thesis suggest two hidden layers may be better than 1:
# http://hdl.handle.net/2429/78331

# Issues to be addressed:

# - Feature engineering:
#   - how to summarize/aggregate hourly forecasts for daily predictions
#   - including lagged values
#   - summarizing/aggregating gridded values
# - Statistical issues:
#   - regularization due to small data size
#   - forecasting for multiple nearby areas (spatially autocorrelated errors)

# A helpful trick: can subtract the lagged portion of TV from the output so it
# doesn't need to be predicted (and add it back afterward)

# To reduce hourly values to daily predictors, can probably get max temperature
# (and other values at the same time?). Eh, I'll just pick a relevant time for
# now

# Fix tensorflow 2.16.1. This must be set in the shell/emacs before running R.
# See
# https://github.com/tensorflow/tensorflow/issues/63362#issuecomment-1988630226.

# LD_LIBRARY_PATH=/home/wmay/.local/lib/python3.10/site-packages/nvidia/cublas/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_cupti/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_nvcc/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_nvrtc/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_runtime/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cudnn/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cufft/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/curand/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cusolver/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cusparse/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/nccl/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/nvjitlink/lib

# set this if you don't want an annoying prompt from reticulate
# Sys.setenv(RETICULATE_PYTHON = '/usr/bin/python3')

# May want to shuffle minibatches. See https://arxiv.org/pdf/1206.5533, p. 6

library(magrittr)
library(stars) # also requires ncmeta
library(keras3)
library(kerastuneR)
library(beepr)
source('R/load_data.R')
source('R/coned_tv.R')

get_effective_temp = function(x) {
  drywet = with(x, dry_wet_ave(tair, relh / 100, pres * 100))
  eff_temp = rolling_mean3(x$time, drywet, 'hours')
  hour_of_day = as.integer(format(x$time, '%H'))
  data.frame(stid = x$stid, time = x$time, eff_temp = eff_temp) %>%
    subset(hour_of_day >= 9 & hour_of_day <= 21)
}

# Get combined TV by averaging site dry+wet temperatures first and *then*
# getting the daily maximum of the averaged values. This ensures that the values
# are all taken from the same time of day
# IMPORTANT: Can't remove holidays prior to this calculation bc they're needed for lagged values
get_combined_tv = function(stids) {
  cols = paste0('eff_temp.', stids)
  combined_temps = rowMeans(et_wide[, cols, drop = FALSE])
  out = data.frame(time = et_wide$time, eff_temp = combined_temps) %>%
    transform(day = as.Date(time, tz = 'EST5EDT')) %>%
    # should `na.rm` be true here? Not sure
    aggregate(eff_temp ~ day, ., max, na.rm = TRUE) %>%
    transform(tv = rolling_mean3(day, eff_temp, 'days', c(.7, .2, .1))) %>%
    subset(select = c(day, tv))
  row.names(out) = as.character(out$day)
  out
}

make_predictor_dataset = function(nc, lon, lat, ltime) {
  if (!is.null(attr(nc, 'dimensions')$refDate$values)) {
    refDates = attr(nc, 'dimensions')$refDate$values
  } else {
    refDates = with(attr(nc, 'dimensions')$refDate, {
      offset + (seq(from, to) - 1) * delta
    })
  }
  times = refDates + as.difftime(3 * ltime, units = 'hours')
  attr(times, 'tzone') = 'America/New_York'
  out = lapply(nc, function(x) x[lon, lat, ltime, ]) %>%
    as.data.frame
  cbind(time = times, out)
}

# based on Rasp+Lerch 2018,
# https://github.com/slerch/ppnn/blob/7af9af3cfb754b8b6da54d2fe5d917cd54e32b9d/nn_postprocessing/nn_src/losses.py#L14
crps_loss = function(y_true, y_pred) {
  # Split y_pred into mean and std dev
  mu = y_pred[, 1]
  sigma = op_abs(y_pred[, 2])
  # Calculate the standardized difference
  z = (y_true[, 1] - mu) / sigma
  # Calculate the CDF and PDF of the standard normal distribution
  pdf = op_exp(-op_square(z) / 2) / sqrt(2 * pi)
  cdf = 0.5 * (1 + op_erf(z / sqrt(2)))
  # Calculate the CRPS components
  term1 = z * (2 * cdf - 1)
  term2 = 2 * pdf - 1 / sqrt(pi)
  crps = sigma * (term1 + term2)
  return(op_mean(crps))
}

# Based on Rasp+Lerch 2018:
# https://github.com/slerch/ppnn/blob/7af9af3cfb754b8b6da54d2fe5d917cd54e32b9d/nn_postprocessing/nn_src/losses.py#L14
# however here I'm using a log link for sigma
crps_loss2 = function(y_true, y_pred) {
  # Split y_pred into mean and std dev
  mu = y_pred[, 1]
  sigma = op_exp(y_pred[, 2])
  # Calculate the standardized difference
  z = (y_true[, 1] - mu) / sigma
  # Calculate the CDF and PDF of the standard normal distribution
  pdf = op_exp(-op_square(z) / 2) / sqrt(2 * pi)
  cdf = 0.5 * (1 + op_erf(z / sqrt(2)))
  # Calculate the CRPS components
  term1 = z * (2 * cdf - 1)
  term2 = 2 * pdf - 1 / sqrt(pi)
  crps = sigma * (term1 + term2)
  return(op_mean(crps))
}

# some neural network helper functions
fit_info = function(x, m, predictors) {
  cur_metrics = tail(as.data.frame(x$metrics), 1)
  min_iter = which.min(x$metrics$val_loss)
  nonzero_weights = lapply(m$weights[c(4, 6)], function(y) {
    table(as.array(y) < .01)
  })
  # how wide is the 95% confidence interval?
  predictions = predict(m, predictors)
  prediction_cis = summary(abs(predictions[, 2]) * 1.96 * 2)
  list(current_fit = cur_metrics, best_iter = min_iter,
       zero_weights = nonzero_weights, prediction_95int = prediction_cis)
}

plot_influences = function(m, predictors) {
  predictor_influence = m$weights[[4]] %>%
    as.array %>%
    abs %>%
    rowSums %>%
    setNames(colnames(predictors)) %>%
    sort(decreasing = T)
  orig_mar = par('mar')
  on.exit(par(mar = orig_mar))
  par(mar = c(5.1, 7.1, 4.1, 2.1))
  barplot(rev(predictor_influence), horiz = T, las = 1,
          main = 'Predictor Influence (sum of weights)')
  grid(ny = NA)
}

plot_influences2 = function(m, predictors) {
  predictor_influence = abs(m$weights[[4]]) %*% abs(m$weights[[6]]) %>%
    as.matrix
  row.names(predictor_influence) = colnames(predictors)
  colnames(predictor_influence) = c('mean', 'sd')
  orig_mar = par('mar')
  orig_mfrow = par('mfrow')
  on.exit(par(mar = orig_mar, mfrow = orig_mfrow))
  par(mar = c(5.1, 7.1, 4.1, 2.1), mfrow = c(1, 2))
  mean_influence = sort(predictor_influence[, 1])
  barplot(mean_influence, horiz = T, las = 1,
          main = 'Point Predictor Influence (sum of weights)')
  grid(ny = NA)
  sd_influence = sort(predictor_influence[, 2])
  barplot(sd_influence, horiz = T, las = 1,
          main = 'Interval Predictor Influence (sum of weights)')
  grid(ny = NA)
}

# needed: TV values, weather predictors

data_wide = readRDS('results/load_vs_weather/tv_and_load.rds')

# set up the effective temp data
tv_vars = c('tair', 'relh', 'pres')
nysm_et = get_combined_nc_data(tv_vars, hour_to_nc_index(7:21)) %>%
  subset(!(is.na(tair) | is.na(relh) | is.na(pres))) %>%
  split(.$stid) %>%
  lapply(get_effective_temp) %>%
  do.call(rbind, .)

asos_et = get_hourly_asos_data(7:21) %>%
  transform(stid = station, tair = as_celsius(tmpf),
            pres = alti * 33.8637526, time = valid_hour) %>%
  subset(!(is.na(tair) | is.na(relh) | is.na(pres))) %>%
  split(.$stid) %>%
  lapply(get_effective_temp) %>%
  do.call(rbind, .)

et_wide = rbind(nysm_et, asos_et) %>%
  reshape(direction = 'wide', idvar = 'time', timevar = 'stid')
# %>%
# subset(as.Date(time, tz = 'EST5EDT') %in% data_wide$day)
et_wide = et_wide[order(et_wide$time), ]


selected_sites = readRDS('results/select_stations/selected_sites.rds')

system_tv_sites = c("BKNYRD", "BKMAPL", "QNDKIL", "JFK", "SIFKIL", "LGA",
                    "QNSOZO", "MHMHIL")
system_tv = get_combined_tv(system_tv_sites)
# only 400 of these, clearly I should collect more weather data for this step


gefs = read_ncdf('data/gefs_nick.nc')
# want to get 9PM utc, which is 4pm EST. Starting from 12pm utc/7am est
ds1 = make_predictor_dataset(gefs, 3, 4, 11)

# subset to match outputs
ds2 = subset(ds1, as.POSIXlt(time)$hour %in% 16:17) %>%
  transform(day = as.Date(time)) %>%
  merge(system_tv) %>%
  subset(select = -c(day, time)) %>%
  na.omit


# before fitting a model, want a baseline-- how accurate is uncorrected ensemble
# forecast?

# ok this is genuinely difficult to calculate. Need to get effective temps from
# forecasts, and 2-day lagged effective temp from observations (for each
# forecast)



# # single hidden layer
# model = keras_model_sequential() %>%
#   layer_dense(16, activation = 'relu') %>%
#   layer_dense(2)

# # add predictor normalization
# layer <- layer_normalization(axis=-1)

# model = keras_model_sequential() %>%
#   layer_normalization(axis=-1) %>%
#   adapt() %>%
#   layer_dense(16, activation = 'relu') %>%
#   layer_dense(2)

# model = keras_model_sequential(input_shape = 10) %>%
#   layer_dense(128, activation = 'relu') %>%
#   layer_dense(1)

# cmodel = compile(model, loss = 'mse', metrics = "accuracy")



# compile(model, loss = crps_loss, metrics = c('mse', 'mae'))

# # predictors = as.matrix(ds2[, -ncol(ds2)])
# predictors = as.matrix(ds2[, c('TMP', 'TMP_sd')])
y = cbind(ds2$tv, 1)


# fit1 = fit(model, predictors, y, epochs = 100, verbose = 0)

# plot(fit1)
# # oh hey honestly that's not bad!

# as_kelvin = function(x) (x - 32) * (5/9) + 273.15
# as_fahrenheit = function(x) (x - 273.15) * (9/5) + 32

# # How do predictions change with temperature forecasts?
# check_data = as.matrix(data.frame(TMP = as_kelvin(seq(50, 85)), TMP_sd = 1))
# predictions = predict(model, check_data)
# plot(as_fahrenheit(check_data[, 'TMP']), predictions[, 1], type = 'b')
# grid()

# # How do predictions change with temperature forecast uncertainty?
# check_data = as.matrix(data.frame(TMP = as_kelvin(75), TMP_sd = seq(.5, 2.5, by = .1)))
# predictions = predict(model, check_data)
# plot(check_data[, 'TMP_sd'], predictions[, 2], type = 'b')
# grid()
# # change is small but in correct direction


# add predictor normalization
predictors = as.matrix(ds2[, c('TMP', 'TMP_sd', 'DPT', 'DPT_sd')])

normalize_predictors = layer_normalization(axis = -1) %>%
  adapt(predictors)

build_model1 = function(hp) {
  learning_rate = hp$Float('learning_rate', min_value = .001, max_value = .9,
                           step = sqrt(10), sampling = 'log')
  layer1_reg = hp$Float('layer1_reg', min_value = .001, max_value = 5,
                        step = sqrt(10), sampling = 'log')
  keras_model_sequential() %>%
    normalize_predictors %>%
    layer_dense(16, activation = 'relu',
                kernel_regularizer = regularizer_l2(layer1_reg)) %>%
    layer_dense(2) %>%
    compile(optimizer = optimizer_adam(learning_rate), loss = crps_loss2,
            metrics = 'mae')
}

# this makes more sense with a higher factor, otherwise it runs tests with only
# 2 epochs! That's nonsensically small
tuner1 = Hyperband(build_model1, 'val_loss', max_epochs = 800, factor = 15,
                   hyperband_iterations = 50, overwrite = T)

search_summary(tuner1)

ind_validate = sample(nrow(ds2), round(nrow(ds2) / 5))
ind_train = setdiff(1:nrow(ds2), ind_validate)
x_train = predictors[ind_train, ]
x_validate = predictors[ind_validate, ]
y_train = y[ind_train, ]
y_validate = y[ind_validate, ]

system.time({
  fit_tuner(tuner1, x_train, y_train,
            validation_data = list(x_validate, y_validate), verbose = 0)
})
#    user  system elapsed 
# 267.449   3.896 272.212
# This seems to take about 5 minutes regardless of the parameters I choose

results_summary(tuner1, num_trials = 8)
# I don't understand this method at all


# compare with grid search

tuner1_rand = RandomSearch(build_model1, 'val_loss', max_trials = 20,
                           overwrite = T)

system.time({
  fit_tuner(tuner1_rand, x_train, y_train, epochs = 800,
            validation_data = list(x_validate, y_validate), verbose = 0)
})
#     user   system  elapsed 
# 1018.875   71.736 1017.587 

results_summary(tuner1_rand, num_trials = 8)
# Results summary
# Results in ./untitled_project
# Showing 8 best trials
# Objective(name="val_loss", direction="min")

# Trial 04 summary
# Hyperparameters:
# learning_rate: 0.0316227766016838
# layer1_reg: 0.010000000000000002
# Score: 0.878132164478302


plot_tuner(tuner1_rand) # bad and ugly









# model2 = keras_model_sequential() %>%
#   normalize_predictors %>%
#   layer_dense(16, activation = 'relu') %>%
#   layer_dense(2)

# compile(model2, loss = crps_loss, metrics = 'mae')

# fit2 = fit(model2, predictors, y, validation_split = 0.2, epochs = 5000,
#            verbose = 0)

# plot(fit2, smooth = F)

# # ok this got a bit weird with all years
# fit_info(fit2, model2, predictors)
# #          loss      mae val_loss val_mae
# # 5000 1.099267 1.363484 1.179444 1.18703
# # great!

# # $prediction_95int
# #     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# #  0.07255  2.09699  4.20486  5.83757  8.39596 21.42089
# # a bit extreme I'd say

# compile(model2, loss = crps_loss2, metrics = 'mae')

# fit2l = fit(model2, predictors, y, validation_split = 0.2, epochs = 5000,
#            verbose = 0)

# plot(fit2l, smooth = F)

# fit_info(fit2l, model2, predictors)
# #          loss      mae  val_loss   val_mae
# # 5000 1.079098 1.015761 0.9262024 0.9026378
# # great


# # Include all predictors
# predictors3 = as.matrix(ds2[, -ncol(ds2)])

# normalize_predictors3 = layer_normalization(axis = -1) %>%
#   adapt(predictors3)

# model3 = keras_model_sequential() %>%
#   normalize_predictors3 %>%
#   layer_dense(128, activation = 'relu') %>%
#   layer_dense(2)

# compile(model3, loss = crps_loss, metrics = 'mae')

# fit3 = fit(model3, predictors3, y, validation_split = 0.2, epochs = 500,
#            verbose = 0)

# plot(fit3, smooth = F)


# Include all predictors-- add regularization
predictors4 = as.matrix(ds2[, -ncol(ds2)])
normalize_predictors4 = layer_normalization(axis = -1) %>%
  adapt(predictors4)

# model4 = keras_model_sequential() %>%
#   normalize_predictors4 %>%
#   layer_dense(16, activation = 'relu',
#               kernel_regularizer = regularizer_l1(.03)) %>%
#   layer_dense(2)
  #layer_dense(2, kernel_regularizer = regularizer_l1(.01))
# L2 regularization seems like a bad idea for the last layer, because the s.d.
# inputs are generally smaller than the point estimate inputs, and we don't want
# to distort that

x_train4 = predictors4[ind_train, ]
x_validate4 = predictors4[ind_validate, ]

build_model4 = function(hp) {
  learning_rate = hp$Float('learning_rate', min_value = .001, max_value = .9,
                           step = sqrt(10), sampling = 'log')
  layer1_reg = hp$Float('layer1_reg', min_value = .001, max_value = 5,
                        step = sqrt(10), sampling = 'log')
  layer1_n = hp$Int('layer1_n', min_value = 16, max_value = 32, step = 8)
  regularizer_type = hp$Int('regularizer_type', min_value = 1, max_value = 2,
                            step = 1)
  regularizer_fun = if (regularizer_type == 1) {
                      regularizer_l1
                    } else {
                      regularizer_l2
                    }
  keras_model_sequential() %>%
    normalize_predictors4 %>%
    layer_dense(layer1_n, activation = 'relu',
                kernel_regularizer = regularizer_fun(layer1_reg)) %>%
    layer_dense(2) %>%
    compile(optimizer = optimizer_adam(learning_rate), loss = crps_loss2,
            metrics = 'mae')
}

tuner4 = RandomSearch(build_model4, 'val_loss', max_trials = 20, overwrite = T)
# have to run `rm(tuner4)` if build function is changed

system.time({
  fit_tuner(tuner4, x_train4, y_train, epochs = 800,
            validation_data = list(x_validate4, y_validate), verbose = 0)
})
beep(3)
#    user  system elapsed 
# 941.582  72.516 959.080

results_summary(tuner4, num_trials = 8)
# Trial 19 summary
# Hyperparameters:
# learning_rate: 0.010000000000000002
# layer1_reg: 0.10000000000000003
# layer1_n: 16
# regularizer_type: 1
# Score: 1.391445279121399
# how is this worse than before. Omfg!!! Do I need more epochs?

m4_1 = get_best_models(tuner4, 1)[[1]]

evaluate(m4_1, x_validate4, y_validate)
# $loss
# [1] 1.319056

# $mae
# [1] 1.029272

# does it just need more training?

trainhist4_1 = fit(m4_1, x_train4, y_train,
                   validation_data = list(x_validate4, y_validate),
                   epochs = 1000, verbose = 0)

plot(trainhist4_1)

plot_influences2(m4_1, predictors4)




# compile(model4, loss = crps_loss2, metrics = 'mae')

# fit4 = fit(model4, predictors4, y, validation_split = 0.2, epochs = 5000,
#            verbose = 0)

# plot(fit4, smooth = T)

fit_info(fit4, model4, predictors4)
#          loss       mae val_loss   val_mae
# 5000 1.425788 0.9682786 1.411143 0.9551771
# current best: l1(.04), l1(.01), 10k epochs

#          loss       mae val_loss  val_mae
# 5000 1.000904 0.9635538 1.128609 1.040081
# l1(.03), 10k epochs

#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.005414 0.547926 1.217413 1.499166 2.097777 5.974933
# this is bonkers high accuracy-- clearly should prioritize CRPS over MAE

# Still slowly decreasing CRPS at 15k epochs, could go longer

# So strange that this model is consistently worse on CRPS than the other, even
# though it beats it on MAE

plot_influences(model4, predictors4)

plot_influences2(model4, predictors4)


as.data.frame(fit4$metrics)[which.min(fit4$metrics$val_mae), ]

# so this can work, just needs lots of regularization for the first layer, and
# lots of training (~20k epochs). Might need a lower learning rate?

# This isn't as good as before with the 2 predictors, but I probably shouldn't
# sweat it with such a small dataset. Let's try later with more data. Update:
# More data fixed it!
