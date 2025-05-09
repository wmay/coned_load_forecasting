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
#   - time-varying load/weather relationships

# Potentially big issue!: Can add April to prepare for a new year, but April and
# May are distorted by cogen maintenance! That could trip up yearly adjustments!

# Probably need regularization due to limited training data. Could use early
# stopping or add regularization to the nodes. Can also try connecting some
# predictors directly to outputs, skipping the hidden layer (implying a linear
# effect only).

# To reduce hourly values to daily predictors, can probably get max temperature
# (and other values at the same time?). Eh, I'll just pick a relevant time for
# now

# Fix tensorflow 2.16.1. This must be set in the shell/emacs before running R.
# See
# https://github.com/tensorflow/tensorflow/issues/63362#issuecomment-1988630226.

# LD_LIBRARY_PATH=/home/wmay/.local/lib/python3.10/site-packages/nvidia/cublas/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_cupti/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_nvcc/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_nvrtc/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_runtime/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cudnn/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cufft/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/curand/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cusolver/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cusparse/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/nccl/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/nvjitlink/lib

library(magrittr)
library(stars) # also requires ncmeta
library(keras3)
library(kerastuneR)
source('R/load_data.R')

# Can use for comparison-- accuracy of polynomial models

# CRPS with a log link function for sigma
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
  prediction_cis = summary(exp(predictions[, 2]) * 1.96 * 2)
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


# needed: load values, weather predictors

loads = read.csv('data/coned/Borough and System Data 2020-2024.csv') %>%
  transform(DT = as.POSIXct(DT, tz = 'EST5EDT', '%m/%d/%Y %H:%M'),
            # some inconsistency, but various forms of 'False' mean data is fine
            BAD = !startsWith(BAD, 'F')) %>%
  subset(!is.na(DT))
names(loads) = tolower(names(loads))

peaks = loads %>%
  subset(!bad) %>%
  transform(day = as.Date(dt, tz = 'EST5EDT')) %>%
  aggregate(reading ~ borough + day, ., max)

system_peaks = peaks %>%
  subset(borough == 'CECONY') %>%
  subset(as.integer(format(day, "%m")) >= 5 &
         as.integer(format(day, "%m")) < 10)


# only 400 of these, clearly I should collect more weather data for this step

gefs = read_ncdf('data/gefs_nick.nc')

# gefs = read_ncdf('data/gefs_nick_old.nc')
# # refDates are wrong due to bug in code
# n_refDates = length(attr(gefs,'dimensions')$refDate$values)
# correct_refDates = as.POSIXct('2023-01-01', tz = 'UTC') +
#   as.difftime((1:n_refDates - 1) * 12, units = 'hours')
# attr(gefs,'dimensions')$refDate$values = correct_refDates

# test run

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


# want to get 9PM utc, which is 4pm EST. Starting from 12pm utc/7am est
ds1 = make_predictor_dataset(gefs, 3, 4, 11)
# ds2021 = make_predictor_dataset(gefs2021, 3, 4, 11)
# ds2023 = make_predictor_dataset(gefs, 3, 4, 11)
# ds1 = rbind(ds2021, ds2023)

# subset to match outputs
ds2 = subset(ds1, as.POSIXlt(time)$hour %in% 16:17) %>%
  transform(day = as.Date(time)) %>%
  merge(system_peaks) %>%
  # add time predictor
  transform(date = as.integer(day - as.Date('2021-01-01'))) %>%
  subset(select = -c(day, time, borough)) %>%
  na.omit


### Baseline models

# before fitting a model, want a baseline-- how accurate is uncorrected ensemble
# forecast?

# what about a GAM?

library(mgcv)

g0 = gam(reading ~ s(TMP, DPT), data = ds2)

mean(abs(resid(g0)))
# [1] 519.2947
sqrt(mean(resid(g0)^2))
# [1] 661.6567

g1 = gam(reading ~ s(TMP, DPT) + s(date), data = ds2)

mean(abs(resid(g1)))
# [1] 495.4976
sqrt(mean(resid(g1)^2))
# [1] 622.7792
# wow this is surprisingly bad compared to our neural networks!

plot(g1) # seems to be working correctly
plot(g1, select = 2) # strong seasonal variation!

# ok let's gam harder
g2 = gam(reading ~ s(TMP, DPT, TMP_sd, DPT_sd) + s(date), data = ds2)

mean(abs(resid(g2))) # only slightly better?
sqrt(mean(resid(g2)^2))



# # single hidden layer

# based on Rasp+Lerch 2018,
# https://github.com/slerch/ppnn/blob/7af9af3cfb754b8b6da54d2fe5d917cd54e32b9d/nn_postprocessing/nn_src/losses.py#L14
# crps_loss = function(y_true, y_pred) {
#   # Split y_pred into mean and std dev
#   mu = y_pred[, 1]
#   sigma = op_abs(y_pred[, 2])
#   # Calculate the standardized difference
#   z = (y_true[, 1] - mu) / sigma
#   # Calculate the CDF and PDF of the standard normal distribution
#   pdf = op_exp(-op_square(z) / 2) / sqrt(2 * pi)
#   cdf = 0.5 * (1 + op_erf(z / sqrt(2)))
#   # Calculate the CRPS components
#   term1 = z * (2 * cdf - 1)
#   term2 = 2 * pdf - 1 / sqrt(pi)
#   crps = sigma * (term1 + term2)
#   return(op_mean(crps))
# }

y = cbind(ds2$reading, 1)

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
model2 = keras_model_sequential() %>%
  normalize_predictors %>%
  layer_dense(16, activation = 'relu') %>%
  layer_dense(2)

compile(model2, optimizer = optimizer_adam(.005), loss = crps_loss2,
        metrics = 'mae')

fit2 = fit(model2, predictors, y, validation_split = 0.2, epochs = 1000,
           verbose = 0)

plot(fit2, smooth = F)

fit_info(fit2, model2, predictors)
#          loss      mae val_loss  val_mae
# 4000 373.2233 264.8913 355.6738 245.9689

# $prediction_95int
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   815.9  2313.8  2749.0  2665.1  3101.7  4043.0 


# Add a time (long-term time trend) predictor
predictors3 = as.matrix(ds2[, c('TMP', 'TMP_sd', 'DPT', 'DPT_sd', 'date')])

normalize_predictors = layer_normalization(axis = -1) %>%
  adapt(predictors3)
model3 = keras_model_sequential() %>%
  normalize_predictors %>%
  layer_dense(16, activation = 'relu') %>%
  layer_dense(2)

compile(model3, loss = crps_loss2, metrics = 'mae')

fit3 = fit(model3, predictors3, y, validation_split = 0.2, epochs = 4000,
           verbose = 0)

plot(fit3, smooth = F)
fit_info(fit3, model3, predictors3)
#          loss      mae val_loss  val_mae
# 4000 389.5543 278.9001 369.7739 257.3088
# Nice. This will definitely overfit though if trained too much


# Include all predictors-- add regularization
predictors4 = as.matrix(subset(ds2, select = -reading))
normalize_predictors4 = layer_normalization(axis = -1) %>%
  adapt(predictors4)
model4 = keras_model_sequential() %>%
  normalize_predictors4 %>%
  layer_dense(16, activation = 'relu',
              kernel_regularizer = regularizer_l1(.5)) %>%
  layer_dense(2)
# L2 regularization seems like a bad idea for the last layer, because the s.d.
# inputs are generally smaller than the point estimate inputs, and we don't want
# to distort that

# set the learning rate higher, following rasp and lerch
compile(model4, optimizer = optimizer_adam(.01), loss = crps_loss2,
        metrics = 'mae')

fit4 = fit(model4, predictors4, y, validation_split = 0.2, epochs = 800,
           verbose = 0)

plot(fit4, smooth = T)

fit_info(fit4, model4, predictors4)
#          loss      mae val_loss  val_mae
# 4000 415.8746 264.3355 410.8071 246.7273
# current best: l2(.5), l1(.01), 16 hidden nodes, 8k epochs
# 159 non-zero weights! That's a lot! should try adding more nodes

# overfitting with l2(.04), l1(0), 32 hidden nodes, 8k epochs

# $prediction_95int
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 170.5  1497.2  2112.6  2199.8  2781.3  8322.4 
# obviously a little too small

# this one is tough nut to crack!

plot_influences(model4, predictors4)

plot_influences2(model4, predictors4)

# so this can work, just needs lots of regularization for the first layer, and
# lots of training (~20k epochs). Might need a lower learning rate?




# ok let's do more formal tuning



build_model4 = function(hp) {
  learning_rate = hp$Float('learning_rate', min_value = .001, max_value = .1, step = 5,
                           sampling = 'log')
  layer1_reg = hp$Float('layer1_reg', min_value = .01, max_value = 5, step = 5,
                        sampling = 'log')
  
  keras_model_sequential() %>%
    normalize_predictors4 %>%
    layer_dense(16, activation = 'relu',
                kernel_regularizer = regularizer_l2(layer1_reg)) %>%
    layer_dense(2) %>%
    compile(optimizer = optimizer_adam(learning_rate), loss = crps_loss2,
            metrics = 'mae')
}

tuner4 = RandomSearch(build_model4, objective = 'val_loss', max_trials = 20,
                      overwrite = T)

search_summary(tuner4)

ind_validate = sample(nrow(ds2), round(nrow(ds2) / 5))
ind_train = setdiff(1:nrow(ds2), ind_validate)
x_train = predictors4[ind_train, ]
x_validate = predictors4[ind_validate, ]
y_train = y[ind_train, ]
y_validate = y[ind_validate, ]

fit_tuner(tuner4, x_train, y_train, epochs = 800,
          validation_data = list(x_validate, y_validate), verbose = 0)

results_summary(tuner4, num_trials = 10)
# Trial 04 summary
# Hyperparameters:
# learning_rate: 0.025
# layer1_reg: 0.05
# Score: 264.57708740234375

# wow! tuning helped!
