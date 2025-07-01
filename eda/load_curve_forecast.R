# Forecast the load curve using a GAM

# This is a state space problem which is technically different than a GAM, but
# the GAM results should be very similar, and it's much easier to run

library(magrittr)
library(timeDate) # holidays
library(ggplot2)
library(mgcv)
library(stars) # also requires ncmeta
# library(mgcViz)
source('R/load_data.R')
source('R/coned_tv.R')
source('R/forecast_dataset.R')

station_obs = get_hourly_asos_data(7:21) %>%
  transform(stid = station, time = valid_hour) %>%
  subset(!(is.na(tmpf) | is.na(dwpf))) %>%
  subset(select = c(time, stid, tmpf, dwpf)) %>%
  reshape(direction = 'wide', idvar = 'time', timevar = 'stid')

# Get combined TV by averaging site dry+wet temperatures first and *then*
# getting the daily maximum of the averaged values. This ensures that the values
# are all taken from the same time of day
# IMPORTANT: Can't remove holidays prior to this calculation bc they're needed for lagged values
get_combined_tv = function(stids) {
  tmp_cols = paste0('tmpf.', stids)
  dwp_cols = paste0('dwpf.', stids)
  tmps = rowMeans(station_obs[, tmp_cols, drop = FALSE], na.rm = TRUE)
  dwps = rowMeans(station_obs[, dwp_cols, drop = FALSE], na.rm = TRUE)
  out = coned_tv(station_obs$time, tmps, dwps)
  row.names(out) = as.character(out$day)
  out
}

# coned's traditional TV
system_tv = get_combined_tv(c('NYC', 'LGA'))

loads = read.csv('data/coned/Borough and System Data 2020-2024.csv') %>%
  transform(DT = as.POSIXct(DT, tz = 'EST5EDT', '%m/%d/%Y %H:%M'),
            # some inconsistency, but various forms of 'False' mean data is fine
            BAD = !startsWith(BAD, 'F')) %>%
  subset(!is.na(DT) & !BAD)
names(loads) = tolower(names(loads))

peaks = loads %>%
  subset(borough == 'CECONY') %>%
  transform(day = as.Date(dt, tz = 'EST5EDT')) %>%
  aggregate(reading ~ day, ., max)

load_curve_dat = merge(peaks, system_tv) %>%
  subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2024))) %>%
  transform(nday = as.integer(day), doy = as.POSIXlt(day)$yday)



# Check the density of observations in terms of doy and TV, to get a sense of
# how our knowledge of the load curve changes over the year

load_curve_dat %>%
  transform(toy = as.Date('2023-01-01') + doy) %>%
  ggplot(aes(x = toy, y = tv)) +
  geom_point(size = 1) +
  # overlay heatmap of point density
  geom_density_2d_filled(aes(fill = ..level..), bins = 5, alpha = 0.5) +
  geom_hline(yintercept = 82, linetype = 'dashed', color = 'red') +
  labs(x = 'Time of Year', y = 'TV', title = 'Observation Density')

get_fit_as_of = function(cutoff) {
  dat = subset(load_curve_dat, day < cutoff)
  # use a simpler model for the first year
  if (cutoff < as.Date('2022-01-01')) {
    model = reading ~ s(tv) + s(nday)
  } else {
    model = reading ~ s(tv) + s(nday) + s(tv, nday)
  }
  # this sometimes fails due to a singular matrix
  fit = try(gamm(model, data = dat))
  if (inherits(fit, 'try-error')) {
    # try a simpler model
    fit = gamm(reading ~ s(tv) + s(tv, nday), data = dat)
  }
  fit
}

first_fit = as.Date('2021-05-01')

year_forecast_dates = function(year) {
  seq(as.Date(paste0(year, '-04-01')), as.Date(paste0(year, '-10-01')),
      by = 'day')
}

forecast_days = lapply(2021:2024, year_forecast_dates) %>%
  unlist %>%
  as.Date %>%
  subset(. >= first_fit & . <= max(load_curve_dat$day) + 1)
forecast_tvs = seq(50, 86)

# forecast_ses = sapply(forecast_days, function(x) {
#   #print(x)
#   fit = try(get_fit_as_of(x))
#   if (inherits(fit, 'try-error')) {
#     return(rep(NA, length(forecast_tvs)))
#   }
#   newdata = data.frame(tv = forecast_tvs, nday = as.integer(x))
#   pred = predict(fit$gam, newdata = newdata, se.fit = T)
#   pred$se.fit
# })

gam_forecasts = sapply(forecast_days, function(x) {
  #print(x)
  fit = try(get_fit_as_of(x))
  if (inherits(fit, 'try-error')) {
    return(rep(NA, length(forecast_tvs)))
  }
  newdata = data.frame(tv = forecast_tvs, nday = as.integer(x))
  pred = predict(fit$gam, newdata = newdata, se.fit = T)
  cbind(fit = pred$fit, se.fit = pred$se.fit)
}, simplify = 'array')

forecast_ses = gam_forecasts[, 2, ]

colnames(forecast_ses) = paste0('se.', as.character(forecast_days))

forecast_ses2 = as.data.frame(forecast_ses) %>%
  transform(tv = forecast_tvs) %>%
  reshape(direction = 'long', varying = 1:(ncol(.) - 1), idvar = 'tv',
          timevar = 'day', v.names = 'se') %>%
  transform(day = forecast_days[day]) %>%
  transform(doy = as.POSIXlt(day)$yday,
            year = as.POSIXlt(day)$year + 1900) %>%
  # adjust day of year for leap years
  transform(toy = as.Date('2023-01-01') + doy - ifelse(year %% 4, 0, 1))

load_curve_dat %>%
  transform(year = as.POSIXlt(day)$year + 1900) %>%
  # adjust day of year for leap years
  transform(toy = as.Date('2023-01-01') + doy - ifelse(year %% 4, 0, 1)) %>%
  ggplot(aes(x = toy, y = tv)) +
  # overlay heatmap from se matrix
  geom_tile(aes(x = toy, y = tv, fill = se), data = forecast_ses2) +
  scale_fill_viridis_c(direction = -1) +
  geom_point(size = 1) +
  geom_hline(yintercept = 82, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = as.Date('2023-05-01'), linetype = 'dashed',
             color = '#555555') +
  labs(x = 'Date', y = 'TV', title = 'Model confidence') +
  facet_wrap(~year, ncol = 1)


f1 = get_fit_as_of(f1_day)

predict(f1, newdata = subset(load_curve_dat, day == f1_day))

predict(f1$gam, newdata = subset(load_curve_dat, day == f1_day), se.fit=T)

  

# remember to use `gamm` to account for serial dependence. Great example here:
# https://petolau.github.io/Analyzing-double-seasonal-time-series-with-GAM-in-R/

g1 = gam(reading ~ s(tv), data = load_curve_dat)

g2 = gamm(reading ~ s(tv) + s(nday), data = load_curve_dat)
# g2 = gam(reading ~ s(tv) + s(doy), data = load_curve_dat)
# g2 = gam(reading ~ s(tv) + s(doy) + s(nday), data = load_curve_dat)

g3 = gam(reading ~ s(tv) + s(nday) + s(tv, nday), data = load_curve_dat)

g4 = gam(reading ~ s(tv) + s(tv, nday), data = load_curve_dat)

plot(g2$gam, residuals = T, scheme = 2)


acf(resid(g2$lme, type = "normalized"), lag.max = 48)
pacf(resid(g2$lme, type = "normalized"), lag.max = 48)
# it appears that an AR(1) model is appropriate

g2 = gamm(reading ~ s(tv) + s(nday), correlation = corAR1(form = ~nday),
          data = load_curve_dat)

g3 = gamm(reading ~ s(tv) + s(nday) + s(tv, nday),
          correlation = corAR1(form = ~nday),
          data = load_curve_dat)

g3 = gamm(reading ~ s(tv) + s(doy) + s(nday) + s(tv, nday),
          correlation = corAR1(form = ~nday),
          data = load_curve_dat)

plot(g2$gam, residuals = T, scheme = 2)


# ok let's make the predictor dataset

# strangely `read_mdim` works fine while `read_ncdf` messes up the coordinates
# gefs_tv = read_ncdf('results/process_nwp_data/gefs_tv2.nc')
gefs_tv = read_mdim('results/process_nwp_data/gefs_tv2.nc')
# tv_fct = make_nwp_tv_dataset(gefs_tv, 3, 4, days_ahead = 2)
tv_fct_list = lapply(2:7, function(x) make_nwp_tv_dataset(gefs_tv, 3, 4, x))

get_gam_forecast = function(model_day, tv, days_ahead = 0) {
  fit = try(get_fit_as_of(model_day))
  n_tv = length(tv)
  n_days = length(days_ahead)
  # if (n_tv > 1 && n_days > 1) stop('Only tv or days_ahead can have length > 1')
  combine_fun = if (n_tv > 1 || n_days > 1) cbind else c
  if (inherits(fit, 'try-error')) {
    return(combine_fun(fit = rep(NA, n_tv), se.fit = rep(NA, n_tv)))
  }
  fct_days = model_day + days_ahead
  newdata = data.frame(tv = tv, nday = as.integer(fct_days))
  pred = predict(fit$gam, newdata = newdata, se.fit = T)
  combine_fun(fit = pred$fit, se.fit = pred$se.fit)
}

gam_forecasts = sapply(forecast_days, function(x) {
  valid_days = x + 2:7
  new_tv = sapply(2:7, function(i) {
    tv_fct = tv_fct_list[[i - 1]]
    tv_fct$tv[match(x + i, tv_fct$day)]
  })
  out = get_gam_forecast(x, new_tv, days_ahead = 2:7)
  # add the forecast error
  obs_load = peaks$reading[match(valid_days, peaks$day)]
  fct_err = obs_load - out[, 'fit']
  cbind(out, err = fct_err, obs = obs_load)
}, simplify = 'array')
gam_forecasts = aperm(gam_forecasts, c(3, 1, 2))

# quick check
plot(gam_forecasts[, 1, 1], gam_forecasts[, 1, 4])
plot(gam_forecasts[, 6, 1], gam_forecasts[, 6, 4])
# strong relationship between GAM predictions and observed loads, as there
# should be

# let's look at the MAPE of the GAM for reference
mae = function(obs, fct) mean(abs(obs - fct), na.rm = TRUE)
mape = function(obs, fct) 100 * mean(abs((obs - fct) / obs), na.rm = TRUE)

sapply(1:6, function(i) mape(gam_forecasts[, i, 1], gam_forecasts[, i, 4]))
# [1] 7.072099 7.368438 7.638298 7.958399 8.517402 9.036653
sapply(1:6, function(i) mae(gam_forecasts[, i, 1], gam_forecasts[, i, 4]))
# [1] 609.0118 635.2708 661.4272 691.1363 741.0981 786.4698

out = list(model_day = forecast_days, forecasts = gam_forecasts)
saveRDS(out, 'results/load_curve_forecast/gam_forecasts.rds')
