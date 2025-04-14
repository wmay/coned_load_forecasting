# Forecast the load curve using a GAM

# This is a state space problem which is technically different than a GAM, but
# the GAM results should be very similar, and it's much easier to run

library(magrittr)
library(timeDate) # holidays
library(mgcv)
library(mgcViz)
source('R/load_data.R')
source('R/coned_tv.R')

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
  subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2023))) %>%
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
  # this sometimes fails due to a singular matrix
  fit = try(gamm(reading ~ s(tv) + s(nday) + s(tv, nday), data = dat))
  if (inherits(fit, 'try-error')) {
    # try a simpler model
    fit = gamm(reading ~ s(tv) + s(tv, nday), data = dat)
  }
  fit
}

data_start = as.Date('2021-04-01')

f1_day = data_start + 90

year_forecast_dates = function(year) {
  seq(as.Date(paste0(year, '-04-01')), as.Date(paste0(year, '-10-01')),
      by = 'day')
}

forecast_days = lapply(2021:2023, year_forecast_dates) %>%
  unlist %>%
  as.Date %>%
  subset(. > data_start + 90)
#forecast_days = load_curve_dat$day[load_curve_dat$day > data_start + 90]
forecast_tvs = seq(50, 86, 2)
forecast_tvs = seq(50, 86)

forecast_ses = sapply(forecast_days, function(x) {
  #print(x)
  fit = try(get_fit_as_of(x))
  if (inherits(fit, 'try-error')) {
    return(rep(NA, length(forecast_tvs)))
  }
  newdata = data.frame(tv = forecast_tvs, nday = as.integer(x))
  pred = predict(fit$gam, newdata = newdata, se.fit = T)
  pred$se.fit
})

colnames(forecast_ses) = paste0('se.', as.character(forecast_days))

forecast_ses2 = as.data.frame(forecast_ses) %>%
  transform(tv = forecast_tvs) %>%
  reshape(direction = 'long', varying = 1:(ncol(.) - 1), idvar = 'tv',
          timevar = 'day', v.names = 'se') %>%
  transform(day = forecast_days[day]) %>%
  transform(toy = as.Date('2023-01-01') + as.POSIXlt(day)$yday,
            year = as.POSIXlt(day)$year + 1900)

load_curve_dat %>%
  transform(toy = as.Date('2023-01-01') + doy,
            year = as.POSIXlt(day)$year + 1900) %>%
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
