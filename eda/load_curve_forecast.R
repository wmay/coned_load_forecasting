# Forecast the load curve using a GAM

# This is a state space problem which is technically different than a GAM, but
# the GAM results should be very similar, and it's much easier to run

library(magrittr)
library(timeDate) # holidays
library(mgcv)
library(mgcViz)
source('../R/load_data.R')
source('../R/coned_tv.R')

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

loads = read.csv('../data/coned/Borough and System Data 2020-2024.csv') %>%
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
