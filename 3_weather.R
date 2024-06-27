# looking at weather patterns

library(sf)
library(psychrolib)
library(ggplot2)
library(smacof) # this overwrites `transform`, ugh
transform = base::transform
source('R/load_data.R')
source('R/plots.R')
source('R/coned_tv.R')

wgs84 = 4326
web_merc = 3857

stations = readRDS('results/station_data/stations.rds') %>%
  st_transform(web_merc)

nj_stids = c('LDJ', 'EWR', 'TEB')


# want to know when peak temps (effective temp?) happen

# for simplicity, just get average across the meso/micronet
mean_temps = get_combined_nc_data('tair', hour_to_nc_index(0:23)) %>%
  aggregate(tair ~ time, ., mean, na.rm = T)

max_times = mean_temps %>%
  transform(day = trunc(time, 'days')) %>%
  by(.$day, function(x) x$time[which.max(x$tair)]) %>%
  as.integer %>%
  as.POSIXct(origin = '1970-01-01')
times_hist = hist(as.POSIXlt(max_times)$hour, 24, right = F)
saveRDS(times_hist, 'results/weather/max_tmp_times.rds')
sort(table(as.POSIXlt(max_times)$hour), decreasing = T)


# how often is there a sea breeze?

# Gedzelman 2003:

# The criterion used to define a day with a strong sea breeze is that at 1600
# local time, temperature at NYC must be at least 2.2°C warmer than at JFK, and
# the southerly component of the wind at JFK must exceed that at NYC by 5 kt.

# drct: Wind Direction in degrees from north

# degrees to radians: (x/360) * 2 * pi

# southerly component is therefore sknt * -cosine of wind direction

asos_temps = get_hourly_asos_data(16) %>%
  subset(station %in% c('JFK', 'NYC')) %>%
  transform(day = as.Date(valid_hour), tmpc = (tmpf - 32) * (5/9),
            drct_rad = (drct / 360) * 2 * pi) %>%
  transform(wind_s = -sknt * cos(drct_rad)) %>%
  subset(select = c(day, station, tmpc, wind_s)) %>%
  reshape(direction = 'wide', timevar = 'station', idvar = 'day') %>%
  # lots of missing wind values due to wind being recorded as "variable"
  transform(wind_s.JFK = replace(wind_s.JFK, is.na(wind_s.JFK), 0),
            wind_s.NYC = replace(wind_s.NYC, is.na(wind_s.NYC), 0)) %>%
  transform(jfk_windier = wind_s.JFK > 5 & wind_s.JFK - wind_s.NYC > 5,
            nyc_warmer = tmpc.NYC - tmpc.JFK >= 2.2)

sea_breeze_freq = table(with(asos_temps, jfk_windier & nyc_warmer))
sum(sea_breeze_freq)
prop.table(sea_breeze_freq)
saveRDS(sea_breeze_freq, 'results/weather/sea_breeze_freq.rds')

sea_breeze_gradients = asos_temps %>%
  subset(jfk_windier & nyc_warmer) %>%
  with((tmpc.NYC - tmpc.JFK) * 9 / 5)
summary(sea_breeze_gradients)
hist(sea_breeze_gradients, main = 'NYC - JFK temps, 4pm sea breeze days',
     freq = F, xlab = '°F', ylab = 'Proportion')
saveRDS(sea_breeze_gradients, 'results/weather/sea_breeze_gradients.rds')


# get matrix of temps at every station, every day at 4pm
nysm_temps = get_combined_nc_data('tair', hour_to_nc_index(16)) %>%
  reshape(direction = 'wide', timevar = 'stid', idvar = 'time')
asos_temps = get_hourly_asos_data(16) %>%
  transform(stid = station, tair = (tmpf - 32) * (5/9)) %>%
  subset(select = c(stid, valid_hour, tair)) %>%
  reshape(direction = 'wide', timevar = 'stid', idvar = 'valid_hour')
tair_wide = merge(nysm_temps, asos_temps, by.x = 'time', by.y = 'valid_hour')

# OK let's plot the means

tair_means = colMeans(tair_wide[, -1], na.rm = T)
names(tair_means) = sub('.*\\.', '', names(tair_means))
hist(tair_means)

stations %>%
  transform(mean_temp = tair_means[stations$stid]) %>%
  plot_station_data(aes(color = mean_temp, shape = network)) +
  scale_shape_manual(values = c(15, 19, 17))
# maybe not the best way to get means due to different missing data

apply(tair_wide[, -1], 2, function(x) sum(is.na(x))) # well, not too bad


# Ok now correlations

tair_cor = cor(as.matrix(tair_wide[, -1]), use = 'pairwise.complete.obs')
colnames(tair_cor) = sub('.*\\.', '', colnames(tair_cor))
row.names(tair_cor) = sub('.*\\.', '', row.names(tair_cor))
hist(tair_cor[lower.tri(tair_cor)])

# also let's try unfolding

# remove NAs, which cause problems
all_missing = apply(tair_cor, 1, function(x) all(is.na(x)))
tair_cor2 = tair_cor[!all_missing, !all_missing]

tair_distances = sim2diss(tair_cor2, to.dist = TRUE)
tair_locations = mds(tair_distances, type = 'interval')
site_boroughs = as.data.frame(stations)[row.names(tair_locations$conf), 'borough']
site_boroughs[is.na(site_boroughs)] = 'NA'
plot(tair_locations, col = factor(site_boroughs), cex = 2)
# that's very reasonable actually

saveRDS(tair_locations, 'results/weather/temperature_mds.rds')
# tair_locations = readRDS('results/temperature_mds.rds')

# does pca reveal any pattern similar to the results from metric unfolding?
pc0 = prcomp(tair_distances)
plot(pc0)

plot(pc0$x[, 1:2])
text(pc0$x[, 1:2], labels = row.names(pc0$x))
# not sure. dim 1 is sea breeze-ish

plot(pc0$x[, 2:3])
text(pc0$x[, 2:3], labels = row.names(pc0$x))
# dim 3 seems like west to east

plot(pc0$x[, 3:4])
text(pc0$x[, 3:4], labels = row.names(pc0$x))

plot(pc0$rotation[, 1:2])
text(pc0$rotation[, 1:2], labels = row.names(pc0$rotation))
# I don't think anything useful is coming out of this

# # what about factor analysis? -- doesn't seem to work with correlations
# fa0 = factanal(tair_distances, 3)
# fa0 = factanal(covmat = tair_cor2, factors = 3)
# plot(fa0)





# use first principal component as sea breeze coordinates
pc1 = prcomp(tair_locations$conf)
summary(pc1)
plot(tair_locations, col = factor(site_boroughs), cex = 2)
abline(0, pc1$rotation[2, 1] / pc1$rotation[1, 1], lty = 'dashed')

png('results/weather/sea_breeze_dimension.png', width = 1100, height = 800, res = 100)
stations %>%
  transform(sea_breeze = pc1$x[stations$stid, 1]) %>%
  plot_station_data(aes(color = sea_breeze, shape = network), size = 7)
dev.off()

# what's the second dimension?
stations %>%
  transform(sb2 = pc1$x[stations$stid, 2]) %>%
  plot_station_data(aes(color = sb2, shape = network)) +
  scale_shape_manual(values = c(15, 19, 17))
# west to east-- cool


# do we get the same result with TV?
SetUnitSystem('SI')
tv_vars = c('tair', 'relh')
nysm_tv = get_combined_nc_data(tv_vars, hour_to_nc_index(7:21)) %>%
  subset(!(is.na(tair) | is.na(relh))) %>%
  transform(dwp = GetTDewPointFromRelHum(tair, relh / 100)) %>%
  split(.$stid) %>%
  lapply(function(x) {
    out = with(x, coned_tv(time, as_fahrenheit(tair), as_fahrenheit(dwp)))
    out$stid = x$stid[1]
    out
  }) %>%
  do.call(rbind, .)

asos_tv = get_hourly_asos_data(7:21) %>%
  transform(time = valid_hour) %>%
  subset(!(is.na(tmpf) | is.na(dwpf))) %>%
  split(.$station) %>%
  lapply(function(x) {
    out = with(x, coned_tv(time, tmpf, dwpf))
    out$stid = x$station[1]
    out
  }) %>%
  do.call(rbind, .)

tv_wide = rbind(nysm_tv, asos_tv) %>%
  reshape(direction = 'wide', idvar = 'day', timevar = 'stid')

saveRDS(tv_wide, 'results/weather/tv.rds')

tv_cor = cor(as.matrix(tv_wide[, -1]), use = 'pairwise.complete.obs')
colnames(tv_cor) = sub('.*\\.', '', colnames(tv_cor))
row.names(tv_cor) = sub('.*\\.', '', row.names(tv_cor))
# hist(tv_cor[lower.tri(tv_cor)])

# remove NAs, which cause problems
all_missing = apply(tv_cor, 1, function(x) all(is.na(x)))
tv_cor2 = tv_cor[!all_missing, !all_missing]
tv_distances = sim2diss(tv_cor2)
tv_locations = mds(tv_distances, type = 'interval')
saveRDS(tv_locations, 'results/weather/tv_mds.rds')
site_boroughs = as.data.frame(stations)[row.names(tv_locations$conf), 'borough']
site_boroughs[is.na(site_boroughs)] = 'NA'
plot(tv_locations, col = factor(site_boroughs), cex = 2)
# add sea breeze dimension line
pc_tv = prcomp(tv_locations$conf)
summary(pc_tv)
# plot(tv_locations, col = factor(site_boroughs), cex = 2)
abline(0, pc_tv$rotation[2, 1] / pc_tv$rotation[1, 1], lty = 'dashed')


stations %>%
  transform(sea_breeze = pc_tv$x[stations$stid, 1]) %>%
  plot_station_data(aes(color = sea_breeze, shape = network), size = 7) +
  scale_shape_manual(values = c(15, 19, 17))

# what's the second dimension?
stations %>%
  transform(dim2 = pc_tv$x[stations$stid, 2]) %>%
  plot_station_data(aes(color = dim2, shape = network)) +
  scale_shape_manual(values = c(15, 19, 17))
# west to east




# # not the best plot
# p1 = stations %>%
#   transform(sea_breeze = pc1$x[stations$stid, 1]) %>%
#   plot_station_data(aes(color = sea_breeze, shape = network), size = 7) +
#   scale_shape_manual(values = c(15, 19, 17))
# p2 = stations %>%
#   transform(dim2 = pc1$x[stations$stid, 2]) %>%
#   plot_station_data(aes(color = dim2, shape = network)) +
#   scale_shape_manual(values = c(15, 19, 17))
# library(gridExtra)
# grid.arrange(p1, p2, ncol = 2)



# ########## OLD ############


# # can I get temperature biases with mgcv?

# tair_long = reshape(tair_wide, direction = 'long', varying = 2:ncol(tair_wide), timevar = 'stid',
#                     idvar = 'time')

# tair_long[, c('x', 'y')] =
#   st_coordinates(st_transform(stations, li_lcc)[tair_long$stid, ])

# library(mgcv)

# res = tair_long %>%
#   subset(time < as.POSIXct('2021-06-01')) %>%
#   transform(stid = factor(stid), time = factor(time)) %>%
#   { gam(tair ~ s(x, y, k = 10) + time + s(stid, bs = 're') + s(x, y, k = 10, by = time), data = .) }
# # maybe not, this is time consuming

# stid_biases = coef(res)[grepl('stid', names(coef(res)))]
# names(stid_biases) = levels(factor(tair_long$stid))[-1]
# stid_biases = c(stid_biases, BKBROW = 0)
# hist(stid_biases)

# stations %>%
#   transform(bias = stid_biases[match(stid, names(stid_biases))]) %>%
#   plot_station_data(aes(color = bias, shape = network)) +
#   scale_shape_manual(values = c(15, 19, 17))

# # try sea breeze coords
# tair_long[, c('d1', 'd2')] =
#   with(site_locations, conf[match(tair_long$stid, row.names(conf)), ])

# res2 = tair_long %>%
#   subset(time < as.POSIXct('2021-06-01') & !is.na(tair)) %>%
#   subset(!stid %in% nj_stids) %>%
#   transform(stid = factor(stid), time = factor(time)) %>%
#   { gam(tair ~ s(d1, d2, k = 10) + time + s(stid, bs = 're') + s(d1, d2, k = 10, by = time), data = .) }
# # maybe not, this is time consuming

# stid_biases = coef(res2)[grepl('stid', names(coef(res2)))]
# bias_stids = levels(factor(tair_long$stid))[-1]
# bias_stids = bias_stids[!bias_stids %in% c('MHLSQR', nj_stids)]
# names(stid_biases) = bias_stids
# stid_biases = c(stid_biases, BKBROW = 0)
# hist(stid_biases)

# stations %>%
#   transform(bias = stid_biases[match(stid, names(stid_biases))]) %>%
#   plot_station_data(aes(color = bias, shape = network)) +
#   scale_shape_manual(values = c(15, 19, 17))

# # just do regression with the means
# res3 = stations %>%
#   transform(stid = factor(stid),
#             d1 = with(site_locations, conf[match(stid, row.names(conf)), 1]),
#             d2 = with(site_locations, conf[match(stid, row.names(conf)), 2])) %>%
#   { gam(mean_temp ~ s(d1, d2, k = 10) + s(stid, bs = 're'), data = .) }

# res3 = stations %>%
#   transform(stid = factor(stid),
#             d1 = with(site_locations, conf[match(stid, row.names(conf)), 1]),
#             d2 = with(site_locations, conf[match(stid, row.names(conf)), 2])) %>%
#   { gam(mean_temp ~ s(d1, d2, k = 10) + stid, data = .) }

# stid_biases = coef(res3)[grepl('stid', names(coef(res3)))]
# names(stid_biases) = levels(factor(stations$stid))[-c(1, 14)]
# stid_biases = c(stid_biases, BKBROW = 0)
# hist(stid_biases)

# stations %>%
#   transform(bias = stid_biases[match(stid, names(stid_biases))]) %>%
#   plot_station_data(aes(color = bias, shape = network)) +
#   scale_shape_manual(values = c(15, 19, 17))


# # Make a temperature variogram to see if we're measuring local differences

# tair_cor # correlations
# colnames(tair_cor) = sub('.*\\.', '', colnames(tair_cor))
# row.names(tair_cor) = sub('.*\\.', '', row.names(tair_cor))

# tair_cov = cov(as.matrix(m3_wide[, -1]), use = 'pairwise.complete.obs')
# colnames(tair_cov) = sub('.*\\.', '', colnames(tair_cov))
# row.names(tair_cov) = sub('.*\\.', '', row.names(tair_cov))

# # distances

# # back to LI coords
# micro_sites = micro_sites %>%
#   st_transform(li_lcc)

# dist_mat = st_distance(micro_sites)
# colnames(dist_mat) = micro_sites$stid
# row.names(dist_mat) = micro_sites$stid

# pairwise_stats = combn(micro_sites$stid, 2) %>%
#   t %>%
#   as.data.frame
# names(pairwise_stats) = c('stid1', 'stid2')

# pairwise_stats$dist = dist_mat[with(pairwise_stats, cbind(stid1, stid2))]
# pairwise_stats$cov = tair_cov[with(pairwise_stats, cbind(stid1, stid2))]
# pairwise_stats$borough1 = substr(pairwise_stats$stid1, 1, 2)
# pairwise_stats$cor = tair_cor[with(pairwise_stats, cbind(stid1, stid2))]

# with(pairwise_stats, plot(dist, cov))

# ggplot(pairwise_stats, aes(x = as.numeric(dist), y = cov, color = borough1)) +
#   geom_point() +
#   geom_smooth() +
#   xlim(0, 1e5)

# ggplot(pairwise_stats, aes(x = as.numeric(dist), y = cov, group = stid1, color = borough1)) +
#   geom_point() +
#   geom_line() + 
#   xlim(0, 1e5)


# # # can't I model covariance as a linear regression?

# # # cov ~ var(stid1) + var(stid2) + var(dist)

# # rvar = lm(cov ~ stid1 + stid2 + dist, pairwise_stats)

# # rvar2 = lm(cov ~ stid1 + stid2 + I(dist^2), pairwise_stats)

# # library(mgcv)

# # rvar3 = gam(cov ~ stid1 + stid2 + s(as.numeric(dist)), data = pairwise_stats)

# # rvar4 = gam(cor ~ stid1 + stid2 + s(as.numeric(dist)), betar(), pairwise_stats)
