# predict load from TV values

# Current approach, used as reference: Average of LGA and NYC

# Simple alternative: stepwise forward selection of sites for averaging

# Can try a locally-weighted averaging or regression approach, like loess or
# GWR. Or fit a 2D spline.

# Or a linear combination approach. Like PLS, projection pursuit regression, or
# single index models. Options are function `stats::ppr`, packages cSEM, simml,
# PLSiMCpp. Also see `?single.index` in mgcv. `?linear.functional.terms` in mgcv
# is similar.

library(magrittr)
library(psychrolib)
library(cv)
library(ggplot2)
library(memoise)
source('R/load_data.R')
source('R/coned_tv.R')
source('R/plots.R')

cv_seed = 468

excluded_sites = c('MHLSQR', 'JRB')

networks = readRDS('results/maps/coned_networks_cleaned.rds')
stations = readRDS('results/station_data/stations.rds')
data_wide = readRDS('results/load_vs_weather/tv_and_load.rds')

# Get combined TV by averaging site temperatures and dewpoints first and then
# deriving TV. `station_obs` needs to be defined for this to work
.get_combined_tv = function(stids) {
  tmp_cols = paste0('tmpf.', stids)
  dwp_cols = paste0('dwpf.', stids)
  tmps = rowMeans(station_obs[, tmp_cols, drop = FALSE], na.rm = TRUE)
  dwps = rowMeans(station_obs[, dwp_cols, drop = FALSE], na.rm = TRUE)
  out = coned_tv(station_obs$time, tmps, dwps)
  row.names(out) = as.character(out$day)
  out
}
# This is the slowest part of the code and needs to be cached to be practical
.get_combined_tv_cached =
  memoise(.get_combined_tv,
          # 2GB cache
          cache = cachem::cache_mem(max_size = 2 * 1024^3))
get_combined_tv = function(stids) {
  # sort so that cachem recognizes the repeat inputs
  .get_combined_tv_cached(sort(stids))
}


station_ave_mse = function(stids, network, dat) {
  # tv = rowMeans(dat[, paste0('tv.', stids), drop = FALSE])
  # To get the best performance out of memoise, need to sort the stids
  # consistently
  tv_dat = get_combined_tv(stids)
  lm_dat = data.frame(load = dat[, paste0('reading.', network)],
                      tv = tv_dat[as.character(dat$day), 'tv'],
                      year = dat$year) %>%
    na.omit
  #f = load ~ factor(year):(tv + I(tv^2) + I(tv^3))
  f = load ~ factor(year) * poly(tv, 3)
  res = lm(f, lm_dat)
  mse = mean(res$residuals^2, na.rm = TRUE)
  mae = mean(abs(res$residuals), na.rm = TRUE)
  list(fit = res, mse = mse, mae = mae)
}

stepwise_forward_ave = function(id, dat) {
  current_stids = character()
  candidate_stids = setdiff(stations$stid, excluded_sites)
  n_mse = vector('numeric', length(candidate_stids))
  res_list = list()
  i = 1
  while(length(candidate_stids)) {
    res = lapply(candidate_stids, function(x) {
      station_ave_mse(c(current_stids, x), id, dat)
    })
    mses = sapply(res, getElement, name = 'mse')
    min_mse = min(mses)
    # stop if mse didn't improve
    if (i > 1 && min_mse >= n_mse[i - 1]) break
    n_mse[i] = min_mse
    min_idx = which.min(mses)
    winner = candidate_stids[min_idx]
    current_stids = c(current_stids, winner)
    candidate_stids = setdiff(candidate_stids, winner)
    res_list[[i]] = res[[min_idx]]
    i = i + 1
  }
  list(mse = n_mse, stids = current_stids, fits = res_list)
}

cv_stepwise_forward = function(id, folds = 5) {
  load_col = paste0('reading.', id)
  row_groups = sample(folds, nrow(data_wide), replace = TRUE)
  mses = numeric(folds)
  maes = numeric(folds)
  for (i in 1:folds) {
    training_data = data_wide[row_groups != i, ]
    test_data = data_wide[row_groups == i, ]
    train_res = stepwise_forward_ave(id, training_data)
    # test predictions with test data
    stids = train_res$stids
    stid_cols = paste0('tv.', stids)
    load = test_data[, load_col]
    fit_j = tail(train_res$fits, 1)[[1]]$fit
    tv_dat = get_combined_tv(stids)
    test_data$tv = tv_dat[as.character(test_data$day), 'tv']
    predicted_load = predict(fit_j, test_data)
    mses[i] = mean((load - predicted_load)^2, na.rm = TRUE)
    maes[i] = mean(abs(load - predicted_load), na.rm = TRUE)
  }
  list(mse = mean(mses), mae = mean(maes))
}

mae = function(y, yhat) mean(abs(yhat - y))



# first we need to set up the effective temp data
SetUnitSystem('SI')
tv_vars = c('tair', 'relh')
nysm_obs = get_combined_nc_data(tv_vars, hour_to_nc_index(7:21)) %>%
  subset(!(is.na(tair) | is.na(relh))) %>%
  transform(dwp = GetTDewPointFromRelHum(tair, relh / 100)) %>%
  transform(tmpf = as_fahrenheit(tair), dwpf = as_fahrenheit(dwp)) %>%
  subset(select = c(time, stid, tmpf, dwpf)) %>%
  reshape(direction = 'wide', idvar = 'time', timevar = 'stid')

asos_obs = get_hourly_asos_data(7:21) %>%
  transform(stid = station, time = valid_hour) %>%
  subset(!(is.na(tmpf) | is.na(dwpf))) %>%
  subset(select = c(time, stid, tmpf, dwpf)) %>%
  reshape(direction = 'wide', idvar = 'time', timevar = 'stid')

station_obs = merge(nysm_obs, asos_obs, all = T) %>%
  subset(as.Date(time, tz = 'EST5EDT') %in% data_wide$day)

# saveRDS(et_wide, 'results/select_stations/et.rds')


### System level analysis ###

# get combined peak loads, ignoring outliers for now
peaks = read.csv('data/coned/Borough and System Data 2020-2024.csv') %>%
  transform(DT = as.POSIXct(DT, tz = 'EST5EDT', '%m/%d/%Y %H:%M'),
            # some inconsistency, but various forms of 'False' mean data is fine
            BAD = !startsWith(BAD, 'F')) %>%
  subset(Borough == 'CECONY') %>%
  subset(as.POSIXlt(DT)$mon >= 4 & as.POSIXlt(DT)$mon < 9) %>%
  subset(as.POSIXlt(DT)$year < 2024) %>% # since we don't have all of 2024 yet
  subset(!BAD & !is.na(DT)) %>%
  subset(as.POSIXlt(DT)$hour %in% 9:21) %>%
  # daily peak loads
  transform(day = as.Date(DT, tz = 'EST5EDT')) %>%
  aggregate(READING ~ day, ., max)
names(peaks) = tolower(names(peaks))

cv_mse = peaks %>%
  transform(year = trunc(day, 'years'),
            tv = get_combined_tv(c('LGA', 'NYC'))[as.character(day), 'tv']) %>%
  na.omit %>%
  lm(reading ~ factor(year) * poly(tv, 3), .) %>%
  cv(k = 10, seed = cv_seed)
cv_mse$`CV crit` # 40737.07

cv_mae = peaks %>%
  transform(year = trunc(day, 'years'),
            tv = get_combined_tv(c('LGA', 'NYC'))[as.character(day), 'tv']) %>%
  na.omit %>%
  lm(reading ~ factor(year) * poly(tv, 3), .) %>%
  cv(criterion = mae, k = 10, seed = cv_seed)
cv_mae$`CV crit` # 162.6026

data_wide$reading.total = peaks$reading[match(data_wide$day, peaks$day)]

set.seed(cv_seed)
combined_stfor_cv = cv_stepwise_forward('total', 10)
# $mse
# [1] 36363.01
# $mae
# [1] 144.297
# good improvement

1 - (combined_stfor_cv$mae / cv_mae$`CV crit`) # 11% improvement

combined_stfor = stepwise_forward_ave('total', data_wide)
combined_stfor$stids
# [1] "BKNYRD" "BKMAPL" "SIFKIL" "QNSOZO" "BRON"  

# save this
system_results = list(
    errors = rbind(c(cv_mse$`CV crit`, cv_mae$`CV crit`),
                   unlist(combined_stfor_cv)),
    stids = combined_stfor$stids
)
row.names(system_results$errors) = c('NYC+LGA', 'Stepwise forward selection')
saveRDS(system_results, 'results/select_stations/system_results.rds')




### Network level analysis ###

# cv = cv::cv # this is needed if raster package is loaded for plots

current_mse = sapply(networks$id, function(x) {
  tv_dat = get_combined_tv(c('LGA', 'NYC'))
  tv = tv_dat[as.character(data_wide$day), 'tv']
  lm_dat = data.frame(load = data_wide[, paste0('reading.', x)], tv = tv,
                      year = data_wide$year) %>%
    na.omit
  f = load ~ factor(year) * poly(tv, 3)
  cv(lm(f, lm_dat), k = 10, seed = cv_seed)$`CV crit`
})

current_mae = sapply(networks$id, function(x) {
  tv_dat = get_combined_tv(c('LGA', 'NYC'))
  tv = tv_dat[as.character(data_wide$day), 'tv']
  lm_dat = data.frame(load = data_wide[, paste0('reading.', x)], tv = tv,
                      year = data_wide$year) %>%
    na.omit
  f = load ~ factor(year) * poly(tv, 3)
  cv(lm(f, lm_dat), criterion = mae, k = 10, seed = cv_seed)$`CV crit`
})


# sort(current_mse, decreasing = T) # seems to be related to network size
# networks %>%
#   transform(current_mse = current_mse[id]) %>%
#   plot_network_data(aes(fill = current_mse))
# with(data_wide, plot(day, reading.6B, type = 'b'))


# # first of all let's do some profiling

# # without memoise
# system.time(cv_stepwise_forward(networks$id[1], 10))
# #   user  system elapsed 
# # 27.311   0.054  27.389 

# # with new instance of memoized function
# system.time(cv_stepwise_forward(networks$id[1], 10))
# #   user  system elapsed 
# # 15.333   0.019  15.367
# # wow big help already

# factorial(28)
# [1] 3.048883e+29
# ok so it's simply not feasible to store all of these results


# this takes ~10 minutes
message('Running network-level stepwise forward selection (takes around 10 minutes)')
set.seed(cv_seed)
system.time({
  cv_stforw_list <- lapply(networks$id, function(id) {
    n_iter = which(networks$id == id)
    message(n_iter, '/', nrow(networks), ' - ', id)
    cv_stepwise_forward(id, 10)
  })
})
#   user  system elapsed 
# 642.63    1.50  647.86 
saveRDS(cv_stforw_list, 'results/select_stations/cv_stforw.rds')

# cv.select <- cvSelect(selectStepAIC, data=D, seed=3791,
#                       model=m.null, direction="forward",
#                       scope=list(lower=~1, 
#                                  upper=formula(m.full)))

# forw_list = lapply(networks$id, stepwise_forward_ave, dat = data_wide)
# mse_mat = t(sapply(forw_list, getElement, name = 'mse'))
# row.names(mse_mat) = networks$id
# hist(apply(mse_mat, 1, which.min)) # nice
# plot(n_mse)



# # compare to forward selection
# cbind(current_mse, new_mse = apply(mse_mat, 1, min))

# plot(current_mse, apply(mse_mat, 1, min), xlab = 'JFK+LGA MSE',
#      ylab = 'Stepwise forward min MSE', asp = 1)
# abline(0, 1)

# plot(current_mse, mse_mat[, 4], xlab = 'JFK+LGA MSE',
#      ylab = 'Stepwise forward 4 station MSE', asp = 1)
# abline(0, 1)

# plot(current_mse, mse_mat[, 2], xlab = 'JFK+LGA MSE',
#      ylab = 'Stepwise forward 2 station MSE', asp = 1)
# abline(0, 1)


# compare to forward selection (CV)

stfor_best_mses = sapply(cv_stforw_list, getElement, name = 'mse')
plot(current_mse, stfor_best_mses, xlab = 'LGA+NYC MSE',
     ylab = 'Stepwise forward MSE', asp = 1)
abline(0, 1)

# cv_mae_mat = t(sapply(cv_stforw_list, getElement, name = 'mae'))
stfor_best_maes = sapply(cv_stforw_list, getElement, name = 'mae')
mae_comparison = data.frame(orig = current_mae,
                            stfor = stfor_best_maes)
saveRDS(mae_comparison, 'results/select_stations/mae_comparison.rds')
# mae_comparison = readRDS('results/mae_comparison.rds')
png('results/select_stations/mae_comparison.png', res = 100)
with(mae_comparison, plot(orig, stfor, xlab = 'LGA+NYC MAE (MW)',
                          ylab = 'Stepwise forward MAE (MW)', asp = 1,
                          xlim = c(0, max(orig)), ylim = c(0, max(stfor))))
grid()
abline(0, 1)
mae_lm = lm(stfor ~ orig, mae_comparison)
abline(mae_lm, lty = 'dashed')
dev.off()

with(mae_comparison, sum(stfor) / sum(orig))
# ~8% reduction in MAE

# is the improvement related to location?
png('results/select_stations/mae_improvement.png', width = 1000, height = 800, res = 100)
networks %>%
  transform(mae_pdiff = with(mae_comparison, 100 * (1 - (stfor / orig)))) %>%
  plot_network_data(aes(fill = mae_pdiff)) +
  scale_fill_viridis_c('MAE % reduction')
dev.off()


# compare % change to unfolded coordinates
unf = readRDS('results/load_vs_weather/tv_load_unfolding.rds')
pc_unf = prcomp(unf$conf.row)
networks %>%
  transform(sea_breeze = pc_unf$x[, 1],
            mae_pdiff = with(mae_comparison, 1 - (stfor / orig))) %>%
  ggplot(aes(sea_breeze, mae_pdiff)) +
  geom_point()

with(mae_comparison, t.test(log(orig), log(stfor)))

lm_identity = lm(stfor ~ -1 + offset(orig), mae_comparison)

anova(lm_identity, mae_lm)

lm_test = lm(stfor ~ orig + offset(orig), mae_comparison)
# slope is significantly less than 1

mae_p_values = list(
    p_45deg = anova(lm_identity, mae_lm)[2, 6],
    p_slope1 = summary(lm_test)$coefficients['orig', 4]
)

saveRDS(mae_p_values, 'results/select_stations/mae_p_values.rds')
     




# get the selected sites

names(cv_stforw_list) = networks$id
stforw_sites_list = lapply(networks$id, function(id) {
  # get the number of sites to include from CV results
  #n_sites = which.min(cv_stforw_list[[id]]$mse)
  res = stepwise_forward_ave(id, data_wide)
  #head(res$stids, n_sites)
  res$stids
})
names(stforw_sites_list) = networks$id
saveRDS(stforw_sites_list, 'results/select_stations/selected_sites.rds')

table(lengths(stforw_sites_list))
 # 1  2  3  4  5  6  8  9 10 11 13 
 # 4  5 20 10 11 12  2  3  1  1  1 

# need to unlist this somehow
stforw_sites = unlist(stforw_sites_list) %>%
  { data.frame(network = names(.), stid = .) } %>%
  transform(network = sub('[0-9]*$', '', network))


tv_mds = readRDS('results/weather/tv_mds.rds')
pc_tv = prcomp(tv_mds$conf)
# tv_load_unf = readRDS('results/load_vs_weather/tv_load_unfolding.rds')
# # get sea breeze dimension
# pc1_points = with(tv_load_unf, {
#   rbind(conf.row[row.names(conf.row) != '14M', ],
#         conf.col[stations[row.names(conf.col), ]$network != 'ASOS', ])
# })
# unf_pc1 = prcomp(pc1_points)
# network_sb = predict(unf_pc1, tv_load_unf$conf.row)
# network_sb_order = row.names(network_sb)[order(network_sb[, 1])]
tv_load_pc = readRDS('results/load_vs_weather/tv_load_pca.rds')
network_sb_order = with(tv_load_pc, row.names(x)[order(x[, 2])])
station_sb_order = with(pc_tv, row.names(x)[order(x[, 1])])
station_sb_order = setdiff(station_sb_order, 'JRB')

# need to order by sea breeze dimension
png('results/select_stations/selected_stations.png', width = 1500, height = 800, res = 100)
stforw_sites %>%
  transform(network = factor(network, network_sb_order),
            station = factor(stid, station_sb_order),
            weather_network = stations[stid, ]$network) %>%
  ggplot(aes(network, station, fill = weather_network)) +
  geom_tile() +
  ggtitle('Stepwise forward selected stations') +
  xlab('← More sea breeze') +
  ylab('← More sea breeze') + 
  scale_y_discrete(drop = FALSE) +
  scale_fill_discrete(type = hcl.colors(3, 'Dark 3')) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

png('results/select_stations/selected_stations_barplot.png', width = 800, height = 900, res = 100)
stforw_sites %>%
  aggregate(network ~ stid, ., length) %>%
  transform(count = network) %>%
  subset(select = -network) %>%
  merge(as.data.frame(stations), by = 'stid', all.y = T) %>%
  transform(count = replace(count, is.na(count), 0)) %>%
  subset(!stid %in% excluded_sites) %>%
  ggplot(aes(count, reorder(stid, count), fill = network)) +
  ggtitle('Stepwise forward selected stations') +
  geom_bar(stat = 'identity', orientation = 'y') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  ylab('Station') +
  # ggplot(aes(reorder(stid, -count), count, fill = network)) +
  # geom_bar(stat = 'identity') +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  # xlab('Station') +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_discrete(type = hcl.colors(3, 'Dark 3'))
dev.off()

# # same but with only 5 sites per network
# stforw_sites_list5 = lapply(networks$id, function(id) {
#   res = stepwise_forward_ave(id, data_wide)
#   head(res$stids, 5)
# })
# names(stforw_sites_list5) = networks$id

# stforw_sites5 = unlist(stforw_sites_list5) %>%
#   { data.frame(network = names(.), stid = .) } %>%
#   transform(network = sub('[0-9]*$', '', network))

# stforw_sites5 %>%
#   transform(network = factor(network, network_sb_order),
#             station = factor(stid, station_sb_order),
#             weather_network = stations[stid, ]$network) %>%
#   ggplot(aes(network, station, fill = weather_network)) +
#   geom_tile() +
#   xlab('← More sea breeze / less sea breeze →') +
#   ylab('← More sea breeze / less sea breeze →') + 
#   scale_y_discrete(drop = FALSE) +
#   scale_fill_discrete(type = hcl.colors(3, 'Dark 3')) +
#   coord_fixed()
#   #scale_fill_brewer(type = 'qual', palette = 'Dark2')

# stforw_sites5 %>%
#   transform(network_coord = -network_sb[network, 1],
#             station_coord = -station_sb[stid, 1],
#             weather_network = stations[stid, ]$network) %>%
#   ggplot(aes(network_coord, station_coord, color = weather_network)) +
#   geom_point() # super ugly

# stforw_sites5 %>%
#   transform(weather_network = stations[stid, ]$network) %>%
#   ggplot(aes(stid, fill = weather_network)) +
#   geom_bar() +
#   scale_fill_discrete(type = hcl.colors(3, 'Dark 3'))

# stforw_sites5 %>%
#   aggregate(network ~ stid, ., length) %>%
#   transform(count = network) %>%
#   subset(select = -network) %>%
#   merge(as.data.frame(stations), by = 'stid', all.y = T) %>%
#   transform(count = replace(count, is.na(count), 0)) %>%
#   ggplot(aes(count, reorder(stid, count), fill = network)) +
#   geom_bar(stat = 'identity', orientation = 'y') +
#   scale_fill_discrete(type = hcl.colors(3, 'Dark 3'))




##### Not sure about this part

# library(simml)

# # sim_data = subset(peaks_wide, select = c(day, reading.10B)) %>%
# #   merge(tv_wide) %>%
# #   subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2023)),
# #          select = -tv.KJRB) %>%
# #   # na values break ppr
# #   na.omit
# sim_data = subset(data_wide, select = -tv.JRB) %>%
#   # na values break ppr
#   na.omit
# #y = sim_data$reading.10B
# x = as.matrix(sim_data[, grepl('tv', names(sim_data))])
# colnames(x) = sub('tv\\.', '', colnames(x))

# sim1 = simml(sim_data$reading.10B, sim_data$year, x, beta.ini = rep(.5, ncol(x)))

# # visualization of the estimated link functions of SIMML
# sim1$beta.coef # estimated single-index coefficients
# g.fit <- simml.obj1$g.fit # estimated trt-specific link functions; "g.fit" is a mgcv::gam object.
# #plot(g.fit)
# plot(sim1$g.fit, pages = 1, scale = 0) # eh

# plot(sim1$g.fit, xlim = c(-.1, .1)) # hmm
# # this must be terribly overfit somehow

# sim_data = subset(data_wide, select = -tv.JRB) %>%
#   # na values break ppr
#   na.omit

# x2 = as.matrix(sim_data[, paste0('tv.', c('BKBROW', 'BKMAPL', 'BKLN', 'QNHBEA'))])
# colnames(x2) = sub('tv\\.', '', colnames(x2))
# sim2 = simml(sim_data$reading.10B, sim_data$year, x2, beta.ini = rep(1, ncol(x2)))

# plot(sim2$g.fit, pages = 1, scale = 0) # what am I looking at?
# sim2$beta.coef

# x3 = as.matrix(sim_data[, paste0('tv.', c('BKMAPL', 'BKLN'))])
# colnames(x3) = sub('tv\\.', '', colnames(x3))
# sim3 = simml(sim_data$reading.10B, factor(sim_data$year), x3,
#              beta.ini = rep(1, ncol(x3)), lambda = 1, k = 4, scale.X = F,
#              center.X = F)

# plot(sim3$g.fit, pages = 1, scale = 0) # what am I looking at?
# sim3$beta.coef

# sim4 = simml(sim_data$reading.10B, factor(sim_data$year), x3, beta.ini = rep(1, ncol(x3)), k = 4)

# plot(sim4$g.fit, pages = 1, scale = 0) # what am I looking at?
# sim4$beta.coef


# with(sim_data, plot(tv.BKBROW, reading.1B, col = factor(year)))

# data_wide$year = factor(trunc(data_wide$day, 'years'))

# with(sim_data, plot(tv.BKBROW, reading.10B, col = year))

# with(sim_data, plot(tv.BKMAPL, reading.10B, col = year))
# with(sim_data, plot(tv.BKLN, reading.10B, col = year))


# # # ok now do ppr

# # ppr_data = subset(peaks_wide, select = c(day, reading.1B)) %>%
# #   merge(tv_wide) %>%
# #   subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2023)),
# #          select = -tv.KJRB) %>%
# #   # na values break ppr
# #   na.omit
# # y = ppr_data[, 2]
# # x = as.matrix(ppr_data[, 3:29])
# # colnames(x) = sub('tv\\.', '', colnames(x))

# ppr1 = ppr(x, sim_data$reading.10B, nterms = 1)
# plot(ppr1)
# ppr1$alpha[order(abs(ppr1$alpha), decreasing = T)]
# # this seems great. Why is simml so bad in comparison?

# # ppr_data %>%
# #   transform(year = trunc(day, 'years')) %>%
# #   with(plot(tv.BKNYRD, reading.1B, col = factor(year)))

# # stations = readRDS('results/stations.rds')

# # stations$ppr_alpha = ppr1$alpha[match(stations$stid, names(ppr1$alpha))]

# # plot_station_data(stations, aes(color = abs(ppr_alpha)))
# # plot_station_data(stations, aes(color = log(abs(ppr_alpha))))

# # # let's try 5Q, around JFK

# # ppr_data = subset(peaks_wide, select = c(day, reading.1X)) %>%
# #   merge(tv_wide) %>%
# #   subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2023)),
# #          select = -tv.KJRB) %>%
# #   # na values break ppr
# #   na.omit
# # y = ppr_data[, 2]
# # x = as.matrix(ppr_data[, 3:29])
# # colnames(x) = sub('tv\\.', '', colnames(x))
# # ppr2 = ppr(x, y, nterms = 1)
# # plot(ppr2)
# # ppr2$alpha[order(abs(ppr2$alpha), decreasing = T)]

# # stations$ppr2_alpha = ppr2$alpha[match(stations$stid, names(ppr2$alpha))]
# # plot_station_data(stations, aes(color = abs(ppr2_alpha)))

# # plot_station_data(stations, aes(color = log(abs(ppr2_alpha))))


# library(PLSiMCpp)

# simml(sim_data$reading.10B, factor(sim_data$year), x3,
#       beta.ini = rep(1, ncol(x3)), lambda = 1, k = 4, scale.X = F,
#       center.X = F)

# z = x3
# xdat = model.matrix(~ -1 + factor(year), data = sim_data)[, -1]

# p1 = plsim.est(xdat, x3, as.matrix(sim_data$reading.10B))
# with(p1, plot(Z_alpha, eta))

# p2 = plsim.est(xdat, x2, as.matrix(sim_data$reading.10B))
# with(p2, plot(Z_alpha, eta))
# # definitely better than simml at least

# with(p2, plot(Z_alpha, y_hat))

# p2$zeta

# with(p2, plot(sim_data$reading.10B, y_hat))
# abline(0, 1)

# with(sim_data, plot(tv.BKLN, reading.10B, col = year))










# # we can also try an mgcv version

# library(mgcv)

# # linear functional-- nah, doesn't work

# gam(sim_data$reading.10B ~ )
