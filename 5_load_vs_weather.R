# examining the relationship between load and weather, especially sea breeze
# patterns

library(timeDate) # holidays
library(MASS) # robust regression
library(magrittr)
library(sf)
library(smacof) # this overwrites `transform`, ugh
transform = base::transform
source('R/load_data.R')
source('R/plots.R')

nj_sites = c('TEB', 'EWR', 'LDJ')

# detect outliers using Tukey's fences
detect_outliers = function(x, k = 3) {
  q13 = quantile(x, c(.25, .75), na.rm = TRUE)
  tf_low = q13[1] - k * diff(q13)
  tf_high = q13[2] + k * diff(q13)
  # (x < tf_low) | (x > tf_high)
  list(low = x < tf_low, high = x > tf_high)
}

get_true_clusters = function(x, n = 5) {
  out = rep(FALSE, length(x))
  cluster_info0 = data.frame(cluster_id = integer(), idx.start = integer(),
                             idx.end = integer())
  empty_res = list(outlier = out, clusters = cluster_info0)
  if (!any(na.omit(x))) return(empty_res)
  clusters = data.frame(idx = 1:length(x), x = x) %>%
    subset(!is.na(x)) %>%
    transform(cluster_id = c(0, cumsum(diff(x) != 0))) %>%
    subset(x) %>%
    transform(cluster_len = ave(x, cluster_id, FUN = length)) %>%
    subset(cluster_len >= n)
  if (!nrow(clusters)) return(empty_res)
  out[clusters$idx] = TRUE
  cluster_starts = aggregate(idx ~ cluster_id, clusters, head, n = 1)
  cluster_ends = aggregate(idx ~ cluster_id, clusters, tail, n = 1)
  cluster_info = merge(cluster_starts, cluster_ends, by = 'cluster_id',
                       suffixes = c('.start', '.end'))
  list(outlier = out, clusters = cluster_info)
}

detect_clustered_outliers = function(x, k = 1.5, n = 5) {
  tf_out = detect_outliers(x, k)
  low_cluster = get_true_clusters(tf_out$low, n)
  high_cluster = get_true_clusters(tf_out$high, n)
  if (nrow(low_cluster$clusters)) {
    low_cluster$clusters$type = 'low'
  } else {
    low_cluster$clusters$type = character()
  }
  if (nrow(high_cluster$clusters)) {
    high_cluster$clusters$type = 'high'
  } else {
    high_cluster$clusters$type = character()
  }
  cluster_info = rbind(low_cluster$clusters, high_cluster$clusters)
  list(high = high_cluster$outlier, low = low_cluster$outlier,
       clusters = cluster_info)
}

get_robust_residuals = function(tv, load, year) {
  # f = load ~ factor(year) + tv + I(tv^2) + I(tv^3)
  f = load ~ factor(year) * poly(tv, 3)
  # `na.exclude` is required here to make the indexing consistent!
  res = rlm(f, na.action = na.exclude)
  resid(res)
}

get_outlier_scores = function(load, k, cluster = FALSE, cluster_n = 5) {
  year = reg_data$year
  if (cluster) {
    fun = function(x) {
      detect_clustered_outliers(x, k, cluster_n)
    }
  } else {
    fun = function(x) {
      detect_outliers(x, k)
    }
  }
  out_list = lapply(test_sites, function(stid) {
    tv = reg_data[, paste0('tv.', stid)]
    resid = get_robust_residuals(tv, load, year)
    fun(resid)
  })
  high_scores = Reduce(`+`, lapply(out_list, getElement, name = 'high'))
  low_scores = Reduce(`+`, lapply(out_list, getElement, name = 'low'))
  list(high = high_scores, low = low_scores)
}

# get matrix of confirmed outliers, and list of outlier clusters
get_network_outliers = function(k = 3, cluster_k = 1.5, cluster_n = 5) {
  cluster_list = list()
  outliers = matrix(F, nrow = nrow(reg_data), ncol = nrow(networks),
                    dimnames = list(as.character(reg_data$day), networks$id))
  for (i in 1:nrow(networks)) {
    id = networks$id[i]
    #message(paste0('reading.', id))
    load = reg_data[, paste0('reading.', id)]
    # first the clusters
    cluster_scores = get_outlier_scores(load, cluster_k, TRUE, cluster_n)
    high_clusters = get_true_clusters(cluster_scores$high == 3, cluster_n)
    low_clusters = get_true_clusters(cluster_scores$low == 3, cluster_n)
    # save cluster info
    if (nrow(low_clusters$clusters)) {
      low_clusters$clusters$type = 'low'
    } else {
      low_clusters$clusters$type = character()
    }
    if (nrow(high_clusters$clusters)) {
      high_clusters$clusters$type = 'high'
    } else {
      high_clusters$clusters$type = character()
    }
    cluster_list[[i]] = rbind(high_clusters$clusters, low_clusters$clusters)
    cluster_list[[i]]$cluster_id = NULL
    if (nrow(cluster_list[[i]])) {
      cluster_list[[i]]$network = id
    } else {
      cluster_list[[i]]$network = character()
    }
    outlier_cluster = cluster_scores$high == 3 | cluster_scores$low == 3
    # now the isolated outliers (removing clustered outliers)
    load[outlier_cluster] = NA
    regular_scores = get_outlier_scores(load, k, cluster = FALSE)
    regular_outlier = regular_scores$high == 3 | regular_scores$low == 3
    outliers[, id] = outlier_cluster | regular_outlier
  }
  outlier_clusters = do.call(rbind, cluster_list) %>%
    transform(start = reg_data$day[idx.start], end = reg_data$day[idx.end],
              days = idx.end - idx.start + 1) %>%
    subset(select = -c(idx.start, idx.end))
  outliers[is.na(outliers)] = FALSE
  list(outliers = outliers, clusters = outlier_clusters)
}

# plot outliers by network, subsetting by both year and weather station
plot_outliers = function(id, outliers) {
  par(mfrow = c(3, 3))
  id_col = paste0('reading.', id)
  stid_cols = paste0('tv.', test_sites)
  x_range = range(as.matrix(reg_data[!outliers[, id], stid_cols]), na.rm = T)
  for (year in 2021:2023) {
    year_idx = reg_data$year == year
    y_lab = paste0(id, ' load (', year, ')')
    for (stid in test_sites) {
      stid_col = paste0('tv.', stid)
      x_lab = paste(stid, 'TV')
      is_outlier = outliers[year_idx, id]
      plot(reg_data[year_idx, stid_col], reg_data[year_idx, id_col],
           col = factor(is_outlier), xlim = x_range, xlab = x_lab,
           ylab = y_lab)
    }
  }
  par(mfrow = c(1, 1))
}

# plot outliers for each network, ordered by number of outliers
plot_outliers_sequence = function(x) {
  n_outliers = colSums(x$outliers, na.rm = TRUE)
  plot_order = order(n_outliers, decreasing = TRUE)
  ids = colnames(x$outliers)
  oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
  for (i in plot_order) plot_outliers(ids[i], x$outliers)
}

get_tv_load_r2_yearly2 = function(tv, load, year) {
  f = load ~ factor(year) * (tv + I(tv^2) + I(tv^3))
  res = lm(f)
  summary(res)$r.squared
}

get_tv_network_distances2 = function(reg_data,
                                     reg_fun = get_tv_load_r2_yearly2) {
  year = trunc(reg_data$day, 'years')
  tv_cols = colnames(reg_data)[startsWith(colnames(reg_data), 'tv')]
  stids = sub('tv\\.', '', tv_cols)
  load_cols = colnames(reg_data)[startsWith(colnames(reg_data), 'reading')]
  net_ids = sub('reading\\.', '', load_cols)
  r2_mat = matrix(NA, length(net_ids), length(stids),
                  dimnames = list(net_ids, stids))
  for (x in stids) {
    tv_col = paste0('tv.', x)
    for (y in net_ids) {
      reading_col = paste0('reading.', y)
      load = reg_data[, reading_col]
      if (any(o4$outliers[, y])) load[o4$outliers[, y]] = NA
      r2_mat[y, x] = reg_fun(reg_data[, tv_col], load, year)
    }
  }
  sim2diss(r2_mat)
}

get_unfolding = function(dist, fix_sites = NULL, init_sites = NULL, ...) {
  if (is.null(fix_sites) && is.null(init_sites)) {
    unfolding(dist, type = 'interval', conditionality = 'row', ...)
  } else if (!is.null(fix_sites)) {
    unfolding(dist, type = 'interval', conditionality = 'row', fixed = 'column',
              fixed.coord = fix_sites, ...)
  } else if (!is.null(init_sites)) {
    res0 = unfolding(dist, type = 'interval', conditionality = 'row',
                     fixed = 'column', fixed.coord = fix_sites)
    unfolding(dist, type = 'interval', conditionality = 'row',
              init = res0[c('conf.row', 'conf.col')], ...)
  }
}

extract_sea_breeze_coords = function(unf_results, site_reference = NULL,
                                     use_networks = FALSE) {
  # remove irrelevant NJ sites
  nj_row = row.names(unf_results$conf.col) %in% nj_sites
  coords0 = with(unf_results, conf.col[!nj_row, ])
  if (use_networks) {
    pc_networks = prcomp(unf_results$conf.row)
    station_coords = predict(pc_networks, unf_results$conf.col)
    network_coords = pc_networks$x
    rotation = pc_networks$rotation
  } else if (is.null(site_reference)) {
    pc_stations = prcomp(coords0)
    station_coords = predict(pc_stations, unf_results$conf.col)
    network_coords = predict(pc_stations, unf_results$conf.row)
    rotation = pc_stations$rotation
  } else {
    # match to reference locations using procrustes
    procr_stations = Procrustes(site_locations[row.names(coords0), ], coords0)
    station_coords0 = procr_stations$Yhat
    network_coords0 =
      with(procr_stations, dilation * unf_results$conf.row %*% rotation +
                           rep(1, unf_results$nind) %*% t(translation))
    # use reference sea breeze dimension
    ref_pc = prcomp(site_locations)
    station_coords = predict(ref_pc, station_coords0)
    network_coords = predict(ref_pc, network_coords0)
    rotation = procr_stations$rotation
  }
  list(stations = station_coords, networks = network_coords,
       rotation = rotation)
}

plot_sea_breeze_dim = function(st_coords, net_coords, sbdim = 1,
                               include_nj = FALSE, max_sb = NULL) {
  networks2 = networks %>%
    transform(type = 'Networks', sea_breeze = net_coords[id, sbdim])
  stations2 = stations %>%
    transform(type = 'Stations',
              sea_breeze = st_coords[match(stid, row.names(st_coords)), sbdim])
  if (include_nj) {
    sb_range = c(stations2$sea_breeze, networks2$sea_breeze) %>%
      range(na.rm = T)
  } else {
    # exclude NJ sites from plotted sb range
    sb_range = stations2$sea_breeze[!stations2$stid %in% nj_sites] %>%
      c(networks2$sea_breeze) %>%
      range(na.rm = T)
  }
  if (!is.null(max_sb)) sb_range[2] = min(sb_range[2], max_sb)
  # x_range = range(c(networks$sea_breeze, stations$sea_breeze), na.rm = T)
  # sb_range[2] = min(c(sb_range[2], max_sb))
  plot_station_data(networks2, aes(fill = sea_breeze), nyc_base, alpha = .9) +
    geom_sf(aes(color = sea_breeze, shape = network), stations2, size = 6) +
    scale_color_viridis_c('Sea breeze', limits = sb_range) +
    scale_fill_viridis_c('Sea breeze', limits = sb_range) +
    scale_shape_manual('Weather network', values = c(15, 19, 17)) +
    facet_wrap(~ factor(type, levels=c('Stations', 'Networks')))
}

get_factor_coefs = function(fit, name) {
  name_regex = paste0('^', name)
  coefs = coef(fit)[grepl(name_regex, names(coef(fit)))]
  names(coefs) = sub(name_regex, '', names(coefs))
  if (inherits(fit, 'lm')) {
    ref_level = fit$xlevels[[name]][1]
  } else if (inherits(fit, 'betareg')) {
    ref_level = fit$levels$mean[[name]][1]
  } else {
    stop('Model type "', class(fit), '" not implemented')
  }
  c(setNames(0, ref_level), coefs)
}


networks = readRDS('results/maps/coned_networks_cleaned.rds')
stations = readRDS('results/station_data/stations.rds')

# these sites are used for outlier detection
test_sites = c('SIFKIL', 'JFK', 'BRON')

peaks = list.files('data/coned', 'Network.*\\.csv', full.names = T) %>%
  lapply(read.csv) %>%
  do.call(rbind, .) %>%
  transform(DT = as.POSIXct(DT, tz = 'EST5EDT', '%m/%d/%Y %H:%M:%S'),
            # some inconsistency, but various forms of 'False' mean data is fine
            BAD = !startsWith(BAD, 'F')) %>%
  subset(!BAD & !is.na(DT)) %>%
  subset(as.POSIXlt(DT)$hour %in% 9:21) %>%
  # get daily load peaks
  transform(day = as.Date(DT, tz = 'EST5EDT')) %>%
  aggregate(READING ~ Network + day, ., max)
names(peaks) = tolower(names(peaks))
peaks_wide = reshape(peaks, direction = 'wide', idvar = 'day',
                     timevar = 'network')
# B=Brooklyn, M=Manhattan, Q=Queens, R=Staten Island (Richmond County), X=Bronx

# OK, compare with loads
tv_wide = readRDS('results/weather/tv.rds')
data_wide = merge(tv_wide, peaks_wide)
stids = sub('tv\\.', '', colnames(tv_wide)[-1])
network_ids = unique(peaks$network)
bus_days = isBizday(as.timeDate(data_wide$day), holidays = holidayNYSE(2021:2023))

# this is the data that will actually be used for regressions-- remove outliers
# first!!
reg_data = data_wide[bus_days, ]
reg_data$year = as.integer(format(trunc(reg_data$day, 'years'), '%Y'))
saveRDS(reg_data, 'results/load_vs_weather/tv_and_load_uncleaned.rds')
# reg_data = readRDS('results/load_vs_weather/tv_and_load_uncleaned.rds')

# but first, need to carefully identify outliers
o1 = get_network_outliers(k = 3, cluster_k = 1.5, cluster_n = 5)
o1$clusters

# now we want to check out the results to verify it looks right
colSums(o1$outliers, na.rm = T)

plot_outliers_sequence(o1)
# potential problems:

# 2R 2021, ~3
# 32M 2022, 8
# 1B 2023, 4
# 2X 2023, 5
# 24M 2022, 4
# 14M 2021, everything

plot_outliers('14M', o1$outliers) # huge mess
plot_outliers('32M', o1$outliers) # maybe missed some

subset(reg_data, year == 2022) %>%
  with(plot(tv.BRON, reading.32M, type = 'b'))
# hmm. Let's try to fix that

o2 = get_network_outliers(k = 2, cluster_k = 1.5, cluster_n = 5)
plot_outliers('32M', o2$outliers) # better

# could just manually fix 32M if it's the only problem

# plot_outliers('2R', o3$outliers)

o3 = get_network_outliers(k = 2, cluster_k = 1.5, cluster_n = 3)
plot_outliers('32M', o3$outliers) # better

plot_outliers('39M', o1$outliers)
plot_outliers('39M', o2$outliers)

# # I think there's a good case for using lower k, maybe 2.5

# o2 = get_network_outliers(k = 2.5, cluster_k = 1.5, cluster_n = 5)
# table(o1$outliers) # 135
# table(o2$outliers) # 157
# plot_outliers_sequence(o2)
# # o2 is better

# plot_outliers('24M', o2$outliers)

# subset(reg_data, year == 2022) %>%
#   with(plot(tv.BRON, reading.32M, type = 'b'))

# subset(reg_data, year == 2021) %>%
#   with(plot(tv.LGA, reading.14M, type = 'b'))
# # I don't see any rhyme or reason to this

# subset(reg_data, year == 2021) %>%
#   with(plot(tv.SIFKIL, reading.2R, type = 'b'))
# # final days of Sept. meh

# subset(reg_data, year == 2023) %>%
#   with(plot(tv.BKNYRD, reading.1B, type = 'b'))
# # also end of Sept.

# subset(reg_data, year == 2023) %>%
#   with(plot(tv.SIFKIL, reading.2X, type = 'b'))

# subset(reg_data, year == 2023) %>%
#   with(plot(day, reading.2X, type = 'b'))
# # also end of Sept. Should be flagged

# subset(reg_data, year == 2022) %>%
#   with(plot(tv.QNDKIL, reading.24M, type = 'b'))

# subset(reg_data, year == 2022) %>%
#   with(plot(day, reading.24M, type = 'b'))
# # early july. Should be flagged


# old list:
# 2R 2021, ~3
# 32M 2022, 8
# 1B 2023, 4
# 2X 2023, 5
# 24M 2022, 4
# 14M 2021, everything

# new:
# 6B, huh
# we got 2X, sweet
# we got 24M, sweet
# 2M, not sure
# 8B is questionable
# 1R is questionable

# I think the answer is lowering cluster k rather than allowing disagreement
o4 = get_network_outliers(k = 2.5, cluster_k = 1, cluster_n = 5)
table(o1$outliers) # 135
table(o2$outliers) # 157
table(o3$outliers) # 211 -- wow
table(o4$outliers) # 178 -- 181 after fix -- 242 after peak fix, 251
o4$clusters[order(o4$clusters$days, decreasing = T), ]
plot_outliers_sequence(o4)

# new:
# 7B, huh
# we got 2X, sweet
# we got 24M, sweet
# 6B, seems fine

# still need:
# 2R 2021, ~3 -- not too bad
# 32M 2022, 8
# 1B 2023, 4
# 14M 2021, everything

plot_outliers('1B', o4$outliers)
plot_outliers('32M', o4$outliers)

subset(reg_data, year == 2022) %>%
  with(plot(tv.MHALPH, reading.32M, type = 'b'))
subset(reg_data, year == 2022 & day <= '2022-05-26') %>%
  with(points(tv.MHALPH, reading.32M, col = 'red'))
# obviously outliers

subset(reg_data, year == 2022) %>%
  with(plot(day, reading.32M, type = 'b'))
# remove up through May 26th

subset(reg_data, year == 2021) %>%
  with(plot(tv.LGA, reading.14M, type = 'b'))

subset(reg_data, year == 2021) %>%
  with(plot(day, reading.14M, type = 'b'))
# I don't see any rhyme or reason to this

subset(reg_data, year == 2021) %>%
  with(plot(tv.SIFKIL, reading.2R, type = 'b'))
# final days of Sept. meh

subset(reg_data, year == 2023) %>%
  with(plot(tv.BKNYRD, reading.1B, type = 'b'))
subset(reg_data, year == 2023 & day < '2023-05-20') %>%
  with(points(tv.BKNYRD, reading.1B, col = 'red'))
# May-- remove through May 15th

subset(reg_data, year == 2023) %>%
  with(plot(day, reading.1B, type = 'b'))


# update `o4` with these new periods removed
n_clust = nrow(o4$clusters)
o4$clusters[n_clust + 1, ] =
  list(type = 'manual', network = '32M', start = as.Date('2022-05-01'),
       end = as.Date('2022-05-26'), days = 19)
o4$clusters[n_clust + 2, ] =
  list(type = 'manual', network = '1B', start = as.Date('2023-05-01'),
       end = as.Date('2023-05-15'), days = 11)

# What about 14M 2021? Maybe I should remove it altogether, but it stays in for
# now.

idx = with(o4$clusters[n_clust + 1, ],
           which(reg_data$day >= start & reg_data$day <= end))
o4$outliers[idx, '32M'] = T
idx = with(o4$clusters[n_clust + 2, ],
           which(reg_data$day >= start & reg_data$day <= end))
o4$outliers[idx, '1B'] = T

table(o4$outliers) # 271

# use `pairs` with neighboring networks to doublecheck these manual flags
reg_data %>%
  subset(year == 2022,
         select = paste0('reading.', c('32M', '10M', '7M', '22M', '6M'))) %>%
  pairs(col = factor(o4$outliers[reg_data$year == 2022, '32M']))

reg_data %>%
  subset(year == 2023,
         select = paste0('reading.', c('32M', '10M', '7M', '22M', '6M'))) %>%
  pairs(col = factor(o4$outliers[reg_data$year == 2023, '32M']))

reg_data %>%
  subset(year == 2021,
         select = paste0('reading.', c('32M', '10M', '7M', '22M', '6M'))) %>%
  pairs(col = factor(o4$outliers[reg_data$year == 2021, '32M']))

subset(reg_data, year == 2021) %>%
  with(plot(tv.BKNYRD, reading.32M, type = 'b'))

# what's up with 22M?


saveRDS(o4, 'results/load_vs_weather/load_outliers.rds')
#o4 = readRDS('results/load_vs_weather/load_outliers.rds')


cleaned_reg_data = reg_data
for (id in colnames(o4$outliers)) {
  network_col = paste0('reading.', id)
  is_outlier = o4$outliers[, id]
  cleaned_reg_data[is_outlier, network_col] = NA
}

saveRDS(cleaned_reg_data, 'results/load_vs_weather/tv_and_load.rds')
#reg_data = readRDS('results/load_vs_weather/tv_and_load.rds')



# We also need to see what regression is most appropriate for this data. How to
# account for yearly differences?

# Options: varying constant only, varying all coefficients, and separate
# regressions

# to identify networks to use for testing, run a regression with yearly-varying
# constant, and look for biggest changes

network_year_offsets = sapply(networks$id, function(x) {
  network_col = paste0('reading.', x)
  cleaned_reg_data %>%
    transform(load = cleaned_reg_data[, network_col], year = factor(year),
              tv = tv.QNMASP) %>%
    lm(load ~ year + tv + I(tv^2) + I(tv^3), .) %>%
    get_factor_coefs('year')
}) %>%
  t
# how did these change from year to year?

yearly_diffs = sapply(networks$id, function(x) {
  network_col = paste0('reading.', x)
  cleaned_reg_data %>%
    transform(load = cleaned_reg_data[, network_col], year = factor(year),
              tv = tv.QNMASP) %>%
    lm(load ~ year + tv + I(tv^2) + I(tv^3), .) %>%
    get_factor_coefs('year') %>%
    range %>%
    diff
})

test_networks = names(head(sort(yearly_diffs, decreasing = T), 5))

f1 = load ~ year + poly(tv.QNMASP, 3)
f2 = load ~ year * poly(tv.QNMASP, 3)

# what are the p-values?
sapply(test_networks, function(x) {
  network_col = paste0('reading.', x)
  test_data = cleaned_reg_data %>%
    transform(load = cleaned_reg_data[, network_col], year = factor(year),
              tv = tv.QNMASP)
  lm1 = lm(f1, test_data)
  lm2 = lm(f2, test_data)
  anova(lm1, lm2)$`Pr(>F)`[2]
}) # no demonstrated advantage of f2

network_anovas = sapply(networks$id, function(x) {
  network_col = paste0('reading.', x)
  test_data = cleaned_reg_data %>%
    transform(load = cleaned_reg_data[, network_col], year = factor(year),
              tv = tv.QNMASP)
  lm1 = lm(f1, test_data)
  lm2 = lm(f2, test_data)
  anova(lm1, lm2)$`Pr(>F)`[2]
})

hist(network_anovas)
# ok this actually really does look like we should use f2

# combine p-values using Fisher's method
x2 = -2 * sum(log(network_anovas))
pooled_pval = pchisq(x2, df = 2 * length(network_anovas), lower.tail = FALSE)
# 2.90887e-27
saveRDS(pooled_pval, 'results/load_vs_weather/pooled_pval.rds')

adj_p = p.adjust(network_anovas)
adj_p[adj_p < .05]
#           3B          11B          12B          29M           5Q          53M          21M 
# 0.0441906013 0.0148952658 0.0122657404 0.0072943553 0.0383273847 0.0006276768 0.0001077183

# so what's up with those?
changing_cubic_networks = names(head(sort(adj_p[adj_p < .05])))
cleaned_reg_data %>%
  subset(select = c('day', 'tv.QNMASP',
                    paste0('reading.', changing_cubic_networks))) %>%
  reshape(direction = 'long', timevar = 'network', idvar = 'day',
          varying = 2 + 1:length(changing_cubic_networks)) %>%
  transform(year = factor(format(day, '%Y'))) %>%
  ggplot(aes(tv.QNMASP, reading, group = year, color = year)) +
  geom_point() +
  geom_smooth(se = F) +
  facet_wrap(~ network, scales = 'free_y')
# ok then, I'm convinced that the cubic polynomial can change

# how did these change from year to year?
tv82_by_year = sapply(networks$id, function(x) {
  network_col = paste0('reading.', x)
  fit = cleaned_reg_data %>%
    transform(load = cleaned_reg_data[, network_col], year = factor(year),
              tv = tv.QNMASP) %>%
    lm(load ~ year + tv + I(tv^2) + I(tv^3), .)
  predictions = data.frame(year = factor(2021:2023), tv = 82) %>%
    predict(fit, .)
}) %>%
  t

tv82_pdiffs = 100 * (tv82_by_year[, 2:3] - tv82_by_year[, 1:2]) /
  tv82_by_year[, 1:2]
colnames(tv82_pdiffs) = paste0(2021:2022, '-', 2022:2023)

tv82_pdiffs %>%
  as.data.frame %>%
  transform(id = row.names(.)) %>%
  reshape(direction = 'long', timevar = 'year', idvar = 'id', varying = 1:2,
          v.names = 'pdiff', times = colnames(tv82_pdiffs)) %>%
  merge(networks, .) %>%
  plot_network_data(aes(fill = pdiff)) +
  scale_fill_gradient2() +
  facet_wrap(~ year)
# general increase from 2021 to 2022



# # first let's see what this network dimension is about

# peak_cors = cor(peaks_wide[bus_days, -1], use = 'pairwise.complete.obs')
# colnames(peak_cors) = sub('.*\\.', '', colnames(peak_cors))
# row.names(peak_cors) = sub('.*\\.', '', row.names(peak_cors))
# # hist(peak_cors[lower.tri(peak_cors)])
# network_coords = cmdscale(sim2diss(peak_cors))
# plot(network_coords)

# # get TV coefficients for each network
# coef_mat = sapply(network_ids, function(x) {
#   cubic_poly = paste0('tv.KLGA + I(tv.KLGA^2) + I(tv.KLGA^3)')
#   reading_col = paste0('reading.', x)
#   f = as.formula(paste(reading_col, '~', cubic_poly))
#   res = lm(f, data_wide, subset = bus_days)
#   coef(res)
# })

# tv_range = range(tv_wide$tv.KLGA, na.rm = T)
# tv_range = c(floor(tv_range[1]), ceiling(tv_range[2]))

# load_range = range(as.matrix(peaks_wide[, -1]), na.rm = T)

# poly3 = function(tv, coefs) {
#   coefs[1] + coefs[2] * tv + coefs[3] * tv^2 + coefs[4] * tv^3
# }

# f1 = function(tv) poly3(tv, coef_mat[, 1])

# curve(poly3, tv_range[1], tv_range[2], coefs = coef_mat[, 1])

# plot(poly3, tv_range[1], tv_range[2], coefs = coef_mat[, 1])

# # use viridis colors
# cols = viridis::plasma(256)
# network_col_ind = round(256 * (network_coords[, 1] - min(network_coords[, 1])) / diff(range(network_coords[, 1])))
# network_cols = cols[network_col_ind]

# f1 = function(tv) poly3(tv, coef_mat[, 1])
# curve(f1, tv_range[1], tv_range[2], ylim = load_range, col = network_cols[1])
# for (i in 2:ncol(coef_mat)) {
#   f1 = function(tv) poly3(tv, coef_mat[, i])
#   curve(f1, tv_range[1], tv_range[2], add = T, col = network_cols[i])
# }

# f2 = function(tv) poly3(tv, coef_mat[, 1]) / poly3(tv_range[2], coef_mat[, 1])
# curve(f2, tv_range[1], tv_range[2], ylim = c(0, 1), col = network_cols[1])
# for (i in 2:ncol(coef_mat)) {
#   f2 = function(tv) poly3(tv, coef_mat[, i]) / poly3(tv_range[2], coef_mat[, i])
#   curve(f2, tv_range[1], tv_range[2], add = T, col = network_cols[i])
# } # hmm

# plot(apply(peaks_wide[, -1], 2, max), network_coords[, 1]) # nope not it

# pairs(cbind(network_coords[, 1], t(coef_mat))) # ?

# networks %>%
#   transform(coord = network_coords[id, 1]) %>%
#   plot_station_data(aes(fill = coord), nyc_base) +
#   scale_fill_viridis_c()
# # it seems to be closely related to peak load density

# tv_r2 = sapply(network_ids, function(x) {
#   cubic_poly = paste0('tv.KLGA + I(tv.KLGA^2) + I(tv.KLGA^3)')
#   reading_col = paste0('reading.', x)
#   f = as.formula(paste(reading_col, '~', cubic_poly))
#   res = lm(f, data_wide, subset = bus_days)
#   summary(res)$r.squared
# })

# networks %>%
#   transform(r2 = tv_r2[id]) %>%
#   plot_station_data(aes(fill = r2), nyc_base) +
#   scale_fill_viridis_c()

# plot(network_coords[, 1], tv_r2)

# # 2023 only-- to see if it's covid related
# network_coords2023 = peaks_wide %>%
#   subset(day > as.Date('2023-01-01'), select = -day) %>%
#   cor(use = 'pairwise.complete.obs') %>%
#   sim2diss %>%
#   cmdscale
# row.names(network_coords2023) = sub('.*\\.', '', row.names(network_coords2023))
# plot(network_coords2023)

# networks %>%
#   transform(coord = network_coords2023[id, 1]) %>%
#   plot_station_data(aes(fill = coord), nyc_base) +
#   scale_fill_viridis_c()
# # nope, not covid related

# # ok, so it's related to density, but can we figure out what it is?

# resid_mat = sapply(network_ids, function(x) {
#   cubic_poly = paste0('tv.KLGA + I(tv.KLGA^2) + I(tv.KLGA^3)')
#   reading_col = paste0('reading.', x)
#   f = as.formula(paste(reading_col, '~', cubic_poly))
#   res = lm(f, data_wide, subset = bus_days, na.action = na.exclude)
#   resid(res)
# })
# row.names(resid_mat) = as.character(data_wide$day[bus_days])

# pc1 = princomp(na.omit(resid_mat))

# plot(as.Date(row.names(pc1$scores)), pc1$scores[, 'Comp.1'])
# # low outliers, what could they be?

# pc1$scores[head(order(pc1$scores[, 'Comp.1'])), 'Comp.1']

# # 2021 had some weirdness
# resid_pc1 = data.frame(day = as.Date(row.names(pc1$scores)),
#                        pc1 = pc1$scores[, 'Comp.1'])
# resid_pc1 %>%
#   subset(day < as.Date('2022-01-01')) %>%
#   with(plot(day, pc1, type = 'l'))

# resid_pc1 %>%
#   subset(day > as.Date('2021-06-15') & day < as.Date('2021-07-15')) %>%
#   with(plot(day, pc1, type = 'l'))






# OK, what we're going to do is take the weather coordinates, and map the
# networks onto that space

tv_locations = readRDS('results/weather/tv_mds.rds')
pc_tv = prcomp(tv_locations$conf)



# dist1 = get_tv_network_distances(reg_data)

# unf1 = get_unfolding(dist1$distance)

# plot(unf1)

# sb1 = extract_sea_breeze_coords(unf1)
# plot_sea_breeze_dim(sb1$stations, sb1$networks, sbdim = 1)
# plot_sea_breeze_dim(sb1$stations, sb1$networks, sbdim = 2)
# # This seems like the most accurate way to capture the sea breeze dimension

# sb2 = extract_sea_breeze_coords(unf1, site_locations)
# plot_sea_breeze_dim(sb2$stations, sb2$networks, sbdim = 1)
# # ok this makes a really big difference. This is actually kind of crappy as a
# # sea breeze though

# plot_sea_breeze_dim(sb2$stations, sb2$networks, sbdim = 2)
# # oh look it's the sea breeze again

dist3 = get_tv_network_distances2(reg_data)

# quick check to see what's better/worse

networks %>%
  transform(mean_dist = rowMeans(dist3[id, ])) %>%
  plot_station_data(aes(fill = mean_dist)) +
  scale_fill_viridis_c()

stations %>%
  transform(mean_dist = colMeans(dist3[, stid])) %>%
  plot_station_data(aes(color = mean_dist, shape = network)) +
  scale_color_viridis_c(direction = -1) + 
  scale_shape_manual('Weather network', values = c(15, 19, 17))

# # mean R^2 by network
# networks$mean_r2 = rowMeans(r2_mat)[match(networks$id, row.names(r2_mat))]
# plot(networks[, 'mean_r2'])

# So the unfolding results capture sea breeze inconsistently. How about
# PCA/factor analysis?

pc_dist = prcomp(dist3)
summary(pc_dist)
plot(pc_dist) # looks like 2 notable components
saveRDS(pc_dist, 'results/load_vs_weather/tv_load_pca.rds')
# pc_dist = readRDS('results/5/tv_load_pca.rds')

plot(pc_dist$x[, 1:2])
text(pc_dist$x[, 1:2], labels = row.names(pc_dist$x))
# dim 2 could be sea breeze

plot(pc_dist$rotation[, 1:2])
text(pc_dist$rotation[, 1:2], labels = row.names(pc_dist$rotation))

stations %>%
  transform(sea_breeze = pc_dist$rotation[stid, 2]) %>%
  plot_station_data(aes(color = -sea_breeze))
# nailed it, wow

networks %>%
  transform(sea_breeze = pc_dist$x[id, 2]) %>%
  plot_network_data(aes(fill = sea_breeze))
# same as unfolding pattern

# what is dim 1?
stations %>%
  transform(dim1 = pc_dist$rotation[stid, 1]) %>%
  plot_station_data(aes(color = dim1))
# not sure this is meaningful

networks %>%
  transform(dim1 = pc_dist$x[id, 1]) %>%
  transform(logdim1 = log(dim1 - min(dim1) + .2)) %>%
  plot_network_data(aes(fill = logdim1))

# # ok let's get super clever-- reconstruct dist3 without 1st component, then get
# # distances
# dist3_filtered = pc_dist$x[, -1] %*% t(pc_dist$rotation[, -1])
# dist3_filtered = dist3_filtered - min(dist3_filtered) + .1

# unf_filtered = get_unfolding(dist3_filtered)
# plot(unf_filtered)
# # the stations seem basically good, but networks need more help

# unf_filtered2 = get_unfolding(dist3_filtered, lambda = .1, omega = 5)
# plot(unf_filtered2)

# unf_filtered3 = get_unfolding(dist3_filtered, lambda = .6, omega = .6)
# plot(unf_filtered3) # closer
# # anyway, yeah, this screwed up the network locations too much to be used


# # mean R^2 by station
# nycm_sites$mean_r2 = colMeans(r2_mat)[match(nycm_sites$stid, colnames(r2_mat))]
# plot(nycm_sites[, 'mean_r2'], pch = 19)
# hist(nycm_sites$mean_r2)

unf3 = get_unfolding(dist3)
saveRDS(unf3, 'results/load_vs_weather/tv_load_unfolding.rds')
# unf3 = readRDS('results/tv_load_unfolding.rds')

plot(unf3, plot.type = "Shepard")

plot(unf3)
# add sea breeze dimension line
sb5 = extract_sea_breeze_coords(unf3)
abline(0, sb5$rotation[2, 1] / sb5$rotation[1, 1], lty = 'dashed')

# not the best, try again
plot(unf3)
# add sea breeze dimension line
sb6 = extract_sea_breeze_coords(unf3, use_networks = T)
abline(0, sb6$rotation[2, 1] / sb6$rotation[1, 1], lty = 'dashed')

plot_sea_breeze_dim(sb6$stations, sb6$networks, sbdim = 1, include_nj = T)
plot_sea_breeze_dim(sb6$stations, sb6$networks, sbdim = 1)
plot_sea_breeze_dim(sb6$stations, sb6$networks, sbdim = 2, include_nj = T)
# subtract median network score to align zero with networks
dim2_network_med = median(sb6$networks[, 2])
plot_sea_breeze_dim(sb6$stations - dim2_network_med,
                    sb6$networks - dim2_network_med, sbdim = 2,
                    include_nj = T) +
  scale_color_continuous_diverging() +
  scale_fill_continuous_diverging()
# clear geographical pattern, seems related to population density

stations %>%
  transform(dim2 = sb6$stations[stid, 2] - dim2_network_med) %>%
  plot_station_data(aes(color = dim2, shape = network)) +
  scale_color_gradient2()

# just networks
networks %>%
  transform(sea_breeze = sb6$networks[id, 1]) %>%
  plot_network_data(aes(fill = sea_breeze))

# very high penalties
x8 = get_unfolding(dist3, lambda = .1, omega = 5)
plot(x8, main = 'Network/Station relationships')

# lower penalties
x12 = get_unfolding(dist3, lambda = .9, omega = .2)
plot(x12, main = 'Network/Station relationships')
sb5 = extract_sea_breeze_coords(x12, use_networks = T)
plot_sea_breeze_dim(sb5$stations, sb5$networks, sbdim = 1)

#sb3 = extract_sea_breeze_coords(unf2)
# plot_sea_breeze_dim(sb3$stations, sb3$networks, sbdim = 1, include_nj = T)
# plot_sea_breeze_dim(sb3$stations, sb3$networks, sbdim = 2, include_nj = T)

# what if we exclude weird stuff?
unf4 = get_unfolding(dist3[!row.names(dist3) %in% c('14M', '32M'), ])
plot(unf4)
unf4 = get_unfolding(dist3[!row.names(dist3) %in% c('14M', '32M'), ])
plot(unf4)
# add sea breeze dimension line
sb6 = extract_sea_breeze_coords(unf3, use_networks = T)
abline(0, sb6$rotation[2, 1] / sb6$rotation[1, 1], lty = 'dashed')

# what if we exclude sites in the middle?
unf4 = get_unfolding(dist3[, !colnames(dist3) %in% c('JFK', 'JRB', 'QNSOZO', 'QNCORO')])
plot(unf4) # more reasonable again

# # exclude the worst sites
# sort(get_factor_coefs(dist_lm, 'stid'))
# worst_5 = names(tail(sort(get_factor_coefs(dist_lm, 'stid')), 5))
# unf5 = get_unfolding(dist3[, !colnames(dist3) %in% worst_5])
# plot(unf5) # didn't help

unf4 = get_unfolding(dist3[, !colnames(dist3) %in% c('JRB')])
plot(unf4) # more reasonable again


unf5 = get_unfolding(dist3, fix_sites = tv_locations$conf)
plot(tv_locations)
plot(unf5)

# very high penalties
unf6 = get_unfolding(dist3, fix_sites = tv_locations$conf, lambda = .1, omega = 5)
plot(unf6, main = 'Network/Station relationships')

x12 = get_unfolding(dist3, fix_sites = tv_locations$conf, lambda = .9, omega = .2)
plot(x12, main = 'Network/Station relationships')


# # ok so now I just want to see if I can get this stretched out a bit to match up
# # site/network locations

# # tweak penalties-- low penalties
# x7 = get_unfolding(dist2$distance, lambda = 1, omega = .1)
# plot(x7, main = 'Network/Station relationships')
# # hahaha this is so awful

# # very high penalties
# x8 = get_unfolding(dist2$distance, lambda = .1, omega = 5)
# plot(x8, main = 'Network/Station relationships')
# sb4 = extract_sea_breeze_coords(x8)
# plot_sea_breeze_dim(sb4$stations, sb4$networks, sbdim = 1)
# # it's good actually?

# # lower penalties
# x12 = get_unfolding(dist2$distance, lambda = .6, omega = .6)
# plot(x12, main = 'Network/Station relationships') # similar
# sb5 = extract_sea_breeze_coords(x12)
# plot_sea_breeze_dim(sb5$stations, sb5$networks, sbdim = 1)

# x9 = get_unfolding(dist2$distance, lambda = .7, omega = .5)
# plot(x9, main = 'Network/Station relationships')
# # kind of broken, but in an interesting way

# x10 = get_unfolding(dist2$distance, lambda = .8, omega = .3)
# plot(x10, main = 'Network/Station relationships')
# # broken




# next step: regression to identify site quality and load TV predictability
# (just like earlier)

# get data frame with network, station, R2 distance, and unfolding distance

tv_load_dist = as.data.frame(dist3) %>%
  reshape(direction = 'long', varying = 1:ncol(dist3), v.names = 'r2_dist',
          timevar = 'stid', idvar = 'network', ids = row.names(dist3),
          times = colnames(dist3))

tv_load_dist$sea_breeze_dist =
  with(tv_load_dist, abs(sb6$stations[stid, 1] - sb6$networks[network, 1]))

dist_lm = lm(r2_dist ~ network + stid + sea_breeze_dist, tv_load_dist)

tail(summary(dist_lm)$coefficients, 1)
#                   Estimate   Std. Error  t value      Pr(>|t|)
# sea_breeze_dist 0.02954703 0.0007444187 39.69141 3.597263e-250

stations %>%
  transform(val = get_factor_coefs(dist_lm, 'stid')[stid]) %>%
  plot_station_data(aes(color = val, shape = network))

networks %>%
  subset(id != '14M') %>% # 14M is a high outlier
  transform(val = get_factor_coefs(dist_lm, 'network')[id]) %>%
  plot_network_data(aes(fill = val))
# I think this is associated with peak loads
