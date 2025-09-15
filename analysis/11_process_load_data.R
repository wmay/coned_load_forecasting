# examining the relationship between load and weather, especially sea breeze
# patterns

library(timeDate) # holidays
library(MASS) # robust regression
library(magrittr)

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

get_robust_residuals = function(load, nload, year) {
  # Predict network load using the load from a neighboring network
  # Some networks are missing years, making `rlm` complain. Just remove the
  # missing year levels
  data_years = unique(year[!is.na(load) & !is.na(nload)])
  year[!year %in% data_years] = data_years[1]
  if (length(data_years) == 1) {
    f = load ~ nload
  } else {
    f = load ~ factor(year) * nload
  }
  # `na.exclude` is required here to make the indexing consistent!
  res = try(rlm(f, na.action = na.exclude))
  if (inherits(res, 'try-error')) return(rep(NA, length(load)))
  resid(res)
}

get_outlier_scores = function(network, k, cluster = FALSE, cluster_n = 5) {
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
  load = reg_data[, paste0('reading.', network)]
  borough = sub('[0-9]*', '', network)
  neighbors = with(networks, id[endsWith(id, borough) & id != network])
  # 14M is unusable due to noise
  neighbors = neighbors[neighbors != '14M']
  out_list = lapply(neighbors, function(nid) {
    nload = reg_data[, paste0('reading.', nid)]
    resid = get_robust_residuals(load, nload, year)
    c(fun(resid), list(resid = resid))
  })
  # make a matrix for each output
  out = list()
  for (v in c('high', 'low', 'resid')) {
    out[[v]] = do.call(cbind, lapply(out_list, getElement, name = v))
    colnames(out[[v]]) = neighbors
  }
  out
}

plot_outlier_scores = function(network, k, cluster = FALSE, cluster_n = 5) {
  x = get_outlier_scores(network, k, cluster = FALSE, cluster_n = 5)
  orig_mfrow = par('mfrow')
  on.exit(par(mfrow = orig_mfrow))
  n_neighbors = ncol(x$resid)
  n_neighbors = min(n_neighbors, 6)
  par(mfrow = c(n_neighbors, 1))
  for (i in seq_len(n_neighbors)) {
    yvals = x$resid[, i]
    outlier = x$low[, i] | x$high[, i]
    # calculate the fence boundaries to plot
    q13 = quantile(x$resid[, i], c(.25, .75), na.rm = TRUE)
    tf_low = q13[1] - k * diff(q13)
    tf_high = q13[2] + k * diff(q13)
    yrange = range(yvals, na.rm = T)
    ylim = c(max(tf_high, yrange[2]), min(tf_low, yrange[1]))
    color = ifelse(outlier, 'red', 'black')
    plot(reg_data$day, yvals, col = color,
         ylim = ylim, type = 'b')
    abline(h = 0, col = '#66666688')
    abline(h = c(tf_low, tf_high), lty = 'dashed')
  }
}

# get matrix of confirmed outliers, and list of outlier clusters
get_network_outliers = function(k = 3, cluster_k = 1.5, cluster_n = 5) {
  cluster_list = list()
  outliers = matrix(F, nrow = nrow(reg_data), ncol = nrow(networks),
                    dimnames = list(as.character(reg_data$day), networks$id))
  single_outliers = matrix(F, nrow = nrow(reg_data), ncol = nrow(networks),
                           dimnames = list(as.character(reg_data$day), networks$id))
  vote_res = function(x) rowMeans(x, na.rm = TRUE) >= .6
  for (i in 1:nrow(networks)) {
    id = networks$id[i]
    #message(paste0('reading.', id))
    load = reg_data[, paste0('reading.', id)]
    # first the clusters
    cluster_scores = get_outlier_scores(id, cluster_k, TRUE, cluster_n)
    high_clusters = get_true_clusters(vote_res(cluster_scores$high), cluster_n)
    low_clusters = get_true_clusters(vote_res(cluster_scores$low), cluster_n)
    # save cluster info
    low_clusters$clusters$type = rep('low', nrow(low_clusters$clusters))
    high_clusters$clusters$type = rep('high', nrow(high_clusters$clusters))
    cluster_list[[i]] = rbind(high_clusters$clusters, low_clusters$clusters)
    cluster_list[[i]]$cluster_id = NULL
    cluster_list[[i]]$network = rep(id, nrow(cluster_list[[i]]))
    outlier_cluster = vote_res(cluster_scores$high) |
      vote_res(cluster_scores$low)
    # now the isolated outliers (removing clustered outliers)
    load[outlier_cluster] = NA
    regular_scores = get_outlier_scores(id, k, cluster = FALSE)
    single_outlier = vote_res(regular_scores$high) | vote_res(regular_scores$low)
    single_outliers[, id] = single_outlier
    outliers[, id] = outlier_cluster | single_outlier
  }
  outlier_clusters = do.call(rbind, cluster_list) %>%
    transform(start = reg_data$day[idx.start], end = reg_data$day[idx.end],
              days = idx.end - idx.start + 1) %>%
    subset(select = -c(idx.start, idx.end))
  outliers[is.na(outliers)] = FALSE
  list(outliers = outliers, clusters = outlier_clusters,
       single_outliers = single_outliers)
}


networks = readRDS('results/maps/coned_networks_cleaned.rds')

peaks = list.files('data/coned', 'Network.*\\.csv', full.names = T) %>%
  lapply(read.csv) %>%
  do.call(rbind, .) %>%
  transform(DT = replace(DT, DT == '', NA),
            # some inconsistency, but various forms of 'False' mean data is fine
            BAD = !startsWith(BAD, 'F'),
            Network = replace(Network, Network == '#N/A', NA)) %>%
  transform(DT = as.POSIXct(DT, tz = 'EST5EDT',
                            tryFormats = c('%m/%d/%Y %H:%M:%S', '%m/%d/%Y %H:%M'))) %>%
  subset(!BAD & !is.na(DT)) %>%
  subset(as.POSIXlt(DT)$hour %in% 9:21) %>%
  # get daily load peaks
  transform(day = as.Date(DT, tz = 'EST5EDT')) %>%
  aggregate(READING ~ Network + day, ., max)
names(peaks) = tolower(names(peaks))
peaks_wide = reshape(peaks, direction = 'wide', idvar = 'day',
                     timevar = 'network')

bus_days = isBizday(as.timeDate(peaks_wide$day), holidays = holidayNYSE(2021:2025))

# is there an extra network?
peak_ids = sub('reading.', '', names(peaks_wide)[-1], fixed = T)
peak_ids[!peak_ids %in% networks$id]
min(peaks_wide$day[!is.na(peaks_wide$reading.47M)]) # 2025-05-01
# hunh

# also need the system data
sys1 = read.csv('data/coned/Network Data - SUNY 2025_05_01-2025_07_10.csv',
                check.names = F) %>%
  subset(`Network Name` == 'CECONY System') %>%
  subset(select = -c(`Network Name`, Network, Borough))
sys2 = read.csv('data/coned/Borough and System Data 2020-2024.csv') %>%
  subset(Borough == 'CECONY') %>%
  subset(select = -Borough)
sys_peaks = rbind(sys1, sys2) %>%
  transform(DT = as.POSIXct(DT, tz = 'EST5EDT', '%m/%d/%Y %H:%M'),
            # some inconsistency, but various forms of 'False' mean data is fine
            BAD = !startsWith(BAD, 'F')) %>%
  subset(!is.na(DT) & !BAD) %>%
  subset(as.POSIXlt(DT)$hour %in% 9:21) %>%
  # get daily load peaks
  transform(day = as.Date(DT, tz = 'EST5EDT')) %>%
  aggregate(READING ~ day, ., max)
names(sys_peaks) = tolower(names(sys_peaks))

peaks_wide$reading.system =
  sys_peaks$reading[match(peaks_wide$day, sys_peaks$day)]


# this is the data that will actually be used for regressions-- remove outliers
# first!!
reg_data = peaks_wide[bus_days, ] %>%
  transform(year = as.integer(substr(trunc(day, 'year'), 1, 4)))

# we're raising `k` compared to the earlier script, because the residuals are
# smaller, now that we're predicting with neighboring loads instead of TV
# o4 = get_network_outliers(k = 2.5, cluster_k = 1, cluster_n = 5)
o5 = get_network_outliers(k = 3.5, cluster_k = 1.4, cluster_n = 5)

o4 = readRDS('results/load_vs_weather/load_outliers.rds')
# # let's compare
# o5$clusters[order(o5$clusters$network), ]
# o4$clusters[order(o4$clusters$network), ]
# table(o4$outliers)
# table(o5$outliers)
# table(o5$single_outliers) # wow
# sum(o4$clusters$days[o4$clusters$type != 'manual'])
# sum(o5$clusters$days)
# sort(colSums(o4$outliers), decreasing = T)
# sort(colSums(o5$outliers), decreasing = T)

# remove the manually flagged outliers from `o4`
idx = with(o4$clusters[15, ],
           which(reg_data$day >= start & reg_data$day <= end))
o5$outliers[idx, '32M'] = T
idx = with(o4$clusters[16, ],
           which(reg_data$day >= start & reg_data$day <= end))
o5$outliers[idx, '1B'] = T
saveRDS(o5, 'results/process_load_data/load_outliers.rds')

# remove outlier values. Note that we're ignoring 47M for now
stopifnot(nrow(o5$outliers) == nrow(reg_data))

table(is.na(reg_data[, -1])) # 783
for (n in colnames(o5$outliers)) {
  reg_col = paste0('reading.', n)
  n_outliers = o5$outliers[, n]
  reg_data[n_outliers, reg_col] = NA
}
table(is.na(reg_data[, -1])) # 987

reg_data %>%
  subset(select = -year) %>%
  write.csv(file = 'results/process_load_data/loads.csv', row.names = FALSE)
