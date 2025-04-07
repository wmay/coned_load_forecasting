# figure out what station attributes are associated with better predictions

source('R/plots.R')

unf = readRDS('results/tv_load_unfolding.rds')

plot(unf)
# add sea breeze dimension line
unf_pc = prcomp(unf$conf.row)
sb6 = list(stations = predict(unf_pc, unf$conf.col),
           networks = unf_pc$x, rotation = unf_pc$rotation)
abline(0, sb6$rotation[2, 1] / sb6$rotation[1, 1], lty = 'dashed')

dim2_network_med = median(sb6$networks[, 2])
stations$dim2 = sb6$stations[stations$stid, 2] - dim2_network_med
plot_sea_breeze_dim(sb6$stations - dim2_network_med,
                    sb6$networks - dim2_network_med, sbdim = 2,
                    include_nj = T) +
  scale_color_continuous_diverging() +
  scale_fill_continuous_diverging()
# *shrug*

stations %>%
  plot_station_data(aes(color = dim2, shape = network)) +
  scale_color_gradient2()

sd1 = read.csv('papers/station_data1.csv', skip = 3, check.names = F)
names(sd1) = sub('\\n', ' ', names(sd1))
sd1$`Station ID` = sub('[0-9]', '', sd1$`Station ID`)
row.names(sd1) = sd1$`Station ID`

stations %>%
  subset(network != 'ASOS') %>%
  transform(Location = sd1[stid, 'Location']) %>%
  plot_station_data(aes(color = dim2, shape = Location)) +
  scale_color_gradient2()

sd1$dim2 = sb6$stations[match(sd1$`Station ID`, row.names(sb6$stations)), 2]

sd1 %>%
  subset(select = -`Station ID`) %>%
  pairs()

sd1['BKBROW', 'Percent Build/Impervious Surface Fraction (25m radius)'] = 87.5
sd1$`Aspect Ratio` = sub(',', '\\.', sd1$`Aspect Ratio`)
for (i in c(8:10, 12)) sd1[, i] = as.numeric(sd1[, i])

pairs(as.matrix(sd1[, c(8:10, 12, 14)]))
# not much help there

surface_tab = sort(table(sd1[, 'Station Surface Type']), decreasing = T)
surface_dict = setNames(ifelse(surface_tab > 1, names(surface_tab), 'Other'),
                        names(surface_tab))

stations %>%
  subset(network != 'ASOS') %>%
  transform(surface = surface_dict[sd1[stid, 'Station Surface Type']]) %>%
  plot_station_data(aes(color = dim2, shape = surface)) +
  scale_color_gradient2()


sd2 = read.csv('papers/station_data2.csv', skip = 3, check.names = F)
sd2$`Station ID` = sub(' [0-9 ]', '', sd2$`Station ID`)
row.names(sd2) = sd2$`Station ID`
sd2$dim2 = sb6$stations[match(sd2$`Station ID`, row.names(sb6$stations)), 2]

pairs(as.matrix(sd2[, c(3, 5:7)]))

stations %>%
  subset(network != 'ASOS') %>%
  transform(surface = surface_dict[sd1[stid, 'Station Surface Type']]) %>%
  plot_station_data(aes(color = dim2, shape = surface)) +
  scale_color_gradient2()

sd2$utz_cat = substr(sd2$UTZ, 1, 2)
sd2$utz_n = as.integer(substr(sd2$UTZ, 3, 3))

ggplot(sd2, aes(UTZ, dim2, color = utz_cat, shape = factor(utz_n))) +
  geom_point()
