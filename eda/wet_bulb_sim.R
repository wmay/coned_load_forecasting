
library(psychrolib)
library(MASS)

SetUnitSystem("SI")

as_celsius = function(x) (x - 32) * (5/9)
as_fahrenheit = function(x) x * (9/5) + 32

get_wet_bulb = function(temp, dwpt) {
  as_fahrenheit(GetTWetBulbFromTDewPoint(as_celsius(temp), as_celsius(dwpt), 101325.0))
}

# what do the error distributions look like?
simulations = 1000
temps = rnorm(simulations, 82, sd = 1)
dwpts = rnorm(simulations, 63, sd = 5)
true_wet_bulb = get_wet_bulb(82, 63)
results = get_wet_bulb(temps, dwpts)

# How skewed are the wet bulb estimates calculated from a single station?
hist(results, main = 'Wet Bulb Estimates (true value is dotted line)')
abline(v = true_wet_bulb, lty = 'dotted') # 69.1
abline(v = mean(results), col = 'red') # result from averaging wet bulbs
abline(v = get_wet_bulb(mean(temps), mean(dwpts)), col = 'blue')
# this is a lot closer than I was expecting, and the distribution is barely
# skewed at all

# is the relationship even curved?
curve(get_wet_bulb(82, x), 50, 80, xlab = 'Dewpoint', ylab = 'Wet Bulb')
# just barely


# simulate weather measurements and see how each calculation compares to the
# true value
simulate_weather = function(stations, temp, dwpt, sd_temp, sd_dwpt, covariance = 0) {
  sigma = matrix(c(sd_temp^2, covariance, covariance, sd_dwpt^2), 2)
  random_vals = mvrnorm(stations, c(temp, dwpt), sigma)
  temps = random_vals[, 1]
  dwpts = random_vals[, 2]
  if (any(dwpts > temps)) {
    warning('Dewpoint greater than temperature. Result will be biased')
    dwpts = min(dwpts, temps)
  }
  true_wet_bulb = get_wet_bulb(temp, dwpt)
  ave_wet_bulb = mean(get_wet_bulb(temps, dwpts))
  ave_dwpt = get_wet_bulb(mean(temps), mean(dwpts))
  c(ave_wet_bulb = ave_wet_bulb - true_wet_bulb,
    ave_dwpt = ave_dwpt - true_wet_bulb)
}

# try to see which is a better estimate of homogeneous wet bulb
compare_calcs = function(n, stations, temp, dwpt, sd_temp, sd_dwpt,
                      covariance = 0) {
  res = replicate(n, {
    simulate_weather(stations, temp, dwpt, sd_temp, sd_dwpt, covariance)
  })
  res2 = rbind(res, colMeans(res))
  row.names(res2) = c('ave_wet_bulb', 'ave_dwpt', 'ave_both')
  biases = rowMeans(res2)
  bias_ses = apply(res2, 1, sd) / sqrt(n)
  maes = rowMeans(abs(res2))
  mses = rowMeans(res2^2)
  list(mae = maes, mse = mses, biases = biases, bias_se = bias_ses,
       corr = cor(res[1, ], res[2, ]), diffs = summary(res[1, ] - res[2, ]))
}


# just checking that my covariance matrix seems reasonable
sigma = matrix(c(1^2, .7, .7, 2^2), 2)
plot(mvrnorm(100, c(82, 63), sigma), xlab = 'Temperature', ylab = 'Dewpoint')
# it looks reasonable

# Which calculation provides a better estimate?
compare_calcs(10000, 2, 82, 60, 2, 4, 1)
# dewpoint averaging is better, but both are (slightly) biased

# Does dewpoint averaging still have a bias if there's no temperature variation?
# (Covariance must also be zero here.)
compare_calcs(10000, 2, 82, 60, 0, 4, 0)
# Yes, it's still biased

# $biases
# ave_wet_bulb     ave_dwpt     ave_both 
#   0.10222805   0.05793747   0.08008276 

# $bias_se
# ave_wet_bulb     ave_dwpt     ave_both 
#   0.01533341   0.01530803   0.01531751


# OK, so far the results of both methods have been very similar. What happens
# when we add a station with more weather variation, like JFK?
sigma = matrix(c(4^2, 3, 3, 4^2), 2)
plot(mvrnorm(100, c(82, 60), sigma), xlab = 'Temperature', ylab = 'Dewpoint')

compare_calcs(10000, 3, 82, 60, 4, 4, 3)
# All the errors got substantially smaller due to the additional station. Not
# clear it's worth caring about at this point

# OK, but if there's a 5 degree F difference between two stations are the two
# calculations different enough to matter? It should be, right?
different_temps = c(82, 77)
different_dwpts = c(63, 62)
mean(get_wet_bulb(different_temps, different_dwpts))
get_wet_bulb(mean(different_temps), mean(different_dwpts))
# this is an absurdly small difference!

# Is it the dewpoints that matter more?
different_temps = c(82, 82)
different_dwpts = c(63, 53)
mean(get_wet_bulb(different_temps, different_dwpts))
get_wet_bulb(mean(different_temps), mean(different_dwpts))
# still an absurdly small difference!

# Is it the dewpoints that matter more?
different_temps = c(77, 77)
different_dwpts = c(63, 53)
mean(get_wet_bulb(different_temps, different_dwpts))
get_wet_bulb(mean(different_temps), mean(different_dwpts))
# still an absurdly small difference!

# Is it the dewpoints that matter more?
different_temps = c(82, 82)
different_dwpts = c(75, 65)
mean(get_wet_bulb(different_temps, different_dwpts))
get_wet_bulb(mean(different_temps), mean(different_dwpts))
# still an absurdly small difference!



# more realistic
compare_calcs(10000, 2, 82, 70, 1, 2.5, .7)
