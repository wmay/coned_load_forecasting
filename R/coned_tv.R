# calculate ConEd's effective temperature and temperature variable (TV)

library(psychrolib)
library(magrittr)

SetUnitSystem('SI')

as_celsius = function(x) (x - 32) * 5/9
as_fahrenheit = function(x) (x * 9/5) + 32

lag_vec = function(x, lag = 1, fill = NA) c(rep(fill, lag), head(x, -lag))

rolling_mean3 = function(time, val, units = c('hours', 'days'),
                         weights = rep(1, 3)) {
  units = match.arg(units)
  # assume times are in order, but we may be missing some, thus making some
  # rolling calculations invalid
  correct_lag = diff(time, 2) == as.difftime(2, units = units)
  valid = c(FALSE, FALSE, correct_lag)
  val3 = weights[1] * val + weights[2] * lag_vec(val) +
    weights[3] * lag_vec(val, 2)
  val3[!valid] = NA
  val3 / sum(weights)
}

# Round .5 up, unlike R's rounding. This is relevant due to the prevalance of .5
# in averaged station temps. Using machine epsilon for a poor man's version of
# `nextafter` in C/Python
coned_round = function(x) round(x + sqrt(.Machine$double.eps))

# Get wet bulb temperature from temperature (`dbc`) and dewpoint (`dpc`).
# Pressure is assumed to be 1013.25
coned_wet_bulb = function(DBC, DPC) {
  # I'm vaguely curious about where this calculation came from, but not going to
  # look it up for now
  PR = 1013.25
  # my $DB= $ARGV[1];    # Dry Bulb
  # my $DP= $ARGV[2];    # Dew Point
  # $DBC=(($DB-32)*5/9);
  # $DPC=(($DP-32)*5/9);
  # my $DDEPC=($DBC-$DPC)*0.18;
  DDEPC = (DBC - DPC) * 0.18
  # $WBC=${DBC}-(0.035*${DDEPC}-0.00072*${DDEPC}*(${DDEPC}-1))*(${DBC}+${DPC}+95.556-(${PR}/30.474));
  WBC = DBC - (0.035 * DDEPC - 0.00072 * DDEPC * (DDEPC - 1)) * (DBC + DPC + 95.556 - (PR / 30.474))
  # $WBF=($WBC*9/5)+32;
  coned_round(as_fahrenheit(WBC))
  # we're supposed to round this in fahrenheit here, but can't because our temps
  # are still in celsius
}

dry_wet_ave = function(temp, rh, pres) {
  # Note: `CalcPsychrometricsFromRelHum` fails with missing data
  missing = is.na(temp) | is.na(rh) | is.na(pres)
  wetbulb = rep(NA, length(temp))
  if (all(missing)) return(wetbulb)
  wetbulb[!missing] =
    CalcPsychrometricsFromRelHum(temp[!missing], rh[!missing],
                                 pres[!missing])$TWetBulb
  as_fahrenheit((temp + wetbulb) / 2)
}

# Calculate ConEd's "effective temperature": Effective temperature = (ET)
# Highest consecutive three (3) hour average of the dry and wet-bulb
# temperatures (hour ending 9AM – 9PM, so 7-9am to 7-9pm)
coned_effective_temp = function(t, temp, rh, pres) {
  out = data.frame(time = t, temp = temp, rh = rh, pres = pres) %>%
    subset(!(is.na(temp) | is.na(rh) | is.na(pres))) %>%
    transform(drywet = dry_wet_ave(temp, rh, pres)) %>%
    transform(ave_3hr = rolling_mean3(time, drywet, 'hours')) %>%
    # remove night/morning
    transform(hour_of_day = as.integer(format(time, '%H'))) %>%
    subset(hour_of_day >= 9 & hour_of_day <= 21) %>%
    # get max by day
    transform(day = as.Date(time, tz = 'EST5EDT'), eff_temp = ave_3hr) %>%
    aggregate(eff_temp ~ day, ., max, na.rm = T)
}

# ConEd's "temperature variable"
coned_tv = function(t, temp, rh, pres) {
  coned_effective_temp(t, temp, rh, pres) %>%
    transform(tv = rolling_mean3(day, eff_temp, 'days', c(.7, .2, .1))) %>%
    subset(select = c(day, tv))
}
