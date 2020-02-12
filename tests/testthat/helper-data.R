library(lubridate, warn.conflicts = FALSE)

month_starts <- as.Date(c('2016-01-02', '2016-02-02', '2016-03-02'))
regions <- 1 : 2
model_ids <- 1 : 3

emissions <- expand.grid(
  month_start = month_starts,
  region = regions
) %>%
  mutate(flux_density = 0)

control <- data.frame(
  model_id = model_ids,
  time = ymd_hm(c('2016-01-03 10:00', '2016-02-03 10:00', '2016-03-03 10:00'))
) %>%
  mutate(
    xco2 = 0,
    latitude = 0
  )

sensitivities <- expand.grid(
  from_month_start = month_starts,
  model_id = model_ids,
  region = regions
) %>%
  mutate(xco2_sensitivity = 1)

soundings <- data.frame(
  sounding_id = 1 : 2,
  instrument_mode = c('LN', 'LG'),
  xco2_error = 1
)
