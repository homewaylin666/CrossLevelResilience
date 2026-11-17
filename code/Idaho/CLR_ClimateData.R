################################################################################
## Build the complete monthly climate record and the seasonal climate variables (Idaho)
##
## Written by Hung-wei Lin Nov 2026
##
## Station: U.S. Sheep Experiment Station, Dubois (GHCN-Daily USC00102707, 1661 m), the station of the
## monthly records in the data set (Zachmann et al. 2010; Dalgleish et al. 2011).
## (1) Monthly mean temperature and total precipitation from the data set ('monthly_mean_temp.csv',
##     'total_monthly_ppt.csv'), 1925-2006.
## (2) Months missing there are filled from the daily records of the same station
##     ('USC00102707.dly', NOAA GHCN-Daily): a day with only its minimum (or maximum) temperature
##     gets the missing value from the other one, by a regression of daily maximum on daily minimum
##     temperature (and the reverse) fitted to every day of the same calendar month in the station
##     record; the monthly mean is then the mean of (maximum + minimum) / 2. Monthly precipitation
##     is filled only from months with every day recorded. Months of the data set and of the daily
##     records are compared where both exist (same station).
## (3) Seasonal variables of each climate year (July to June, named after the year it ends in;
##     the census is in June, so climate year t+1 runs between the censuses of t and t+1):
##       Pp, Tp   pregrowing season (October-March) precipitation total and mean temperature
##       Pg, Tg   growing season (April-June) precipitation total and mean temperature
##                (seasons of Blaisdell 1958, at this station; October start as the fall-spring
##                window of Adler et al. 2012)
##       Tmax, Tmin  mean temperature of the warmest and of the coldest month (July-June)
##       P_year, T_year  July-June precipitation total and mean temperature
##       P_pg, T_pg      October-June precipitation total and mean temperature
##       Pg_share        Pg / (Pp + Pg);  Tg_Tp  Tg - Tp
## Snowfall is left out (records start in 1930). Uses 'site_dir' and 'site_tag'.


# Settings ---------------------------------------------------------------------
  # Years of the station record used (the monthly records of the data set)
clim_years <- 1925:2006
  # Least number of recorded days for a monthly mean temperature from the daily records
min_days_t <- 25

src <- function(f) here("data", "data.source", site_dir, f)
days_in <- function(y, m) as.integer(format(as.Date(sprintf("%d-%02d-01", y + (m == 12), m %% 12 + 1)) - 1, "%d"))
der_dir <- here("data", "data.derived", site_tag)
tab_dir <- here("result", "tables", site_tag)

# Monthly records of the data set ----------------------------------------------
read_month <- function(f) {
  read.csv(src(f), na.strings = c("NA", ""), check.names = FALSE) %>%
    rename_with(~ trimws(toupper(.x))) %>%
    select(-ANNUAL) %>%
    pivot_longer(-YEAR, names_to = "month_name", values_to = "value") %>%
    transmute(cal_year = YEAR, month = match(month_name, toupper(month.abb)), value)
}
monthly <- expand.grid(month = 1:12, cal_year = clim_years)[, c("cal_year", "month")] %>%
  left_join(read_month("monthly_mean_temp.csv") %>% transmute(cal_year, month, tmean_c = (value - 32) * 5 / 9),
            by = c("cal_year", "month")) %>%
  left_join(read_month("total_monthly_ppt.csv") %>% transmute(cal_year, month, ppt_mm = value * 25.4),
            by = c("cal_year", "month")) %>%
  mutate(source_t = ifelse(is.na(tmean_c), NA, "data set"), source_p = ifelse(is.na(ppt_mm), NA, "data set"))

# Daily records of the station (GHCN-Daily) -------------------------------------
  # Fixed-width file: one line per element and month, 31 days of value (tenths of degC or mm) and
  # three flags; values with a quality flag (failed a GHCN check) are left out
dly <- readLines(src("USC00102707.dly"))
dly <- dly[substr(dly, 18, 21) %in% c("TMAX", "TMIN", "PRCP")]
daily <- bind_rows(lapply(1:31, function(d) {
  pos <- 22 + 8 * (d - 1)
  data.frame(cal_year = as.integer(substr(dly, 12, 15)), month = as.integer(substr(dly, 16, 17)), day = d,
             element = substr(dly, 18, 21), value = as.numeric(substr(dly, pos, pos + 4)),
             qflag = substr(dly, pos + 6, pos + 6))
})) %>%
  filter(value != -9999, qflag == " ", day <= days_in(cal_year, month)) %>%
  mutate(value = value / 10) %>%
  select(-qflag) %>%
  pivot_wider(names_from = element, values_from = value)

  # Fill a missing daily maximum (minimum) from the minimum (maximum) of the same day:
  # linear regression within each calendar month, every day of the station record
both_t <- daily %>% filter(!is.na(TMAX), !is.na(TMIN))
fit_max <- lapply(1:12, function(m) lm(TMAX ~ TMIN, data = both_t[both_t$month == m, ]))
fit_min <- lapply(1:12, function(m) lm(TMIN ~ TMAX, data = both_t[both_t$month == m, ]))
fill_r2 <- data.frame(month = 1:12, r2_max_on_min = sapply(fit_max, function(f) summary(f)$r.squared),
                      r2_min_on_max = sapply(fit_min, function(f) summary(f)$r.squared))
cf_max <- t(sapply(fit_max, coef)); cf_min <- t(sapply(fit_min, coef))
daily <- daily %>%
  mutate(TMAX_f = ifelse(is.na(TMAX), cf_max[month, 1] + cf_max[month, 2] * TMIN, TMAX),
         TMIN_f = ifelse(is.na(TMIN), cf_min[month, 1] + cf_min[month, 2] * TMAX, TMIN))
daily_month <- daily %>%
  group_by(cal_year, month) %>%
  summarise(n_t = sum(!is.na(TMAX_f) & !is.na(TMIN_f)), n_filled = sum(xor(is.na(TMAX), is.na(TMIN))),
            tmean_daily = mean((TMAX_f + TMIN_f) / 2, na.rm = TRUE),
            n_p = sum(!is.na(PRCP)), ppt_daily = sum(PRCP, na.rm = TRUE), .groups = "drop") %>%
  mutate(n_days = days_in(cal_year, month))

  # Agreement of the two sources where both exist (complete months of the daily records)
check <- monthly %>%
  inner_join(daily_month, by = c("cal_year", "month")) %>%
  summarise(n_months_t = sum(!is.na(tmean_c) & n_t == n_days & n_filled == 0),
            mean_diff_t = mean((tmean_daily - tmean_c)[!is.na(tmean_c) & n_t == n_days & n_filled == 0]),
            max_abs_diff_t = max(abs(tmean_daily - tmean_c)[!is.na(tmean_c) & n_t == n_days & n_filled == 0]),
            n_months_p = sum(!is.na(ppt_mm) & n_p == n_days),
            mean_diff_p = mean((ppt_daily - ppt_mm)[!is.na(ppt_mm) & n_p == n_days]),
            cor_p = cor(ppt_daily[!is.na(ppt_mm) & n_p == n_days], ppt_mm[!is.na(ppt_mm) & n_p == n_days]))
cat(sprintf("Data set vs daily records: temperature %d months, mean difference %.2f degC (largest %.2f); precipitation %d months, mean difference %.2f mm (r = %.3f).\n",
            check$n_months_t, check$mean_diff_t, check$max_abs_diff_t, check$n_months_p, check$mean_diff_p, check$cor_p))

# Fill the missing months ------------------------------------------------------
monthly <- monthly %>%
  left_join(daily_month, by = c("cal_year", "month")) %>%
  mutate(fill_t = is.na(tmean_c) & !is.na(n_t) & n_t >= min_days_t,
         fill_p = is.na(ppt_mm) & !is.na(n_p) & n_p == n_days,
         source_t = ifelse(fill_t, sprintf("daily records (%d days, %d with one value estimated)", n_t, n_filled), source_t),
         source_p = ifelse(fill_p, "daily records", source_p),
         tmean_c = ifelse(fill_t, tmean_daily, tmean_c),
         ppt_mm = ifelse(fill_p, ppt_daily, ppt_mm))
  # Long-term mean of each month, to show how far a filled month is from it
lt_mean <- monthly %>% group_by(month) %>% summarise(lt_tmean = mean(tmean_c, na.rm = TRUE), lt_ppt = mean(ppt_mm, na.rm = TRUE))
filled <- monthly %>%
  filter(fill_t | fill_p | is.na(tmean_c) | is.na(ppt_mm)) %>%
  left_join(lt_mean, by = "month") %>%
  transmute(cal_year, month, tmean_c, lt_tmean, source_t, ppt_mm, lt_ppt, source_p)
cat(sprintf("Months filled from the daily records: %d (temperature), %d (precipitation); months still missing: %d, %d.\n",
            sum(monthly$fill_t), sum(monthly$fill_p), sum(is.na(monthly$tmean_c)), sum(is.na(monthly$ppt_mm))))
print(as.data.frame(mutate(filled, across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
monthly <- select(monthly, cal_year, month, tmean_c, ppt_mm, source_t, source_p)

# Seasonal variables of each climate year --------------------------------------
  # Climate year cy: July (cy - 1) to June (cy); pregrowing season October (cy - 1) to March (cy);
  # growing season April to June (cy). A variable needs every month of its window.
clim_year <- monthly %>%
  mutate(cy = ifelse(month >= 7, cal_year + 1L, cal_year),
         season = case_when(month %in% 4:6 ~ "g", month %in% c(10:12, 1:3) ~ "p", TRUE ~ "summer")) %>%
  group_by(cy) %>%
  filter(n() == 12) %>%
  summarise(Pp = sum(ppt_mm[season == "p"]), Pg = sum(ppt_mm[season == "g"]),
            Tp = mean(tmean_c[season == "p"]), Tg = mean(tmean_c[season == "g"]),
            Tmax = max(tmean_c), Tmin = min(tmean_c),
            P_year = sum(ppt_mm), T_year = mean(tmean_c),
            P_pg = sum(ppt_mm[season != "summer"]), T_pg = mean(tmean_c[season != "summer"]),
            .groups = "drop") %>%
  mutate(Pg_share = Pg / P_pg, Tg_Tp = Tg - Tp) %>%
  filter(cy %in% (min(clim_years) + 1):max(clim_years))   # complete July-June years
cat(sprintf("Climate years with every variable: %d of %d (%d to %d).\n",
            sum(complete.cases(clim_year)), nrow(clim_year), min(clim_year$cy), max(clim_year$cy)))

# Save -------------------------------------------------------------------------
write.csv(monthly, file.path(der_dir, "CLR_Climate_monthly.csv"), row.names = FALSE)
write.csv(clim_year, file.path(der_dir, "CLR_Climate_year.csv"), row.names = FALSE)
write.csv(filled, file.path(tab_dir, "CLR_Climate_filled.csv"), row.names = FALSE)
write.csv(cbind(check, n_fill_models_month = 12, mean_r2_max_on_min = mean(fill_r2$r2_max_on_min)),
          file.path(tab_dir, "CLR_Climate_check.csv"), row.names = FALSE)
cat("Saved:\n",
    " data/data.derived/", site_tag, ": CLR_Climate_monthly.csv (complete monthly record), CLR_Climate_year.csv (seasonal variables)\n",
    " result/tables/", site_tag, ": CLR_Climate_filled.csv (filled months), CLR_Climate_check.csv (data set vs daily records)\n", sep = "")
rm(list = intersect(c("read_month", "src", "days_in", "dly", "daily", "both_t", "fit_max", "fit_min", "cf_max", "cf_min", "daily_month", "check",
                      "lt_mean", "filled", "fill_r2", "min_days_t"), ls()))
