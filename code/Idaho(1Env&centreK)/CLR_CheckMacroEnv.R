################################################################################
## Check the raw monthly macro-environment data of the Idaho site
## (US Sheep Experiment Station, Dubois; station records shipped with the data)
##
## Written by Hung-wei Lin Oct 2026
##
## Reads the three monthly climate files (mean temperature, total precipitation,
## total snowfall), puts them in one long table and checks each variable for
## coverage, broad range, monthly outliers and the ANNUAL column. Raw values stay
## as they are (flags only); cleaning decisions belong to the next step.
## Uses 'site_dir' and 'site_tag' from the runscript.

# Read the raw csv files -------------------------------------------------------
  # Units come from the data paper metadata (Zachmann et al. 2010):
  # temperature in degF, precipitation and snowfall in inches.
var_info <- data.frame(
  code  = c("TMEAN", "PPT", "SNOW"),
  file  = c("monthly_mean_temp.csv", "total_monthly_ppt.csv", "total_monthly_sno.csv"),
  name  = c("monthly mean temperature", "monthly total precipitation", "monthly total snowfall"),
  unit  = c("degF", "in", "in"),
  # lower/upper = broad plausibility screens in raw units (sensor or typing errors only)
  lower = c(-30, 0, 0),
  upper = c(100, 15, 80),
  # ANNUAL should equal the mean (temperature) or the sum (precipitation, snowfall) of the months
  annual_fun = c("mean", "sum", "sum"),
  annual_tol = c(0.5, 0.1, 0.5),
  stringsAsFactors = FALSE
)
month_names <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                 "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

climate_wide <- list()
for (i in seq_len(nrow(var_info))) {
  d <- read.csv(here("data", "data.source", site_dir, var_info$file[i]),
                na.strings = c("NA", ""), check.names = FALSE)
  names(d) <- trimws(toupper(names(d)))
  # Guard the column layout before using it (a changed layout would shift months silently)
  stopifnot(identical(names(d), c("YEAR", month_names, "ANNUAL")))
  d$code <- var_info$code[i]
  climate_wide[[var_info$code[i]]] <- d
}
climate_wide <- bind_rows(climate_wide)
dim(climate_wide) # one row per variable and year

  # Long version: one row per variable, year and month
climate_raw <- climate_wide %>%
  select(-ANNUAL) %>%
  pivot_longer(all_of(month_names), names_to = "month_name", values_to = "value") %>%
  mutate(month = match(month_name, month_names),
         code  = factor(code, levels = var_info$code)) %>%
  rename(Year = YEAR) %>%
  arrange(code, Year, month)

# Check the year sequence ------------------------------------------------------
  # One row per year and variable, no duplicates; years missing inside the record are listed.
for (v in var_info$code) {
  yrs  <- climate_wide$YEAR[climate_wide$code == v]
  gaps <- setdiff(seq(min(yrs), max(yrs)), yrs)
  cat(sprintf("%-5s %d to %d, %d years, %d duplicated, missing years: %s\n",
              v, min(yrs), max(yrs), length(yrs), sum(duplicated(yrs)),
              ifelse(length(gaps) == 0, "none", paste(gaps, collapse = ", "))))
}

# Single-variable checks -------------------------------------------------------
  # For each variable:
  # (1) coverage = months with a value within the record span,
  # (2) broad range,
  # (3) monthly outliers = values far from the usual value of the same calendar month
  #     (temperature: more than 4 robust SD from the monthly median, both sides;
  #      precipitation and snowfall are skewed with many zeros, so only the high side
  #      is screened, above Q3 + 3 IQR of the non-zero values of the same calendar month;
  #      months with fewer than 10 non-zero values, e.g. summer snowfall, skip this check),
  # (4) the ANNUAL column against the mean or sum of the twelve months.
  # Every flagged month goes into flag_list with the check name and a note.
flag_list <- list()
summary_list <- list()

for (i in seq_len(nrow(var_info))) {
  v  <- var_info$code[i]
  dv <- filter(climate_raw, code == v)

  # (1) Coverage
  n_missing <- sum(is.na(dv$value))

  # (2) Range
  out_range <- filter(dv, value < var_info$lower[i] | value > var_info$upper[i])
  if (nrow(out_range) > 0) {
    flag_list[[paste(v, "range")]] <- transmute(out_range, code = v, Year, month, value,
      check = "range",
      note = sprintf("outside %g to %g %s", var_info$lower[i], var_info$upper[i], var_info$unit[i]))
  }

  # (3) Monthly outliers
  if (v == "TMEAN") {
    out_month <- dv %>%
      group_by(month) %>%
      mutate(centre = median(value, na.rm = TRUE), spread = mad(value, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(abs(value - centre) > 4 * spread) %>%
      mutate(note = sprintf("%.1f degF from the monthly median (4 robust SD = %.1f)",
                            value - centre, 4 * spread))
  } else {
    out_month <- dv %>%
      group_by(month) %>%
      mutate(n_wet = sum(value > 0, na.rm = TRUE),
             cut = quantile(value[value > 0], 0.75, na.rm = TRUE) +
               3 * IQR(value[value > 0], na.rm = TRUE)) %>%
      ungroup() %>%
      filter(n_wet >= 10, value > cut) %>%
      mutate(note = sprintf("above the monthly Q3 + 3 IQR (%.2f %s)", cut, var_info$unit[i]))
  }
  if (nrow(out_month) > 0) {
    flag_list[[paste(v, "month")]] <- transmute(out_month, code = v, Year, month, value,
                                                check = "monthly_outlier", note)
  }

  # (4) ANNUAL column (checked only in years with all twelve months present)
  annual_check <- climate_wide %>%
    filter(code == v) %>%
    mutate(n_month = rowSums(!is.na(across(all_of(month_names)))),
           from_months = if (var_info$annual_fun[i] == "mean") {
             rowMeans(across(all_of(month_names)))
           } else {
             rowSums(across(all_of(month_names)))
           }) %>%
    filter(!is.na(ANNUAL), n_month == 12,
           abs(ANNUAL - from_months) > var_info$annual_tol[i])
  if (nrow(annual_check) > 0) {
    flag_list[[paste(v, "annual")]] <- transmute(annual_check, code = v, Year = YEAR,
      month = NA_integer_, value = ANNUAL, check = "annual_mismatch",
      note = sprintf("ANNUAL %.2f vs %s of months %.2f", ANNUAL, var_info$annual_fun[i], from_months))
  }

  summary_list[[v]] <- data.frame(
    code = v, name = var_info$name[i], unit = var_info$unit[i],
    first_year = min(dv$Year), last_year = max(dv$Year),
    n_values = sum(!is.na(dv$value)), n_missing_months = n_missing,
    pct_missing = round(100 * n_missing / nrow(dv), 1),
    min = min(dv$value, na.rm = TRUE), median = median(dv$value, na.rm = TRUE),
    max = max(dv$value, na.rm = TRUE),
    n_out_range = nrow(out_range), n_monthly_outlier = nrow(out_month),
    n_annual_filled = sum(!is.na(climate_wide$ANNUAL[climate_wide$code == v])),
    n_annual_mismatch = nrow(annual_check))
}
macroenv_summary <- bind_rows(summary_list)

# Cross-variable consistency check ---------------------------------------------
  # Snowfall in a month with a warm monthly mean is possible (one cold spell is enough),
  # so this is a soft flag: snowfall > 1 in while the monthly mean is above 60 degF.
snow_warm <- climate_raw %>%
  select(code, Year, month, value) %>%
  pivot_wider(names_from = code, values_from = value) %>%
  filter(SNOW > 1, TMEAN > 60)
if (nrow(snow_warm) > 0) {
  flag_list[["SNOW warm"]] <- transmute(snow_warm, code = "SNOW", Year, month, value = SNOW,
    check = "snow_in_warm_month",
    note = sprintf("snowfall %.1f in with monthly mean %.1f degF", SNOW, TMEAN))
}

macroenv_flags <- bind_rows(flag_list) %>% arrange(code, Year, month)
cat(sprintf("Flagged %d variable-months across %d checks.\n",
            nrow(macroenv_flags), length(unique(macroenv_flags$check))))
table(macroenv_flags$code, macroenv_flags$check)

# Yearly coverage --------------------------------------------------------------
  # Months with a value (0 to 12), per variable and year; years absent from a file count as 0
macroenv_coverage <- climate_raw %>%
  group_by(code, Year) %>%
  summarise(n_months = sum(!is.na(value)), .groups = "drop") %>%
  complete(code, Year = seq(min(Year), max(Year)), fill = list(n_months = 0))

  # Census years of the mapped quadrats (from the inventory), to see which censuses
  # have climate for the preceding twelve months
census_years <- read.csv(here("data", "data.source", site_dir, "quad_inventory.csv"))$year + 1900

# Save tables ------------------------------------------------------------------
tab_dir <- here("result", "tables", site_tag)
fig_dir <- here("result", "figs", site_tag)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(macroenv_summary, file.path(tab_dir, "CLR_CheckMacroEnv_summary.csv"), row.names = FALSE)
write.csv(macroenv_flags, file.path(tab_dir, "CLR_CheckMacroEnv_flags.csv"), row.names = FALSE)
write.csv(macroenv_coverage, file.path(tab_dir, "CLR_CheckMacroEnv_coverage.csv"), row.names = FALSE)

# Figures ----------------------------------------------------------------------
  # Coverage heat map; census years marked with ticks along the top
p_cov <- ggplot(macroenv_coverage, aes(x = Year, y = code, fill = n_months)) +
  geom_tile() +
  geom_point(data = data.frame(Year = census_years), aes(x = Year, y = 3.65),
             inherit.aes = FALSE, shape = 25, size = 1.6, fill = "black") +
  scale_fill_gradient(low = "grey95", high = "steelblue4", limits = c(0, 12),
                      breaks = c(0, 6, 12), name = "Months\nwith value") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL,
       title = paste0(site_tag, " station climate: yearly coverage (triangles = census years)")) +
  theme_bw(base_size = 9)
ggsave(file.path(fig_dir, "CLR_CheckMacroEnv_coverage.png"), p_cov, width = 10, height = 3, dpi = 300)

  # Monthly series, one panel per variable; flagged months in red,
  # census period shaded (one page, so a png)
series_data <- climate_raw %>%
  filter(!is.na(value)) %>%
  mutate(date = Year + (month - 0.5) / 12)
flag_points <- macroenv_flags %>%
  filter(!is.na(month)) %>%
  mutate(code = factor(code, levels = var_info$code),
         date = Year + (month - 0.5) / 12)
p_series <- ggplot(series_data, aes(x = date, y = value)) +
  annotate("rect", xmin = min(census_years), xmax = 1958, ymin = -Inf, ymax = Inf,
           fill = "khaki", alpha = 0.35) +
  geom_line(linewidth = 0.2, colour = "grey30") +
  geom_point(data = flag_points, colour = "red", size = 0.9) +
  facet_wrap(~ code, scales = "free_y", ncol = 1,
             labeller = as_labeller(setNames(paste0(var_info$name, " (", var_info$unit, ")"),
                                             var_info$code))) +
  labs(x = NULL, y = NULL,
       title = paste0(site_tag, " station climate (shaded = main census period; red = flagged)")) +
  theme_bw(base_size = 9)
ggsave(file.path(fig_dir, "CLR_CheckMacroEnv_series.png"), p_series, width = 11, height = 7, dpi = 300)

# Report saved files -----------------------------------------------------------
saved_files <- data.frame(
  folder = c(rep(file.path("result/tables", site_tag), 3), rep(file.path("result/figs", site_tag), 2)),
  file   = c("CLR_CheckMacroEnv_summary.csv", "CLR_CheckMacroEnv_flags.csv",
             "CLR_CheckMacroEnv_coverage.csv", "CLR_CheckMacroEnv_coverage.png",
             "CLR_CheckMacroEnv_series.png"),
  content = c("one row per variable: record span, gaps, min/median/max, counts per check",
              "one row per flagged variable-month: year, month, raw value, check, note",
              "months with a value, per variable and year",
              "heat map of the coverage table, census years marked",
              "monthly series per variable, flagged months in red")
)
stopifnot(all(file.exists(here(saved_files$folder, saved_files$file))))
cat("\nCLR_CheckMacroEnv.R (", site_tag, ") saved these files:\n", sep = "")
for (fd in unique(saved_files$folder)) {
  cat(sprintf("  %s/\n", fd))
  sf <- saved_files[saved_files$folder == fd, ]
  cat(sprintf("    %-34s %s\n", sf$file, sf$content), sep = "")
}
cat("Workspace objects kept: climate_raw, climate_wide, var_info, macroenv_summary, macroenv_flags, macroenv_coverage\n")

  # Tidy up temporary objects (keep the raw tables, var_info and the three result tables)
rm(list = intersect(c("i", "v", "d", "dv", "yrs", "gaps", "n_missing", "out_range", "out_month",
                      "annual_check", "snow_warm", "summary_list", "flag_list", "month_names",
                      "census_years", "tab_dir", "fig_dir", "p_cov", "p_series", "series_data",
                      "flag_points", "saved_files", "fd", "sf"), ls()))
