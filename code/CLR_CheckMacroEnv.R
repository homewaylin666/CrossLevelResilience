################################################################################
## Check the raw daily macro-environment data (Schofield Pass SNOTEL 737)
## 
## Written by Hung-wei Lin Oct 2026
##
## Reads the raw NRCS csv, gives every column a short code, and checks every
## numeric variable for coverage, broad range, QC flags, spikes, flat runs and
## cross-variable consistency. Raw values stay as they are (flags only);
## cleaning decisions belong to the next step.

# Read the raw csv -------------------------------------------------------------
macroenv_file <- here("data", "data.source", "SNOTEL737, macroenv, 1985-2026", "sp737_all.csv")

  # Lines starting with "#" are NRCS header notes.
  # Every element comes as a value column followed by its QC flag column.
n_elem <- 14
macroenv_raw <- read.csv(macroenv_file, comment.char = "#", check.names = FALSE,
                         colClasses = c("character", rep(c("numeric", "character"), n_elem)))
dim(macroenv_raw) # one row per day, 1 + 2 * 14 columns

# Variable look-up table -------------------------------------------------------
  # Order follows the download order in README_macroenvironment.md.
  # lower/upper = broad plausibility screens in raw units (in, degF, pct),
  # meant to catch sensor failures only (tune before using them for cleaning).
  # flat_run = number of identical consecutive values that counts as a flat run;
  # flat_skip = a plateau value that is physical rather than a frozen sensor
  # (soil under snow sits at the freezing point, 32 degF, for months).
  # Soil temperature comes in whole degF and deep soil changes slowly, hence its longer run.
soil_codes <- paste0(rep(c("SMS", "STO"), each = 3), "_", c(2, 8, 20))
var_info <- data.frame(
  code    = c("PREC", "PRCP", "PRCPSA", "TMAX", "TMIN", "TAVG", "WTEQ", "SNWD", soil_codes),
  keyword = c("Precipitation Accumulation", "Precipitation Increment (", "Snow-adj",
              "Air Temperature Maximum", "Air Temperature Minimum", "Air Temperature Average",
              "Snow Water Equivalent", "Snow Depth",
              rep("Soil Moisture Percent", 3), rep("Soil Temperature Observed", 3)),
  group   = c(rep("precipitation", 3), rep("air temperature", 3), rep("snow", 2),
              rep("soil moisture", 3), rep("soil temperature", 3)),
  unit    = c("in", "in", "in", rep("degF", 3), "in", "in", rep("pct", 3), rep("degF", 3)),
  lower   = c(0, 0, 0, rep(-50, 3), 0, 0, rep(0, 3), rep(-20, 3)),
  upper   = c(100, 5, 5, rep(100, 3), 80, 250, rep(60, 3), rep(100, 3)),
  flat_run  = c(NA, NA, NA, rep(5, 3), NA, NA, rep(30, 3), rep(60, 3)),
  flat_skip = c(rep(NA, 11), rep(32, 3)),
  stringsAsFactors = FALSE
)
var_info$depth_in <- ifelse(grepl("_", var_info$code), sub(".*_", "", var_info$code), NA)
stopifnot(nrow(var_info) == n_elem)

# Rename columns to short codes ------------------------------------------------
  # Check that each column really is what we think it is before renaming
  # (a changed download order would otherwise shift every name silently).
hdr <- sub("^Schofield Pass \\(737\\) ", "", names(macroenv_raw))
stopifnot(ncol(macroenv_raw) == 1 + 2 * n_elem)
stopifnot(hdr[1] == "Date")
val_hdr <- hdr[seq(2, ncol(macroenv_raw), by = 2)]
qc_hdr  <- hdr[seq(3, ncol(macroenv_raw), by = 2)]
stopifnot(all(grepl("QC Flag", qc_hdr)))
stopifnot(all(!grepl("QC Flag", val_hdr)))
stopifnot(all(mapply(function(k, h) grepl(k, h, fixed = TRUE), var_info$keyword, val_hdr)))
  # Soil columns also need the right depth
soil_rows <- !is.na(var_info$depth_in)
stopifnot(all(mapply(function(d, h) grepl(paste0("-", d, "in"), h, fixed = TRUE),
                     var_info$depth_in[soil_rows], val_hdr[soil_rows])))

names(macroenv_raw) <- c("Date", as.vector(rbind(var_info$code, paste0(var_info$code, "_qc"))))
macroenv_raw$Date <- as.Date(macroenv_raw$Date)
macroenv_raw$Year <- as.integer(format(macroenv_raw$Date, "%Y"))
  # Water year starts on 1 Oct (PREC resets then)
macroenv_raw$WaterYear <- macroenv_raw$Year + (as.integer(format(macroenv_raw$Date, "%m")) >= 10)

# Check the date sequence ------------------------------------------------------
  # One row per calendar day, no duplicates, no skipped days.
n_dup  <- sum(duplicated(macroenv_raw$Date))
n_skip <- length(seq(min(macroenv_raw$Date), max(macroenv_raw$Date), by = "day")) -
  length(unique(macroenv_raw$Date))
cat(sprintf("Dates: %s to %s, %d rows, %d duplicated, %d skipped.\n",
            min(macroenv_raw$Date), max(macroenv_raw$Date), nrow(macroenv_raw), n_dup, n_skip))

# Single-variable checks -------------------------------------------------------
  # For each numeric variable (QC flag columns are tallied, never checked as numbers):
  # (1) coverage within its own record, (2) broad range, (3) NRCS QC flags,
  # (4) spikes = day-to-day jumps far beyond the variable's usual jumps,
  # (5) flat runs = the same value repeated for many days (frozen sensor).
  # Every flagged day goes into flag_list with the check name and a note.
flag_list <- list()
summary_list <- list()

for (i in seq_len(nrow(var_info))) {
  v   <- var_info$code[i]
  x   <- macroenv_raw[[v]]
  qc  <- macroenv_raw[[paste0(v, "_qc")]]
  has <- which(!is.na(x))
  
  # (1) Coverage: count gaps between the first and last recorded day only
  span      <- has[1]:has[length(has)]
  n_missing <- sum(is.na(x[span]))
  
  # (2) Range
  out_range <- which(x < var_info$lower[i] | x > var_info$upper[i])
  if (length(out_range) > 0) {
    flag_list[[paste(v, "range")]] <- data.frame(
      Date = macroenv_raw$Date[out_range], code = v, value = x[out_range], check = "range",
      note = sprintf("outside %g to %g %s", var_info$lower[i], var_info$upper[i], var_info$unit[i]))
  }
  
  # (3) QC flags: tally every letter; "S" (suspect) days also go to the flag list
  qc_tab <- table(qc[has], useNA = "no")
  qc_txt <- paste(names(qc_tab), qc_tab, sep = "=", collapse = "; ")
  suspect <- which(qc %in% "S" & !is.na(x))
  if (length(suspect) > 0) {
    flag_list[[paste(v, "qc")]] <- data.frame(
      Date = macroenv_raw$Date[suspect], code = v, value = x[suspect], check = "qc_suspect",
      note = "NRCS QC flag S")
  }
  
  # (4) Spikes: jump larger than 8 robust SDs of the day-to-day changes
  # (robust SD from non-zero changes only, since snow and soil sit still for
  # months; precipitation increments and the precipitation running total jump on
  # every storm by nature, so they skip this check)
  n_spike <- NA
  if (!v %in% c("PREC", "PRCP", "PRCPSA")) {
    dx <- c(NA, diff(x))
    dx_move <- dx[!is.na(dx) & dx != 0]
    robust_sd <- if (length(dx_move) > 10) mad(dx_move) else NA
    if (is.finite(robust_sd) && robust_sd > 0) {
      spike <- which(abs(dx) > 8 * robust_sd)
      n_spike <- length(spike)
      if (n_spike > 0) {
        flag_list[[paste(v, "spike")]] <- data.frame(
          Date = macroenv_raw$Date[spike], code = v, value = x[spike], check = "spike",
          note = sprintf("day-to-day change %.2f (8 robust SD = %.2f)", dx[spike], 8 * robust_sd))
      }
    }
  }
  
  # (5) Flat runs: identical values repeated >= flat_run days
  # (zeros are skipped: zero snow or zero rain for weeks is normal;
  # flat_skip +- 1 is skipped too, see var_info)
  n_flat <- NA
  if (!is.na(var_info$flat_run[i])) {
    r       <- rle(ifelse(is.na(x), NA_real_, x))
    r_end   <- cumsum(r$lengths)
    r_beg   <- r_end - r$lengths + 1
    plateau <- !is.na(var_info$flat_skip[i]) & abs(r$values - var_info$flat_skip[i]) <= 1
    long    <- which(r$lengths >= var_info$flat_run[i] & !is.na(r$values) & r$values != 0 &
                       !plateau %in% TRUE)
    n_flat  <- length(long)
    if (n_flat > 0) {
      days <- unlist(mapply(seq, r_beg[long], r_end[long], SIMPLIFY = FALSE))
      flag_list[[paste(v, "flat")]] <- data.frame(
        Date = macroenv_raw$Date[days], code = v, value = x[days], check = "flat_run",
        note = sprintf("same value for >= %d days", var_info$flat_run[i]))
    }
  }
  
  summary_list[[v]] <- data.frame(
    code = v, group = var_info$group[i], unit = var_info$unit[i],
    first_date = macroenv_raw$Date[has[1]], last_date = macroenv_raw$Date[has[length(has)]],
    n_values = length(has), n_missing_in_span = n_missing,
    pct_missing_in_span = round(100 * n_missing / length(span), 1),
    min = min(x, na.rm = TRUE), median = median(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
    n_out_range = length(out_range), n_spike = n_spike, n_flat_run = n_flat,
    qc_flags = qc_txt)
}
macroenv_summary <- bind_rows(summary_list)

# Cross-variable consistency checks --------------------------------------------
add_flag <- function(rows, code, value, check, note) {
  if (length(rows) == 0) return(NULL)
  data.frame(Date = macroenv_raw$Date[rows], code = code, value = value[rows],
             check = check, note = note)
}

  # (a) Daily air temperature: minimum <= average <= maximum
bad <- with(macroenv_raw, which(TMIN > TAVG | TAVG > TMAX | TMIN > TMAX))
flag_list[["TAVG order"]] <- add_flag(bad, "TAVG", macroenv_raw$TAVG, "min_avg_max_order",
                                      "TMIN <= TAVG <= TMAX broken")

  # (b) PREC is a running total within each water year, so it should only go up.
  # A drop of more than 0.1 in (outside the 1 Oct reset) is suspicious.
d_prec <- c(NA, diff(macroenv_raw$PREC))
d_prec[c(FALSE, diff(macroenv_raw$WaterYear) != 0)] <- NA
bad <- which(d_prec < -0.1)
flag_list[["PREC drop"]] <- add_flag(bad, "PREC", macroenv_raw$PREC, "prec_drop",
                                     "accumulation drops within a water year")

  # (c) PRCP of day t should match PREC(t+1) - PREC(t) (PREC is start-of-day)
d_next <- c(diff(macroenv_raw$PREC), NA)
d_next[c(diff(macroenv_raw$WaterYear) != 0, FALSE)] <- NA
bad <- which(abs(macroenv_raw$PRCP - d_next) > 0.05)
flag_list[["PRCP mismatch"]] <- add_flag(bad, "PRCP", macroenv_raw$PRCP, "prcp_vs_prec",
                                         "increment differs from change in accumulation by > 0.05 in")

  # (d) Snow-adjusted increment should stay at or above the raw increment
bad <- with(macroenv_raw, which(PRCPSA < PRCP - 0.05))
flag_list[["PRCPSA below"]] <- add_flag(bad, "PRCPSA", macroenv_raw$PRCPSA, "prcpsa_vs_prcp",
                                        "snow-adjusted increment below raw increment")

  # (e) Snow depth and snow water equivalent should agree on snow / no snow
bad <- with(macroenv_raw, which((SNWD == 0 & WTEQ > 1) | (SNWD > 12 & WTEQ == 0)))
flag_list[["SNWD WTEQ"]] <- add_flag(bad, "SNWD", macroenv_raw$SNWD, "snwd_vs_wteq",
                                     "snow depth and SWE disagree on snow presence")

macroenv_flags <- bind_rows(flag_list) %>% arrange(code, Date)
cat(sprintf("Flagged %d variable-days across %d checks.\n",
            nrow(macroenv_flags), length(unique(macroenv_flags$check))))
table(macroenv_flags$code, macroenv_flags$check)

# Yearly coverage --------------------------------------------------------------
  # Percent of days with a value, per variable and calendar year
macroenv_coverage <- macroenv_raw %>%
  select(Year, all_of(var_info$code)) %>%
  pivot_longer(-Year, names_to = "code", values_to = "value") %>%
  group_by(code, Year) %>%
  summarise(n_days = n(), n_values = sum(!is.na(value)),
            pct_present = round(100 * n_values / n_days, 1), .groups = "drop") %>%
  mutate(code = factor(code, levels = var_info$code))

# Save tables ------------------------------------------------------------------
dir.create(here("result", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(macroenv_summary, here("result", "tables", "CLR_CheckMacroEnv_summary.csv"), row.names = FALSE)
write.csv(macroenv_flags, here("result", "tables", "CLR_CheckMacroEnv_flags.csv"), row.names = FALSE)
write.csv(macroenv_coverage, here("result", "tables", "CLR_CheckMacroEnv_coverage.csv"), row.names = FALSE)

# Figures ----------------------------------------------------------------------
dir.create(here("result", "figs"), recursive = TRUE, showWarnings = FALSE)

  # Coverage heat map: which variable has data in which year
p_cov <- ggplot(macroenv_coverage, aes(x = Year, y = code, fill = pct_present)) +
  geom_tile() +
  scale_fill_gradient(low = "grey95", high = "steelblue4", limits = c(0, 100),
                      name = "% days\nwith value") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, title = "Schofield Pass SNOTEL 737: yearly coverage") +
  theme_bw(base_size = 9)
ggsave(here("result", "figs", "CLR_CheckMacroEnv_coverage.png"), p_cov, width = 10, height = 5, dpi = 300)

  # Daily series, one page per variable group; flagged days in red
  # (several pages, so this one stays a pdf)
macroenv_long <- macroenv_raw %>%
  select(Date, all_of(var_info$code)) %>%
  pivot_longer(-Date, names_to = "code", values_to = "value") %>%
  filter(!is.na(value)) %>%
  left_join(var_info[, c("code", "group")], by = "code") %>%
  mutate(code = factor(code, levels = var_info$code))
flag_points <- macroenv_flags %>%
  left_join(var_info[, c("code", "group")], by = "code") %>%
  mutate(code = factor(code, levels = var_info$code))

pdf(here("result", "figs", "CLR_CheckMacroEnv_series.pdf"), width = 11, height = 8)
for (g in unique(var_info$group)) {
  fp <- filter(flag_points, group == g)
  p <- ggplot(filter(macroenv_long, group == g), aes(x = Date, y = value)) +
    geom_line(linewidth = 0.2, colour = "grey30") +
    {if (nrow(fp) > 0) geom_point(data = fp, colour = "red", size = 0.6)} +
    facet_wrap(~ code, scales = "free_y", ncol = 1) +
    labs(x = NULL, y = NULL, title = paste0("Schofield Pass SNOTEL 737: ", g,
                                           " (red = flagged by any check)")) +
    theme_bw(base_size = 9)
  print(p)
}
dev.off()

# Report saved files -----------------------------------------------------------
saved_files <- data.frame(
  folder = c(rep("result/tables", 3), rep("result/figs", 2)),
  file   = c("CLR_CheckMacroEnv_summary.csv", "CLR_CheckMacroEnv_flags.csv",
             "CLR_CheckMacroEnv_coverage.csv", "CLR_CheckMacroEnv_coverage.png",
             "CLR_CheckMacroEnv_series.pdf"),
  content = c("one row per variable: record span, gaps, min/median/max, counts per check, QC flags",
              "one row per flagged variable-day: date, variable, raw value, check, note",
              "percent of days with a value, per variable and year",
              "heat map of the coverage table",
              "daily series per variable group, flagged days in red")
)
stopifnot(all(file.exists(here(saved_files$folder, saved_files$file))))
cat("\nCLR_CheckMacroEnv.R saved these files:\n")
for (fd in unique(saved_files$folder)) {
  cat(sprintf("  %s/\n", fd))
  sf <- saved_files[saved_files$folder == fd, ]
  cat(sprintf("    %-34s %s\n", sf$file, sf$content), sep = "")
}
cat("Workspace objects kept: macroenv_raw, var_info, macroenv_summary, macroenv_flags, macroenv_coverage\n")

  # Tidy up temporary objects (keep macroenv_raw, var_info and the three result tables)
rm(list = intersect(c("i", "v", "x", "qc", "has", "span", "n_missing", "out_range", "qc_tab",
                      "qc_txt", "suspect", "n_spike", "n_flat", "r", "r_end", "r_beg", "plateau",
                      "long", "dx", "dx_move", "robust_sd", "spike", "days", "bad", "d_prec",
                      "d_next", "hdr", "val_hdr", "qc_hdr", "soil_rows", "soil_codes",
                      "summary_list", "flag_list", "n_dup", "n_skip", "n_elem", "g", "p", "fp",
                      "p_cov", "macroenv_long", "flag_points", "add_flag", "macroenv_file",
                      "saved_files", "fd", "sf"), ls()))
