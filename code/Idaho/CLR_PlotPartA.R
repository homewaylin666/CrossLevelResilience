################################################################################
## Visualisation of the Part A model results (Idaho): forest plots and convergence
##
## Written by Hung-wei Lin Oct 2026
##
## Reads the posterior samples of the main model and of the two submodels, and draws forest
## plots of the estimates: posterior median (point), 89% interval (thick bar) and 95% interval
## (thin bar). One panel per vital rate. Colours: individual level (log size, grazing) orange and
## precipitation blue, as in the model diagram (C3_model_diagram_Adler); each crowding class has
## its own colour (same species, other grasses, other shrubs, forbs), used for its main effect and
## for every interaction with it. The rows come in blocks: main effects (filled points), then the
## interactions with precipitation (delta), log size (theta) and grazing (omega) (open points);
## within an interaction block the colour gives the term that the block's covariate modifies.
## Variance components: year effects blue, group, quadrat and plant effects orange, variation
## among species red.
## Every estimate is per SD of its covariate, on the scale of each vital rate: logit (survival),
## log size at t+1 (growth) and log recruits per unit of effective cover (recruitment). Covariates:
## log size (log area); precipitation anomaly in SD of the long-term record (+ wetter, - drier
## than the long-term mean); crowding by other grasses, other shrubs and forbs per SD over all
## focal plants, same-species crowding per SD within each species (CLR_InputPartA.R); recruitment
## uses the square root of quadrat cover in place of size and crowding.
## A convergence table (Rhat, effective sample size) comes with the figures.
## Uses 'site_tag'.
################################################################################

# Settings ---------------------------------------------------------------------
  # TRUE: plot the test runs ('_test' files) to check this script
plot_test <- FALSE
  # Stage of the main model to plot (1 = no trait)
plot_stage <- 1

der_dir <- here("data", "data.derived", site_tag)
fig_dir <- here("result", "figs", site_tag)
tab_dir <- here("result", "tables", site_tag)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
tag <- if (plot_test) "_test" else ""
fit_file <- function(name) file.path(der_dir, paste0(name, tag, ".rds"))

# Colours: blocks of the model diagram, and one colour per crowding class
col_ind   <- "#b0641a"   # individual level: size, grazing; group, quadrat and plant effects
col_clim  <- "#1560b8"   # community level: precipitation, year effects
col_pop   <- "#b03a3a"   # population level: variation among species
col_W     <- c("#1e6b2a",  # crowding by the same species (green of the diagram)
               "#a08a00",  # crowding by other grasses
               "#6a3d9a",  # crowding by other shrubs
               "#c2457f")  # crowding by forbs
vr_names  <- c("Survival", "Growth", "Recruitment")
  # Crowding classes (W columns: same species, other grasses, other shrubs, forbs). In the main
  # model every forb is another species (forbs form their own pooled population)
crowd_lab <- c("Crowding by same species", "Crowding by other grasses", "Crowding by other shrubs",
               "Crowding by forbs")
crowd_by  <- tolower(crowd_lab)
names(col_W) <- crowd_lab
size_lab  <- "Log size"
clim_lab  <- "Precip. anomaly"
cov_note  <- paste("Estimates per SD of each covariate. Log size = log area; precip. anomaly: + wetter, \u2212 drier",
                   "than the long-term mean; same-species crowding scaled within each species.")

# ====================================================================================
# Helper functions
# ====================================================================================

# --- Posterior draws of the parameters matching a name, with their indices ----------
#   Returns one row per parameter: name, index columns i1, i2, ..., and the summary
#   (median, 95% and 89% intervals)
post_summary <- function(draws, stem) {
  cols <- grep(paste0("^", gsub(".", "\\.", stem, fixed = TRUE), "(\\[|$)"), colnames(draws), value = TRUE)
  stopifnot(length(cols) > 0)
  idx <- lapply(strsplit(sub("^[^[]*\\[?([^]]*)\\]?$", "\\1", cols), ",\\s*"), as.integer)
  n_idx <- max(lengths(idx))
  out <- data.frame(name = cols)
  for (j in seq_len(n_idx)) out[[paste0("i", j)]] <- sapply(idx, `[`, j)
  q <- apply(draws[, cols, drop = FALSE], 2, quantile, probs = c(0.5, 0.025, 0.975, 0.055, 0.945))
  out$med  <- q[1, ]; out$lo95 <- q[2, ]; out$hi95 <- q[3, ]; out$lo89 <- q[4, ]; out$hi89 <- q[5, ]
  out
}

# --- draw_forest_panel: one forest plot ----------------------------------------------
#   dat: one row per estimate, with columns row (label) and the summary.
#   Rows run from the top in the order of 'row_levels'; 'row_cols' gives the colour of each
#   row (named by row). 'show_labels = FALSE' leaves the row labels out (panels to the right
#   of the first one). 'sp_dat' (optional): species medians (columns row, med), drawn as small
#   faded points in the colour of their row, spread vertically in species order.
#   'headers' (optional): rows of 'row_levels' that head a block of rows; drawn as bold labels
#   with a darker line above them, and no estimate. 'open_rows' (optional): rows drawn with
#   open points (interaction terms).
draw_forest_panel <- function(dat, row_levels, row_cols, main, xlab = "", panel_label = NULL,
                              xlim = NULL, show_labels = TRUE, sp_dat = NULL, headers = character(0),
                              open_rows = character(0)) {
  n_r <- length(row_levels)
  centre <- rev(seq(1, n_r * 3, by = 3))
  if (is.null(xlim)) {
    xlim <- range(c(dat$lo95, dat$hi95, sp_dat$med, 0), na.rm = TRUE)
    xlim <- xlim + diff(xlim) * c(-0.08, 0.08)
  }
  plot(NA, xlim = xlim, ylim = c(-0.5, max(centre) + 1.5), xlab = xlab, ylab = "", main = main,
       yaxt = "n", font.main = 1)
  abline(v = 0, lty = 2, col = "grey50")
  is_head <- row_levels %in% headers
  if (n_r > 1) abline(h = centre[-n_r] - 1.5, col = ifelse(is_head[-1], "grey55", "grey90"))
  if (show_labels) {
    axis(2, at = centre[!is_head], labels = row_levels[!is_head], las = 2, tick = FALSE)
    if (any(is_head)) axis(2, at = centre[is_head], labels = row_levels[is_head], las = 2, tick = FALSE, font = 2)
  }
  for (i in seq_len(n_r)) {
    if (is_head[i]) next
    col_i <- row_cols[[row_levels[i]]]
    if (!is.null(sp_dat)) {
      d_sp <- sp_dat[sp_dat$row == row_levels[i], ]
      if (nrow(d_sp) > 0)
        points(d_sp$med, centre[i] + seq(1.1, -1.1, length.out = nrow(d_sp)),
               pch = 16, cex = 0.7, col = adjustcolor(col_i, alpha.f = 0.45))
    }
    d <- dat[dat$row == row_levels[i], ]
    if (nrow(d) == 0) next
    segments(d$lo95, centre[i], d$hi95, centre[i], col = col_i, lwd = 1.5)
    segments(d$lo89, centre[i], d$hi89, centre[i], col = col_i, lwd = 3.5)
    if (row_levels[i] %in% open_rows) {
      points(d$med, centre[i], pch = 21, bg = "white", col = col_i, cex = 1.25, lwd = 2)
    } else {
      points(d$med, centre[i], pch = 21, bg = col_i, col = "white", cex = 1.2)
    }
  }
  if (!is.null(panel_label)) mtext(panel_label, side = 3, line = 0.2, adj = 0, font = 2, cex = 1.1, xpd = NA)
}

# --- Panels for the vital rates side by side, with the row labels on the first ----------
#   dat: estimates with columns row, series (vital rate) and the summary
draw_vr_panels <- function(dat, rows, row_cols, vr = vr_names, sp_dat = NULL, title_prefix = "",
                           letters_from = 1, headers = character(0), open_rows = character(0)) {
  for (v in seq_along(vr)) {
    par(mar = c(4, if (v == 1) 20 else 1, 2.5, 1), mgp = c(2.0, 0.8, 0), cex.lab = 1.1, cex.axis = 0.95)
    draw_forest_panel(dat[dat$series == vr[v], ], rows, row_cols, main = paste0(title_prefix, vr[v]),
                      xlab = "Posterior estimate (per SD)", panel_label = panel_letters[letters_from - 1 + v],
                      show_labels = v == 1, headers = headers, open_rows = open_rows,
                      sp_dat = if (is.null(sp_dat)) NULL else sp_dat[sp_dat$series == vr[v], ])
  }
}

# --- Colour key along the bottom of the figure (needs par(oma) with a bottom margin) -----
#   'open': labels drawn with an open point (interaction terms); 'note': a line of text below the key
add_type_legend <- function(labels, cols, open = rep(FALSE, length(labels)), note = NULL) {
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", legend = labels, col = cols, pch = 21, pt.bg = ifelse(open, "white", cols), pt.cex = 1.2,
         pt.lwd = ifelse(open, 2, 1), lwd = 2, bty = "n", ncol = min(length(labels), 4), cex = 0.95,
         inset = if (is.null(note)) 0.005 else 0.025, x.intersp = 0.6,
         text.width = max(strwidth(labels, cex = 0.95)) + 0.01)
  if (!is.null(note)) mtext(note, side = 1, line = -1, cex = 0.8, col = "grey30")
}

# --- Convergence: Rhat and effective sample size of every monitored parameter -----------
convergence_table <- function(samples, model_name) {
  m <- as.matrix(samples)
  keep <- colnames(m)[apply(m, 2, sd) > 0]   # constant nodes (e.g. diagonal of a correlation matrix) left aside
  rhat <- coda::gelman.diag(samples[, keep], autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  data.frame(model = model_name, parameter = keep, rhat = unname(rhat),
             ess = unname(coda::effectiveSize(samples[, keep])))
}

panel_letters <- paste0("(", letters, ")")

# ====================================================================================
# Read the posterior samples
# ====================================================================================
fit_main <- readRDS(fit_file(sprintf("CLR_PartA_stage%d", plot_stage)))
draws    <- as.matrix(fit_main$samples)
sp_lab   <- fit_main$species$code
cat(sprintf("Main model: %d posterior draws, %d species.\n", nrow(draws), length(sp_lab)))

# ====================================================================================
# FIG 1: All effects of the main model, one panel per vital rate
#   Main effects (filled points): log size and grazing (individual), precipitation anomaly and
#   the four crowding classes, as the means over species of the population-level coefficients
#   (pop.*.mu); the small faded points are the species medians (pop.*[k]), whose spread shows
#   response diversity.
#   Interaction terms (open points), shared by all species, in three blocks: precipitation x log
#   size and x crowding (delta), log size x crowding (theta), grazing x crowding (omega).
#   Recruitment has no size terms (it uses quadrat cover in place of size and crowding).
# File saved: CLR_PartA_main_effects.png
# ====================================================================================
head_main  <- "Main effects"
head_delta <- "\u00d7 Precip. anomaly (\u03b4)"
head_theta <- "\u00d7 Log size (\u03b8)"
head_omega <- "\u00d7 Grazing (\u03c9)"
main_rows  <- c(paste0(size_lab, " (\u03b11)"), paste0(size_lab, "\u00b2 (\u03b12)"), "Grazing (\u03b3)",
                paste0(clim_lab, ", + wetter (\u03b2 climate)"), paste0(crowd_lab, " (\u03b2)"))
delta_rows <- c(paste0(clim_lab, " \u00d7 log size"), paste0(clim_lab, " \u00d7 ", crowd_by))
theta_rows <- paste0("Log size \u00d7 ", crowd_by)
omega_rows <- paste0("Grazing \u00d7 ", crowd_by)
int_rows   <- c(delta_rows, theta_rows, omega_rows)
headers    <- c(head_main, head_delta, head_theta, head_omega)
  # Colour of a row: its covariate (main effects) or the term the block's covariate modifies
main_cols  <- setNames(c(col_ind, col_ind, col_ind, col_clim, col_W), main_rows)
int_cols   <- setNames(c(col_ind, col_W, col_W, col_W), int_rows)
all_rows   <- c(head_main, main_rows, head_delta, delta_rows, head_theta, theta_rows, head_omega, omega_rows)
all_cols   <- c(main_cols, int_cols)
type_lab   <- c(paste(size_lab, "and grazing"), "Precipitation", sub("Crowding by ", "Crowding: ", crowd_lab))
type_col   <- c(col_ind, col_clim, col_W)

d_main <- bind_rows(
  post_summary(draws, "pop.slope.size.mu")  %>% mutate(row = main_rows[1], series = vr_names[i1]),
  post_summary(draws, "pop.slope.size2.mu") %>% mutate(row = main_rows[2], series = vr_names[i1]),
  post_summary(draws, "pop.slope.graze.mu") %>% mutate(row = main_rows[3], series = vr_names[i1]),
  post_summary(draws, "pop.slope.clim.mu")  %>% mutate(row = main_rows[4], series = vr_names[i1]),
  post_summary(draws, "pop.slope.W.mu")     %>% mutate(row = main_rows[4 + i2], series = vr_names[i1]),
  post_summary(draws, "com.size.x.clim")    %>% mutate(row = delta_rows[1], series = vr_names[i1]),
  post_summary(draws, "com.W.x.clim")       %>% mutate(row = delta_rows[1 + i2], series = vr_names[i1]),
  post_summary(draws, "com.W.x.size")       %>% mutate(row = theta_rows[i2], series = vr_names[i1]),
  post_summary(draws, "com.W.x.graze")      %>% mutate(row = omega_rows[i2], series = vr_names[i1]))
  # Species medians of the population-level coefficients (one point per species and row)
d_main_sp <- bind_rows(
  post_summary(draws, "pop.slope.size")  %>% mutate(row = main_rows[1], series = vr_names[i2]),
  post_summary(draws, "pop.slope.size2") %>% mutate(row = main_rows[2], series = vr_names[i2]),
  post_summary(draws, "pop.slope.graze") %>% mutate(row = main_rows[3], series = vr_names[i2]),
  post_summary(draws, "pop.slope.clim")  %>% mutate(row = main_rows[4], series = vr_names[i2]),
  post_summary(draws, "pop.slope.W")     %>% mutate(row = main_rows[4 + i3], series = vr_names[i2])) %>%
  arrange(i1)

png(file.path(fig_dir, sprintf("CLR_PartA_main_effects%s.png", tag)), width = 15.5, height = 11, units = "in", res = 300)
layout(matrix(1:3, 1, 3), widths = c(1.65, 1, 1))
par(oma = c(6.5, 0, 0, 0))
draw_vr_panels(d_main, all_rows, all_cols, sp_dat = d_main_sp, headers = headers, open_rows = int_rows)
add_type_legend(c(type_lab, "Interaction (open)", "Species"),
                c(type_col, "grey30", adjustcolor("grey40", alpha.f = 0.45)),
                open = c(rep(FALSE, length(type_lab)), TRUE, FALSE), note = cov_note)
dev.off()
cat("Fig 1 saved: CLR_PartA_main_effects.png\n")

# ====================================================================================
# FIG 1b: Interaction terms of the main model only (own x scale), in the blocks and colours of
#   Fig 1
# File saved: CLR_PartA_interactions.png
# ====================================================================================
int_all_rows <- c(head_delta, delta_rows, head_theta, theta_rows, head_omega, omega_rows)
png(file.path(fig_dir, sprintf("CLR_PartA_interactions%s.png", tag)), width = 15.5, height = 8.5, units = "in", res = 300)
layout(matrix(1:3, 1, 3), widths = c(1.65, 1, 1))
par(oma = c(6.5, 0, 0, 0))
draw_vr_panels(d_main[d_main$row %in% int_rows, ], int_all_rows, int_cols, headers = headers, open_rows = int_rows)
add_type_legend(c(size_lab, type_lab[3:6]), c(col_ind, col_W), open = rep(TRUE, 5), note = cov_note)
dev.off()
cat("Fig 1b saved: CLR_PartA_interactions.png\n")

# ====================================================================================
# FIG 2: Crowding effects per species (beta_x[k, v]) and their community mean
#   Rows of panels: vital rates; columns: crowding classes
# File saved: CLR_PartA_crowding_species.png
# ====================================================================================
sp_rows <- c(sp_lab, "Mean")
d_beta <- bind_rows(
  post_summary(draws, "pop.slope.W")    %>% mutate(row = sp_lab[i1], series = vr_names[i2], class = crowd_lab[i3]),
  post_summary(draws, "pop.slope.W.mu") %>% mutate(row = "Mean", series = vr_names[i1], class = crowd_lab[i2]))

png(file.path(fig_dir, sprintf("CLR_PartA_crowding_species%s.png", tag)), width = 14, height = 11, units = "in", res = 300)
par(mfrow = c(3, 4), mar = c(4, 5, 2.5, 1), mgp = c(2.0, 0.8, 0), cex.lab = 1.1, cex.axis = 0.95)
p <- 0
for (v in vr_names) {
  for (cl in crowd_lab) {
    p <- p + 1
    draw_forest_panel(d_beta[d_beta$series == v & d_beta$class == cl, ], sp_rows,
                      setNames(rep(col_W[[cl]], length(sp_rows)), sp_rows),
                      main = paste0(v, ": ", tolower(cl)), xlab = if (v == vr_names[3]) "Posterior estimate (per SD)" else "",
                      panel_label = panel_letters[p])
  }
}
dev.off()
cat("Fig 2 saved: CLR_PartA_crowding_species.png\n")

# ====================================================================================
# FIG 3: Precipitation anomaly (beta_climate[k, v], blue) and grazing (gamma[k, v], orange)
#   effects per species. Rows of panels: precipitation, grazing; columns: vital rates
# File saved: CLR_PartA_climate_grazing_species.png
# ====================================================================================
d_bclim <- bind_rows(post_summary(draws, "pop.slope.clim") %>% mutate(row = sp_lab[i1], series = vr_names[i2]),
                     post_summary(draws, "pop.slope.clim.mu") %>% mutate(row = "Mean", series = vr_names[i1]))
d_graze <- bind_rows(post_summary(draws, "pop.slope.graze") %>% mutate(row = sp_lab[i1], series = vr_names[i2]),
                     post_summary(draws, "pop.slope.graze.mu") %>% mutate(row = "Mean", series = vr_names[i1]))

png(file.path(fig_dir, sprintf("CLR_PartA_climate_grazing_species%s.png", tag)), width = 13, height = 9, units = "in", res = 300)
par(mfrow = c(2, 3), mar = c(4, 5, 2.5, 1), mgp = c(2.0, 0.8, 0), cex.lab = 1.1, cex.axis = 0.95)
for (v in seq_along(vr_names)) {
  draw_forest_panel(d_bclim[d_bclim$series == vr_names[v], ], sp_rows, setNames(rep(col_clim, length(sp_rows)), sp_rows),
                    main = paste0(vr_names[v], ": precip. anomaly, + wetter (\u03b2 climate)"), panel_label = panel_letters[v])
}
for (v in seq_along(vr_names)) {
  draw_forest_panel(d_graze[d_graze$series == vr_names[v], ], sp_rows, setNames(rep(col_ind, length(sp_rows)), sp_rows),
                    main = paste0(vr_names[v], ": grazing (\u03b3)"), xlab = "Posterior estimate (per SD; grazing 0/1)",
                    panel_label = panel_letters[3 + v])
}
dev.off()
cat("Fig 3 saved: CLR_PartA_climate_grazing_species.png\n")

# ====================================================================================
# FIG 4: Variance components (standard deviations), one panel per vital rate
#   Year effects (blue), group, quadrat and plant effects (orange), variation among species
#   of each coefficient (sigma of the population level, red)
# File saved: CLR_PartA_variance.png
# ====================================================================================
yr_rows   <- c("Year", "Year \u00d7 species", "Year: log size slope", paste0("Year: ", crowd_by, " slope"))
  # Plant effects: rows of the plant (individual) effects present in the fit
has_plant <- any(grepl("^ind\\.sd", colnames(draws)))
plot_rows <- c("Group", "Quadrat", if (has_plant) "Plant")
sig_rows  <- c("Intercept (\u03c3 \u03b10)", "Log size (\u03c3 \u03b11)", "Log size\u00b2 (\u03c3 \u03b12)", "Grazing (\u03c3 \u03b3)",
               "Precip. anomaly (\u03c3 \u03b2 climate)", paste0(crowd_lab, " (\u03c3 \u03b2)"))
var_rows <- c(yr_rows, plot_rows, sig_rows)
var_cols <- setNames(c(rep(col_clim, length(yr_rows)), rep(col_ind, length(plot_rows)), rep(col_pop, length(sig_rows))),
                     var_rows)
d_var <- bind_rows(
  post_summary(draws, "yr.int.sd")          %>% mutate(row = yr_rows[1], series = vr_names[i1]),
  post_summary(draws, "yrsp.int.sd")        %>% mutate(row = yr_rows[2], series = vr_names[i1]),
  post_summary(draws, "yr.slope.size.sd")   %>% mutate(row = yr_rows[3], series = vr_names[i1]),
  post_summary(draws, "yr.slope.W.sd")      %>% mutate(row = yr_rows[3 + i2], series = vr_names[i1]),
  post_summary(draws, "grp.sd")             %>% mutate(row = "Group", series = vr_names[i1]),
  post_summary(draws, "quad.sd")            %>% mutate(row = "Quadrat", series = vr_names[i1]),
  if (has_plant) post_summary(draws, "ind.sd") %>% mutate(row = "Plant", series = vr_names[i1]),
  post_summary(draws, "pop.int.sd")         %>% mutate(row = sig_rows[1], series = vr_names[i1]),
  post_summary(draws, "pop.slope.size.sd")  %>% mutate(row = sig_rows[2], series = vr_names[i1]),
  post_summary(draws, "pop.slope.size2.sd") %>% mutate(row = sig_rows[3], series = vr_names[i1]),
  post_summary(draws, "pop.slope.graze.sd") %>% mutate(row = sig_rows[4], series = vr_names[i1]),
  post_summary(draws, "pop.slope.clim.sd")  %>% mutate(row = sig_rows[5], series = vr_names[i1]),
  post_summary(draws, "pop.slope.W.sd")     %>% mutate(row = sig_rows[5 + i2], series = vr_names[i1]))

png(file.path(fig_dir, sprintf("CLR_PartA_variance%s.png", tag)), width = 15, height = 8.5, units = "in", res = 300)
layout(matrix(1:3, 1, 3), widths = c(1.55, 1, 1))
par(oma = c(2.5, 0, 0, 0))
draw_vr_panels(d_var, var_rows, var_cols)
add_type_legend(c("Year effects", paste(plot_rows, collapse = ", "), "Variation among species (\u03c3)"),
                c(col_clim, col_ind, col_pop))
dev.off()
cat("Fig 4 saved: CLR_PartA_variance.png\n")

# ====================================================================================
# FIG 5: Submodels (pooled forb; pooled other shrubs), laid out as Fig 1
#   Blocks and colours as Fig 1 (main effects filled, interactions with precipitation open).
#   The submodels have one pooled population each, so these are their only coefficients. In the
#   forb submodel, forb crowding is crowding by the other forb plants (same pooled population).
# File saved: CLR_PartA_submodels.png (when both submodel fits exist)
# ====================================================================================
conv <- convergence_table(fit_main$samples, sprintf("main stage %d", plot_stage))
sub_files <- c(forb = fit_file("CLR_PartA_forb"), shrub = fit_file("CLR_PartA_shrub"))
if (all(file.exists(sub_files))) {
  fit_forb  <- readRDS(sub_files[["forb"]])
  fit_shrub <- readRDS(sub_files[["shrub"]])
  conv <- bind_rows(conv, convergence_table(fit_forb$samples, "forb"),
                    convergence_table(fit_shrub$samples, "other shrubs"))
  draws_forb  <- as.matrix(fit_forb$samples)
  draws_shrub <- as.matrix(fit_shrub$samples)

    # Forb: vital rates survival and recruitment; crowding classes forb, grass, shrub
  vr_forb <- c("Survival", "Recruitment")
  forb_lab <- "Crowding by other forbs"
  rows_forb_main <- c(main_rows[3:4], paste0(forb_lab, " (\u03b2)"), main_rows[4 + c(2, 3)])
  rows_forb_int  <- c(paste0(clim_lab, " \u00d7 ", tolower(forb_lab)), delta_rows[1 + c(2, 3)])
  d_forb <- bind_rows(
    post_summary(draws_forb, "slope.graze")  %>% mutate(row = rows_forb_main[1], series = vr_forb[i1]),
    post_summary(draws_forb, "slope.clim")   %>% mutate(row = rows_forb_main[2], series = vr_forb[i1]),
    post_summary(draws_forb, "slope.W")      %>% mutate(row = rows_forb_main[2 + i2], series = vr_forb[i1]),
    post_summary(draws_forb, "com.W.x.clim") %>% mutate(row = rows_forb_int[i2], series = vr_forb[i1]))

    # Other shrubs: vital rates survival, growth, recruitment; the four crowding classes
  rows_shrub_main <- main_rows[-2]   # (no size² in the other-shrub submodel)
  rows_shrub_int  <- delta_rows
  d_shrub <- bind_rows(
    post_summary(draws_shrub, "slope.size")      %>% mutate(row = main_rows[1], series = vr_names[i1]),
    post_summary(draws_shrub, "slope.graze")     %>% mutate(row = main_rows[3], series = vr_names[i1]),
    post_summary(draws_shrub, "slope.clim")      %>% mutate(row = main_rows[4], series = vr_names[i1]),
    post_summary(draws_shrub, "slope.W")         %>% mutate(row = main_rows[4 + i2], series = vr_names[i1]),
    post_summary(draws_shrub, "com.size.x.clim") %>% mutate(row = delta_rows[1], series = vr_names[i1]),
    post_summary(draws_shrub, "com.W.x.clim")    %>% mutate(row = delta_rows[1 + i2], series = vr_names[i1]))
  sub_cols <- c(all_cols, setNames(col_W[4], rows_forb_main[3]), setNames(col_W[c(4, 2, 3)], rows_forb_int))

  png(file.path(fig_dir, sprintf("CLR_PartA_submodels%s.png", tag)), width = 15, height = 11.5, units = "in", res = 300)
  layout(matrix(1:6, 2, 3, byrow = TRUE), widths = c(1.55, 1, 1), heights = c(8, 13))
  par(oma = c(6.5, 0, 0, 0))
  draw_vr_panels(d_forb, c(head_main, rows_forb_main, head_delta, rows_forb_int), sub_cols, vr = vr_forb,
                 title_prefix = "Forb (pooled): ", headers = headers, open_rows = rows_forb_int)
  plot.new()
  draw_vr_panels(d_shrub, c(head_main, rows_shrub_main, head_delta, rows_shrub_int), sub_cols,
                 title_prefix = "Other shrubs (pooled): ", letters_from = 3, headers = headers,
                 open_rows = rows_shrub_int)
  add_type_legend(c(type_lab, "Interaction (open)"), c(type_col, "grey30"),
                  open = c(rep(FALSE, length(type_lab)), TRUE), note = cov_note)
  dev.off()
  cat("Fig 5 saved: CLR_PartA_submodels.png\n")
} else {
  cat("Fig 5 left out: the submodel fits are missing (", paste(basename(sub_files[!file.exists(sub_files)]), collapse = ", "), ").\n", sep = "")
}

# ====================================================================================
# Convergence table
# ====================================================================================
write.csv(conv, file.path(tab_dir, sprintf("CLR_PartA_convergence%s.csv", tag)), row.names = FALSE)
conv_summary <- conv %>% group_by(model) %>%
  summarise(n_par = n(), max_rhat = max(rhat, na.rm = TRUE), n_rhat_over_1.05 = sum(rhat > 1.05, na.rm = TRUE),
            min_ess = min(ess), .groups = "drop")
print(as.data.frame(conv_summary))

cat("Saved:\n",
    " result/figs/", site_tag, ": CLR_PartA_main_effects, _interactions, _crowding_species, _climate_grazing_species, _variance",
    if (all(file.exists(sub_files))) ", _submodels" else "", " (", tag, ".png)\n",
    " result/tables/", site_tag, ": CLR_PartA_convergence", tag, ".csv\n", sep = "")
rm(list = intersect(c("draws", "d_main", "d_main_sp", "d_beta", "d_bclim", "d_graze", "d_var", "d_forb", "d_shrub",
                      "draws_forb", "draws_shrub", "vr_forb", "rows_forb_main", "rows_forb_int", "rows_shrub_main",
                      "rows_shrub_int", "main_rows", "delta_rows", "theta_rows", "omega_rows", "int_rows",
                      "main_cols", "int_cols", "sub_cols", "all_rows", "all_cols", "int_all_rows", "headers",
                      "head_main", "head_delta", "head_theta", "head_omega", "type_lab", "type_col", "has_plant",
                      "forb_lab", "sp_rows", "yr_rows", "plot_rows", "sig_rows", "var_rows",
                      "var_cols", "p", "v", "cl", "conv_summary"), ls()))
