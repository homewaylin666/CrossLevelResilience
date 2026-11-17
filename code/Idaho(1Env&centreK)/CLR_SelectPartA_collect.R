################################################################################
## Model selection of Part A (Idaho): rank the models of one phase
##
## Written by Hung-wei Lin Oct 2026
##
## Phase 'select_phase' (set in the runscript):
##   1: climate variable. Current precipitation stays unless another variable predicts left-out years
##      better by more than 2 SE (leave-one-year-out; the year is the replicate of climate).
##   2: interaction terms. The model with the lowest AICc (Gupta et al. 2025), BIC alongside.
##   3: model_after (exploratory): AICc ranking of the trait terms; model_after keeps the structure
##      of model_before.
## Tables in 'result/tables/Idaho/PartA_model.selection', figures in 'result/figs/...'; the choices
## that later steps read ('CLR_select_climate_chosen.csv', 'CLR_select_selected.csv') in the
## selection folder.


# Settings ---------------------------------------------------------------------
  # Models left out of phase 2 after a look at their estimates (biologically implausible), e.g. "surv_00123"
exclude_models <- character(0)

tab_dir <- here("result", "tables", site_tag, "PartA_model.selection")
fig_dir <- here("result", "figs", site_tag, "PartA_model.selection")
dir.create(tab_dir, showWarnings = FALSE); dir.create(fig_dir, showWarnings = FALSE)

  # All summaries of a phase in one table
read_phase <- function(phase) {
  tk <- read.csv(file.path(sel_dir(), sprintf("CLR_select_tasks_phase%d.csv", phase)))
  res <- lapply(file.path(sel_dir(), tk$file), readRDS)
  get <- function(x) sapply(res, function(r) if (is.null(r[[x]])) NA else r[[x]])
  data.frame(tk[c("vr", "model", "terms", "tterms", "clim_use", "year_out")],
             ok = get("ok"), df = get("df"), logLik = get("logLik"), AICc = get("AICc"), BIC = get("BIC"),
             lpd = get("lpd"))
}

  # Rank figure: terms of each model (filled = in), and delta AICc or BIC (square-root axis)
draw_rank <- function(tabs, crit, file, mark, n_show = 20) {
  png(file, width = 17, height = 8, units = "in", res = 300)
  par(mfrow = c(1, 3), oma = c(4, 0, 1, 0), mar = c(4.5, 1, 3, 1))
  for (v in vrs) {
    x <- tabs[[v]]; dx <- x[[crit]] - min(x[[crit]])
    o <- order(dx); o <- o[seq_along(o) <= n_show | x$mark[o]]
    x <- x[o, ]; dx <- dx[o]; cds <- attr(tabs[[v]], "codes")
    n <- nrow(x); yy <- n:1; xr <- sqrt(max(dx, 4)); step <- 0.11 * xr
    plot(NA, xlim = c(-length(cds) * step, xr * 1.12), ylim = c(0.5, n + 2), axes = FALSE, xlab = "", ylab = "",
         main = tools::toTitleCase(v), font.main = 1, cex.main = 1.3)
    for (j in seq_along(cds)) {
      on <- x[[cds[j]]] == 1; cl <- c(C = "#1560b8", S = "#c9822e", G = "#8c6d46", T = "#b03a3a")[substr(cds[j], 1, 1)]
      points(rep(-(length(cds) - j + 1) * step, n), yy, pch = ifelse(on, 15, 0), cex = 1.2, col = ifelse(on, cl, "grey80"))
      text(-(length(cds) - j + 1) * step, n + 1, cds[j], srt = 90, adj = 0, font = 2, col = cl, cex = 0.8)
    }
    rect(0, yy - 0.35, sqrt(dx), yy + 0.35, col = ifelse(dx <= 2, "grey35", "grey75"), border = NA)
    text(sqrt(dx[x$mark]) + 0.02 * xr, yy[x$mark], mark, adj = 0, font = 2, cex = 0.85)
    at <- c(0, 2, 10, 50, 100, 500, 1000); at <- at[sqrt(at) <= xr]
    axis(1, at = sqrt(at), labels = at); mtext(paste0("Δ", crit), side = 1, line = 2.6, cex = 0.8)
  }
  lab <- rbind(ix_terms[c("code", "label")], trait_terms[c("code", "label")])
  lab <- lab[lab$code %in% unlist(lapply(tabs, attr, "codes")), ]
  lines <- split(paste(lab$code, lab$label), ceiling(seq_len(nrow(lab)) / 6))   # six terms per line
  for (i in seq_along(lines)) mtext(paste(lines[[i]], collapse = ";  "), side = 1, outer = TRUE, line = i - 0.5, cex = 0.7, col = "grey25")
  dev.off()
}


# Phase 1: climate variable ----------------------------------------------------
if (select_phase == 1) {
  r <- read_phase(1)
  clim_tab <- do.call(rbind, lapply(vrs, function(v) {
    cur <- r$lpd[r$vr == v & r$clim_use == "current" & !is.na(r$year_out)]
    x <- do.call(rbind, lapply(clim_choices, function(u) {
      a <- r$lpd[r$vr == v & r$clim_use == u & !is.na(r$year_out)]
      f <- r[r$vr == v & r$clim_use == u & is.na(r$year_out), ]
      data.frame(vr = v, clim_use = u, AICc = f$AICc, BIC = f$BIC, lpd = sum(a), dlpd = sum(a - cur),
                 se = sqrt(length(a)) * sd(a - cur))
    }))
    x$dAICc <- x$AICc - x$AICc[1]; x$dBIC <- x$BIC - x$BIC[1]
    better <- x$dlpd > 2 * x$se & x$clim_use != "current"
    x$chosen <- if (any(better)) x$dlpd == max(x$dlpd[better]) else x$clim_use == "current"
    x
  }))
  write.csv(clim_tab, file.path(tab_dir, "CLR_PartA_selection_climate.csv"), row.names = FALSE)
  write.csv(clim_tab[clim_tab$chosen, c("vr", "clim_use")], file.path(sel_dir(), "CLR_select_climate_chosen.csv"), row.names = FALSE)
  print(clim_tab, digits = 4)

    # Figure: difference from current precipitation in leave-one-year-out log predictive density (± 2 SE)
  lab <- c(current = "Precipitation (current year)", lag = "Precipitation (previous year)",
           both = "Precipitation (both years)", temp = "Summer temperature")
  png(file.path(fig_dir, "CLR_PartA_selection_climate.png"), width = 14, height = 4.2, units = "in", res = 300)
  par(mfrow = c(1, 3), mar = c(4.5, 13, 3, 6.5))
  for (v in vrs) {
    x <- clim_tab[clim_tab$vr == v, ]; se <- ifelse(is.na(x$se), 0, x$se)
    plot(x$dlpd, 4:1, xlim = range(x$dlpd - 2 * se, x$dlpd + 2 * se, 0) + c(-1, 1), ylim = c(0.5, 4.5), yaxt = "n",
         pch = 21, bg = ifelse(x$chosen, "#1560b8", "white"), col = "#1560b8", cex = 1.6, ylab = "",
         xlab = "Δ log predictive density (leave one year out)", main = tools::toTitleCase(v), font.main = 1)
    segments(x$dlpd - 2 * se, 4:1, x$dlpd + 2 * se, 4:1, col = "#1560b8")
    axis(2, at = 4:1, labels = lab[x$clim_use], las = 1, tick = FALSE)
    axis(4, at = 4:1, labels = sprintf("ΔAICc %.1f\nΔBIC %.1f", x$dAICc, x$dBIC), las = 1, tick = FALSE, cex.axis = 0.75)
    abline(v = 0, lty = 2, col = "grey50")
  }
  dev.off()
}


# Phase 2: interaction terms ---------------------------------------------------
if (select_phase == 2) {
  r <- read_phase(2)
  clim <- read.csv(file.path(sel_dir(), "CLR_select_climate_chosen.csv"))
  tabs <- list(); check <- NULL
  for (v in vrs) {
    cds <- terms_of(v)$code
    x <- r[r$vr == v & r$ok & !r$model %in% exclude_models, ]
    for (cd in cds) x[[cd]] <- as.integer(grepl(paste0("\\b", cd, "\\b"), x$terms))
    x$subset_id <- as.integer(sub(".*_", "", x$model))
    x$dAICc <- x$AICc - min(x$AICc); x$weight <- exp(-x$dAICc / 2) / sum(exp(-x$dAICc / 2))
    x$dBIC <- x$BIC - min(x$BIC)
      # Wald approximation of each subset from the full model (check of the exact fits)
    full <- readRDS(file.path(sel_dir(), sprintf("phase1_%s_full_%s.rds", vr_short[v], clim$clim_use[clim$vr == v])))
    x <- merge(x, wald_aicc(full), by = "subset_id")
    x <- x[order(x$AICc), ]
    x$mark <- seq_len(nrow(x)) == 1
    attr(x, "codes") <- cds
    tabs[[v]] <- x
    top <- 1:10
    check <- rbind(check, data.frame(vr = v, n_fits = nrow(x),
                                     spearman = cor(x$AICc, x$AICc_wald, method = "spearman"),
                                     max_diff_top10 = max(abs(x$dAICc[top] - (x$AICc_wald[top] - min(x$AICc_wald)))),
                                     same_best = x$subset_id[1] == x$subset_id[which.min(x$AICc_wald)]))
    write.csv(x[c("model", cds, "df", "logLik", "AICc", "dAICc", "weight", "BIC", "dBIC", "AICc_wald")],
              file.path(tab_dir, sprintf("CLR_PartA_selection_%s.csv", v)), row.names = FALSE)
  }
    # Importance of each term: sum of the AICc weights of the models that hold it (Burnham & Anderson 2002)
  imp <- do.call(rbind, lapply(vrs, function(v) {
    x <- tabs[[v]]; cds <- attr(x, "codes")
    data.frame(vr = v, code = cds, term = ix_terms$label[match(cds, ix_terms$code)],
               importance = sapply(cds, function(cd) sum(x$weight[x[[cd]] == 1])))
  }))
  write.csv(imp, file.path(tab_dir, "CLR_PartA_selection_importance.csv"), row.names = FALSE)
  write.csv(check, file.path(tab_dir, "CLR_PartA_selection_wald_check.csv"), row.names = FALSE)
    # Selected model of each vital rate (read by CLR_ModelPartA.R with structure "selected")
  sel <- do.call(rbind, lapply(vrs, function(v) {
    x <- tabs[[v]][1, ]
    o <- data.frame(vr = v, model = x$model, terms = x$terms, clim_use = clim$clim_use[clim$vr == v], AICc = x$AICc)
    for (cd in ix_terms$code) o[[cd]] <- if (cd %in% names(x)) x[[cd]] else NA
    o
  }))
  write.csv(sel, file.path(sel_dir(), "CLR_select_selected.csv"), row.names = FALSE)
  write.csv(sel, file.path(tab_dir, "CLR_PartA_selection_selected.csv"), row.names = FALSE)
  draw_rank(tabs, "AICc", file.path(fig_dir, "CLR_PartA_selection_AICcrank.png"), "selected")
  draw_rank(tabs, "BIC", file.path(fig_dir, "CLR_PartA_selection_BICrank.png"), "selected")
  print(sel[c("vr", "model", "terms", "clim_use")]); print(check, digits = 3)
}


# Phase 3: model_after (exploratory) -------------------------------------------
if (select_phase == 3) {
  r <- read_phase(3)
  tabs <- list()
  for (v in vrs) {
    x <- r[r$vr == v & r$ok, ]
    for (cd in trait_terms$code) x[[cd]] <- as.integer(grepl(paste0("\\b", cd, "\\b"), x$tterms))
    x$dAICc <- x$AICc - min(x$AICc); x$dBIC <- x$BIC - min(x$BIC)
    x <- x[order(x$AICc), ]
    x$mark <- x$tterms == ""
    attr(x, "codes") <- trait_terms$code
    tabs[[v]] <- x
    write.csv(x[c("model", trait_terms$code, "df", "AICc", "dAICc", "BIC", "dBIC")],
              file.path(tab_dir, sprintf("CLR_PartA_selection_after_%s.csv", v)), row.names = FALSE)
  }
  draw_rank(tabs, "AICc", file.path(fig_dir, "CLR_PartA_selection_after_AICcrank.png"), "model_before")
}

cat("Saved:\n  result/tables/", site_tag, "/PartA_model.selection: CLR_PartA_selection_* (phase ", select_phase, ")\n",
    "  result/figs/", site_tag, "/PartA_model.selection: CLR_PartA_selection_*.png\n", sep = "")
rm(list = intersect(c("r", "clim", "clim_tab", "tabs", "check", "imp", "sel", "full", "x", "v", "cd", "cds", "top",
                      "lab", "se", "lines", "i", "read_phase", "draw_rank", "tab_dir", "fig_dir"), ls()))
