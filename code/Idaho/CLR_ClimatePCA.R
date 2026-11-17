################################################################################
## Principal component analysis of the climate variables (Idaho), for exploration
##
## Written by Hung-wei Lin Nov 2026
##
## Three sets of six seasonal variables (CLR_ClimateData.R), each taken for two climate years:
## the climate year of the transition (suffix 0: July t to June t+1, between the censuses of t and t+1)
## and the year before it (suffix 1), so 12 variables per PCA:
##   (1) Pp, Pg, Tp, Tg, Tmax, Tmin;
##   (2) P_year, Pg_share, T_year, Tg_Tp, Tmax, Tmin (amount and seasonal share);
##   (3) P_pg, Pg_share, T_pg, Tg_Tp, Tmax, Tmin (as 2, with October-June in place of July-June).
## The PCA uses the climate years of the station record (centred and scaled variables, so the
## components describe the long-term climate, as the anomalies of the Part A models); the years of
## the Part A transitions are marked. Figure 'CLR_ClimatePCA.png' (PC1 against PC2 of the three sets);
## loadings and the variance explained in 'CLR_ClimatePCA_loadings.csv'.
## Uses 'clim_year' (CLR_ClimateData.R), 'genet' and 'site_tag'.


# Settings ---------------------------------------------------------------------
  # Climate years of each variable: 0 = the year of the transition, 1 = the year before
pca_periods <- c(0, 1)
var_sets <- list(
  `Set 1: seasonal totals and means` = c("Pp", "Pg", "Tp", "Tg", "Tmax", "Tmin"),
  `Set 2: July-June amount and seasonal share` = c("P_year", "Pg_share", "T_year", "Tg_Tp", "Tmax", "Tmin"),
  `Set 3: October-June amount and seasonal share` = c("P_pg", "Pg_share", "T_pg", "Tg_Tp", "Tmax", "Tmin"))
  # Components kept in the tables
n_pc <- 4

tab_dir <- here("result", "tables", site_tag)
fig_dir <- here("result", "figs", site_tag)
stopifnot(exists("clim_year"), exists("genet"))

# Variables of each year -------------------------------------------------------
  # Row = climate year cy of the transition (period 0); period k takes the values of year cy - k
lagged <- function(vars) {
  out <- data.frame(cy = clim_year$cy)
  for (k in pca_periods) for (v in vars)
    out[[paste0(v, "_", k)]] <- clim_year[[v]][match(out$cy - k, clim_year$cy)]
  out[complete.cases(out), ]
}
  # Climate years of the Part A transitions: transition t -> t+1 has climate year t + 1
trans_cy <- sort(unique(genet$Year[!is.na(genet$survives_tplus1)])) + 1

# PCA of each set --------------------------------------------------------------
pcs <- lapply(var_sets, function(vars) {
  d <- lagged(vars)
  p <- prcomp(d[-1], center = TRUE, scale. = TRUE)
  list(cy = d$cy, pca = p, var_exp = p$sdev^2 / sum(p$sdev^2))
})
  # Loadings of the first components, with the variance each explains on top
loadings <- bind_rows(lapply(names(pcs), function(s) {
  r <- pcs[[s]]$pca$rotation[, 1:n_pc]
  v <- pcs[[s]]$var_exp[1:n_pc]
  data.frame(set = s, variable = c("variance explained", "cumulative variance", rownames(r)),
             rbind(v, cumsum(v), r), row.names = NULL)
}))
for (s in names(pcs)) {
  cat(sprintf("%s: %d climate years; PC1 %.0f%%, PC2 %.0f%%, PC3 %.0f%% of the variance.\n", s, length(pcs[[s]]$cy),
              100 * pcs[[s]]$var_exp[1], 100 * pcs[[s]]$var_exp[2], 100 * pcs[[s]]$var_exp[3]))
  print(round(pcs[[s]]$pca$rotation[, 1:3], 2))
}

# Figure: PC1 against PC2 of the three sets ------------------------------------
  # Years as points (Part A transitions in black); arrows = loadings, stretched to the spread of the
  # scores, red for temperature and blue for precipitation, solid for the year of the transition (0)
  # and dashed for the year before (1)
scores <- bind_rows(lapply(names(pcs), function(s) {
  data.frame(set = s, cy = pcs[[s]]$cy, transition_year = pcs[[s]]$cy %in% trans_cy, pcs[[s]]$pca$x[, 1:2])
}))
arrows <- bind_rows(lapply(names(pcs), function(s) {
  p <- pcs[[s]]$pca
  k <- 0.8 * max(abs(p$x[, 1:2])) / max(abs(p$rotation[, 1:2]))
  data.frame(set = s, variable = rownames(p$rotation), x = k * p$rotation[, 1], y = k * p$rotation[, 2])
})) %>%
  mutate(type = ifelse(substr(variable, 1, 1) == "T", "Temperature", "Precipitation"),
         period = ifelse(grepl("_1$", variable), "Year before (1)", "Year of the transition (0)"))
set_lab <- sapply(names(pcs), function(s) sprintf("%s\nPC1 %.0f%%, PC2 %.0f%%", s, 100 * pcs[[s]]$var_exp[1], 100 * pcs[[s]]$var_exp[2]))
scores$lab <- set_lab[scores$set]; arrows$lab <- set_lab[arrows$set]
p_pca <- ggplot(scores, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(shape = transition_year), colour = "grey55", size = 1.3) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = x, yend = y, colour = type, linetype = period),
               arrow = arrow(length = unit(0.15, "cm"))) +
  geom_text(data = arrows, aes(x = 1.1 * x, y = 1.1 * y, label = variable, colour = type,
                               hjust = ifelse(x > 0, 0, 1)), size = 2.8, show.legend = FALSE) +
  facet_wrap(~ lab, nrow = 1) +
  scale_colour_manual(values = c(Temperature = "#c0392b", Precipitation = "#1f5fbf")) +
  scale_linetype_manual(values = c(`Year of the transition (0)` = "solid", `Year before (1)` = "22")) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c(`FALSE` = "Other years of the station record", `TRUE` = "Part A transitions")) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  labs(colour = NULL, linetype = NULL, shape = NULL, title = "Climate PCA, Dubois (U.S. Sheep Experiment Station)",
       subtitle = paste("Suffix 0: climate year of the transition (July-June, between the censuses); 1: the year before.",
                        "Pp, Tp: October-March; Pg, Tg: April-June; Tmax, Tmin: warmest and coldest month")) +
  coord_equal() + theme_bw() + theme(legend.position = "bottom", panel.grid = element_blank())
ggsave(file.path(fig_dir, "CLR_ClimatePCA.png"), p_pca, width = 16, height = 6.5, dpi = 300)

# Save -------------------------------------------------------------------------
write.csv(loadings, file.path(tab_dir, "CLR_ClimatePCA_loadings.csv"), row.names = FALSE)
cat("Saved:\n",
    " result/tables/", site_tag, ": CLR_ClimatePCA_loadings.csv (with the variance explained)\n",
    " result/figs/", site_tag, ": CLR_ClimatePCA.png\n", sep = "")
rm(list = intersect(c("lagged", "trans_cy", "pcs", "loadings", "variance", "scores", "arrows", "axis_lab",
                      "set_lab", "p_pca", "s", "n_pc"), ls()))
