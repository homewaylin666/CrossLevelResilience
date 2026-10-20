################################################################################
## How much do the tracking choices change the vital rates? (Idaho)
##
## Written by Hung-wei Lin Oct 2026
##
## First calibrates the first-year grouping distance (buffGenet) against the hand
## assignment of Chu & Adler (2015). Then re-tracks the main size species under
## alternative tracking choices, fits simple
## vital-rate regressions per species (survival and growth on log size, recruits per
## cm2 of parent cover), and expresses every change against the main tracking in units
## of the standard error of the main estimate. A change below 1 SE sits inside sampling
## noise; a change above 2 SE would move the Part A estimates by more than their own
## uncertainty. Also checks the year-to-year pattern of survival and recruitment (the
## signal the climate terms use).
## Uses objects of CLR_TrackGenet.R ('records', 'tracked', 'track_records', 'flag_fates',
## settings and lookups), so it runs right after that module.

# Settings ---------------------------------------------------------------------
  # Species: polygon species with at least this many annual survival records (main tracking)
sens_min_survival <- 500
  # Re-tracking takes a few minutes per alternative; saved and read again on later runs
resens <- FALSE

# Calibrate buffGenet ----------------------------------------------------------
  # Agreement between the distance rule and the hand assignment, over the first-year
  # polygons of each quadrat (adjusted Rand index: 1 = same grouping, 0 = chance level)
ari <- function(a, b) {
  t <- table(a, b); s2 <- function(x) sum(choose(x, 2))
  ea <- s2(rowSums(t)); eb <- s2(colSums(t)); ex <- ea * eb / choose(sum(t), 2)
  if ((ea + eb) / 2 == ex) return(NA)
  (s2(t) - ex) / ((ea + eb) / 2 - ex)
}
calib <- expand.grid(Species = unname(chu_code), buffGenet = c(0.005, 0.01, 0.02, 0.03, 0.05),
                     stringsAsFactors = FALSE)
calib$ari <- mapply(function(sp, bg) {
  d <- records[!is.na(records$hand_genet) & records$Species == sp, ]
  v <- sapply(split(d, paste(d$Quad, d$Year)), function(x) {
    if (nrow(x) < 2) return(NA)
    ari(x$hand_genet, groupByGenet(x, buffGenet = bg))
  })
  mean(v, na.rm = TRUE)
}, calib$Species, calib$buffGenet)
calib$code <- names(chu_code)[match(calib$Species, chu_code)]
write.csv(calib, file.path(tab_dir, "CLR_TrackSensitivity_buffGenet.csv"), row.names = FALSE)
print(pivot_wider(calib[, c("code", "buffGenet", "ari")], names_from = buffGenet, values_from = ari), digits = 2)

# Species to check -------------------------------------------------------------
annual_next <- function(q, y) mapply(function(qq, yy) (yy + 1L) %in% inv[[qq]], q, y)
main_df <- st_drop_geometry(tracked)
n_surv <- main_df %>%
  filter(record_type == "polygon") %>%
  distinct(Species, Quad, Year, trackID, .keep_all = TRUE) %>%
  filter(!is.na(survives_tplus1), annual_next(Quad, Year)) %>%
  count(Species)
sens_species <- n_surv$Species[n_surv$n >= sens_min_survival]
cat("Species checked:", paste(sens_species, collapse = ", "), "\n")

# Alternative trackings --------------------------------------------------------
  # Each alternative changes one choice of the main tracking:
  #   'grouping 0.5 cm', 'grouping 3 cm' - first-year pieces of every species grouped within
  #   the narrow or the wide end of the calibration range, without the hand genets
  #   'no seedling rule' - seedlings may link to a plant of the year before
  #   'merged = death'   - coalescence counts as death (as plantTracker and Chu & Adler 2015)
  #   grass seedling flags of Chu & Adler (2015) for HECO, POSE, PSSP: the opposite of the
  #   main setting ('Chu seedling flags' or 'public seedling flags only')
sens_file <- file.path(der_dir, "CLR_TrackSensitivity_ramet.rds")
if (file.exists(sens_file) && !resens) {
  alt <- readRDS(sens_file)
} else {
  recs <- records[records$Species %in% sens_species, ]
  keep <- list(use_hand_genets = use_hand_genets, buffGenet = buffGenet)
  alt <- list()

  use_hand_genets <- FALSE
  for (bg in c(0.005, 0.03)) {
    buffGenet$buffGenet <- bg
    alt[[sprintf("grouping %g cm", bg * 100)]] <- st_drop_geometry(flag_fates(track_records(recs)))
  }
  use_hand_genets <- keep$use_hand_genets
  buffGenet <- keep$buffGenet

  r2 <- recs; r2$seedling <- NA
  alt[["no seedling rule"]] <- st_drop_geometry(flag_fates(track_records(r2)))

  r3 <- recs
  if (seedling_from_chu) {
    r3$seedling[r3$record_type == "polygon"] <- NA
    lab_seed <- "public seedling flags only"
  } else {
    r3$seedling[r3$record_id %in% chu_seed] <- "Y"
    lab_seed <- "Chu seedling flags"
  }
  alt[[lab_seed]] <- st_drop_geometry(flag_fates(track_records(r3)))
  saveRDS(alt, sens_file)
}
m <- main_df[main_df$Species %in% sens_species, ]
alt[["main"]] <- m
alt[["merged = death"]] <- mutate(m, survives_tplus1 = survives_raw)

  # Chu & Adler (2015) tracked data, four species (area in cm2)
alt[["Chu & Adler 2015"]] <- bind_rows(lapply(names(chu_code), function(k) {
  read.csv(file.path(chu_dir, k, paste0(k, "_buf5_dorm2.csv"))) %>%
    transmute(Quad = quad, Year = year + 1900L, Species = chu_code[k],
              trackID = paste(quad, trackID), basalArea_ramet = area / 1e4,
              survives_tplus1 = survives, recruit = as.integer(age == 1))
}))

# Genet-year tables ------------------------------------------------------------
to_genet <- function(d) {
  g <- d %>%
    group_by(Species, Quad, Year, trackID) %>%
    summarise(area = sum(basalArea_ramet) * 1e4,
              surv = suppressWarnings(max(survives_tplus1, na.rm = TRUE)),
              recruit = suppressWarnings(max(recruit, na.rm = TRUE)), .groups = "drop") %>%
    mutate(surv = ifelse(is.finite(surv), surv, NA), recruit = ifelse(is.finite(recruit), recruit, NA),
           logsize = log(area), annual = annual_next(Quad, Year))
  nxt <- transmute(g, Species, Quad, trackID, Year = Year - 1L, logsize_t1 = logsize)
  g %>% left_join(nxt, by = c("Species", "Quad", "trackID", "Year")) %>%
    mutate(surv = ifelse(annual, surv, NA), logsize_t1 = ifelse(annual & surv %in% 1, logsize_t1, NA))
}
genets <- lapply(alt, to_genet)

# Vital-rate estimates ---------------------------------------------------------
  # Log size centred on each species' median under the main tracking, so the survival and
  # growth intercepts describe a plant of typical size
centre <- genets[["main"]] %>% group_by(Species) %>% summarise(mid = median(logsize))
fit_one <- function(g, sp) {
  g <- g[g$Species == sp, ]
  if (nrow(g) == 0) return(NULL)
  g$z <- g$logsize - centre$mid[centre$Species == sp]
  s <- glm(surv ~ z, binomial, data = g[!is.na(g$surv), ])
  gr <- lm(logsize_t1 ~ z, data = g[!is.na(g$logsize_t1), ])
  sg <- sigma(gr)
    # Recruits at t+1 per cm2 of parent cover at t, over annual quadrat transitions
  cov <- g %>% filter(annual) %>% group_by(Quad, Year) %>% summarise(cover = sum(area), .groups = "drop")
  rec <- g %>% filter(recruit %in% 1) %>% count(Quad, Year = Year - 1L, name = "rec")
  r <- left_join(cov, rec, by = c("Quad", "Year")) %>% mutate(rec = replace_na(rec, 0))
  n_rec <- sum(r$rec)
  data.frame(Species = sp,
             parameter = c("survival: intercept", "survival: size slope",
                           "growth: intercept", "growth: size slope", "growth: SD",
                           "recruits per cover (log)"),
             estimate = c(coef(s), coef(gr), sg, log(n_rec / sum(r$cover))),
             se = c(sqrt(diag(vcov(s))), sqrt(diag(vcov(gr))), sg / sqrt(2 * (nobs(gr) - 2)),
                    1 / sqrt(n_rec)),
             n = c(nobs(s), nobs(s), nobs(gr), nobs(gr), nobs(gr), n_rec))
}
est <- bind_rows(lapply(names(genets), function(v) {
  bind_rows(lapply(sens_species, function(sp) fit_one(genets[[v]], sp))) %>% mutate(variant = v)
}))
est <- est %>%
  left_join(est %>% filter(variant == "main") %>% select(Species, parameter, est_main = estimate, se_main = se),
            by = c("Species", "parameter")) %>%
  mutate(change_in_se = (estimate - est_main) / se_main)

  # Year-to-year pattern: correlation of yearly survival and yearly recruits with the main tracking
yearly <- lapply(genets, function(g) g %>% filter(annual) %>% group_by(Species, Year) %>%
                   summarise(surv = mean(surv, na.rm = TRUE), .groups = "drop"))
yearly_rec <- lapply(genets, function(g) g %>% filter(recruit %in% 1) %>% count(Species, Year))
year_cor <- bind_rows(lapply(setdiff(names(genets), "main"), function(v) {
  bind_rows(lapply(sens_species, function(sp) {
    a <- inner_join(filter(yearly[["main"]], Species == sp), filter(yearly[[v]], Species == sp), by = "Year")
    b <- inner_join(filter(yearly_rec[["main"]], Species == sp), filter(yearly_rec[[v]], Species == sp), by = "Year")
    if (nrow(a) < 5) return(NULL)
    data.frame(Species = sp, variant = v, cor_yearly_survival = cor(a$surv.x, a$surv.y, use = "complete"),
               cor_yearly_recruits = if (nrow(b) >= 5) cor(b$n.x, b$n.y) else NA)
  }))
}))

write.csv(est, file.path(tab_dir, "CLR_TrackSensitivity_estimates.csv"), row.names = FALSE)
write.csv(year_cor, file.path(tab_dir, "CLR_TrackSensitivity_years.csv"), row.names = FALSE)

  # Summary: largest change (in SE) per alternative and parameter, over species
sens_summary <- est %>% filter(variant != "main") %>%
  group_by(variant, parameter) %>%
  summarise(max_abs_change_se = max(abs(change_in_se)), species_above_2se = sum(abs(change_in_se) > 2),
            n_species = n(), .groups = "drop")
print(as.data.frame(sens_summary), digits = 2)
print(year_cor %>% group_by(variant) %>%
        summarise(min_cor_survival = min(cor_yearly_survival), min_cor_recruits = min(cor_yearly_recruits, na.rm = TRUE)),
      digits = 3)

# Figure -----------------------------------------------------------------------
est$code <- toupper(paste0(substr(est$Species, 1, 2), substr(sub("^\\S+\\s*", "", est$Species), 1, 2)))
p_sens <- ggplot(filter(est, variant != "main"), aes(change_in_se, parameter)) +
  annotate("rect", xmin = -2, xmax = 2, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_point(aes(colour = code), size = 2, alpha = 0.8, position = position_jitter(height = 0.15, width = 0)) +
  facet_wrap(~ variant, ncol = 2) +
  labs(x = "Change against the main tracking (in SE of the main estimate)", y = NULL, colour = "Species",
       title = "Sensitivity of vital rates to tracking choices, Idaho",
       subtitle = "Grey band = within 2 SE of the main estimate") +
  theme_bw()
ggsave(file.path(fig_dir, "CLR_TrackSensitivity.png"), p_sens, width = 10, height = 8, dpi = 200)
