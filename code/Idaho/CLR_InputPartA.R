################################################################################
## Build the input file of the Part A models (Idaho)
##
## Written by Hung-wei Lin Oct 2026
##
## Collects everything the Part A models read into one file, 'CLR_PartA_input.rds' (a list):
##   (1) main: the main model (species with at least 'main_min_survival' survival records):
##       survival and growth records of focal plants, recruitment per species x quadrat x
##       transition, sizes of recruits;
##   (2) forb: the forb submodel (all point species pooled as one population): survival records
##       and recruitment per quadrat x transition;
##   (3) shrub: the other-shrub submodel (shrubs below the threshold, pooled as one population):
##       survival and growth records, recruitment per quadrat x transition, sizes of recruits;
##   (4) index tables (species, years, quadrats, groups) and the scaling constants used for size,
##       crowding, quadrat cover and climate (kept for the IPMs of Part B).
## Covariates: size z = log area (cm2), centred and scaled; crowding W (CLR_CrowdIndex.R): crowding by
## other grasses, other shrubs and forbs centred and scaled per class over all focal plants (same
## units for every focal species), same-species crowding centred and scaled within each species
## (its units belong to the species: basal area of a grass, canopy area of a shrub); grazing E = 1 when the quadrat's paddock was grazed in the transition (0 in
## exclosures); climate = water-year precipitation anomaly: the difference from the long-term mean of the
## station record, in units of its long-term SD (sign kept: below 0 = drier than the long-term mean;
## as the rainfall anomalies of Calleja-Solanas et al. 2025).
## Recruitment uses quadrat cover in place of z and W (square root of cover, as Chu & Adler 2015;
## own cover scaled within each species, as same-species crowding).
## The other non-shrub species below the threshold enter only as neighbours (listed in the file).
## Uses 'genet_W', 'recruit_W', 'recruit_forb' (CLR_CrowdIndex.R), 'species', 'sp_class',
## 'census' and 'quad_info'.


# Settings ---------------------------------------------------------------------
  # Main model species: polygon species (with a size) with at least this many annual survival
  # records (the focal species of the kernel fit; next species below has 344)
main_min_survival <- 500
  # Climate covariate of the transition (water year July(t) to June(t+1))
clim_var <- "ppt_wy"
  # Reference of the climate anomaly: "record" = mean and SD of every complete water year of the
  # station record (1925-2006 at Dubois); "transitions" = the transition years of the study
clim_ref <- "record"

der_dir <- here("data", "data.derived", site_tag)
tab_dir <- here("result", "tables", site_tag)
stopifnot(exists("genet_W"), exists("recruit_W"), exists("recruit_forb"), exists("species"),
          exists("sp_class"), exists("census"), exists("quad_info"))
W_cols <- c("W_intra", "W_grass", "W_shrub", "W_forb")
C_cols <- c("C_intra", "C_grass", "C_shrub", "C_forb")

# Species sets -----------------------------------------------------------------
polygon_sp <- species$Species[sp_class[species$Species] != "forb"]
main_sp  <- species$Species[species$Species %in% polygon_sp & species$n_survival >= main_min_survival]
shrub_sp <- species$Species[species$Species %in% polygon_sp & !species$Species %in% main_sp &
                              species$growth_form == "shrub"]
fixed_sp <- setdiff(polygon_sp, c(main_sp, shrub_sp))
sp_tab <- species %>%
  filter(Species %in% main_sp) %>%
  arrange(desc(n_survival)) %>%
  transmute(k = row_number(), Species, code, growth_form, crowd_class = unname(sp_class[Species]),
            n_survival, seedling_info)
cat(sprintf("Main model: %d species (%s); other-shrub submodel: %d species (%s); %d species as neighbours only.\n",
            nrow(sp_tab), paste(sp_tab$code, collapse = ", "), length(shrub_sp),
            paste(species$code[match(shrub_sp, species$Species)], collapse = ", "), length(fixed_sp)))

# Index tables: years, quadrats, groups ----------------------------------------
  # Years: every annual transition t -> t+1 (recruitment table holds all of them)
yr_tab <- recruit_W %>%
  distinct(Year, clim = .data[[clim_var]]) %>%
  arrange(Year) %>%
  mutate(t = row_number(), .before = Year)
stopifnot(!any(duplicated(yr_tab$Year)), !anyNA(yr_tab$clim))
  # Long-term reference from the water-year table of CLR_ProcessDataPartA.R (complete water years)
stopifnot(exists("water_year"), clim_ref %in% c("record", "transitions"), clim_var == "ppt_wy")
ref_values <- if (clim_ref == "record") water_year$ppt_mm[!is.na(water_year$ppt_mm)] else yr_tab$clim
clim_scale <- c(mean = mean(ref_values), sd = sd(ref_values), n_years = length(ref_values))
yr_tab$clim_anom <- (yr_tab$clim - clim_scale[["mean"]]) / clim_scale[["sd"]]
cat(sprintf("Climate anomaly: reference %s (%d water years, mean %.0f mm, SD %.0f mm); transition years range %.2f to %.2f SD, mean %.2f.\n",
            clim_ref, length(ref_values), clim_scale[["mean"]], clim_scale[["sd"]],
            min(yr_tab$clim_anom), max(yr_tab$clim_anom), mean(yr_tab$clim_anom)))

grp_tab  <- distinct(recruit_W, Group) %>% arrange(Group) %>% mutate(g = row_number())
quad_tab <- distinct(recruit_W, Quad, Group) %>% arrange(Quad) %>%
  mutate(q = row_number()) %>% left_join(grp_tab, by = "Group")
stopifnot(!any(duplicated(quad_tab$Quad)))

add_index <- function(d) {
  d %>%
    left_join(select(yr_tab, Year, t, clim_anom), by = "Year") %>%
    left_join(select(quad_tab, Quad, q, g), by = "Quad") %>%
    mutate(E = as.integer(graze_total > 0))
}

# Size at t+1 and its measure (point records carry a fixed 0.25 cm2) -------------
next_rec <- genet %>%
  transmute(Quad, Species, trackID, Year = Year - 1L, size_measure_next = size_measure)

  # Survival and growth records of focal plants of a species set
vital_records <- function(sp_set) {
  genet_W %>%
    filter(Species %in% sp_set, focal_ok, !is.na(survives_tplus1)) %>%
    left_join(next_rec, by = c("Quad", "Species", "trackID", "Year")) %>%
    add_index() %>%
      # Growth: survivors with a mapped size at t and t+1 (point records have a fixed area)
    mutate(grow_ok = survives_tplus1 == 1 & !is.na(area_tplus1_cm2) &
             size_measure != "point" & size_measure_next %in% c("basal", "canopy"))
}

# Main model: survival and growth ----------------------------------------------
main_rec <- vital_records(main_sp)
  # Size: log area, centred and scaled over the focal survival records (same scale at t and t+1)
z_scale <- c(mean = mean(log(main_rec$area_cm2)), sd = sd(log(main_rec$area_cm2)))
scale_z <- function(a) (log(a) - z_scale[["mean"]]) / z_scale[["sd"]]
  # Crowding by other grasses, other shrubs and forbs: centred and scaled per class over the focal
  # survival records (the same neighbour area means the same for every focal species)
W_scale <- sapply(W_cols[-1], function(w) c(mean = mean(main_rec[[w]]), sd = sd(main_rec[[w]])))
  # Same-species crowding: centred and scaled within each species. Its units belong to the
  # species (conspecific basal area of a grass, conspecific canopy area of a shrub), so a pooled
  # scale gives one SD that the shrubs set: the grasses with few conspecific neighbours then span
  # a small part of a unit and the shrubs many units. Within each species, beta_intra and its
  # shared interactions (theta, omega, delta) are effects per SD of the crowding each species meets.
main_rec <- main_rec %>%
  left_join(select(sp_tab, Species, k), by = "Species")
W_intra_scale <- main_rec %>%
  group_by(k) %>%
  summarise(mean = mean(W_intra), sd = sd(W_intra), .groups = "drop") %>%
  mutate(code = sp_tab$code[k], .after = k)
main_rec <- main_rec %>%
  mutate(i = match(trackID, unique(trackID)),
         z = scale_z(area_cm2), z_next = ifelse(grow_ok, scale_z(area_tplus1_cm2), NA_real_),
         W_intra_z = (W_intra - W_intra_scale$mean[k]) / W_intra_scale$sd[k],
         across(all_of(W_cols[-1]), ~ (.x - W_scale["mean", cur_column()]) / W_scale["sd", cur_column()],
                .names = "{.col}_z")) %>%
  arrange(k, t, i)
cat("Same-species crowding (cm2) per species, mean and SD before scaling:\n")
print(as.data.frame(mutate(W_intra_scale, across(c(mean, sd), ~ round(.x, 2)))), row.names = FALSE)
stopifnot(!anyNA(main_rec[c("k", "t", "q", "g", "z", "E", paste0(W_cols, "_z"))]),
          all(is.finite(main_rec$z_next[main_rec$grow_ok])))
cat(sprintf("Main model: %d survival records, %d growth records, %d plants, %d years.\n",
            nrow(main_rec), sum(main_rec$grow_ok), max(main_rec$i), n_distinct(main_rec$t)))

# Main model: recruitment ------------------------------------------------------
  # One row per species x quadrat x transition (zeros included). Parent cover: quadrat cover and
  # group mean cover at t (the effective cover mixes them inside the model). Density dependence:
  # square root of the quadrat cover of other grasses and other shrubs and of the forb count,
  # centred and scaled (the intra term uses the effective cover inside the model, scaled within
  # each species with 'C_intra_scale', as same-species crowding).
sqrt_scale <- function(x) c(mean = mean(sqrt(x)), sd = sd(sqrt(x)))
main_recruit <- recruit_W %>%
  filter(Species %in% main_sp) %>%
  add_index() %>%
  left_join(select(sp_tab, Species, k), by = "Species")
C_intra_scale <- main_recruit %>%
  group_by(k) %>%
  summarise(mean = mean(sqrt(cover_t_cm2)), sd = sd(sqrt(cover_t_cm2)), .groups = "drop") %>%
  mutate(code = sp_tab$code[k], .after = k)
C_scale <- cbind(C_grass = sqrt_scale(main_recruit$other_grass_t_cm2),
                 C_shrub = sqrt_scale(main_recruit$other_shrub_t_cm2),
                 C_forb  = sqrt_scale(main_recruit$n_forb_t))
main_recruit <- main_recruit %>%
  mutate(C_grass_z = (sqrt(other_grass_t_cm2) - C_scale["mean", "C_grass"]) / C_scale["sd", "C_grass"],
         C_shrub_z = (sqrt(other_shrub_t_cm2) - C_scale["mean", "C_shrub"]) / C_scale["sd", "C_shrub"],
         C_forb_z  = (sqrt(n_forb_t)          - C_scale["mean", "C_forb"])  / C_scale["sd", "C_forb"]) %>%
  select(k, Species, Quad, q, g, Year, t, clim_anom, E, recruits = recruits_tplus1,
         cover_q = cover_t_cm2, cover_g = group_cover_t_cm2, C_grass_z, C_shrub_z, C_forb_z) %>%
  arrange(k, t, q)
stopifnot(nrow(main_recruit) == nrow(sp_tab) * nrow(filter(recruit_W, Species == main_sp[1])),
          !anyNA(main_recruit))
  # Recruits without parents in the whole group (seed bank or arrival from outside the quadrats)
no_parent <- main_recruit$cover_g == 0 & main_recruit$recruits > 0
cat(sprintf("Main recruitment: %d rows, %d recruits; %d rows (%d recruits) with recruits and no parent cover in the group.\n",
            nrow(main_recruit), sum(main_recruit$recruits), sum(no_parent), sum(main_recruit$recruits[no_parent])))

  # Sizes of recruits at t+1 (recruits of annual transitions)
recruit_size <- function(sp_set) {
  genet %>%
    filter(Species %in% sp_set, recruit %in% 1) %>%
    mutate(Year = Year - 1L) %>%
    semi_join(recruit_W, by = c("Quad", "Year")) %>%
    transmute(Species, Quad, Year, z_rec = scale_z(area_cm2), point = size_measure == "point")
}
main_recsize <- recruit_size(main_sp) %>%
  left_join(select(sp_tab, Species, k), by = "Species") %>%
  arrange(k)
stopifnot(nrow(main_recsize) == sum(main_recruit$recruits))

# Forb submodel ----------------------------------------------------------------
  # Survival of forb plants (no size); W_forb = other forb plants (forb kernel), W_grass and
  # W_shrub from the main kernel; scaled on their own over the forb focal records
forb_rec <- genet_W %>%
  filter(crowd_class == "forb", focal_ok, !is.na(survives_tplus1)) %>%
  add_index()
W_scale_forb <- sapply(W_cols[-1], function(w) c(mean = mean(forb_rec[[w]]), sd = sd(forb_rec[[w]])))
forb_rec <- forb_rec %>%
  mutate(across(all_of(W_cols[-1]), ~ (.x - W_scale_forb["mean", cur_column()]) / W_scale_forb["sd", cur_column()],
                .names = "{.col}_z")) %>%
  select(Species, Quad, q, g, Year, t, clim_anom, E, trackID, survives = survives_tplus1,
         W_forb, W_grass, W_shrub, W_forb_z, W_grass_z, W_shrub_z)
stopifnot(!anyNA(forb_rec))
forb_recruit <- recruit_forb %>%
  add_index() %>%
  mutate(C_forb_z  = (sqrt(n_forb_t)          - C_scale["mean", "C_forb"])  / C_scale["sd", "C_forb"],
         C_grass_z = (sqrt(other_grass_t_cm2) - C_scale["mean", "C_grass"]) / C_scale["sd", "C_grass"],
         C_shrub_z = (sqrt(other_shrub_t_cm2) - C_scale["mean", "C_shrub"]) / C_scale["sd", "C_shrub"]) %>%
  select(Quad, q, g, Year, t, clim_anom, E, recruits = recruits_tplus1, n_q = n_forb_t, n_g = group_n_forb_t,
         C_forb_z, C_grass_z, C_shrub_z)
stopifnot(!anyNA(forb_recruit))
cat(sprintf("Forb submodel: %d survival records, %d recruitment rows, %d recruits.\n",
            nrow(forb_rec), nrow(forb_recruit), sum(forb_recruit$recruits)))

# Other-shrub submodel ---------------------------------------------------------
  # Shrubs below the threshold pooled as one population; size and crowding on the main scales,
  # same-species crowding scaled within the group (as within each main species)
shrub_rec <- vital_records(shrub_sp)
W_intra_scale <- bind_rows(W_intra_scale, tibble(k = NA_integer_, code = "other shrubs",
                                                 mean = mean(shrub_rec$W_intra), sd = sd(shrub_rec$W_intra)))
shrub_rec <- shrub_rec %>%
  mutate(z = scale_z(area_cm2), z_next = ifelse(grow_ok, scale_z(area_tplus1_cm2), NA_real_),
         W_intra_z = (W_intra - W_intra_scale$mean[nrow(W_intra_scale)]) / W_intra_scale$sd[nrow(W_intra_scale)],
         across(all_of(W_cols[-1]), ~ (.x - W_scale["mean", cur_column()]) / W_scale["sd", cur_column()],
                .names = "{.col}_z"),
         i = match(trackID, unique(trackID))) %>%
  select(Species, Quad, q, g, Year, t, clim_anom, E, trackID, i, z, survives = survives_tplus1, grow_ok, z_next,
         all_of(W_cols), all_of(paste0(W_cols, "_z")))
  # Recruitment: pooled cover of the group's species (quadrat and group mean) as parents;
  # the group's own cover leaves the shrub total
shrub_q <- genet %>% filter(Species %in% shrub_sp) %>%
  group_by(Quad, Year) %>% summarise(cover_q = sum(area_cm2), .groups = "drop")
shrub_g <- census %>%
  select(Quad, Year) %>%
  left_join(select(quad_info, Quad, Group), by = "Quad") %>%
  left_join(shrub_q, by = c("Quad", "Year")) %>%
  mutate(cover_q = replace_na(cover_q, 0)) %>%
  group_by(Group, Year) %>% summarise(cover_g = mean(cover_q), .groups = "drop")
shrub_recruit <- recruit_W %>%
  filter(Species %in% shrub_sp) %>%
  group_by(Quad, Year) %>%
  summarise(recruits = sum(recruits_tplus1), other_grass_t_cm2 = first(other_grass_t_cm2),
            shrub_all = first(other_shrub_t_cm2) + first(cover_t_cm2), n_forb_t = first(n_forb_t),
            graze_total = first(graze_total), .groups = "drop") %>%
  left_join(shrub_q, by = c("Quad", "Year")) %>%
  mutate(cover_q = replace_na(cover_q, 0)) %>%
  left_join(select(quad_tab, Quad, Group), by = "Quad") %>%
  left_join(shrub_g, by = c("Group", "Year")) %>%
  add_index() %>%
  mutate(other_shrub = shrub_all - cover_q,
         C_grass_z = (sqrt(other_grass_t_cm2) - C_scale["mean", "C_grass"]) / C_scale["sd", "C_grass"],
         C_shrub_z = (sqrt(pmax(other_shrub, 0)) - C_scale["mean", "C_shrub"]) / C_scale["sd", "C_shrub"],
         C_forb_z  = (sqrt(n_forb_t)          - C_scale["mean", "C_forb"])  / C_scale["sd", "C_forb"]) %>%
  select(Quad, q, g, Year, t, clim_anom, E, recruits, cover_q, cover_g, C_grass_z, C_shrub_z, C_forb_z) %>%
  arrange(t, q)
stopifnot(!anyNA(shrub_recruit), all(shrub_recruit$cover_g > 0 | shrub_recruit$cover_q == 0))
C_intra_scale <- bind_rows(C_intra_scale, tibble(k = NA_integer_, code = "other shrubs",
                                                 mean = mean(sqrt(shrub_recruit$cover_q)), sd = sd(sqrt(shrub_recruit$cover_q))))
shrub_recsize <- recruit_size(shrub_sp)
stopifnot(nrow(shrub_recsize) == sum(shrub_recruit$recruits))
cat(sprintf("Other-shrub submodel: %d survival records, %d growth records, %d recruits.\n",
            nrow(shrub_rec), sum(shrub_rec$grow_ok), sum(shrub_recruit$recruits)))

# Save -------------------------------------------------------------------------
partA_input <- list(
  main  = list(rec = main_rec, recruit = main_recruit, recsize = main_recsize),
  forb  = list(rec = forb_rec, recruit = forb_recruit),
  shrub = list(rec = shrub_rec, recruit = shrub_recruit, recsize = shrub_recsize),
  species = sp_tab, years = yr_tab, quadrats = quad_tab, groups = grp_tab,
  species_sets = list(main = main_sp, shrub = shrub_sp, forb = species$Species[sp_class[species$Species] == "forb"],
                      neighbour_only = fixed_sp),
  scaling = list(z = z_scale, W = W_scale, W_intra = W_intra_scale, W_forb_model = W_scale_forb, C = C_scale,
                 C_intra = C_intra_scale, clim = clim_scale,
                 clim_var = clim_var, clim_ref = clim_ref),
  settings = list(main_min_survival = main_min_survival, R_main = R_main, R_forb = R_forb)
)
saveRDS(partA_input, file.path(der_dir, "CLR_PartA_input.rds"))

  # Sample sizes per species (main model) and per submodel
input_n <- bind_rows(
  main_rec %>% group_by(Species) %>%
    summarise(n_survival = n(), n_growth = sum(grow_ok), n_plant = n_distinct(trackID), .groups = "drop") %>%
    left_join(main_recruit %>% group_by(Species) %>% summarise(n_recruit = sum(recruits), .groups = "drop"), by = "Species") %>%
    mutate(model = "main"),
  tibble(Species = "forb (pooled)", n_survival = nrow(forb_rec), n_growth = 0L,
         n_plant = n_distinct(forb_rec$trackID), n_recruit = sum(forb_recruit$recruits), model = "forb"),
  tibble(Species = "other shrubs (pooled)", n_survival = nrow(shrub_rec), n_growth = sum(shrub_rec$grow_ok),
         n_plant = n_distinct(shrub_rec$trackID), n_recruit = sum(shrub_recruit$recruits), model = "shrub")) %>%
  mutate(code = species$code[match(Species, species$Species)], .after = Species)
write.csv(input_n, file.path(tab_dir, "CLR_PartA_input_n.csv"), row.names = FALSE)
print(as.data.frame(input_n))
cat("Saved:\n",
    " data/data.derived/", site_tag, ": CLR_PartA_input.rds\n",
    " result/tables/", site_tag, ": CLR_PartA_input_n.csv\n", sep = "")
rm(list = intersect(c("polygon_sp", "add_index", "next_rec", "vital_records", "scale_z", "sqrt_scale",
                      "no_parent", "recruit_size", "shrub_q", "shrub_g", "W_cols", "C_cols", "ref_values"), ls()))
