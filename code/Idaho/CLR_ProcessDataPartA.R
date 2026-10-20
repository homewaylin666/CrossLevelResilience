################################################################################
## Build the Part A data of the Idaho site (vital-rate model)
##
## Written by Hung-wei Lin Oct 2026
##
## Starts from the tracked records (CLR_TrackGenet.R) and builds three tables
## (log size, circle radius and climate z-scores are left to the model step):
##   (1) genet: one row per plant (genet) per census year, with size, fate at the
##       next census, location, quadrat structure, grazing and climate;
##   (2) recruit: one row per species, quadrat and annual transition, with the
##       recruits at t+1 and the parent cover at t (quadrat and group);
##   (3) species: one row per species, with growth form, size measure and sample sizes.
## Crowding (INTRA, INTER_grass, INTER_shrub) is added in a later step.
## Uses 'site_dir' and 'site_tag' from the runscript.


# Settings ---------------------------------------------------------------------
  # Focal plants need their neighbourhood inside the quadrat. A plant is a focal plant
  # when its outline stays at least this far from the quadrat border and its circle
  # (centroid, radius from area) stays inside. Revisit once the kernel range R is set.
focal_edge_cm <- 5
  # Climate window of a transition t -> t+1 (census in late May to early July):
  # water year July(t) to June(t+1); the lagged window is July(t-1) to June(t).
wy_start_month <- 7

der_dir <- here("data", "data.derived", site_tag)
tab_dir <- here("result", "tables", site_tag)
src <- function(f) here("data", "data.source", site_dir, f)

# Read the tracked records -----------------------------------------------------
  # One row per polygon or point (ramet), with its trackID (genet)
ramet <- readRDS(file.path(der_dir, "CLR_TrackGenet_ramet.rds"))
cent  <- st_coordinates(st_centroid(st_geometry(ramet)))
  # Outline extent of every ramet (for its distance to the quadrat border)
bb <- as.data.frame(st_coordinates(st_geometry(ramet))) %>%
  mutate(id = .data[[tail(grep("^L", names(.), value = TRUE), 1)]]) %>%
  group_by(id) %>%
  summarise(xmin = min(X), xmax = max(X), ymin = min(Y), ymax = max(Y), .groups = "drop") %>%
  arrange(id)
stopifnot(nrow(bb) == nrow(ramet))
ramet_df <- st_drop_geometry(ramet) %>%
  mutate(cx = cent[, 1], cy = cent[, 2],
         edge_dist = pmin(bb$xmin, 1 - bb$xmax, bb$ymin, 1 - bb$ymax))

species_raw <- read.csv(src("species_list.csv"), na.strings = c("NA", ""))
growth_form <- setNames(species_raw$growthForm, species_raw$species)

# Genet table: one row per genet per census year -------------------------------
  # Size = sum of the ramet areas; location = area-weighted centroid of the ramets
genet <- ramet_df %>%
  group_by(Site, Quad, Year, Species, trackID) %>%
  summarise(n_ramet = n(),
            n_point = sum(record_type == "point"),
            seedling = any(seedling %in% "Y"),
            area_m2 = sum(basalArea_ramet),
            x = sum(cx * basalArea_ramet) / area_m2,
            y = sum(cy * basalArea_ramet) / area_m2,
            edge_dist = min(edge_dist),
            age = max(age),
            recruit = max(recruit),
            survives_tplus1 = max(survives_tplus1),
            fate_contested = any(fate_contested),
            fate_merged = any(fate_merged),
            .groups = "drop")

  # Size measure: shrub polygons = canopy, other polygons = basal,
  # genets made of points only = assigned point area (0.25 cm2, no real size)
genet <- genet %>%
  mutate(growth_form = growth_form[Species],
         size_measure = case_when(n_point == n_ramet ~ "point",
                                  growth_form == "shrub" ~ "canopy",
                                  TRUE ~ "basal"),
         area_cm2 = area_m2 * 1e4,
         x_cm = x * 100, y_cm = y * 100,
         edge_dist_cm = edge_dist * 100,
           # Circle of the genet (radius from area) must stay inside the quadrat as well
         focal_ok = edge_dist_cm >= focal_edge_cm &
           sqrt(area_cm2 / pi) <= pmin(x_cm, 100 - x_cm, y_cm, 100 - y_cm)) %>%
  select(-x, -y, -edge_dist)

sp <- sort(unique(genet$Species))

# Fate at the next census -------------------------------------------------------
  # A transition t -> t+1 counts only when the quadrat was mapped in both years.
  # plantTracker links plants across census gaps up to 'dorm' years, so survival and
  # size at t+1 are kept for annual transitions only.
inventory_raw <- read.csv(src("quad_inventory.csv"))
census <- inventory_raw %>%
  pivot_longer(starts_with("Q"), names_to = "Quad", values_to = "yy") %>%
  filter(!is.na(yy)) %>%
  transmute(Quad, Year = yy + 1900L) %>%
  arrange(Quad, Year) %>%
  group_by(Quad) %>%
  mutate(next_census = lead(Year), first_census = Year == min(Year)) %>%
  ungroup() %>%
  mutate(annual_transition = !is.na(next_census) & next_census == Year + 1L)

next_size <- genet %>% transmute(Quad, Species, trackID, Year = Year - 1L, area_tplus1_cm2 = area_cm2)
genet <- genet %>%
  left_join(census, by = c("Quad", "Year")) %>%
  left_join(next_size, by = c("Quad", "Species", "trackID", "Year")) %>%
  mutate(survives_tplus1 = ifelse(annual_transition, survives_tplus1, NA),
         area_tplus1_cm2 = ifelse(annual_transition, area_tplus1_cm2, NA),
           # Flags from tracking (annual transitions only): death decided against a neighbour
           # within 5 cm, and coalescence (fate censored for basal plants)
         fate_contested = annual_transition & fate_contested,
         fate_merged = annual_transition & fate_merged)
stopifnot(!any(genet$survives_tplus1 %in% 0 & !is.na(genet$area_tplus1_cm2)))

# Quadrat structure and grazing ------------------------------------------------
  # Group = the quadrat name without its last part ("P1 E1 Q2" -> "P1 E1");
  # grazing yes/no and paddock as given in quad_info
quad_info <- read.csv(src("quad_info.csv"), na.strings = c("NA", "")) %>%
  transmute(Quad = shapefiles,
            Group = sub("\\s*Q\\d+$", "", quadrat),
            grazed_quad = grazing, paddock)

  # Stocking rates (sheep-days per acre, metadata) per paddock and year;
  # 'yearRange' holds two 2-digit years with a Windows-1252 dash between them (read as latin1)
grazing_raw <- read.csv(src("grazing_info.csv"), na.strings = c("NA", ""), fileEncoding = "latin1")
grazing_raw$from <- as.integer(sapply(regmatches(grazing_raw$yearRange, gregexpr("\\d+", grazing_raw$yearRange)), `[`, 1)) + 1900L
grazing_raw$to   <- as.integer(sapply(regmatches(grazing_raw$yearRange, gregexpr("\\d+", grazing_raw$yearRange)), `[`, 2)) + 1900L
graze_year <- grazing_raw %>%
  rowwise() %>%
  mutate(Year = list(from:to)) %>%
  unnest(Year) %>%
  transmute(paddock, Year, timing,
            spring = replace_na(springRates, 0), fall = replace_na(fallRates, 0))
stopifnot(!any(duplicated(graze_year[c("paddock", "Year")])))

  # Grazing of a transition t -> t+1: fall grazing after the census in year t,
  # and spring grazing before the census in year t+1. Exclosure quadrats get 0.
graze_trans <- quad_info %>%
  select(Quad, paddock) %>%
  crossing(Year = sort(unique(census$Year))) %>%
  left_join(select(graze_year, paddock, Year, timing_t = timing, graze_fall_t = fall),
            by = c("paddock", "Year")) %>%
  left_join(transmute(graze_year, paddock, Year = Year - 1L, graze_spring_t1 = spring),
            by = c("paddock", "Year")) %>%
  mutate(timing_t = ifelse(is.na(paddock), "Exclosure", timing_t),
         across(c(graze_fall_t, graze_spring_t1), ~ replace_na(.x, 0)),
         graze_total = graze_fall_t + graze_spring_t1) %>%
  select(-paddock)

genet <- genet %>%
  left_join(quad_info, by = "Quad") %>%
  left_join(graze_trans, by = c("Quad", "Year"))

# Climate ----------------------------------------------------------------------
  # Water-year precipitation (mm) and mean temperature (degC), July to June
read_month <- function(f) {
  read.csv(src(f), na.strings = c("NA", ""), check.names = FALSE) %>%
    rename_with(~ trimws(toupper(.x))) %>%
    select(-ANNUAL) %>%
    pivot_longer(-YEAR, names_to = "month_name", values_to = "value") %>%
    mutate(month = match(month_name, toupper(month.abb))) %>%
    rename(cal_year = YEAR)
}
monthly <- full_join(read_month("total_monthly_ppt.csv") %>% transmute(cal_year, month, ppt_mm = value * 25.4),
                     read_month("monthly_mean_temp.csv") %>% transmute(cal_year, month, tmean_c = (value - 32) * 5 / 9),
                     by = c("cal_year", "month")) %>%
    # Water year named after the year it ends in (July 1944 to June 1945 = 1945)
  mutate(wy = ifelse(month >= wy_start_month, cal_year + 1L, cal_year))
  # Complete water years only (12 months with values)
water_year <- monthly %>%
  group_by(wy) %>%
  summarise(ppt_mm = if (sum(!is.na(ppt_mm)) == 12) sum(ppt_mm) else NA,
            tmean_c = if (sum(!is.na(tmean_c)) == 12) mean(tmean_c) else NA,
            .groups = "drop")

  # Transition t -> t+1: main window = water year t+1, lagged window = water year t
climate <- data.frame(Year = sort(unique(census$Year))) %>%
  left_join(rename(water_year, ppt_wy = ppt_mm, tmean_wy = tmean_c) %>% mutate(Year = wy - 1L) %>% select(-wy), by = "Year") %>%
  left_join(rename(water_year, ppt_wy_lag = ppt_mm, tmean_wy_lag = tmean_c) %>% mutate(Year = wy) %>% select(-wy), by = "Year")
genet <- left_join(genet, climate, by = "Year")

# Final genet table ------------------------------------------------------------
genet <- genet %>%
  select(Site, Group, Quad, grazed_quad, paddock, Year, next_census,
         Species, growth_form, size_measure, trackID, n_ramet, n_point, seedling, age, recruit,
         area_cm2, survives_tplus1, fate_contested, fate_merged, area_tplus1_cm2,
         x_cm, y_cm, edge_dist_cm, focal_ok,
         timing_t, graze_fall_t, graze_spring_t1, graze_total,
         ppt_wy, tmean_wy, ppt_wy_lag, tmean_wy_lag) %>%
  arrange(Quad, Species, Year, trackID)
stopifnot(!any(duplicated(genet[c("Quad", "Species", "Year", "trackID")])))
cat(sprintf("Genet table: %d genet-years, %d species; %d with an annual survival record, %d with size at t+1.\n",
            nrow(genet), n_distinct(genet$Species), sum(!is.na(genet$survives_tplus1)),
            sum(!is.na(genet$area_tplus1_cm2))))

# Recruitment table: species x quadrat x annual transition ---------------------
  # Recruits at t+1 (new genets) against the parent cover at t in the quadrat and
  # in its group (mean over the group's quadrats mapped in year t). Zeros included.
trans <- census %>% filter(annual_transition) %>% select(Quad, Year) %>% left_join(quad_info, by = "Quad")
parent <- genet %>%
  group_by(Species, Quad, Year) %>%
  summarise(cover_t_cm2 = sum(area_cm2), n_genet_t = n(), .groups = "drop")
recruits <- genet %>%
  filter(recruit %in% 1) %>%
  group_by(Species, Quad, Year = Year - 1L) %>%
  summarise(recruits_tplus1 = n(), recruit_area_cm2 = sum(area_cm2), .groups = "drop")
  # Group mean cover at t uses every quadrat of the group mapped in year t
group_cover <- census %>%
  select(Quad, Year) %>%
  left_join(select(quad_info, Quad, Group), by = "Quad") %>%
  crossing(Species = sp) %>%
  left_join(parent, by = c("Species", "Quad", "Year")) %>%
  mutate(cover_t_cm2 = replace_na(cover_t_cm2, 0)) %>%
  group_by(Species, Group, Year) %>%
  summarise(group_cover_t_cm2 = mean(cover_t_cm2), n_quad_group_t = n(), .groups = "drop")

recruit <- trans %>%
  crossing(Species = sp) %>%
  left_join(parent, by = c("Species", "Quad", "Year")) %>%
  left_join(recruits, by = c("Species", "Quad", "Year")) %>%
  left_join(group_cover, by = c("Species", "Group", "Year")) %>%
  mutate(across(c(cover_t_cm2, n_genet_t, recruits_tplus1, recruit_area_cm2), ~ replace_na(.x, 0)),
         growth_form = growth_form[Species]) %>%
  left_join(graze_trans, by = c("Quad", "Year")) %>%
  left_join(climate, by = "Year") %>%
  relocate(Species, growth_form, .before = Quad) %>%
  arrange(Species, Quad, Year)
  # Every recruit of an annual transition should land in this table
stopifnot(sum(recruit$recruits_tplus1) ==
            sum(genet$recruit %in% 1 & paste(genet$Quad, genet$Year - 1L) %in% paste(trans$Quad, trans$Year)))
cat(sprintf("Recruitment table: %d rows (%d species x %d quadrat transitions), %d recruits.\n",
            nrow(recruit), length(sp), nrow(trans), sum(recruit$recruits_tplus1)))

# Species table ----------------------------------------------------------------
  # Short species codes: two letters of the genus + two of the epithet (PSSP, ARTR, ...),
  # with a number added where two species share a code
code <- toupper(paste0(substr(sp, 1, 2), substr(sub("^\\S+\\s*", "", sp), 1, 2)))
code <- setNames(make.unique(gsub("[^A-Z]", "X", code), sep = ""), sp)
species <- genet %>%
  mutate(code = code[Species]) %>%
  group_by(Species, code, growth_form) %>%
  summarise(size_measure = paste(sort(unique(size_measure)), collapse = "+"),
            n_genet_year = n(),
            n_survival = sum(!is.na(survives_tplus1)),
            n_survival_focal = sum(!is.na(survives_tplus1) & focal_ok),
            n_growth = sum(!is.na(area_tplus1_cm2) & size_measure != "point"),
            n_recruit = sum(recruit %in% 1),
            n_quad = n_distinct(Quad),
              # Seedling flags exist for this species (point records, or Chu's grass flags):
              # recruits near a parent are recognised; species without them (sensitivity test)
            seedling_info = any(seedling),
            .groups = "drop") %>%
  arrange(desc(n_survival))

# Save -------------------------------------------------------------------------
write.csv(genet,   file.path(der_dir, "CLR_PartA_genet.csv"), row.names = FALSE)
write.csv(recruit, file.path(der_dir, "CLR_PartA_recruit.csv"), row.names = FALSE)
write.csv(species, file.path(tab_dir, "CLR_PartA_species.csv"), row.names = FALSE)
write.csv(climate, file.path(tab_dir, "CLR_PartA_climate.csv"), row.names = FALSE)
print(as.data.frame(head(species, 20)))
