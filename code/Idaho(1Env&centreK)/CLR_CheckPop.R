################################################################################
## Check the raw population data of the Idaho site (mapped quadrats)
##
## Written by Hung-wei Lin Oct 2026
##
## Reads the polygon records (plants with an area = size), the point records
## (plants with a location only = no size) and the counts of unmapped plants,
## checks the records against the sampling inventory, and shows the abundance
## of every species in every year. Also compares polygon cover with the number
## of point plants, to judge how much of the community the size species carry.
## Raw records stay as they are (flags only).
## Uses 'site_dir' and 'site_tag' from the runscript.

# Read the raw csv files -------------------------------------------------------
pop_path <- function(f) here("data", "data.source", site_dir, f)

  # Polygons: basal cover of grasses, canopy cover of shrubs; area in m2
cover_raw   <- read.csv(pop_path("allrecords_cover.csv"), na.strings = c("NA", "null"))
  # Points: forbs, single grass tillers and small young shrubs (without canopy);
  # shrubs with a canopy polygon also have a stem point (stem = "Y") linked to it
density_raw <- read.csv(pop_path("allrecords_density.csv"), na.strings = c("NA", "null"))
  # Counts per quadrat and year (plants counted, without a map)
counts_raw  <- read.csv(pop_path("annuals_counts_v3.csv"), na.strings = c("NA", ""))
species_raw <- read.csv(pop_path("species_list.csv"), na.strings = c("NA", ""))
  # Inventory: one row per census year, one column per quadrat (year if sampled, else NA)
inventory_raw <- read.csv(pop_path("quad_inventory.csv"), na.strings = c("NA", ""))

stopifnot(all(c("quad", "year", "species", "area", "x", "y") %in% names(cover_raw)))
stopifnot(all(c("quad", "year", "species", "seedling", "stem", "x", "y") %in% names(density_raw)))
stopifnot(all(c("quad", "year", "species", "count") %in% names(counts_raw)))
dim(cover_raw)   # one row per polygon
dim(density_raw) # one row per point

  # Stem points duplicate the shrub canopy polygons they belong to, so they leave the
  # point records here (as the data paper metadata advises: remove stem = "Y" before analysis)
n_stem <- sum(density_raw$stem %in% "Y")
density_raw <- density_raw[!density_raw$stem %in% "Y", ]
cat(sprintf("Removed %d shrub stem points (stem = Y); %d point records left.\n", n_stem, nrow(density_raw)))

  # Years come as two digits (23 = 1923, 73 = 1973)
cover_raw$Year     <- cover_raw$year + 1900
density_raw$Year   <- density_raw$year + 1900
counts_raw$Year    <- counts_raw$year + 1900
inventory_raw$Year <- inventory_raw$year + 1900

# Sampling effort --------------------------------------------------------------
  # Abundance per year is divided by the number of quadrats sampled that year,
  # since the number of sampled quadrats changes a lot (4 in the first years, 26 later).
sampled <- inventory_raw %>%
  select(-year) %>%
  pivot_longer(-Year, names_to = "quad", values_to = "flag") %>%
  filter(!is.na(flag)) %>%
  select(quad, Year)
effort <- count(sampled, Year, name = "n_quad")
cat(sprintf("Census years: %d (%d to %d); quadrats per year: %d to %d.\n",
            nrow(effort), min(effort$Year), max(effort$Year), min(effort$n_quad), max(effort$n_quad)))

# Record checks ----------------------------------------------------------------
  # (1) every quadrat-year with records should appear in the inventory,
  # (2) areas should be positive, (3) coordinates should sit inside the 1 x 1 m quadrat,
  # (4) exact duplicates (same quadrat, year, species, position and area),
  # (5) species missing from the species list.
flag_list <- list()
add_flag <- function(d, source, check, note) {
  if (nrow(d) == 0) return(NULL)
  data.frame(source = source, quad = d$quad, Year = d$Year, species = d$species,
             check = check, note = note)
}
recorded <- bind_rows(distinct(cover_raw, quad, Year), distinct(density_raw, quad, Year),
                      distinct(counts_raw, quad, Year)) %>% distinct()
off_inventory <- anti_join(recorded, sampled, by = c("quad", "Year"))
if (nrow(off_inventory) > 0) {
  flag_list[["inventory"]] <- data.frame(source = "all", quad = off_inventory$quad,
    Year = off_inventory$Year, species = NA, check = "outside_inventory",
    note = "quadrat-year with records, absent from quad_inventory")
}
flag_list[["area"]] <- add_flag(filter(cover_raw, !(area > 0)), "cover", "area", "area missing or <= 0")
flag_list[["xy cover"]] <- add_flag(filter(cover_raw, x < 0 | x > 1 | y < 0 | y > 1),
                                    "cover", "position", "centroid outside the quadrat")
flag_list[["xy density"]] <- add_flag(filter(density_raw, x < 0 | x > 1 | y < 0 | y > 1),
                                      "density", "position", "point outside the quadrat")
flag_list[["dup cover"]] <- add_flag(
  cover_raw[duplicated(cover_raw[, c("quad", "Year", "species", "x", "y", "area")]), ],
  "cover", "duplicate", "same quadrat, year, species, position and area")
flag_list[["dup density"]] <- add_flag(
  density_raw[duplicated(density_raw[, c("quad", "Year", "species", "x", "y")]), ],
  "density", "duplicate", "same quadrat, year, species and position")
unlisted <- setdiff(unique(c(cover_raw$species, density_raw$species)), species_raw$species)
if (length(unlisted) > 0) {
  flag_list[["unlisted"]] <- data.frame(source = "cover/density", quad = NA, Year = NA,
    species = unlisted, check = "unlisted_species", note = "absent from species_list.csv")
}
pop_flags <- bind_rows(flag_list)
cat(sprintf("Flagged %d records across %d checks.\n", nrow(pop_flags), length(unique(pop_flags$check))))
if (nrow(pop_flags) > 0) table(pop_flags$check)

# Data type of every species ---------------------------------------------------
  # size (basal)  = polygon of a grass or forb (basal area)
  # size (canopy) = polygon of a shrub (canopy area, a different size measure)
  # no size       = point record (location only)
  # count only    = plants counted per quadrat, without a map
  # A species can carry two types (shrubs: canopy polygons + stem points), so the
  # unit of analysis below is species x type.
type_levels <- c("size (basal)", "size (canopy)", "no size", "count only")
type_colour <- c("size (basal)" = "#1e6b2a", "size (canopy)" = "#b0641a",
                 "no size" = "#1560b8", "count only" = "grey45")
growth_form <- setNames(species_raw$growthForm, species_raw$species)

  # Polygons: cover in percent of the sampled area (area m2 / n quadrats x 100)
abund_cover <- cover_raw %>%
  mutate(type = ifelse(growth_form[species] %in% "shrub", "size (canopy)", "size (basal)")) %>%
  group_by(species, type, Year) %>%
  summarise(n_records = n(), total = sum(area), .groups = "drop") %>%
  left_join(effort, by = "Year") %>%
  mutate(abundance = 100 * total / n_quad, unit = "cover (%)")

  # Points: plants per m2 (shrub stem points already removed above)
abund_point <- density_raw %>%
  group_by(species, Year) %>%
  summarise(n_records = n(), total = n(), .groups = "drop") %>%
  left_join(effort, by = "Year") %>%
  mutate(type = "no size", abundance = total / n_quad, unit = "plants per m2")

  # Counts: plants per m2
abund_count <- counts_raw %>%
  group_by(species, Year) %>%
  summarise(n_records = n(), total = sum(count, na.rm = TRUE), .groups = "drop") %>%
  left_join(effort, by = "Year") %>%
  mutate(type = "count only", abundance = total / n_quad, unit = "plants per m2")

  # Fill in zeros for census years where a species x type has no record,
  # so a disappearance shows as a drop to zero (instead of a line bridging the gap)
pop_abundance <- bind_rows(abund_cover, abund_point, abund_count) %>%
  select(species, type, unit, Year, n_records, total, abundance) %>%
  group_by(species, type, unit) %>%
  complete(Year = effort$Year, fill = list(n_records = 0, total = 0, abundance = 0)) %>%
  ungroup() %>%
  mutate(type = factor(type, levels = type_levels))
stopifnot(!anyNA(pop_abundance$abundance))

  # One row per species x type: how much and how often
pop_species <- pop_abundance %>%
  group_by(species, type, unit) %>%
  summarise(growth_form = growth_form[species[1]],
            n_years_present = sum(n_records > 0),
            first_year = min(Year[n_records > 0]), last_year = max(Year[n_records > 0]),
            n_records = sum(n_records),
            mean_abundance = mean(abundance), max_abundance = max(abundance),
            .groups = "drop") %>%
  left_join(distinct(cover_raw, species, quad) %>% count(species, name = "n_quad_polygon"),
            by = "species") %>%
  left_join(distinct(density_raw, species, quad) %>% count(species, name = "n_quad_point"),
            by = "species") %>%
  mutate(across(c(n_quad_polygon, n_quad_point), ~ replace_na(.x, 0L))) %>%
  arrange(type, desc(n_records))
table(pop_species$type) # species x type per data type

# Polygon cover against point plants -------------------------------------------
  # Point plants have no area, so their share of the total cover is unknown.
  # Here every point plant (forbs and grass tillers, shrub stems left out) gets an assumed
  # basal area, and the share is computed per quadrat-year; the share changes with
  # the assumption, so several areas are shown side by side.
  # The empirical anchors: single tillers mapped as points were given 0.25 cm2 in the
  # data paper; Antennaria polygons in this data set have a median of a few cm2.
assumed_cm2 <- c(0.25, 1, 2, 5)
basal_qy <- cover_raw %>%
  filter(!growth_form[species] %in% "shrub") %>%
  group_by(quad, Year) %>%
  summarise(basal_cm2 = sum(area) * 1e4, .groups = "drop")
point_qy <- density_raw %>%
  filter(!growth_form[species] %in% "shrub") %>%
  count(quad, Year, name = "n_point")
pointshare <- sampled %>%
  left_join(basal_qy, by = c("quad", "Year")) %>%
  left_join(point_qy, by = c("quad", "Year")) %>%
  mutate(across(c(basal_cm2, n_point), ~ replace_na(.x, 0))) %>%
  crossing(assumed_cm2 = assumed_cm2) %>%
  mutate(point_cm2 = n_point * assumed_cm2,
         point_share = ifelse(point_cm2 + basal_cm2 > 0, point_cm2 / (point_cm2 + basal_cm2), NA))
pointshare_summary <- pointshare %>%
  group_by(assumed_cm2) %>%
  summarise(median_share = median(point_share, na.rm = TRUE),
            q90_share = quantile(point_share, 0.9, na.rm = TRUE), .groups = "drop")
print(pointshare_summary)

# Save tables ------------------------------------------------------------------
tab_dir <- here("result", "tables", site_tag)
fig_dir <- here("result", "figs", site_tag)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(pop_species, file.path(tab_dir, "CLR_CheckPop_species.csv"), row.names = FALSE)
write.csv(pop_abundance, file.path(tab_dir, "CLR_CheckPop_abundance.csv"), row.names = FALSE)
write.csv(effort, file.path(tab_dir, "CLR_CheckPop_effort.csv"), row.names = FALSE)
write.csv(pop_flags, file.path(tab_dir, "CLR_CheckPop_flags.csv"), row.names = FALSE)
write.csv(pointshare_summary, file.path(tab_dir, "CLR_CheckPop_pointshare.csv"), row.names = FALSE)

# Figures ----------------------------------------------------------------------
  # (a) Overview heat map: every species x type (rows) by census year (columns).
  # Colour = data type; shade = abundance relative to the species' own maximum,
  # so rare and common species are equally readable. Rows ordered by total records.
heat_data <- pop_abundance %>%
  group_by(species, type) %>%
  mutate(rel = if (max(abundance) > 0) abundance / max(abundance) else 0,
         n_total = sum(n_records)) %>%
  ungroup() %>%
  mutate(row = paste0(species, ifelse(type == "size (canopy)" | (type == "no size" &
                        growth_form[species] %in% "shrub"), paste0(" [", type, "]"), "")),
         # year labels carry the number of sampled quadrats, e.g. "1923 (4)"
         Year_f = factor(Year, levels = effort$Year,
                         labels = paste0(effort$Year, " (", effort$n_quad, ")")))
row_order <- heat_data %>% distinct(row, type, n_total) %>% arrange(type, n_total) %>% pull(row)
heat_data$row <- factor(heat_data$row, levels = unique(row_order))
p_heat <- ggplot(heat_data, aes(x = Year_f, y = row, fill = type, alpha = rel)) +
  geom_tile(colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = type_colour, name = "Data type") +
  scale_alpha_continuous(range = c(0, 1), name = "Abundance\n(relative to\nspecies max)") +
  facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Census year (number of quadrats sampled)", y = NULL,
       title = paste0(site_tag, ": abundance of every species in every census year")) +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text.y = element_text(angle = 0), panel.grid = element_blank())
ggsave(file.path(fig_dir, "CLR_CheckPop_heatmap.png"), p_heat,
       width = 9, height = 0.11 * length(levels(heat_data$row)) + 2, dpi = 300, limitsize = FALSE)

  # (b) Abundance series, one panel per species x type, one page per data type
  # (several pages, so this one stays a pdf). Lines break across census gaps
  # longer than five years, so the lone 1973 census stands as a separate point.
heat_data <- heat_data %>%
  group_by(row, type) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(segment = cumsum(c(1, diff(Year) > 5))) %>%
  ungroup()
pdf(file.path(fig_dir, "CLR_CheckPop_abundance.pdf"), width = 11, height = 8.5)
for (ty in type_levels) {
  d_ty <- filter(heat_data, type == ty)
  if (nrow(d_ty) == 0) next
  rows_ty <- rev(levels(droplevels(d_ty$row)))
  for (pg in split(rows_ty, ceiling(seq_along(rows_ty) / 30))) {
    p <- ggplot(filter(d_ty, row %in% pg), aes(x = Year, y = abundance)) +
      geom_line(aes(group = segment), colour = type_colour[ty], linewidth = 0.4) +
      geom_point(colour = type_colour[ty], size = 0.7) +
      facet_wrap(~ factor(row, levels = pg), scales = "free_y", ncol = 5) +
      labs(x = NULL, y = unique(d_ty$unit),
           title = paste0(site_tag, ": ", ty, " (", unique(d_ty$unit), ", per sampled m2)")) +
      theme_bw(base_size = 7) +
      theme(strip.text = element_text(size = 6))
    print(p)
  }
}
dev.off()

  # (c) Point plants against polygon cover: share of total basal cover that point
  # plants would carry under each assumed area (one box per assumption)
p_share <- ggplot(filter(pointshare, !is.na(point_share)),
                  aes(x = factor(assumed_cm2), y = 100 * point_share)) +
  geom_boxplot(outlier.size = 0.4, fill = "grey90") +
  labs(x = expression("Assumed basal area per point plant (cm"^2*")"),
       y = "Share of total basal cover (%)",
       title = paste0(site_tag, ": share of point plants in basal cover, per quadrat-year"),
       subtitle = "Polygons of grasses and forbs vs point plants (shrubs left out)") +
  theme_bw(base_size = 9)
ggsave(file.path(fig_dir, "CLR_CheckPop_pointshare.png"), p_share, width = 6, height = 4.5, dpi = 300)

# Report saved files -----------------------------------------------------------
saved_files <- data.frame(
  folder = c(rep(file.path("result/tables", site_tag), 5), rep(file.path("result/figs", site_tag), 3)),
  file   = c("CLR_CheckPop_species.csv", "CLR_CheckPop_abundance.csv", "CLR_CheckPop_effort.csv",
             "CLR_CheckPop_flags.csv", "CLR_CheckPop_pointshare.csv",
             "CLR_CheckPop_heatmap.png", "CLR_CheckPop_abundance.pdf", "CLR_CheckPop_pointshare.png"),
  content = c("one row per species x data type: records, years present, quadrats, abundance",
              "one row per species x data type x census year: abundance per sampled m2",
              "number of quadrats sampled per census year",
              "one row per flagged record: source, quadrat, year, check, note",
              "median and 90th percentile share of point plants, per assumed area",
              "species x year overview, coloured by data type",
              "abundance series per species, one page per data type",
              "share of point plants in basal cover under assumed areas")
)
stopifnot(all(file.exists(here(saved_files$folder, saved_files$file))))
cat("\nCLR_CheckPop.R (", site_tag, ") saved these files:\n", sep = "")
for (fd in unique(saved_files$folder)) {
  cat(sprintf("  %s/\n", fd))
  sf <- saved_files[saved_files$folder == fd, ]
  cat(sprintf("    %-30s %s\n", sf$file, sf$content), sep = "")
}
cat("Workspace objects kept: cover_raw, density_raw, counts_raw, species_raw, sampled, effort,",
    "pop_abundance, pop_species, pop_flags, pointshare\n")

  # Tidy up temporary objects
rm(list = intersect(c("pop_path", "n_stem", "inventory_raw", "flag_list", "add_flag", "recorded",
                      "off_inventory", "unlisted", "type_levels", "type_colour", "growth_form",
                      "abund_cover", "abund_point", "abund_count", "assumed_cm2", "basal_qy",
                      "point_qy", "pointshare_summary", "tab_dir", "fig_dir", "heat_data",
                      "row_order", "p_heat", "ty", "d_ty", "rows_ty", "pg", "p", "p_share",
                      "saved_files", "fd", "sf"), ls()))
