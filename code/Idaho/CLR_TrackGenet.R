################################################################################
## Track individual plants (genets) of the Idaho mapped quadrats
##
## Written by Hung-wei Lin Oct 2026
##
## Reads the polygon and point shapefiles of every quadrat-year, turns every point
## into a small circle, and gives every plant an individual ID (trackID) through time
## with plantTracker (Stears et al. 2022), which follows the tracking rules of
## Lauenroth & Adler (2008), with three additions: seedlings start new genets, the
## first-year genets of ARTR, HECO, POSE and PSSP come from the hand assignment of
## Chu & Adler (2015), and plants that merge into a neighbour get a censored fate.
## Raw files stay as they are; the tracked records go to 'data/data.derived/<site>'.
## Uses 'site_dir' and 'site_tag' from the runscript.


# Settings ---------------------------------------------------------------------
  # Same settings as the tracked Idaho data of Chu & Adler (2015) ('buf5_dorm2'):
  # a plant may move up to 5 cm between years, and may stay invisible for up to 2 years.
track_buff  <- 0.05   # m, buffer for overlap between year t and t+1
track_dorm  <- 2      # years of dormancy allowed
  # First year of a species in a quadrat: polygons closer than 'buffGenet' are taken as
  # pieces (ramets) of one plant; later years link pieces through overlap with last year.
  # Values chosen by agreement with the hand assignment of Chu & Adler (2015), see
  # CLR_TrackSensitivity.R and 'CLR_TrackSensitivity_buffGenet.csv'.
track_buffGenet_shrub <- 0.02  # m, shrub canopy pieces
track_buffGenet_other <- 0.01  # m, grass clumps and forb mats
  # TRUE: the hand-assigned first-year genets of Chu & Adler (2015) replace the distance
  # rule for ARTR, HECO, POSE and PSSP (assigned by eye from the first two maps of each quadrat).
use_hand_genets <- TRUE
  # TRUE: a point mapped as a seedling (seedling = "Y") starts a new genet.
seedling_rule <- TRUE
  # Seedling flags exist in the public data for point plants only (shrubs, forbs), from 1935 on;
  # grass seedlings sit among the polygons without a flag. The tracked files of Chu & Adler (2015)
  # carry seedling flags for the polygons of HECO, POSE and PSSP (same SP_ID as here).
  # TRUE adds those flags (three grasses only); FALSE keeps every species on the public data.
  # TRUE here: the flags change recruitment of these grasses by 3-14 SE (CLR_TrackSensitivity.R).
seedling_from_chu <- TRUE
  # TRUE: a basal plant whose outline overlaps a plant of the same species at the next census
  # that carries another trackID merged into it (coalescence), so its fate is censored (NA)
  # in place of death (Lauenroth & Adler 2008). Shrub canopies overlap without merging, so
  # shrubs keep the death.
censor_merged <- TRUE
  # Point plants get this area (as in Chu & Adler 2015, where every point record is 0.25 cm2)
point_area <- 0.25e-4     # m2
  # Records left out of tracking (non-vascular)
drop_forms <- c("moss", "lichen")
  # Tracking takes about 10 min; the result is saved and read again on later runs.
  # Set TRUE after changing a setting above.
retrack <- FALSE

  # Output folders of this site
der_dir <- here("data", "data.derived", site_tag)
tab_dir <- here("result", "tables", site_tag)
fig_dir <- here("result", "figs", site_tag)
for (d in c(der_dir, tab_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Read the shapefiles ----------------------------------------------------------
  # One pair of layers per quadrat-year: Qn_yy_C (polygons) and Qn_yy_D (points).
  # Units are metres inside a 1 x 1 m quadrat (no CRS).
shp_dir <- here("data", "data.source", site_dir, "shapefiles")
read_layers <- function(kind) {
  files <- list.files(shp_dir, pattern = paste0("_", kind, "\\.shp$"))
  layers <- lapply(files, function(f) {
    s <- st_read(file.path(shp_dir, f), quiet = TRUE)
    s$Quad <- sub("_.*", "", f)
    s$Year <- as.integer(sub("^Q\\d+_(\\d+)_.*", "\\1", f)) + 1900L
    s
  })
    # Join the layers through their attribute tables and geometries (much faster than rbind)
  geom <- do.call(c, lapply(layers, st_geometry))
  attr <- bind_rows(lapply(layers, st_drop_geometry))
  st_sf(attr, geometry = geom)
}
poly_raw  <- read_layers("C")
point_raw <- read_layers("D")
cat(sprintf("Read %d polygons and %d points from the shapefiles.\n", nrow(poly_raw), nrow(point_raw)))

# Clean the records ------------------------------------------------------------
species_raw <- read.csv(here("data", "data.source", site_dir, "species_list.csv"), na.strings = c("NA", ""))
growth_form <- setNames(species_raw$growthForm, species_raw$species)
stopifnot(all(c(poly_raw$species, point_raw$species) %in% names(growth_form)))

  # Polygons: repair self-crossing outlines (they keep their area)
n_bad <- sum(!st_is_valid(poly_raw))
poly_raw <- st_make_valid(poly_raw)
poly_raw <- poly_raw[!st_is_empty(poly_raw), ]
cat(sprintf("Repaired %d polygons with invalid outlines.\n", n_bad))
  # Repairs can turn a polygon into a collection (polygon + stray line); keep its polygon part
mixed <- st_is(poly_raw, "GEOMETRYCOLLECTION")
if (any(mixed)) {
  st_geometry(poly_raw)[mixed] <- st_collection_extract(st_geometry(poly_raw)[mixed], "POLYGON")
}
poly <- data.frame(record_id = paste(poly_raw$Quad, poly_raw$Year, "C", poly_raw$SP_ID, sep = "_"),
                   Quad = poly_raw$Quad, Year = poly_raw$Year, Species = poly_raw$species,
                   record_type = "polygon", seedling = NA)
poly <- st_sf(poly, geometry = st_cast(st_geometry(poly_raw), "MULTIPOLYGON"))

  # Points: stem points belong to shrub canopy polygons (metadata: remove before analysis);
  # the rest become circles with an area of 'point_area'
point_raw <- point_raw[!point_raw$stem %in% "Y", ]
  # A buffered circle is a 32-sided polygon, slightly smaller than the true circle;
  # the radius is scaled so that the polygon area equals 'point_area'
r0 <- sqrt(point_area / pi)
r_adj <- r0 * sqrt(point_area / as.numeric(st_area(st_buffer(st_point(c(0, 0)), r0, nQuadSegs = 8))))
circle <- st_buffer(st_geometry(point_raw), r_adj, nQuadSegs = 8)
point <- data.frame(record_id = paste(point_raw$Quad, point_raw$Year, "D", point_raw$OBJECTID, sep = "_"),
                    Quad = point_raw$Quad, Year = point_raw$Year, Species = point_raw$species,
                    record_type = "point", seedling = point_raw$seedling)
point <- st_sf(point, geometry = st_cast(circle, "MULTIPOLYGON"))

records <- rbind(poly, point)
chu_seed <- bind_rows(lapply(c("HECO", "POSE", "PSSP"), function(k) {
  read.csv(here("data", "data.source", "Chu, 2015", "wrapper-ID", "speciesData", k,
                paste0(k, "_buf5_dorm2.csv")))
}))
chu_seed <- with(chu_seed[chu_seed$seedling == "Y", ], paste(quad, year + 1900, "C", SP_ID, sep = "_"))
if (seedling_from_chu) {
  records$seedling[records$record_id %in% chu_seed] <- "Y"
  cat(sprintf("Seedling flags from Chu & Adler (2015): %d grass polygons.\n", sum(records$record_id %in% chu_seed)))
}
records <- records[!growth_form[records$Species] %in% drop_forms, ]
if (!seedling_rule) records$seedling <- NA
records$Site <- site_tag
stopifnot(!any(duplicated(records$record_id)))
cat(sprintf("Records to track: %d (%d polygons, %d points), %d species.\n", nrow(records),
            sum(records$record_type == "polygon"), sum(records$record_type == "point"),
            n_distinct(records$Species)))

# Sampling inventory -----------------------------------------------------------
  # Named list: for every quadrat, the years it was mapped (also years without any plant)
inventory_raw <- read.csv(here("data", "data.source", site_dir, "quad_inventory.csv"))
quads <- grep("^Q", names(inventory_raw), value = TRUE)
inv <- lapply(setNames(quads, quads), function(q) sort(inventory_raw[[q]][!is.na(inventory_raw[[q]])]) + 1900L)
  # Shapefile years should all be inventory years
stopifnot(all(mapply(function(q, y) y %in% inv[[q]], records$Quad, records$Year)))

# Hand-assigned first-year genets (Chu & Adler 2015) ----------------------------
  # One row per polygon of the first year of a species in a quadrat, with the genet number
  # given by eye; polygons join through SP_ID (some ARTR rows point to other species: left out)
chu_dir <- here("data", "data.source", "Chu, 2015", "wrapper-ID", "speciesData")
chu_code <- c(ARTR = "Artemisia tripartita", HECO = "Hesperostipa comata",
              POSE = "Poa secunda", PSSP = "Pseudoroegneria spicata")
hand_files <- c(ARTR = "ARTR/assign_genetsARTR.csv", HECO = "HECO/assign_genetsHECO.csv",
                POSE = "POSE/POSEassign_genets.csv", PSSP = "PSSP/assign_genets_PSSP.csv")
hand <- bind_rows(lapply(names(hand_files), function(k) {
  read.csv(file.path(chu_dir, hand_files[k])) %>%
    transmute(record_id = paste(quad, year + 1900L, "C", SP_ID, sep = "_"),
              Species = chu_code[k], GenID = paste(k, GenID))
}))
hand <- semi_join(hand, st_drop_geometry(records), by = c("record_id", "Species"))
records$hand_genet <- hand$GenID[match(records$record_id, hand$record_id)]
cat(sprintf("Hand-assigned first-year genets: %d polygons of %d species.\n",
            sum(!is.na(records$hand_genet)), n_distinct(hand$Species)))

# Track the plants -------------------------------------------------------------
  # Species with polygons can have several pieces per plant (grass clumps, shrub canopy
  # lobes, Antennaria mats), so they are clonal here; point-only species have one point per plant.
clonal <- records %>% st_drop_geometry() %>% group_by(Species) %>%
  summarise(clonal = any(record_type == "polygon")) %>% as.data.frame()
buffGenet <- data.frame(Species = clonal$Species,
                        buffGenet = ifelse(growth_form[clonal$Species] %in% "shrub",
                                           track_buffGenet_shrub, track_buffGenet_other))

  # Seedling rule: a point mapped as a seedling (seedling = "Y") is a new plant, so it starts
  # a new genet and never links to a plant of the year before (in the tracked data of
  # Chu & Adler 2015 every seedling record has age 1). plantTracker has no option for this;
  # the line of its internal 'assign' function that overlaps the plants of year t (buffered)
  # with those of year t+1 gets the seedlings of year t+1 left out, so they become recruits.
  # trackSpp passes only its own columns to 'assign' (plus a row index 'indexStore'), so the
  # seedling flags of the data in hand sit in 'patch_env$is_seedling', in the same row order.
assign_src <- deparse(plantTracker:::assign, width.cutoff = 500L)
target <- grep("st_intersection(tempPreviousBuff, tempCurrentYear)", assign_src, fixed = TRUE)
stopifnot(length(target) == 1)
assign_src[target] <- sub("tempCurrentYear))",
                          "tempCurrentYear[!is_seedling[tempCurrentYear$indexStore], ]))",
                          assign_src[target], fixed = TRUE)
  # Hand genets: inside 'ifClonal' (first-year grouping), polygons with a hand-assigned
  # genet take it; the rest keep the distance rule. Lookup in 'patch_env$hand_genet'.
target <- grep("cloneDat$genetID <- groupByGenet(cloneDat, buffGenet)", assign_src, fixed = TRUE)
stopifnot(length(target) == 1)
assign_src[target] <- paste0(assign_src[target], "; hg <- hand_genet[cloneDat$indexStore]; ",
  "if (any(!is.na(hg))) cloneDat$genetID <- as.integer(factor(ifelse(is.na(hg), ",
  "paste0(\"rule \", cloneDat$genetID), hg)))")
assign_new <- eval(parse(text = assign_src))
patch_env <- new.env(parent = asNamespace("plantTracker"))
environment(assign_new) <- patch_env
assignInNamespace("assign", assign_new, ns = "plantTracker")

  # plantTracker 1.2.1 leaves out the records of a species in the last census year of a
  # quadrat when that year follows a gap longer than 'dorm' (1957 -> 1973), or follows a
  # census year without that species after such a gap. These plants cannot link to any
  # earlier plant, so they enter here as new genets (pieces grouped with 'buffGenet');
  # they are recruits when the quadrat was mapped the year before, otherwise unknown.
set_patch <- function(d) {
  patch_env$is_seedling <- d$seedling %in% "Y"
  patch_env$hand_genet  <- if (use_hand_genets) d$hand_genet else rep(NA, nrow(d))
}
add_left_out <- function(d, out, q) {
  left <- d[!d$record_id %in% out$record_id, ]
  if (nrow(left) == 0) return(out)
  prev_census <- function(y) { p <- inv[[q]][inv[[q]] < y]; if (length(p) == 0) NA else max(p) }
  added <- lapply(split(left, paste(left$Species, left$Year)), function(s) {
    set_patch(s)
    a <- trackSpp(dat = s, inv = list(unique(s$Year)) %>% setNames(q), dorm = track_dorm,
                  buff = track_buff, buffGenet = buffGenet[buffGenet$Species %in% s$Species, ],
                  clonal = clonal[clonal$Species %in% s$Species, ],
                  aggByGenet = FALSE, printMessages = FALSE)
    new_recruit <- prev_census(unique(s$Year)) %in% (unique(s$Year) - 1L)
    a$recruit <- ifelse(new_recruit, 1L, NA_integer_)
    a$age     <- ifelse(new_recruit, 0L, NA_integer_)
    a
  })
  added <- do.call(rbind, added)[, names(out)]
  stopifnot(!any(added$trackID %in% out$trackID[out$Species %in% added$Species]))
  rbind(out, added)
}

  # Track a set of records quadrat by quadrat (prints progress; plantTracker loops over
  # species inside); trackIDs become unique across quadrats
track_records <- function(recs) {
  t0 <- Sys.time()
  out <- lapply(intersect(quads, unique(recs$Quad)), function(q) {
    d <- recs[recs$Quad == q, ]
    set_patch(d)
    o <- trackSpp(dat = d, inv = inv[q], dorm = track_dorm, buff = track_buff,
                  buffGenet = buffGenet[buffGenet$Species %in% d$Species, ],
                  clonal = clonal[clonal$Species %in% d$Species, ],
                  aggByGenet = FALSE, printMessages = FALSE)
    n_left <- sum(!d$record_id %in% o$record_id)
    o <- add_left_out(d, o, q)
    cat(sprintf("  %s: %d records (%d added after tracking), %d genet-years\n", q, nrow(o),
                n_left, n_distinct(paste(o$Species, o$Year, o$trackID))))
    o
  })
  out <- do.call(rbind, out)
  cat(sprintf("Tracking took %.1f min.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  stopifnot(nrow(out) == nrow(recs), setequal(out$record_id, recs$record_id))
  out$trackID <- paste(out$Quad, out$trackID, sep = "_")
  out[order(out$Quad, out$Species, out$Year, out$trackID), ]
}

  # Fates decided by a neighbour: for every genet that dies at the next census of its
  # quadrat, 'fate_contested' = a plant of the same species with another trackID sits within
  # 'track_buff' of it at that census (the overlap rule gave that plant to a neighbour);
  # 'fate_merged' = that plant overlaps its outline itself (coalescence).
  # 'survives_raw' keeps the fate before censoring.
flag_fates <- function(tr) {
  tr$survives_raw <- tr$survives_tplus1
  tr$fate_contested <- FALSE
  tr$fate_merged <- FALSE
  dead <- which(tr$survives_tplus1 %in% 0)
  next_year <- mapply(function(q, y) { n <- inv[[q]][inv[[q]] > y]; if (length(n)) n[1] else NA },
                      tr$Quad[dead], tr$Year[dead])
  key_dead <- paste(tr$Quad[dead], tr$Species[dead], next_year)
  key_all  <- paste(tr$Quad, tr$Species, tr$Year)
  for (k in unique(key_dead[!is.na(next_year)])) {
    i <- dead[key_dead == k]
    nx <- tr[key_all == k, ]
    if (nrow(nx) == 0) next
    tr$fate_contested[i] <- lengths(st_intersects(st_buffer(tr[i, ], track_buff), nx)) > 0
    tr$fate_merged[i]    <- lengths(st_intersects(tr[i, ], nx)) > 0
  }
    # Genet level: a genet is contested or merged when any of its pieces is
  g <- paste(tr$trackID, tr$Year)
  tr$fate_contested <- ave(tr$fate_contested, g, FUN = any)
  tr$fate_merged    <- ave(tr$fate_merged, g, FUN = any) & !growth_form[tr$Species] %in% "shrub"
  if (censor_merged) tr$survives_tplus1[tr$fate_merged] <- NA
  tr
}

if (file.exists(file.path(der_dir, "CLR_TrackGenet_ramet.rds")) && !retrack) {
  tracked <- readRDS(file.path(der_dir, "CLR_TrackGenet_ramet.rds"))
  cat("Read the tracked records saved before (set 'retrack <- TRUE' to track again).\n")
} else {
  tracked <- flag_fates(track_records(records))
  saveRDS(tracked, file.path(der_dir, "CLR_TrackGenet_ramet.rds"))
}
fates <- st_drop_geometry(tracked) %>% filter(survives_raw %in% 0) %>% distinct(trackID, Year, .keep_all = TRUE)
cat(sprintf("Deaths: %d genet-years; %.0f%% with a neighbour of the same species within %g cm at the next census, %.0f%% merged (censored: %s).\n",
            nrow(fates), 100 * mean(fates$fate_contested), track_buff * 100, 100 * mean(fates$fate_merged), censor_merged))
