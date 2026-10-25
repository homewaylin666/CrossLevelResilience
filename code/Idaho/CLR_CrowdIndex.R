################################################################################
## Build the crowding indices of the Part A model (Idaho)
##
## Written by Hung-wei Lin Oct 2026
##
## Turns the neighbour table of CLR_CrowdKernel.R into the crowding of every plant with an
## annual transition, with the two fitted kernels cut at their range R (Teller et al. 2016):
##   W_intra, W_grass, W_shrub = sum over annuli of F_main(d) x neighbour area (cm2)
##     (intra: same species; grass: other non-shrub species, basal area; shrub: canopy area);
##   W_forb = sum over annuli of F_forb(d) x forb plants (count).
## A forb plant belongs to the pooled forb population, so its own class is forb: W_forb holds
## the other forb plants around it, and W_intra stays empty (NA).
## Adds the quadrat totals at t for the recruitment model (other grass cover, other shrub cover,
## forb plants), and builds the recruitment table of the pooled forb population (forb submodel).
## Neighbours outside the quadrat are unmapped; 'focal_ok' (5-cm rule) marks the plants for the
## survival and growth models (edge check in CLR_CrowdKernel.R).
## Uses 'genet', 'recruit', 'species', 'census', 'quad_info' (CLR_ProcessDataPartA.R) and
## 'rings', 'w_main', 'w_forb', 'R_main', 'R_forb', 'sp_class', 'K' (CLR_CrowdKernel.R).


# Settings ---------------------------------------------------------------------
  # Correlations among the crowding indices: species with at least this many focal
  # survival records (fewer give unstable correlations)
cor_min_n <- 100
  # Correlations above this value (absolute) are listed in the console report
cor_flag <- 0.7

der_dir <- here("data", "data.derived", site_tag)
tab_dir <- here("result", "tables", site_tag)
stopifnot(exists("genet"), exists("recruit"), exists("species"), exists("census"), exists("quad_info"),
          exists("rings"), exists("w_main"), exists("w_forb"), exists("sp_class"))
stopifnot(length(w_main) == K, length(w_forb) == K, all(species$Species %in% names(sp_class)))
W_cols <- c("W_intra", "W_grass", "W_shrub", "W_forb")

# Crowding of every plant ------------------------------------------------------
ring_mat <- function(cl) as.matrix(rings[paste(cl, seq_len(K), sep = "_")])
crowd <- rings %>%
  select(Quad, Year, Species, trackID) %>%
  mutate(W_intra = drop(ring_mat("intra") %*% w_main),
         W_grass = drop(ring_mat("grass") %*% w_main),
         W_shrub = drop(ring_mat("shrub") %*% w_main),
         W_forb  = drop(ring_mat("forb")  %*% w_forb))
  # Forb plants: the other forb plants sit in their 'intra' annuli (same pooled population),
  # so their forb crowding uses those annuli with the forb kernel
is_forb <- sp_class[crowd$Species] == "forb"
stopifnot(all(ring_mat("forb")[is_forb, ] == 0))
crowd$W_forb[is_forb]  <- drop(ring_mat("intra")[is_forb, , drop = FALSE] %*% w_forb)
crowd$W_intra[is_forb] <- NA

genet_W <- genet %>%
  filter(!is.na(survives_tplus1)) %>%
  mutate(crowd_class = unname(sp_class[Species])) %>%
  left_join(crowd, by = c("Quad", "Year", "Species", "trackID"))
  # Every plant with an annual transition has its neighbour table
stopifnot(nrow(genet_W) == nrow(crowd), !anyNA(genet_W[c("W_grass", "W_shrub", "W_forb")]),
          all(is.na(genet_W$W_intra) == (genet_W$crowd_class == "forb")))
cat(sprintf("Crowding of %d plants (main kernel to %g cm, forb kernel to %g cm); %d of them focal plants.\n",
            nrow(genet_W), R_main, R_forb, sum(genet_W$focal_ok)))

# Quadrat totals for recruitment -----------------------------------------------
  # Cover of grass and shrub and the number of forb plants in each quadrat at t
  # (every plant mapped at t). Recruitment uses these totals in place of W.
quad_total <- genet %>%
  mutate(cl = sp_class[Species]) %>%
  group_by(Quad, Year) %>%
  summarise(grass_cover_t_cm2 = sum(area_cm2[cl == "grass"]),
            shrub_cover_t_cm2 = sum(area_cm2[cl == "shrub"]),
            n_forb_t = sum(cl == "forb"), .groups = "drop")

recruit_W <- recruit %>%
  mutate(crowd_class = unname(sp_class[Species])) %>%
  left_join(quad_total, by = c("Quad", "Year")) %>%
  mutate(across(c(grass_cover_t_cm2, shrub_cover_t_cm2, n_forb_t), ~ replace_na(.x, 0)),
           # Other grass and other shrub: the species' own cover leaves the total of its class
         other_grass_t_cm2 = grass_cover_t_cm2 - ifelse(crowd_class == "grass", cover_t_cm2, 0),
         other_shrub_t_cm2 = shrub_cover_t_cm2 - ifelse(crowd_class == "shrub", cover_t_cm2, 0))
stopifnot(nrow(recruit_W) == nrow(recruit),
          all(recruit_W$other_grass_t_cm2 > -1e-6), all(recruit_W$other_shrub_t_cm2 > -1e-6))
  # Rounding of the sums leaves tiny negative values where the species holds all the cover
recruit_W <- recruit_W %>%
  mutate(across(c(other_grass_t_cm2, other_shrub_t_cm2), ~ pmax(.x, 0))) %>%
  select(-grass_cover_t_cm2, -shrub_cover_t_cm2)

# Recruitment table of the pooled forb population ------------------------------
  # One row per quadrat and annual transition: forb plants at t (parents), forb recruits at
  # t+1, and the group mean of forb plants at t (every quadrat of the group mapped in year t),
  # as the parent cover and group mean cover of the other species
qy_cols <- c("Quad", "Year", "Group", "grazed_quad", "paddock", "timing_t", "graze_fall_t",
             "graze_spring_t1", "graze_total", "ppt_wy", "tmean_wy", "ppt_wy_lag", "tmean_wy_lag")
forb_group <- census %>%
  select(Quad, Year) %>%
  left_join(select(quad_info, Quad, Group), by = "Quad") %>%
  left_join(select(quad_total, Quad, Year, n_forb_t), by = c("Quad", "Year")) %>%
  mutate(n_forb_t = replace_na(n_forb_t, 0)) %>%
  group_by(Group, Year) %>%
  summarise(group_n_forb_t = mean(n_forb_t), n_quad_group_t = n(), .groups = "drop")

recruit_forb <- recruit_W %>%
  filter(crowd_class == "forb") %>%
  group_by(across(all_of(qy_cols))) %>%
  summarise(n_forb_check = first(n_forb_t), n_forb_t = sum(n_genet_t), recruits_tplus1 = sum(recruits_tplus1),
            other_grass_t_cm2 = first(other_grass_t_cm2), other_shrub_t_cm2 = first(other_shrub_t_cm2),
            .groups = "drop") %>%
  left_join(forb_group, by = c("Group", "Year")) %>%
  arrange(Quad, Year)
  # Parents counted from the species rows match the quadrat total; one row per quadrat transition
stopifnot(all(recruit_forb$n_forb_t == recruit_forb$n_forb_check),
          nrow(recruit_forb) == n_distinct(paste(recruit$Quad, recruit$Year)),
          !anyNA(recruit_forb$group_n_forb_t))
recruit_forb <- select(recruit_forb, -n_forb_check)
cat(sprintf("Forb recruitment table: %d quadrat transitions, %d forb plants at t, %d forb recruits.\n",
            nrow(recruit_forb), sum(recruit_forb$n_forb_t), sum(recruit_forb$recruits_tplus1)))

# Distribution and correlation of the crowding indices -------------------------
  # Focal plants only (the plants of the survival and growth models)
focal <- filter(genet_W, focal_ok)
crowd_summary <- focal %>%
  pivot_longer(all_of(W_cols), names_to = "index", values_to = "W") %>%
  filter(!is.na(W)) %>%
  group_by(Species, crowd_class, index) %>%
  summarise(n = n(), mean = mean(W), sd = sd(W), median = median(W), share_zero = mean(W == 0),
            .groups = "drop") %>%
  mutate(code = species$code[match(Species, species$Species)], .after = Species)

  # Pairwise correlation of the indices within each species (indices that vary)
crowd_cor <- bind_rows(lapply(split(focal, focal$Species), function(d) {
  if (nrow(d) < cor_min_n) return(NULL)
  W <- as.matrix(d[W_cols])
  W <- W[, colSums(is.na(W)) == 0, drop = FALSE]
  W <- W[, apply(W, 2, sd) > 0, drop = FALSE]
  if (ncol(W) < 2) return(NULL)
  r <- cor(W)
  pr <- which(upper.tri(r), arr.ind = TRUE)
  data.frame(Species = d$Species[1], code = species$code[match(d$Species[1], species$Species)],
             crowd_class = d$crowd_class[1], n = nrow(d),
             index_1 = colnames(r)[pr[, 1]], index_2 = colnames(r)[pr[, 2]], r = r[pr])
}))
high_cor <- filter(crowd_cor, abs(r) > cor_flag)
cat(sprintf("Correlations among crowding indices (%d species with at least %d focal plants): median |r| %.2f, largest |r| %.2f; %d pairs above %.1f.\n",
            n_distinct(crowd_cor$Species), cor_min_n, median(abs(crowd_cor$r)), max(abs(crowd_cor$r)),
            nrow(high_cor), cor_flag))
if (nrow(high_cor) > 0) print(as.data.frame(mutate(high_cor, r = round(r, 2))))
  # Forb crowding against the other indices
forb_cor <- filter(crowd_cor, index_1 == "W_forb" | index_2 == "W_forb")
cat(sprintf("W_forb against the other indices: median |r| %.2f, largest |r| %.2f.\n",
            median(abs(forb_cor$r)), max(abs(forb_cor$r))))

# Save -------------------------------------------------------------------------
write.csv(genet_W, file.path(der_dir, "CLR_CrowdIndex_genet.csv"), row.names = FALSE)
write.csv(recruit_W, file.path(der_dir, "CLR_CrowdIndex_recruit.csv"), row.names = FALSE)
write.csv(recruit_forb, file.path(der_dir, "CLR_CrowdIndex_forbrecruit.csv"), row.names = FALSE)
write.csv(crowd_summary, file.path(tab_dir, "CLR_CrowdIndex_summary.csv"), row.names = FALSE)
write.csv(crowd_cor, file.path(tab_dir, "CLR_CrowdIndex_cor.csv"), row.names = FALSE)
cat("Saved:\n",
    " data/data.derived/", site_tag, ": CLR_CrowdIndex_genet.csv, CLR_CrowdIndex_recruit.csv, CLR_CrowdIndex_forbrecruit.csv\n",
    " result/tables/", site_tag, ": CLR_CrowdIndex_summary.csv, CLR_CrowdIndex_cor.csv\n", sep = "")
rm(list = intersect(c("ring_mat", "crowd", "is_forb", "quad_total", "qy_cols", "forb_group", "focal",
                      "high_cor", "forb_cor", "cor_min_n", "cor_flag"), ls()))
  # The ramet polygons and the neighbour table stay on disk ('CLR_TrackGenet_ramet.rds',
  # 'CLR_CrowdKernel_rings.rds'); later steps use the crowding indices, so they leave the workspace
  # (it is saved after every step, and these two are its largest objects)
rm(list = intersect(c("ramet", "rings"), ls()))
