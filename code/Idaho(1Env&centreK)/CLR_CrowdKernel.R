################################################################################
## Estimate the crowding kernel and its range (Idaho)
##
## Written by Hung-wei Lin Oct 2026
##
## Follows Teller et al. (2016, Methods Ecol Evol 7: 171-183; Appendix S1 and Data S1).
## (1) Neighbour table: for every plant with an annual transition, the neighbours in
##     concentric annuli around its centroid (rings of Teller et al. 2016), by class:
##     conspecific, heterospecific grass (basal area), heterospecific shrub (canopy area),
##     and forb (all point species pooled as one population, counted as plants).
##     Area inside an annulus is the part of each neighbour polygon that lies in it.
## (2) Kernels: monotone decreasing spline kernels fitted by maximising the survival + growth
##     likelihood of the main size species, with a roughness penalty chosen by approximate AIC.
##     The main kernel is shared by the three area classes (same species, grass, shrub); the
##     forb class (counted as plants) has its own kernel, fitted with the main kernel fixed.
## (3) Range R of each kernel: the kernel is refitted with its reach cut at a set of distances,
##     and R is the shortest reach at which every crowding coefficient (per SD of crowding) and
##     its interaction with water-year precipitation stay within 'R_se' standard errors of their
##     values at the reference reach 'R_ref'.
## Uses 'genet', 'species' and 'ramet' of CLR_ProcessDataPartA.R, and 'site_tag'.


# Settings ---------------------------------------------------------------------
  # Annulus boundaries (cm), as Teller et al. (2016): 2 cm wide to 20 cm, 5 cm to 50 cm,
  # 10 cm to 150 cm (the quadrat diagonal is 141 cm). Neighbours count at the mid-radius.
ring_edges <- c(seq(0, 20, 2), seq(25, 50, 5), seq(60, 150, 10))
  # Polygon area is spread over a regular grid of this spacing (cm) to split it among annuli
  # (polygons smaller than a grid cell keep their full area at the centroid).
grid_cm <- 0.5
  # Focal species for the kernel fit: polygon species with at least this many annual
  # survival records (same rule as CLR_TrackSensitivity.R)
kernel_min_survival <- 500
  # Spline knots (cm) for the kernel, dense at short distances where the kernel changes
kernel_breaks <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50, 80, 145)
  # Roughness penalties tried (Teller et al. 2016 used 10^-6 to 10^6). Three values: the
  # first fits (Oct 2026) gave nearly the same kernel and AIC from 10^-4 to 10^0.
kernel_lambda <- 10^c(-3, 0, 3)
  # Multistart effort per penalty (Teller et al. 2016: 249 random starts, 5 rough fits)
kernel_starts <- 249
kernel_rough  <- 5
  # Random seed of the random starts (2026 throughout the pipeline)
kernel_seed <- 2026
  # TRUE: year as a factor in the survival and growth regressions (Teller et al. 2016,
  # Figure S1.2). TRUE here: Part A has year effects, and without them the kernel keeps a long
  # tail (to about 100 cm) that picks up differences between years (first fit, Oct 2026:
  # neighbours beyond 25 cm add 20 log-likelihood units without year, 1 unit with year).
kernel_year <- TRUE
  # TRUE: the kernel regressions also hold the quadrat-year totals of grass cover, shrub cover
  # and forb plants (standardised). These totals carry how good a quadrat is in a year, which
  # would otherwise load on the far annuli (first fits: a tail to 100 cm, mostly from
  # heterospecific and forb neighbours). Used only to fit the kernels; Part A stays as it is.
kernel_controls <- TRUE
  # TRUE: the kernel regressions also hold the distance of the plant to the quadrat border
  # (natural spline, capped at 50 cm). Survival changes with this distance, and plants in the
  # middle of a quadrat see more of their far annuli, so position would otherwise load on the
  # far annuli. Used only to fit the kernels.
kernel_border <- FALSE
  # TRUE: neighbour amounts divided by the share of each annulus inside the quadrat
  # (edge correction); FALSE: neighbours outside the quadrat count as absent (Teller et al. 2016)
edge_correct <- FALSE
  # Range R. Distances run from the centre of the focal plant (most focal plants have an
  # equivalent radius below 1 cm). The fit keeps improving with reach (checks, Oct 2026), so
  # the deviance cannot set R. R is set instead by the estimates used later: the crowding
  # coefficients (beta, per SD of crowding) and their interaction with standardised water-year
  # precipitation (delta). R = shortest candidate reach at which every beta and delta of the
  # kernel's classes stays within 'R_se' standard errors of its value at the reference reach
  # 'R_ref' (in the spirit of the one-standard-error rule: the simplest model within one SE of
  # the best; Breiman et al. 1984, Hastie et al. 2009). R_ref = 25 cm, the next annulus edge beyond
  # the candidates: the check asks whether neighbours one annulus further away change the estimates.
  # Both kernels end in a flat tail before 25 cm, and a full reach stretches the main kernel's tail
  # over the quadrat, a measure of quadrat-wide abundance (check, Oct 2026: 0.59 from 17 to 145 cm).
R_candidates <- c(10, 12, 14, 16, 18, 20)
R_ref <- 25
R_se <- 1
  # Full reach (cm), refitted for the main kernel and reported for the record
R_full <- 150
  # The neighbour table takes a few minutes and the kernel fits 3-4 h.
  # Both are saved and read again on later runs while the settings above stay the same;
  # a changed setting builds them again.

der_dir <- here("data", "data.derived", site_tag)
tab_dir <- here("result", "tables", site_tag)
fig_dir <- here("result", "figs", site_tag)
stopifnot(exists("genet"), exists("species"), exists("ramet"))

K <- length(ring_edges) - 1
ring_mid <- (head(ring_edges, -1) + tail(ring_edges, -1)) / 2
classes <- c("intra", "grass", "shrub", "forb")
classes_main <- c("intra", "grass", "shrub")   # area classes, main kernel; 'forb' has its own

  # Neighbour class of every species: point species (no size) pooled as 'forb',
  # shrubs (canopy), every other polygon species (basal, mostly grasses)
sp_class <- with(species, setNames(case_when(size_measure == "point" ~ "forb",
                                             growth_form == "shrub" ~ "shrub",
                                             TRUE ~ "grass"), Species))
  # Population identity used for 'intra': the pooled forb counts as one population
crowd_sp <- function(s) ifelse(sp_class[s] == "forb", "forb", s)
cat(sprintf("Neighbour classes: %d grass (basal), %d shrub (canopy), %d point species pooled as forb.\n",
            sum(sp_class == "grass"), sum(sp_class == "shrub"), sum(sp_class == "forb")))

# Neighbour table --------------------------------------------------------------
  # Share of the disk of radius r around (cx, cy) inside the 100 x 100 cm quadrat,
  # by numerical integration across x (vectorised over plants)
disk_in_quad <- function(cx, cy, r, n = 400) {
  t <- (seq_len(n) - 0.5) / n * 2 * r - r
  x <- outer(cx, t, "+")
  half <- matrix(sqrt(r^2 - t^2), length(cx), n, byrow = TRUE)
  len <- pmax(0, pmin(100, cy + half) - pmax(0, cy - half)) * (x >= 0 & x <= 100)
  rowSums(len) * 2 * r / n / (pi * r^2)
}

  # Neighbours of every focal plant of one quadrat-year, as a matrix
  # (rows = focal plants, columns = class x annulus)
ring_one <- function(foc, pieces) {
  n_f <- nrow(foc); n_p <- nrow(pieces)
  out <- matrix(0, n_f, length(classes) * K)
  if (n_p == 0) return(out)
  d <- sqrt(outer(foc$x_cm, pieces$x, "-")^2 + outer(foc$y_cm, pieces$y, "-")^2)
  k <- findInterval(d, ring_edges)
  cls <- matrix(match(pieces$class, classes), n_f, n_p, byrow = TRUE)
  same <- outer(foc$crowd_sp, pieces$crowd_sp, "==")
  cls[same] <- 1L
    # The focal plant itself is left out of its own neighbourhood
  own <- same & outer(foc$trackID, pieces$trackID, "==")
  keep <- k <= K & !own
  key <- ((cls[keep] - 1L) * K + k[keep] - 1L) * n_f + row(d)[keep]
  s <- rowsum(rep(pieces$w, each = n_f)[keep], key)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

ring_file <- file.path(der_dir, "CLR_CrowdKernel_rings.rds")
ring_set <- list(ring_edges = ring_edges, grid_cm = grid_cm, sp_class = sp_class)
rings <- if (file.exists(ring_file)) readRDS(ring_file) else NULL
if (identical(attr(rings, "settings"), ring_set)) {
  cat("Read the neighbour table saved before (same settings).\n")
} else {
  t0 <- Sys.time()
    # Focal plants: every genet-year with an annual transition (all species)
  foc_all <- genet %>%
    filter(!is.na(survives_tplus1)) %>%
    mutate(crowd_sp = crowd_sp(Species)) %>%
    select(Quad, Year, Species, trackID, crowd_sp, x_cm, y_cm)
    # Neighbours with an area: every ramet of the polygon species (their point records
    # are 0.25 cm2 circles); neighbours counted as plants: genets of the forb class
  ramet_area <- ramet[sp_class[ramet$Species] != "forb", ]
  ramet_key <- split(seq_len(nrow(ramet_area)), paste(ramet_area$Quad, ramet_area$Year))
  forb_genet <- genet %>%
    filter(sp_class[Species] == "forb") %>%
    transmute(Quad, Year, trackID, crowd_sp = "forb", class = "forb", x = x_cm, y = y_cm, w = 1)
  grid_xy <- expand.grid(x = seq(grid_cm / 2, 100 - grid_cm / 2, by = grid_cm),
                         y = seq(grid_cm / 2, 100 - grid_cm / 2, by = grid_cm))
  grid_geom <- st_geometry(st_as_sf(grid_xy / 100, coords = c("x", "y")))  # metres, as the ramets

  ring_names <- c(paste(rep(classes, each = K), rep(seq_len(K), length(classes)), sep = "_"),
                  paste("propQ", seq_len(K), sep = "_"))
  qy <- distinct(foc_all, Quad, Year)
  rings <- vector("list", nrow(qy))
  for (i in seq_len(nrow(qy))) {
    q <- qy$Quad[i]; y <- qy$Year[i]
    foc <- foc_all[foc_all$Quad == q & foc_all$Year == y, ]
    ri <- ramet_key[[paste(q, y)]]
    ra <- ramet_area[if (is.null(ri)) integer(0) else ri, ]
    pieces <- NULL
    if (nrow(ra) > 0) {
        # Grid points inside each polygon share its area; polygons without a grid point
        # keep their full area at their centroid
      hit <- st_intersects(st_geometry(ra), grid_geom)
      n_hit <- lengths(hit)
      a_cm2 <- ra$basalArea_ramet * 1e4
      cen <- st_coordinates(st_centroid(st_geometry(ra))) * 100
      j <- rep(seq_len(nrow(ra)), pmax(n_hit, 1))
      g <- unlist(lapply(seq_len(nrow(ra)), function(r) if (n_hit[r] > 0) hit[[r]] else NA_integer_))
      pieces <- data.frame(trackID = ra$trackID[j],
                           crowd_sp = crowd_sp(ra$Species[j]),
                           class = sp_class[ra$Species[j]],
                           x = ifelse(is.na(g), cen[j, 1], grid_xy$x[g]),
                           y = ifelse(is.na(g), cen[j, 2], grid_xy$y[g]),
                           w = a_cm2[j] / pmax(n_hit[j], 1))
    }
    fb <- forb_genet[forb_genet$Quad == q & forb_genet$Year == y, c("trackID", "crowd_sp", "class", "x", "y", "w")]
    pieces <- rbind(pieces, fb)
    m <- ring_one(foc, pieces)
      # Share of every annulus inside the quadrat
    inside <- sapply(ring_edges[-1], function(r) disk_in_quad(foc$x_cm, foc$y_cm, r) * pi * r^2)
    inside <- matrix(inside, nrow(foc))
    prop_q <- (inside - cbind(0, inside[, -K, drop = FALSE])) / matrix(diff(pi * ring_edges^2), nrow(foc), K, byrow = TRUE)
    rings[[i]] <- cbind(foc[c("Quad", "Year", "Species", "trackID")],
                        setNames(as.data.frame(cbind(m, pmin(pmax(prop_q, 0), 1))), ring_names))
    if (i %% 50 == 0) cat(sprintf("  %d of %d quadrat-years\n", i, nrow(qy)))
  }
  rings <- bind_rows(rings)
  attr(rings, "settings") <- ring_set
  saveRDS(rings, ring_file)
  cat(sprintf("Neighbour table: %d plants in %d quadrat-years (%.1f min).\n",
              nrow(rings), nrow(qy), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(foc_all, ramet_area, ramet_key, forb_genet, grid_xy, grid_geom, qy, ring_names)
}
stopifnot(!any(duplicated(rings[c("Quad", "Year", "Species", "trackID")])))


# Data for the kernel fit ------------------------------------------------------
  # Focal plants of the main size species, away from the quadrat edge (focal_ok, 5 cm as
  # Teller et al. 2016); growth only for plants with a real size at t
kernel_species <- species$Species[species$size_measure != "point" & species$n_survival >= kernel_min_survival]
cat("Focal species for the kernel:", paste(species$code[match(kernel_species, species$Species)], collapse = ", "), "\n")

  # Quadrat-year totals by class (all plants of the quadrat), standardised
quad_total <- genet %>%
  mutate(cl = sp_class[Species]) %>%
  group_by(Quad, Year) %>%
  summarise(Q_grass = sum(area_cm2[cl == "grass"]), Q_shrub = sum(area_cm2[cl == "shrub"]),
            Q_forb = sum(cl == "forb"), .groups = "drop") %>%
  mutate(across(starts_with("Q_"), ~ as.numeric(scale(.x))))

ring_cols <- function(cl) paste(cl, seq_len(K), sep = "_")
fit_base <- genet %>%
  filter(Species %in% kernel_species, focal_ok, !is.na(survives_tplus1)) %>%
  select(Quad, Year, Species, trackID, Group, edge_dist_cm, area_cm2, area_tplus1_cm2, survives_tplus1, size_measure) %>%
  left_join(rings, by = c("Quad", "Year", "Species", "trackID")) %>%
  left_join(quad_total, by = c("Quad", "Year"))
  # Water-year precipitation of the transition, standardised over years (for delta)
ppt_year <- genet %>% distinct(Year, ppt_wy) %>% mutate(ppt_z = as.numeric(scale(ppt_wy)))
fit_base <- left_join(fit_base, select(ppt_year, Year, ppt_z), by = "Year")
stopifnot(!anyNA(fit_base$intra_1), !anyNA(fit_base$Q_grass), !anyNA(fit_base$ppt_z))
if (edge_correct) {
  pq <- as.matrix(fit_base[ring_cols("propQ")])
  for (cl in classes) fit_base[ring_cols(cl)] <- ifelse(pq > 0, as.matrix(fit_base[ring_cols(cl)]) / pq, 0)
}

make_part <- function(d, y) {
  f <- paste("~ factor(Group) + z", if (kernel_year) "+ factor(Year)",
             if (kernel_controls) "+ Q_grass + Q_shrub + Q_forb",
             if (kernel_border) "+ splines::ns(pmin(edge_dist_cm, 50), 3)")
  X0 <- model.matrix(as.formula(f), data = mutate(d, z = log(area_cm2)))
  X0 <- X0[, colSums(X0 != 0) > 0, drop = FALSE]
  C <- lapply(setNames(classes, classes), function(cl) as.matrix(d[ring_cols(cl)]))
  C <- C[sapply(C, function(m) any(m > 0))]   # classes absent around this species
  list(y = y, X0 = X0, C = C)
}
fit_data <- lapply(setNames(kernel_species, kernel_species), function(sp) {
  s <- fit_base[fit_base$Species == sp, ]
  g <- s[s$survives_tplus1 == 1 & !is.na(s$area_tplus1_cm2) & s$size_measure != "point", ]
  list(s = make_part(s, s$survives_tplus1), g = make_part(g, log(g$area_tplus1_cm2)))
})
cat(sprintf("Kernel fit data: %d survival and %d growth records.\n",
            sum(sapply(fit_data, function(d) length(d$s$y))), sum(sapply(fit_data, function(d) length(d$g$y)))))

# Kernel functions (Teller et al. 2016, Appendix S1) ----------------------------
  # Kernel F(d) = exp(-f(d)), f(d) = integral of H(g) from the first mid-radius to d,
  # g = cubic B-spline, H(g) = g+^2 / (0.02 + g+): F is positive, decreasing, F = 1 at
  # the first mid-radius.
d_grid <- seq(min(ring_mid), max(ring_mid), by = 0.02)
knots <- c(rep(min(kernel_breaks), 3), kernel_breaks, rep(max(kernel_breaks), 3))
B1 <- splines::splineDesign(knots, d_grid, ord = 4)
nb <- ncol(B1)
  # Second-derivative roughness penalty matrix, by the midpoint rule
h_pen <- 0.01
d_pen <- seq(min(kernel_breaks) + h_pen / 2, max(kernel_breaks) - h_pen / 2, by = h_pen)
B2 <- splines::splineDesign(knots, d_pen, ord = 4, derivs = rep(2, length(d_pen)))
P <- crossprod(B2) * h_pen

kernel_at <- function(p, x = ring_mid) {
  u <- pmax(B1 %*% p, 0)
  cdf <- cumsum(u^2 / (0.02 + u)) * (d_grid[2] - d_grid[1])
  exp(-approx(d_grid, cdf - cdf[1], xout = x, rule = 2)$y)
}
  # Derivative of the kernel values with respect to each spline coefficient (finite difference)
kernel_grad <- function(p, eps = 0.01) {
  sapply(seq_len(nb), function(j) {
    up <- p; up[j] <- up[j] + eps; dn <- p; dn[j] <- dn[j] - eps
    (kernel_at(up) - kernel_at(dn)) / (2 * eps)
  })
}

  # Survival (logistic) and growth (normal) regressions of every focal species, for a main
  # kernel and a forb kernel (NULL = forb crowding left out). The last survival fit of each
  # species starts the next one (5 times faster), but a start far from the new optimum can make
  # glm.fit stop at a wrong point while reporting convergence; the fit is kept only when the
  # score (gradient of the log-likelihood) is near zero, and is redone from scratch otherwise.
warm <- new.env()
fit_parts <- function(w_main, w_forb = NULL) {
  lapply(names(fit_data), function(sp) {
    lapply(c("s", "g"), function(v) {
      d <- fit_data[[sp]][[v]]
      C <- if (is.null(w_forb)) d$C[names(d$C) != "forb"] else d$C
      W <- sapply(names(C), function(cl) C[[cl]] %*% (if (cl == "forb") w_forb else w_main))
      X <- cbind(d$X0, matrix(W, nrow(d$X0), dimnames = list(NULL, names(C))))
      if (v == "s") {
        key <- paste(sp, ncol(X))
        f <- suppressWarnings(glm.fit(X, d$y, family = binomial(), start = warm[[key]]))
        ok <- !is.na(f$coefficients)
        if (!f$converged || max(abs(crossprod(X[, ok, drop = FALSE], d$y - f$fitted.values))) > 0.01)
          f <- suppressWarnings(glm.fit(X, d$y, family = binomial()))
        warm[[key]] <- ifelse(is.na(f$coefficients), 0, f$coefficients)
        mu <- f$fitted.values
        list(ll = sum(dbinom(d$y, 1, mu, log = TRUE)), beta = f$coefficients, X = X, C = C,
             omega = pmin(1 / (mu * (1 - mu)), 50))
      } else {
        f <- .lm.fit(X, d$y)
        n <- length(d$y); rss <- sum(f$residuals^2)
          # Coefficients back in column order (.lm.fit pivots aliased columns to the end)
        beta <- rep(NA_real_, ncol(X)); beta[f$pivot[seq_len(f$rank)]] <- f$coefficients[seq_len(f$rank)]
        list(ll = -n / 2 * (log(2 * pi * rss / n) + 1), beta = beta, X = X, C = C,
             omega = rep(n / rss, n))
      }
    })
  })
}
total_ll <- function(parts) sum(sapply(parts, function(s) sapply(s, `[[`, "ll")))
  # Kernels of one stage: 'main' fits the main kernel (forb crowding left out); 'forb' fits the
  # forb kernel with the main kernel fixed ('w_fixed'). 'reach' limits the kernel to distances
  # below it (kernel set to zero beyond).
stage_w <- function(p, stage, w_fixed, reach) {
  w <- kernel_at(p) * (ring_mid < reach)
  if (stage == "main") list(w_main = w, w_forb = NULL) else list(w_main = w_fixed, w_forb = w)
}
objective <- function(p, lambda, stage, w_fixed = NULL, reach = Inf) {
  w <- stage_w(p, stage, w_fixed, reach)
  -total_ll(fit_parts(w$w_main, w$w_forb)) + lambda * drop(t(p) %*% P %*% p)
}

  # Approximate AIC: effective degrees of freedom from the linearised penalised regression
  # (Appendix S1), trace((X'WX + lambda M)^-1 X'WX), built species by species
approx_aic <- function(p, lambda, stage, w_fixed = NULL, reach = Inf) {
  w <- stage_w(p, stage, w_fixed, reach); dB <- kernel_grad(p) * (ring_mid < reach)
  fitted_cl <- if (stage == "main") classes_main else "forb"
  parts <- fit_parts(w$w_main, w$w_forb)
  blocks <- unlist(parts, recursive = FALSE)
  q <- sapply(blocks, function(b) ncol(b$X))
  A <- matrix(0, nb + sum(q), nb + sum(q))
  pos <- nb + c(0, cumsum(q))
  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    bw <- b$beta[colnames(b$X) %in% names(b$C)]; names(bw) <- names(b$C); bw[is.na(bw)] <- 0
    S <- matrix(0, nrow(b$X), nb)
    for (cl in intersect(names(b$C), fitted_cl)) S <- S + bw[[cl]] * (b$C[[cl]] %*% dB)
    Xb <- b$X; Xb[, is.na(b$beta)] <- 0
    idx <- (pos[i] + 1):pos[i + 1]
    A[1:nb, 1:nb] <- A[1:nb, 1:nb] + crossprod(S, S * b$omega)
    A[1:nb, idx] <- crossprod(S, Xb * b$omega)
    A[idx, 1:nb] <- t(A[1:nb, idx])
    A[idx, idx] <- crossprod(Xb, Xb * b$omega)
  }
  M <- matrix(0, nrow(A), ncol(A)); M[1:nb, 1:nb] <- P
  dof <- sum(diag(MASS::ginv(A + lambda * M) %*% A))
  c(aic = -2 * total_ll(parts) + 2 * dof, dof = dof)
}

  # One penalty, reach 'R_ref': best of 'kernel_starts' random starts (Latin hypercube on normal
  # quantiles, as Teller et al. 2016), 'kernel_rough' rough NEWUOA fits, then a full fit
fit_full <- function(lambda, stage, w_fixed = NULL, reach = R_ref) {
  q0 <- 0.5 * qnorm(seq_len(kernel_starts) / (kernel_starts + 1))
  best <- list(par = rep(0, nb), fval = Inf)
  for (r in seq_len(kernel_rough)) {
    start <- sapply(seq_len(nb), function(j) sample(q0))
    fv <- apply(start, 1, objective, lambda = lambda, stage = stage, w_fixed = w_fixed, reach = reach)
    f <- suppressWarnings(minqa::newuoa(start[which.min(fv), ], objective, lambda = lambda, stage = stage,
                                        w_fixed = w_fixed, reach = reach,
                                        control = list(rhobeg = 0.5, maxfun = 1000, npt = 2 * nb + 1)))
    if (f$fval < best$fval) best <- f
  }
  f <- minqa::newuoa(best$par, objective, lambda = lambda, stage = stage, w_fixed = w_fixed, reach = reach,
                     control = list(rhobeg = 0.25, maxfun = 5000, npt = 2 * nb + 1))
  f <- minqa::bobyqa(f$par, objective, lambda = lambda, stage = stage, w_fixed = w_fixed, reach = reach,
                     control = list(rhobeg = 0.1, maxfun = 5000, npt = 2 * nb + 1))
  c(list(par = f$par), as.list(approx_aic(f$par, lambda, stage, w_fixed, reach)))
}
  # Limited reach: started from the full-reach kernel and from the best of 'kernel_starts'
  # random starts; deviance (-2 log-likelihood) of the best fit
fit_reach <- function(reach, lambda, stage, start, w_fixed = NULL) {
  q0 <- 0.5 * qnorm(seq_len(kernel_starts) / (kernel_starts + 1))
  rand <- sapply(seq_len(nb), function(j) sample(q0))
  fv <- apply(rand, 1, objective, lambda = lambda, stage = stage, w_fixed = w_fixed, reach = reach)
  fits <- lapply(list(start, rand[which.min(fv), ]), function(s0)
    suppressWarnings(minqa::newuoa(s0, objective, lambda = lambda, stage = stage, w_fixed = w_fixed,
                                   reach = reach, control = list(rhobeg = 0.25, maxfun = 2000, npt = 2 * nb + 1))))
  f <- fits[[which.min(sapply(fits, `[[`, "fval"))]]
  f <- minqa::bobyqa(f$par, objective, lambda = lambda, stage = stage, w_fixed = w_fixed, reach = reach,
                     control = list(rhobeg = 0.1, maxfun = 3000, npt = 2 * nb + 1))
  w <- stage_w(f$par, stage, w_fixed, reach)
  list(par = f$par, deviance = -2 * total_ll(fit_parts(w$w_main, w$w_forb)))
}
  # Crowding effects for a pair of kernels: per species and vital rate, the regressions of the
  # kernel fit (log size, group, year, quadrat-year totals) with every crowding term scaled to
  # SD 1 (beta) and its product with standardised precipitation (delta); rows for 'report' classes
crowd_effects <- function(w_main, w_forb, report) {
  bind_rows(lapply(kernel_species, function(sp) {
    s <- fit_base[fit_base$Species == sp, ]
    g <- s[s$survives_tplus1 == 1 & !is.na(s$area_tplus1_cm2) & s$size_measure != "point", ]
    parts <- list(survival = list(d = s, y = s$survives_tplus1, fam = binomial()),
                  growth = list(d = g, y = log(g$area_tplus1_cm2), fam = gaussian()))
    bind_rows(lapply(names(parts), function(v) {
      d <- parts[[v]]$d
      pt <- make_part(d, parts[[v]]$y)
      C <- if (is.null(w_forb)) pt$C[names(pt$C) != "forb"] else pt$C
      W <- sapply(names(C), function(cl) C[[cl]] %*% (if (cl == "forb") w_forb else w_main))
      W <- matrix(W, nrow(pt$X0), dimnames = list(NULL, names(C)))
      W <- scale(W[, apply(W, 2, sd) > 0, drop = FALSE])
      WP <- W * d$ppt_z; colnames(WP) <- paste0(colnames(W), ":ppt")
      X <- cbind(pt$X0, W, WP)
      ct <- coef(summary(glm(pt$y ~ X - 1, family = parts[[v]]$fam)))
      rn <- sub("^X", "", rownames(ct))
      keep <- rn %in% c(intersect(colnames(W), report), paste0(intersect(colnames(W), report), ":ppt"))
      data.frame(Species = sp, vital_rate = v, term = rn[keep], estimate = ct[keep, 1], se = ct[keep, 2])
    }))
  }))
}
  # Range of one stage: the kernel refitted with every candidate reach (the reference reach reuses
  # the fit at R_ref); crowding effects at each reach compared with those at R_ref, in SE units
  # of the reference estimate
fit_stability <- function(lambda, stage, ref, w_fixed = NULL, t0, extra = numeric(0)) {
  reaches <- c(R_candidates, R_ref, extra)
  report <- if (stage == "main") classes_main else "forb"
  fits <- lapply(reaches, function(r) {
    f <- if (r == R_ref) {
      w <- stage_w(ref$par, stage, w_fixed, R_ref)
      list(par = ref$par, deviance = -2 * total_ll(fit_parts(w$w_main, w$w_forb)))
    } else fit_reach(r, lambda, stage, ref$par, w_fixed)
    cat(sprintf("  %s kernel, reach %g cm: deviance %.1f (%.0f min so far)\n", stage, r, f$deviance,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    f
  })
  kern <- sapply(seq_along(reaches), function(i) kernel_at(fits[[i]]$par) * (ring_mid < reaches[i]))
  eff <- bind_rows(lapply(seq_along(reaches), function(i) {
    w <- stage_w(fits[[i]]$par, stage, w_fixed, reaches[i])
    crowd_effects(w$w_main, w$w_forb, report) %>% mutate(stage = stage, reach_cm = reaches[i])
  }))
  cmp <- eff %>% filter(reach_cm != R_ref) %>%
    inner_join(filter(eff, reach_cm == R_ref) %>% select(Species, vital_rate, term, est_ref = estimate, se_ref = se),
               by = c("Species", "vital_rate", "term")) %>%
    mutate(type = ifelse(grepl(":ppt", term), "delta", "beta"), change_in_se = (estimate - est_ref) / se_ref,
           flip = abs(est_ref / se_ref) > 2 & sign(estimate) != sign(est_ref))
  tab <- cmp %>% group_by(stage, reach_cm) %>%
    summarise(max_change_beta = max(abs(change_in_se[type == "beta"])), max_change_delta = max(abs(change_in_se[type == "delta"])),
              median_change_beta = median(abs(change_in_se[type == "beta"])), median_change_delta = median(abs(change_in_se[type == "delta"])),
              n_over_1se = sum(abs(change_in_se) > 1), n_over_2se = sum(abs(change_in_se) > 2), sign_flips = sum(flip),
              n_coef = n(), .groups = "drop") %>%
    left_join(data.frame(reach_cm = reaches, deviance = sapply(fits, `[[`, "deviance")), by = "reach_cm") %>%
    mutate(delta_deviance_vs_ref = deviance - fits[[match(R_ref, reaches)]]$deviance)
  ok <- tab$reach_cm[tab$reach_cm %in% R_candidates & pmax(tab$max_change_beta, tab$max_change_delta) <= R_se]
  R <- if (length(ok)) min(ok) else R_ref
  tab$chosen <- tab$reach_cm == R
  list(table = tab, effects = eff, R = R, kernels = kern, reaches = reaches,
       par = fits[[match(R, reaches)]]$par)
}

# Fit the kernels --------------------------------------------------------------
kernel_file <- file.path(der_dir, "CLR_CrowdKernel_fit.rds")
fit_set <- list(ring_set = ring_set, kernel_species = kernel_species, kernel_breaks = kernel_breaks,
                kernel_lambda = kernel_lambda, kernel_starts = kernel_starts, kernel_rough = kernel_rough,
                kernel_year = kernel_year, kernel_controls = kernel_controls, kernel_border = kernel_border,
                edge_correct = edge_correct, R_candidates = R_candidates, R_ref = R_ref, R_se = R_se, R_full = R_full,
                kernel_seed = kernel_seed)
kernel_fit <- if (file.exists(kernel_file)) readRDS(kernel_file) else NULL
if (identical(kernel_fit$settings, fit_set)) {
  cat("Read the kernel fits saved before (same settings).\n")
} else {
  t0 <- Sys.time()
  set.seed(kernel_seed)
    # (1) Main kernel, reach R_ref, every penalty
  fits <- lapply(kernel_lambda, function(l) {
    f <- fit_full(l, "main")
    cat(sprintf("  main kernel, lambda 10^%g: AIC %.1f, DoF %.1f (%.0f min so far)\n", log10(l), f$aic, f$dof,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    f
  })
    # Smoothest local minimum of AIC along the effective degrees of freedom
    # (Teller et al. 2016); the end of the curve counts when its AIC is lower than its neighbour
  aic <- sapply(fits, `[[`, "aic"); dof <- sapply(fits, `[[`, "dof")
  o <- order(dof)
  best <- o[which(diff(sign(diff(c(Inf, aic[o], Inf)))) > 0)[1]]
  lambda_best <- kernel_lambda[best]
    # (2) Main kernel, range (and the full reach for the record)
  stab_main <- fit_stability(lambda_best, "main", fits[[best]], t0 = t0, extra = R_full)
  w_main <- kernel_at(stab_main$par) * (ring_mid < stab_main$R)
    # (3) Forb kernel with the main kernel fixed (same penalty), reach R_ref and range
  ref_forb <- fit_full(lambda_best, "forb", w_main)
  stab_forb <- fit_stability(lambda_best, "forb", ref_forb, w_main, t0)
  w_forb <- kernel_at(stab_forb$par) * (ring_mid < stab_forb$R)
  kernel_fit <- list(lambda = kernel_lambda, aic = aic, dof = dof, best = best,
                     kernel_main_lambda = sapply(fits, function(f) kernel_at(f$par) * (ring_mid < R_ref)),
                     stability = bind_rows(stab_main$table, stab_forb$table),
                     effects = bind_rows(stab_main$effects, stab_forb$effects),
                     kernels_main = stab_main$kernels %>% `colnames<-`(paste0("main_reach_", stab_main$reaches)),
                     kernels_forb = stab_forb$kernels %>% `colnames<-`(paste0("forb_reach_", stab_forb$reaches)),
                     R_main = stab_main$R, R_forb = stab_forb$R, w_main = w_main, w_forb = w_forb,
                     ring_edges = ring_edges, ring_mid = ring_mid, settings = fit_set)
  saveRDS(kernel_fit, kernel_file)
  cat(sprintf("Kernel fits: %.0f min.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
w_main <- kernel_fit$w_main; w_forb <- kernel_fit$w_forb
R_main <- kernel_fit$R_main; R_forb <- kernel_fit$R_forb
R_edge <- max(R_main, R_forb)
cat(sprintf("Chosen penalty 10^%g; range R = %g cm (main kernel), %g cm (forb kernel).\n",
            log10(kernel_fit$lambda[kernel_fit$best]), R_main, R_forb))
print(as.data.frame(kernel_fit$stability %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
  # Share of focal plants (size species) whose outline stays at least R from the border
focal_R <- genet %>% filter(Species %in% kernel_species, !is.na(survives_tplus1))
cat(sprintf("Plants of the focal species at least 5 cm from the border: %.0f%%; at least %g cm: %.0f%%.\n",
            100 * mean(focal_R$edge_dist_cm >= 5), R_edge, 100 * mean(focal_R$edge_dist_cm >= R_edge)))

# Check the 5-cm edge rule -----------------------------------------------------
  # The kernels were fitted with focal plants at least 5 cm from the border (R unknown then).
  # With the chosen kernels, the regressions are refitted with the focal plants at least R from
  # the border (the larger R of the two kernels), whose neighbourhood lies inside the quadrat.
  # The crowding coefficients are compared by the z value of their difference: the R-cm plants are
  # a subset of the 5-cm plants, so the variance of the difference is about the variance of the
  # subset estimate minus that of the full estimate (|z| > 2: the 5-cm rule shifts the estimate
  # beyond sampling noise).
crowd_coefs <- function(base) {
  empty <- data.frame(Species = character(), vital_rate = character(), class = character(),
                      estimate = numeric(), se = numeric(), n = integer())
  bind_rows(empty, lapply(kernel_species, function(sp) {
    s <- base[base$Species == sp, ]
    g <- s[s$survives_tplus1 == 1 & !is.na(s$area_tplus1_cm2) & s$size_measure != "point", ]
    parts <- list(survival = list(d = s, y = s$survives_tplus1, fam = binomial()),
                  growth = list(d = g, y = log(g$area_tplus1_cm2), fam = gaussian()))
    bind_rows(lapply(names(parts), function(v) {
      d <- parts[[v]]$d
      if (nrow(d) < 50) return(NULL)   # too few plants to refit
      pt <- make_part(d, parts[[v]]$y)
      W <- sapply(names(pt$C), function(cl) pt$C[[cl]] %*% (if (cl == "forb") w_forb else w_main))
      X <- cbind(pt$X0, matrix(W, nrow(pt$X0), dimnames = list(NULL, names(pt$C))))
      ct <- coef(summary(glm(pt$y ~ X - 1, family = parts[[v]]$fam)))
      ct <- ct[rownames(ct) %in% paste0("X", names(pt$C)), , drop = FALSE]
      data.frame(Species = sp, vital_rate = v, class = sub("^X", "", rownames(ct)),
                 estimate = ct[, 1], se = ct[, 2], n = nrow(d))
    }))
  }))
}
edge_tab <- full_join(crowd_coefs(fit_base),
                      crowd_coefs(fit_base[fit_base$edge_dist_cm >= R_edge, ]),
                      by = c("Species", "vital_rate", "class"), suffix = c("_5cm", "_R")) %>%
  mutate(code = species$code[match(Species, species$Species)],
         z_diff = (estimate_R - estimate_5cm) / sqrt(pmax(se_R^2 - se_5cm^2, 1e-12)))
n_check <- sum(!is.na(edge_tab$z_diff))
cat(sprintf("Edge check (plants at least %g cm from the border): %d of %d crowding coefficients refitted%s.\n",
            R_edge, n_check, nrow(edge_tab),
            if (n_check > 0) sprintf("; %d with |z| > 2 (about %.0f expected by chance)",
                                     sum(abs(edge_tab$z_diff) > 2, na.rm = TRUE), 0.046 * n_check) else ""))

# Save tables and figure -------------------------------------------------------
kernel_tab <- data.frame(ring_inner_cm = head(ring_edges, -1), ring_outer_cm = tail(ring_edges, -1),
                         mid_cm = ring_mid, kernel_main = w_main, kernel_forb = w_forb,
                         kernel_fit$kernels_main, kernel_fit$kernels_forb)
aic_tab <- data.frame(log10_lambda = log10(kernel_fit$lambda), aic = kernel_fit$aic, dof = kernel_fit$dof,
                      chosen = seq_along(kernel_fit$lambda) == kernel_fit$best)
write.csv(kernel_tab, file.path(tab_dir, "CLR_CrowdKernel_kernel.csv"), row.names = FALSE)
write.csv(aic_tab, file.path(tab_dir, "CLR_CrowdKernel_aic.csv"), row.names = FALSE)
write.csv(kernel_fit$stability, file.path(tab_dir, "CLR_CrowdKernel_range.csv"), row.names = FALSE)
write.csv(mutate(kernel_fit$effects, code = species$code[match(Species, species$Species)]),
          file.path(tab_dir, "CLR_CrowdKernel_effects.csv"), row.names = FALSE)
write.csv(edge_tab, file.path(tab_dir, "CLR_CrowdKernel_edge.csv"), row.names = FALSE)

lab_main <- sprintf("main (same species, grass, shrub), R = %g cm", R_main)
lab_forb <- sprintf("forb, R = %g cm", R_forb)
plot_k <- bind_rows(
  data.frame(x = ring_mid, y = w_main, kernel = lab_main, line = "chosen"),
  data.frame(x = ring_mid, y = w_forb, kernel = lab_forb, line = "chosen"),
  data.frame(x = ring_mid, y = kernel_fit$kernels_main[, paste0("main_reach_", R_ref)], kernel = lab_main,
             line = sprintf("reference reach %g cm", R_ref)),
  data.frame(x = ring_mid, y = kernel_fit$kernels_forb[, paste0("forb_reach_", R_ref)], kernel = lab_forb,
             line = sprintf("reference reach %g cm", R_ref))) %>%
  filter(x <= 40) %>%
  mutate(panel = "(a) Kernel (relative interaction strength) by distance from the plant centre (cm)")
plot_s <- kernel_fit$stability %>%
  filter(reach_cm <= R_ref) %>%
  pivot_longer(c(max_change_beta, max_change_delta), names_to = "line", values_to = "y") %>%
  transmute(x = reach_cm, y, kernel = ifelse(stage == "forb", lab_forb, lab_main),
            line = ifelse(line == "max_change_beta", "beta", "delta"), chosen,
            panel = sprintf("(b) Largest change of beta and delta from the %g-cm reach (SE), by reach (cm)", R_ref))
p_kernel <- ggplot(mapping = aes(x, y, colour = kernel)) +
  geom_step(data = plot_k, aes(linetype = line, group = interaction(kernel, line)), direction = "mid") +
  geom_hline(data = data.frame(panel = plot_s$panel[1], h = R_se), aes(yintercept = h), linetype = 2, colour = "grey40") +
  geom_line(data = plot_s, aes(linetype = line, group = interaction(kernel, line))) +
  geom_point(data = plot_s, aes(shape = chosen), size = 2.5) +
  facet_wrap(~ panel, scales = "free") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8), guide = "none") +
  labs(x = NULL, y = NULL, colour = NULL, linetype = NULL, title = "Crowding kernels, Idaho",
       subtitle = sprintf("Star: chosen range R (shortest reach with every beta and delta within %g SE of the %g-cm reach)", R_se, R_ref)) +
  theme_bw() + theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(fig_dir, "CLR_CrowdKernel.png"), p_kernel, width = 13, height = 5.5, dpi = 300)

cat("Saved:\n",
    " data/data.derived/", site_tag, ": CLR_CrowdKernel_rings.rds, CLR_CrowdKernel_fit.rds\n",
    " result/tables/", site_tag, ": CLR_CrowdKernel_kernel.csv, _aic.csv, _range.csv, _effects.csv, _edge.csv\n",
    " result/figs/", site_tag, ": CLR_CrowdKernel.png\n", sep = "")
rm(list = intersect(c("fit_base", "fit_data", "warm", "B1", "B2", "P", "d_grid", "d_pen", "knots", "quad_total",
                      "plot_k", "plot_s", "focal_R", "pq", "t0", "crowd_coefs", "n_check", "fits", "ppt_year",
                      "lab_main", "lab_forb", "stab_main", "stab_forb", "ref_forb"), ls()))
