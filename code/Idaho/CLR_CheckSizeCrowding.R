################################################################################
## Check of the size x crowding terms (Part A, Idaho): grasses only vs all species
##
## Written by Hung-wei Lin Nov 2026
##
## The selected survival and growth models hold size x crowding (same species, grasses) with negative
## estimates. Crowding is measured from the focal centroid, so large shrubs (radius up to ~7 cm) have
## no neighbours inside their own canopy; grasses are small (90% below ~2 cm radius). If the terms keep
## their sign with grasses only, the shrub geometry does not drive them. Selected GLMMs of
## CLR_SelectPartA (same likelihoods and random effects). Table in 'PartA_model.selection'.


# Settings ---------------------------------------------------------------------
check_terms <- c("z:W_intra_z", "z:W_grass_z")
tab_dir <- here("result", "tables", site_tag, "PartA_model.selection")

input <- readRDS(here("data", "data.derived", site_tag, "CLR_PartA_input.rds"))
setup <- readRDS(file.path(sel_dir(), "CLR_select_setup.rds"))
sel <- read.csv(file.path(sel_dir(), "CLR_select_selected.csv"))
grass_sp <- input$species$code[input$species$growth_form == "grass"]

# Fit: all species and grasses only ----------------------------------------------
res <- NULL; pred <- NULL
for (v in c("survival", "growth")) {
  s <- sel[sel$vr == v, ]
  d_all <- vr_data(v, input, setup$p_eff, s$clim_use)
  for (set in c("all species", "grasses")) {
    d <- if (set == "grasses") droplevels(d_all[d_all$sp %in% grass_sp, ]) else d_all
    f <- model_formula(v, to_codes(s$terms), s$clim_use, re_drop = to_codes(setup$re_drop[[v]]))
    fit <- fit_model(v, f, d); m <- fit$model
    b <- fixef(m)$cond; se <- sqrt(diag(vcov(m)$cond))
    res <- rbind(res, data.frame(vr = v, species = set, term = check_terms, estimate = b[check_terms],
                                 se = se[check_terms], n = nrow(d), ok = fit$ok))
      # Survival on the probability scale: grasses at the median and the 90th percentile of size,
      # crowding by the same species from 0 to +2 SD (other covariates 0, mean species intercept and size slope)
    if (v == "survival" && set == "grasses") {
      sps <- levels(d$sp)
      a0 <- mean(b[paste0("sp", sps)]); a1 <- mean(b[paste0("sp", sps, ":z")])
      z <- quantile(d$z, c(0.5, 0.9))
      pred <- expand.grid(size = c("median", "90th percentile"), W_intra_sd = 0:2)
      zz <- z[match(pred$size, c("median", "90th percentile"))]
      eta <- a0 + a1 * zz + b["W_intra_z"] * pred$W_intra_sd + b["z:W_intra_z"] * zz * pred$W_intra_sd
      pred$survival <- plogis(eta)
    }
  }
}
res$z_value <- res$estimate / res$se
write.csv(res, file.path(tab_dir, "CLR_PartA_check_size_crowding.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)
cat("\nPredicted grass survival (probability):\n"); print(pred, digits = 3, row.names = FALSE)
rm(check_terms, tab_dir, input, setup, sel, grass_sp, res, pred, v, s, d_all, set, d, f, fit, m, b, se, sps, a0, a1, z, zz, eta)
