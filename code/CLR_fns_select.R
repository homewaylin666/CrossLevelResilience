################################################################################
## Functions for the model selection of Part A
##
## Written by Hung-wei Lin Oct 2026
##
## glmmTMB versions of the vital-rate models (maximum likelihood; Brooks et al. 2017) with the
## likelihoods of the NIMBLE model: survival binomial, growth Gaussian (residual SD per species),
## recruitment negative binomial (dispersion per species, offset log C').
## The main effects (size, grazing, climate, four crowding classes, each with species-level slopes)
## are in every model. The selection is over the climate variable and the interaction terms; every
## subset of interaction terms is fitted, as the candidate sets of Gupta et al. (2025).


# Interaction terms ------------------------------------------------------------
  # Codes: C climate, G grazing, S size, Wi/Wg/Ws/Wf crowding by the same species, grasses, shrubs, forbs
W_names <- c("W_intra_z", "W_grass_z", "W_shrub_z", "W_forb_z")
ix_terms <- data.frame(
  code  = c("CxWi", "CxWg", "CxWs", "CxWf", "CxS", "SxWi", "SxWg", "SxWs", "SxWf",
            "GxWi", "GxWg", "GxWs", "GxWf", "GxS", "GxC"),
  term  = c(paste0("clim:", W_names), "z:clim", paste0("z:", W_names), paste0("E:", W_names), "z:E", "E:clim"),
  label = c(paste("Climate × crowding", c("(same sp.)", "(grasses)", "(shrubs)", "(forbs)")), "Climate × size",
            paste("Size × crowding", c("(same sp.)", "(grasses)", "(shrubs)", "(forbs)")),
            paste("Grazing × crowding", c("(same sp.)", "(grasses)", "(shrubs)", "(forbs)")),
            "Grazing × size", "Grazing × climate"))
  # Recruitment has no size terms (15 terms for survival and growth, 9 for recruitment)
terms_of <- function(vr) if (vr == "recruitment") ix_terms[!grepl("S", ix_terms$code), ] else ix_terms
  # Trait terms of model_after: the trait moderates the slopes of grazing, climate and crowding
trait_terms <- data.frame(code = c("TxG", "TxC", "TxWi", "TxWg", "TxWs", "TxWf"),
                          term = c("trait:E", "trait:clim", paste0("trait:", W_names)),
                          label = c("Trait × grazing", "Trait × climate",
                                    paste("Trait × crowding", c("(same sp.)", "(grasses)", "(shrubs)", "(forbs)"))))
  # Climate variable: precipitation of the transition water year, of the previous water year, both
  # (previous year as an extra main effect), or summer temperature in place of precipitation
clim_choices <- c("current", "lag", "both", "temp")
vrs <- c("survival", "growth", "recruitment")
vr_short <- c(survival = "surv", growth = "grow", recruitment = "recr")
sel_dir <- function() here("data", "data.derived", site_tag, "selection")
  # Term codes kept as one string ("CxWi+GxS"; "" = none)
to_codes <- function(s) if (is.na(s) || s == "") character(0) else strsplit(s, "+", fixed = TRUE)[[1]]


# Data of one vital rate -------------------------------------------------------
  # Recruitment: C' = p * quadrat cover + (1 - p) * group cover, with p at the posterior median of the
  # full Bayesian model; rows with C' = 0 hold only background recruitment, so they are left out here.
vr_data <- function(vr, input, p_eff, clim_use) {
  yrs <- input$years
  if (vr == "recruitment") {
    d <- input$main$recruit
    d$Ceff <- p_eff[d$k] * d$cover_q + (1 - p_eff[d$k]) * d$cover_g
    d <- d[d$Ceff > 0, ]
    cs <- input$scaling$C_intra[match(input$species$code[d$k], input$scaling$C_intra$code), ]
    d$W_intra_z <- (sqrt(d$Ceff) - cs$mean) / cs$sd
    d$W_grass_z <- d$C_grass_z; d$W_shrub_z <- d$C_shrub_z; d$W_forb_z <- d$C_forb_z
    d$y <- d$recruits
  } else {
    d <- input$main$rec
    if (vr == "growth") d <- d[d$grow_ok, ]
    d$y <- if (vr == "survival") d$survives_tplus1 else d$z_next
    d$z2 <- d$z^2
    d$plant <- factor(d$i)
  }
  d$clim <- switch(clim_use, lag = yrs$clim_lag_anom, temp = yrs$temp_anom, yrs$clim_anom)[d$t]
  d$clim_lag <- yrs$clim_lag_anom[d$t]
  d$sp <- factor(input$species$code[d$k], levels = input$species$code)
  d$yr <- factor(d$t); d$ys <- factor(paste(d$t, d$k)); d$grp <- factor(d$g); d$quad <- factor(d$q)
  d
}


# Model formula ----------------------------------------------------------------
  # Main effects plus the interaction terms 'codes' (and the trait terms 'tcodes' for model_after).
  # 're_drop': random effects left out, e.g. "grp", "sp:clim" (species slope), "yr:W_shrub_z" (yearly slope)
model_formula <- function(vr, codes, clim_use = "current", tcodes = NULL, re_drop = NULL) {
  sz <- vr != "recruitment"
  fixed <- c("0", "sp", if (sz) c("sp:z", "sp:z2"), "E", "clim", if (clim_use == "both") "clim_lag", W_names,
             ix_terms$term[match(codes, ix_terms$code)], trait_terms$term[match(tcodes, trait_terms$code)])
  int <- setdiff(c("yr", "ys", "grp", "quad", if (sz) "plant"), re_drop)
  yr_s <- setdiff(c(if (sz) "z", W_names), sub("yr:", "", re_drop))
  sp_s <- setdiff(c("E", "clim", W_names), sub("sp:", "", re_drop))
  random <- c(sprintf("(1 | %s)", int),
              if (length(yr_s)) sprintf("diag(0 + %s | yr)", paste(yr_s, collapse = " + ")),
              if (length(sp_s)) sprintf("diag(0 + %s | sp)", paste(sp_s, collapse = " + ")))
  as.formula(paste("y ~", paste(c(fixed, random, if (!sz) "offset(log(Ceff))"), collapse = " + ")))
}


# Fit one model ----------------------------------------------------------------
  # glmmTMB gives no likelihood when the Hessian is non-positive-definite; such fits are left out of
  # the ranking. 'keep_vcov': keep the estimates and covariance of the interaction terms (full models).
fit_model <- function(vr, f, d, keep_vcov = FALSE) {
  t0 <- Sys.time()
  m <- suppressWarnings(glmmTMB(f, data = d, dispformula = if (vr == "survival") ~1 else ~0 + sp,
                                family = switch(vr, survival = binomial(), growth = gaussian(), recruitment = nbinom2()),
                                control = glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000), parallel = 1)))
  k <- attr(logLik(m), "df"); n <- nrow(d)
  res <- list(ok = isTRUE(m$sdr$pdHess) && m$fit$convergence == 0, df = k, logLik = as.numeric(logLik(m)), n = n,
              AICc = AIC(m) + 2 * k * (k + 1) / (n - k - 1), BIC = BIC(m),
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")), model = m)
  if (keep_vcov && res$ok) {
    b <- fixef(m)$cond; tm <- terms_of(vr)
      # (glmmTMB may name an interaction either way round, e.g. "clim:z" for "z:clim")
    sorted <- function(x) sapply(strsplit(x, ":"), function(p) paste(sort(p), collapse = ":"))
    idx <- match(sorted(tm$term), sorted(names(b)))
    res$ix_coef <- setNames(b[idx], tm$code)
    res$ix_vcov <- vcov(m)$cond[idx, idx]
  }
  res
}

  # Leave-one-year-out log predictive density: fit without year y, predict year y as a new year
loyo_lpd <- function(vr, f, d, y) {
  fit <- fit_model(vr, f, d[d$t != y, ])
  test <- d[d$t == y, ]
  mu <- predict(fit$model, newdata = test, type = "response", allow.new.levels = TRUE)
  disp <- if (vr != "survival") predict(fit$model, newdata = test, type = "disp", allow.new.levels = TRUE)
  lpd <- switch(vr, survival = dbinom(test$y, 1, mu, log = TRUE), growth = dnorm(test$y, mu, disp, log = TRUE),
                recruitment = dnbinom(test$y, mu = mu, size = disp, log = TRUE))
  list(lpd = sum(lpd), ok = fit$ok, minutes = fit$minutes)
}

  # SD of every random effect, named by the codes of 're_drop'
re_sd <- function(m) {
  vc <- VarCorr(m)$cond
  unlist(lapply(names(vc), function(nm) {
    sds <- attr(vc[[nm]], "stddev"); g <- sub("\\.[0-9]+$", "", nm)
    setNames(sds, ifelse(names(sds) == "(Intercept)", g, paste0(g, ":", names(sds))))
  }))
}


# Set-up of the selection ------------------------------------------------------
  # p of C' (posterior median, full Bayesian model) and the random effects left out of every glmmTMB
  # model of a vital rate: while its full model has a non-positive-definite Hessian, the random effect
  # with the smallest SD (at or near zero) is left out (Bates et al. 2015). Saved in 'CLR_select_setup.rds'.
select_setup <- function() {
  input <- readRDS(here("data", "data.derived", site_tag, "CLR_PartA_input.rds"))
  draws <- as.matrix(readRDS(here("data", "data.derived", site_tag, "CLR_PartA_stage1_full.rds"))$samples)
  p_eff <- apply(draws[, sprintf("p.eff[%d]", seq_len(nrow(input$species)))], 2, median)
  re_drop <- sapply(vrs, function(v) {
    d <- vr_data(v, input, p_eff, "current"); drop <- character(0)
    repeat {
      fit <- fit_model(v, model_formula(v, terms_of(v)$code, re_drop = drop), d)
      if (fit$ok) break
      sds <- re_sd(fit$model); drop <- c(drop, names(which.min(sds)))
      cat(sprintf("%s: left out %s (SD %.2g)\n", v, tail(drop, 1), min(sds)))
    }
    paste(drop, collapse = "+")
  })
  saveRDS(list(p_eff = p_eff, re_drop = re_drop), file.path(sel_dir(), "CLR_select_setup.rds"))
}


# Tasks (one model fit each) ----------------------------------------------------
  # Phase 1: full model with each climate variable, plus leave-one-year-out refits (22 years)
  # Phase 2: every subset of interaction terms, with the chosen climate variable
  # Phase 3: model_after (exploratory): every subset of trait terms on the selected model
make_tasks <- function(phase) {
  if (phase == 1) {
    tk <- expand.grid(year_out = c(NA, 1:22), clim_use = clim_choices, vr = vrs, stringsAsFactors = FALSE)
    tk$model <- paste0(vr_short[tk$vr], "_full")
    tk$terms <- sapply(tk$vr, function(v) paste(terms_of(v)$code, collapse = "+"))
    tk$tterms <- ""
  } else if (phase == 2) {
    clim <- read.csv(file.path(sel_dir(), "CLR_select_climate_chosen.csv"))
    tk <- do.call(rbind, lapply(vrs, function(v) {
      cd <- terms_of(v)$code
      bits <- as.matrix(expand.grid(rep(list(0:1), length(cd))))
      data.frame(vr = v, model = sprintf("%s_%05d", vr_short[v], bits %*% 2^(seq_along(cd) - 1)),
                 terms = apply(bits, 1, function(r) paste(cd[r == 1], collapse = "+")), tterms = "",
                 clim_use = clim$clim_use[clim$vr == v], year_out = NA)
    }))
  } else {
    sel <- read.csv(file.path(sel_dir(), "CLR_select_selected.csv"))
    bits <- as.matrix(expand.grid(rep(list(0:1), nrow(trait_terms))))
    tt <- apply(bits, 1, function(r) paste(trait_terms$code[r == 1], collapse = "+"))
    tk <- do.call(rbind, lapply(vrs, function(v) {
      s <- sel[sel$vr == v, ]
      data.frame(vr = v, model = sprintf("%s_after%02d", vr_short[v], seq_along(tt) - 1), terms = s$terms,
                 tterms = tt, clim_use = s$clim_use, year_out = NA)
    }))
  }
  tk$terms[is.na(tk$terms)] <- ""
  tk$file <- sprintf("phase%d_%s_%s%s.rds", phase, tk$model, tk$clim_use,
                     ifelse(is.na(tk$year_out), "", sprintf("_out%02d", tk$year_out)))
  write.csv(tk, file.path(sel_dir(), sprintf("CLR_select_tasks_phase%d.csv", phase)), row.names = FALSE)
  tk
}

  # Run one task and save its summary (the model object is left out: it is large)
run_task <- function(phase, i, trait_file = NULL) {
  tk <- read.csv(file.path(sel_dir(), sprintf("CLR_select_tasks_phase%d.csv", phase)))[i, ]
  if (file.exists(file.path(sel_dir(), tk$file))) return(invisible(NULL))
  input <- readRDS(here("data", "data.derived", site_tag, "CLR_PartA_input.rds"))
  setup <- readRDS(file.path(sel_dir(), "CLR_select_setup.rds"))
  d <- vr_data(tk$vr, input, setup$p_eff, tk$clim_use)
  if (phase == 3) {
    tr <- read.csv(trait_file)
    d$trait <- as.numeric(scale(tr$trait))[match(as.character(d$sp), tr$code)]
  }
  f <- model_formula(tk$vr, to_codes(tk$terms), tk$clim_use, to_codes(tk$tterms), to_codes(setup$re_drop[[tk$vr]]))
  res <- if (is.na(tk$year_out)) fit_model(tk$vr, f, d, keep_vcov = phase == 1) else loyo_lpd(tk$vr, f, d, tk$year_out)
  res$model <- NULL
  saveRDS(c(as.list(tk), res), file.path(sel_dir(), tk$file))
}

  # TRUE when every task of a phase has its summary
phase_done <- function(phase) {
  f <- file.path(sel_dir(), sprintf("CLR_select_tasks_phase%d.csv", phase))
  file.exists(f) && all(file.exists(file.path(sel_dir(), read.csv(f)$file)))
}


# Wald approximation (check of the exact fits) --------------------------------
  # Dropping the set D of interaction terms from the full model changes -2 log-likelihood by about
  # b_D' V_DD^-1 b_D (Lawless & Singhal 1978); compared with the exact AICc of every subset.
wald_aicc <- function(full) {
  b <- full$ix_coef; V <- full$ix_vcov
  bits <- as.matrix(expand.grid(rep(list(0:1), length(b))))
  W <- apply(bits, 1, function(keep) { dr <- keep == 0; if (any(dr)) sum(b[dr] * solve(V[dr, dr, drop = FALSE], b[dr])) else 0 })
  k <- full$df - rowSums(bits == 0)
  data.frame(subset_id = as.vector(bits %*% 2^(seq_along(b) - 1)),
             AICc_wald = -2 * full$logLik + W + 2 * k + 2 * k * (k + 1) / (full$n - k - 1))
}
