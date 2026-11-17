################################################################################
## Part A forb submodel (Idaho): NIMBLE code, data, inits and MCMC
##
## Written by Hung-wei Lin Oct 2026
##
## All point species pooled as one forb population (point records carry no size), fitted
## separately from the main model (the two share no parameters). Survival of forb plants and
## recruitment per quadrat x transition, with the predictors of the main model except size and
## trait: crowding by forbs (forb kernel), grasses and shrubs (main kernel), grazing E, water-year
## precipitation and its modulation of each crowding effect (delta_x), year effects and quadrat
## (in group) effects. Gives the forb dynamics for Part B (optional there).
## Reads 'CLR_PartA_input.rds' (CLR_InputPartA.R); uses 'site_tag' and 'run_chains' (CLR_fns_nimble.R).


# Settings ---------------------------------------------------------------------
n_chains <- 4
n_iter   <- 20000
n_burn   <- 10000
n_thin   <- 10
  # Cores for the chains: 1 = one chain after another (progress bar in the console); 2-4 = chains
  # in parallel (each core builds and compiles its own copy of the model, so memory grows with
  # it; the progress bars go to a log file in 'data/data.derived/Idaho')
n_cores  <- 2
  # Test run: a short chain to check the code (FALSE for the real run)
test_run <- FALSE

der_dir <- here("data", "data.derived", site_tag)


# ═══════════════════════════════════════════════════════════════════════════
# 模型
# ═══════════════════════════════════════════════════════════════════════════
# 生命率编号：1 存活、2 补充；擁擠度编号 x：1 forb、2 禾草、3 灌木（与 W 的列同序）
# ★ 只有一个（合并的）族群，所以系数没有 [k]，也没有 trait
model_forb <- nimbleCode({

  # ═══════════════════ 固定效应 ═══════════════════
  for (v in 1:2){
    intercept[v]         ~ dnorm(0, sd = 5)      # α0
    slope.graze[v] ~ dnorm(0, sd = 2)      # γ
    slope.clim[v]  ~ dnorm(0, sd = 2)      # β_climate
    for (x in 1:nX){
      slope.W[v,x]     ~ dnorm(0, sd = 2)  # β_x：擁擠度（补充里为样方的 √值）
      com.W.x.clim[v,x] ~ dnorm(0, sd = 2) # δ_x：降水 × 擁擠度
    }
  }

  # ═══════════════════ 年份层与样方层（随机） ═══════════════════
  # ★ 年份效应中心化在截距与降水上（同主模型）：yr.int = 该年的总水平，
  #   intercept 与 β_climate 在它的平均里，不再单独写进线性预测子
  for (v in 1:2){
    yr.int.sd[v] ~ dexp(2)
    grp.sd[v]    ~ dexp(1)
    quad.sd[v]   ~ dexp(1)
    for (t in 1:nYr){   yr.int[t,v]   ~ dnorm(intercept[v] + slope.clim[v] * clim[t], sd = yr.int.sd[v]) }
    for (g in 1:nGrp){  grp.int[g,v]  ~ dnorm(0, sd = grp.sd[v]) }
    for (q in 1:nQuad){ quad.int[q,v] ~ dnorm(0, sd = quad.sd[v]) }
  }

  # ═══════════════════ 观测层 ① 存活（无大小） ═══════════════════
  for (r in 1:nSurv){
    lpS[r] <- slope.graze[1] * graze[r] +
      inprod(slope.W[1,1:nX] + com.W.x.clim[1,1:nX] * clim[idYr[r]], W[r,1:nX]) +
      yr.int[idYr[r],1] + grp.int[idGrp[r],1] + quad.int[idQuad[r],1]
    surv[r] ~ dbern(ilogit(lpS[r]))
  }

  # ═══════════════════ 观测层 ② 补充（样方 × 转换） ═══════════════════
  # 有效亲本数 N′ = p·样方 forb 株数 + (1 − p)·组平均株数（同主模型的 C′）
  # 密度依赖：forb 用 √N′（模型内算，按 scaling$C 的 C_forb 标准化），禾草 / 灌木用样方盖度的 √值
  p.eff ~ dbeta(1, 1)
  nb.size ~ dexp(0.1)
  for (j in 1:nRec){
    Neff[j] <- p.eff * nQ[j] + (1 - p.eff) * nG[j]
    CRf[j]  <- (sqrt(Neff[j]) - cForbMean) / cForbSd
    betaR[j,1:nX] <- slope.W[2,1:nX] + com.W.x.clim[2,1:nX] * clim[idYrR[j]]
    lpR[j] <- slope.graze[2] * grazeR[j] +
      betaR[j,1] * CRf[j] + inprod(betaR[j,2:nX], CRo[j,1:2]) +
      yr.int[idYrR[j],2] + grp.int[idGrpR[j],2] + quad.int[idQuadR[j],2]
    lambdaR[j] <- Neff[j] * exp(lpR[j])
    recruits[j] ~ dnegbin(prob = nb.size / (nb.size + lambdaR[j]), size = nb.size)
  }
  # ★ 组内无 forb 时 N′ = 0，λ = 0；输入资料里没有「无亲本却有新株」的格子（建模前检查）
})


# ═══════════════════════════════════════════════════════════════════════════
# data / constants
# ═══════════════════════════════════════════════════════════════════════════
partA_input <- readRDS(file.path(der_dir, "CLR_PartA_input.rds"))
rec  <- partA_input$forb$rec
recr <- partA_input$forb$recruit
if (test_run) {
  rec <- rec[seq(1, nrow(rec), by = 10), ]   # one record in ten
  n_chains <- 2; n_iter <- 200; n_burn <- 100; n_thin <- 1
}

constList_forb <- list(
  nYr = nrow(partA_input$years), nQuad = nrow(partA_input$quadrats), nGrp = nrow(partA_input$groups), nX = 3,
  nSurv = nrow(rec), nRec = nrow(recr),
  idYr = rec$t, idQuad = rec$q, idGrp = rec$g, graze = rec$E,
  W = as.matrix(rec[c("W_forb_z", "W_grass_z", "W_shrub_z")]),
  clim = partA_input$years$clim_anom,
  idYrR = recr$t, idQuadR = recr$q, idGrpR = recr$g, grazeR = recr$E,
  nQ = recr$n_q, nG = recr$n_g, CRo = as.matrix(recr[c("C_grass_z", "C_shrub_z")]),
  cForbMean = partA_input$scaling$C["mean", "C_forb"], cForbSd = partA_input$scaling$C["sd", "C_forb"]
)
dataList_forb <- list(surv = rec$survives, recruits = recr$recruits)

# ── 建模前的检查 ──────────────────────────────────────────────────────────
with(constList_forb, stopifnot(
  all(idYr %in% 1:nYr), all(idQuad %in% 1:nQuad), all(idGrp %in% 1:nGrp), all(is.finite(W)),
  all(graze %in% 0:1), all(is.finite(CRo)), all(is.finite(clim)),
  !any(nG == 0 & dataList_forb$recruits > 0)))
stopifnot(all(dataList_forb$surv %in% 0:1))
cat(sprintf("Forb submodel: %d survival records, %d recruitment rows, %d years.\n",
            constList_forb$nSurv, constList_forb$nRec, constList_forb$nYr))


# ═══════════════════════════════════════════════════════════════════════════
# inits
# ═══════════════════════════════════════════════════════════════════════════
inits_forb <- function() {
  nYr <- constList_forb$nYr
  l <- list(intercept = c(rnorm(1, 0.5, 0.3), rnorm(1, -1, 0.3)), slope.graze = rnorm(2, 0, 0.1),
       slope.clim = rnorm(2, 0, 0.1), slope.W = matrix(rnorm(6, 0, 0.1), 2, 3),
       com.W.x.clim = matrix(rnorm(6, 0, 0.1), 2, 3),
       yr.int.sd = rexp(2, 4), grp.sd = rexp(2, 4), quad.sd = rexp(2, 4),
       yr.int = matrix(rnorm(nYr * 2, 0, 0.1), nYr, 2),
       grp.int = matrix(rnorm(constList_forb$nGrp * 2, 0, 0.1), constList_forb$nGrp, 2),
       quad.int = matrix(rnorm(constList_forb$nQuad * 2, 0, 0.1), constList_forb$nQuad, 2),
       p.eff = runif(1, 0.3, 0.7), nb.size = rgamma(1, 4, 2))
    # 年份效应从它们的平均（截距）附近出发
  l$yr.int <- sweep(l$yr.int, 2, l$intercept, "+")
  l
}
params_forb <- c("intercept", "slope.graze", "slope.clim", "slope.W", "com.W.x.clim",
                 "yr.int.sd", "grp.sd", "quad.sd", "yr.int", "grp.int", "quad.int", "p.eff", "nb.size")


# ═══════════════════════════════════════════════════════════════════════════
# 建模、编译、MCMC
# ═══════════════════════════════════════════════════════════════════════════
t0 <- Sys.time()
m_forb <- nimbleModel(model_forb, constants = constList_forb, data = dataList_forb, inits = inits_forb(),
                      calculate = FALSE)
ll0 <- m_forb$calculate()
cat(sprintf("Forb submodel built; log probability at the inits: %.1f\n", ll0))
stopifnot(is.finite(ll0))
rm(m_forb); invisible(gc())
samples_forb <- run_chains(model_forb, constList_forb, dataList_forb, inits_forb, monitors = params_forb,
                           thin = n_thin, n_iter = n_iter, n_burn = n_burn, n_chains = n_chains,
                           n_cores = n_cores, seed = 2026,
                           log_file = file.path(der_dir, "CLR_PartA_forb_chains.log"))
cat(sprintf("Forb submodel done (%.1f min).\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Save -------------------------------------------------------------------------
out_file <- file.path(der_dir, sprintf("CLR_PartA_forb%s.rds", if (test_run) "_test" else ""))
saveRDS(list(samples = samples_forb, scaling = partA_input$scaling,
             settings = list(n_iter = n_iter, n_burn = n_burn, n_thin = n_thin, n_chains = n_chains,
                             test_run = test_run)),
        out_file)
cat("Saved:\n  data/data.derived/", site_tag, ": ", basename(out_file), "\n", sep = "")
rm(list = intersect(c("rec", "recr", "t0", "ll0", "out_file"), ls()))
