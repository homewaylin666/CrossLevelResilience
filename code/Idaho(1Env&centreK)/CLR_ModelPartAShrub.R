################################################################################
## Part A other-shrub submodel (Idaho): NIMBLE code, data, inits and MCMC
##
## Written by Hung-wei Lin Oct 2026
##
## Shrubs below the sample-size threshold of the main model (GUSA, ERMI, TECA, PUTR) pooled as one
## population, fitted separately from the main model (the two share no parameters). Survival and
## growth with size, recruitment per quadrat x transition and sizes of recruits, with the predictors
## of the main model except trait: size, crowding (same species, grass, shrub, forb), grazing E,
## water-year precipitation and its modulation of size and crowding (delta), year effects, quadrat
## (in group) effects and plant effects in survival and growth (correlated, non-centred; random effects
## follow the sampling design, as the main model). Size² stays out: no gain in AIC for these few
## records (ΔAIC: survival +2.0, growth +1.2; glmer, 22 Oct 2026). Gives the dynamics of these shrubs for
## Part B (optional there).
## The data are few (about 650 survival and 220 growth records), so the size x crowding (theta) and
## grazing x crowding (omega) terms of the main model stay out.
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
# 生命率编号：1 存活、2 生长、3 补充；擁擠度编号 x：1 种内、2 禾草、3 灌木、4 forb
# ★ 一个（合并的）族群，系数没有 [k]，也没有 trait
model_shrub <- nimbleCode({

  # ═══════════════════ 固定效应 ═══════════════════
  for (v in 1:3){
    intercept[v]         ~ dnorm(0, sd = 5)       # α0
    slope.graze[v] ~ dnorm(0, sd = 2)       # γ
    slope.clim[v]  ~ dnorm(0, sd = 2)       # β_climate
    for (x in 1:nX){
      slope.W[v,x]      ~ dnorm(0, sd = 2)  # β_x
      com.W.x.clim[v,x] ~ dnorm(0, sd = 2)  # δ_x
    }
  }
  for (v in 1:2){
    slope.size[v]      ~ dnorm(0, sd = 2)   # α1（size² 的 AIC 没有改善：存活 +2.0，生长 +1.2）
    com.size.x.clim[v] ~ dnorm(0, sd = 2)   # δ_size
  }

  # ═══════════════════ 年份层与样方层（随机） ═══════════════════
  # ★ 年份效应中心化在截距与降水上（同主模型）：yr.int = 该年的总水平，
  #   intercept 与 β_climate 在它的平均里，不再单独写进线性预测子
  for (v in 1:3){
    yr.int.sd[v] ~ dexp(2)
    grp.sd[v]    ~ dexp(1)
    quad.sd[v]   ~ dexp(1)
    for (t in 1:nYr){   yr.int[t,v]   ~ dnorm(intercept[v] + slope.clim[v] * clim[t], sd = yr.int.sd[v]) }
    for (g in 1:nGrp){  grp.int[g,v]  ~ dnorm(0, sd = grp.sd[v]) }
    for (q in 1:nQuad){ quad.int[q,v] ~ dnorm(0, sd = quad.sd[v]) }
  }

  # ═══════════════════ 个体层（植株随机效应，非中心化） ═══════════════════
  # ★ 存活与生长的个体效应都放，两者相关 ρ（同主模型：随机效应按抽样设计放）
  for (v in 1:2){ ind.sd[v] ~ dexp(1) }
  ind.rho01 ~ dbeta(2, 2)                        # ρ = 2·ρ01 − 1，即 2×2 的 LKJ(η = 2)
  ind.rho <- 2 * ind.rho01 - 1
  for (i in 1:nInd){
    ind.raw[i,1] ~ dnorm(0, sd = 1)
    ind.raw[i,2] ~ dnorm(0, sd = 1)
    ind.eff[i,1] <- ind.sd[1] * ind.raw[i,1]
    ind.eff[i,2] <- ind.sd[2] * (ind.rho * ind.raw[i,1] + sqrt(1 - ind.rho^2) * ind.raw[i,2])
  }

  # ═══════════════════ 观测层 ① 存活 ═══════════════════
  for (r in 1:nSurv){
    lpS[r] <- slope.size[1] * size[r] +
      slope.graze[1] * graze[r] +
      com.size.x.clim[1] * size[r] * clim[idYr[r]] +
      inprod(slope.W[1,1:nX] + com.W.x.clim[1,1:nX] * clim[idYr[r]], W[r,1:nX]) +
      yr.int[idYr[r],1] + grp.int[idGrp[r],1] + quad.int[idQuad[r],1] + ind.eff[idInd[r],1]
    surv[r] ~ dbern(ilogit(lpS[r]))
  }

  # ═══════════════════ 观测层 ② 生长 ═══════════════════
  for (s in 1:nGrow){
    lpG[s] <- slope.size[2] * size[idGrow[s]] +
      slope.graze[2] * graze[idGrow[s]] +
      com.size.x.clim[2] * size[idGrow[s]] * clim[idYr[idGrow[s]]] +
      inprod(slope.W[2,1:nX] + com.W.x.clim[2,1:nX] * clim[idYr[idGrow[s]]], W[idGrow[s],1:nX]) +
      yr.int[idYr[idGrow[s]],2] + grp.int[idGrp[idGrow[s]],2] + quad.int[idQuad[idGrow[s]],2] +
      ind.eff[idInd[idGrow[s]],2]
    sizeNext[s] ~ dnorm(lpG[s], sd = grow.sd)
  }
  grow.sd ~ dexp(1)

  # ═══════════════════ 观测层 ③ 补充（样方 × 转换） ═══════════════════
  # 有效盖度 C′ = p·样方盖度 + (1 − p)·组平均盖度；种内用 √C′，其余用样方的 √值（同主模型）
  p.eff ~ dbeta(1, 1)
  nb.size ~ dexp(0.1)
  for (j in 1:nRec){
    Ceff[j] <- p.eff * covQ[j] + (1 - p.eff) * covG[j]
    CRi[j]  <- (sqrt(Ceff[j]) - cIntraMean) / cIntraSd
    betaR[j,1:nX] <- slope.W[3,1:nX] + com.W.x.clim[3,1:nX] * clim[idYrR[j]]
    lpR[j] <- slope.graze[3] * grazeR[j] +
      betaR[j,1] * CRi[j] + inprod(betaR[j,2:nX], CRo[j,1:3]) +
      yr.int[idYrR[j],3] + grp.int[idGrpR[j],3] + quad.int[idQuadR[j],3]
    lambdaR[j] <- Ceff[j] * exp(lpR[j])
    recruits[j] ~ dnegbin(prob = nb.size / (nb.size + lambdaR[j]), size = nb.size)
  }
  # ★ 组内无亲本时 C′ = 0，λ = 0；输入资料里没有「无亲本却有新株」的格子（建模前检查）

  # ═══════════════════ 新个体大小分布 c(z') ═══════════════════
  # ★ 只用有面积的新株；记为点的新株在 IPM 里按观测比例进入「点」阶段
  sizeRec.mu ~ dnorm(0, sd = 5)
  sizeRec.sd ~ dexp(1)
  for (j in 1:nRecSize){ sizeRec[j] ~ dnorm(sizeRec.mu, sd = sizeRec.sd) }
})


# ═══════════════════════════════════════════════════════════════════════════
# data / constants
# ═══════════════════════════════════════════════════════════════════════════
partA_input <- readRDS(file.path(der_dir, "CLR_PartA_input.rds"))
rec   <- partA_input$shrub$rec
recr  <- partA_input$shrub$recruit
recsz <- filter(partA_input$shrub$recsize, !point)
if (test_run) { n_chains <- 2; n_iter <- 200; n_burn <- 100; n_thin <- 1 }
W_z  <- c("W_intra_z", "W_grass_z", "W_shrub_z", "W_forb_z")
grow <- which(rec$grow_ok)

constList_shrub <- list(
  nYr = nrow(partA_input$years), nQuad = nrow(partA_input$quadrats), nGrp = nrow(partA_input$groups),
  nInd = max(rec$i), nX = 4, nSurv = nrow(rec), nGrow = length(grow), nRec = nrow(recr),
  nRecSize = nrow(recsz),
  idYr = rec$t, idQuad = rec$q, idGrp = rec$g, idInd = rec$i, idGrow = grow,
  size = rec$z, graze = rec$E, W = as.matrix(rec[W_z]),
  clim = partA_input$years$clim_anom,
  idYrR = recr$t, idQuadR = recr$q, idGrpR = recr$g, grazeR = recr$E,
  covQ = recr$cover_q, covG = recr$cover_g, CRo = as.matrix(recr[c("C_grass_z", "C_shrub_z", "C_forb_z")]),
    # 种内盖度的标准化常数：合并族群自己的一组（scaling$C_intra 的 "other shrubs" 行）
  cIntraMean = with(partA_input$scaling$C_intra, mean[code == "other shrubs"]),
  cIntraSd   = with(partA_input$scaling$C_intra, sd[code == "other shrubs"])
)
dataList_shrub <- list(surv = rec$survives, sizeNext = rec$z_next[grow], recruits = recr$recruits,
                       sizeRec = recsz$z_rec)

# ── 建模前的检查 ──────────────────────────────────────────────────────────
with(constList_shrub, stopifnot(
  all(idYr %in% 1:nYr), all(idQuad %in% 1:nQuad), all(idGrp %in% 1:nGrp),
  all(sort(unique(idInd)) == 1:nInd), all(is.finite(size)), all(is.finite(W)), all(graze %in% 0:1),
  all(is.finite(CRo)), all(is.finite(clim)), !any(covG == 0 & dataList_shrub$recruits > 0)))
stopifnot(all(dataList_shrub$surv %in% 0:1), all(is.finite(dataList_shrub$sizeNext)),
          all(rec$survives[grow] == 1))
cat(sprintf("Other-shrub submodel: %d survival, %d growth, %d recruitment rows, %d recruit sizes (%d of %d recruits as points), %d plants.\n",
            constList_shrub$nSurv, constList_shrub$nGrow, constList_shrub$nRec, constList_shrub$nRecSize,
            sum(partA_input$shrub$recsize$point), nrow(partA_input$shrub$recsize), constList_shrub$nInd))


# ═══════════════════════════════════════════════════════════════════════════
# inits
# ═══════════════════════════════════════════════════════════════════════════
inits_shrub <- function() {
  nYr <- constList_shrub$nYr; nGrp <- constList_shrub$nGrp; nQuad <- constList_shrub$nQuad
  l <- list(intercept = c(rnorm(2, 0, 0.3), rnorm(1, -1, 0.3)), slope.graze = rnorm(3, 0, 0.1),
       slope.clim = rnorm(3, 0, 0.1), slope.W = matrix(rnorm(12, 0, 0.1), 3, 4),
       com.W.x.clim = matrix(rnorm(12, 0, 0.1), 3, 4),
       slope.size = rnorm(2, 0.5, 0.2), com.size.x.clim = rnorm(2, 0, 0.1),
       yr.int.sd = rexp(3, 4), grp.sd = rexp(3, 4), quad.sd = rexp(3, 4),
       yr.int = matrix(rnorm(nYr * 3, 0, 0.1), nYr, 3),
       grp.int = matrix(rnorm(nGrp * 3, 0, 0.1), nGrp, 3),
       quad.int = matrix(rnorm(nQuad * 3, 0, 0.1), nQuad, 3),
       ind.sd = runif(2, 0.1, 1.5), ind.rho01 = runif(1, 0.3, 0.7),
       ind.raw = matrix(rnorm(constList_shrub$nInd * 2), constList_shrub$nInd, 2),
       grow.sd = rexp(1, 2) + 0.3, p.eff = runif(1, 0.3, 0.7), nb.size = rgamma(1, 4, 2),
       sizeRec.mu = rnorm(1, -1, 0.3), sizeRec.sd = rexp(1, 2) + 0.3)
    # 年份效应从它们的平均（截距）附近出发
  l$yr.int <- sweep(l$yr.int, 2, l$intercept, "+")
  l
}
params_shrub <- c("intercept", "slope.size", "slope.graze", "slope.clim", "slope.W",
                  "com.size.x.clim", "com.W.x.clim", "yr.int.sd", "grp.sd", "quad.sd",
                  "yr.int", "grp.int", "quad.int", "ind.sd", "ind.rho", "grow.sd", "p.eff", "nb.size",
                  "sizeRec.mu", "sizeRec.sd")


# ═══════════════════════════════════════════════════════════════════════════
# 建模、编译、MCMC
# ═══════════════════════════════════════════════════════════════════════════
t0 <- Sys.time()
m_shrub <- nimbleModel(model_shrub, constants = constList_shrub, data = dataList_shrub, inits = inits_shrub(),
                       calculate = FALSE)
ll0 <- m_shrub$calculate()
cat(sprintf("Other-shrub submodel built; log probability at the inits: %.1f\n", ll0))
stopifnot(is.finite(ll0))
rm(m_shrub); invisible(gc())
samples_shrub <- run_chains(model_shrub, constList_shrub, dataList_shrub, inits_shrub, monitors = params_shrub,
                            thin = n_thin, n_iter = n_iter, n_burn = n_burn, n_chains = n_chains,
                            n_cores = n_cores, seed = 2026,
                            log_file = file.path(der_dir, "CLR_PartA_shrub_chains.log"))
cat(sprintf("Other-shrub submodel done (%.1f min).\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Save -------------------------------------------------------------------------
out_file <- file.path(der_dir, sprintf("CLR_PartA_shrub%s.rds", if (test_run) "_test" else ""))
saveRDS(list(samples = samples_shrub, scaling = partA_input$scaling,
             point_share_recruits = mean(partA_input$shrub$recsize$point),
             settings = list(n_iter = n_iter, n_burn = n_burn, n_thin = n_thin, n_chains = n_chains,
                             test_run = test_run)),
        out_file)
cat("Saved:\n  data/data.derived/", site_tag, ": ", basename(out_file), "\n", sep = "")
rm(list = intersect(c("rec", "recr", "recsz", "grow", "W_z", "t0", "ll0", "out_file"), ls()))
