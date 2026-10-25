################################################################################
## Part A multilevel vital-rate model (Idaho): NIMBLE code, data, inits and MCMC
##
## Written by Hung-wei Lin Oct 2026
##
## Joint model of survival, growth (focal plants) and recruitment (species x quadrat x
## transition) of the main species, with four blocks (Fig. C3_model_diagram_Adler):
## individual (size, grazing), population (life-history trait), community (water-year
## precipitation, year effects) and interactions (crowding W by same species, grass, shrub, forb).
## Stage 1 ('stage <- 1'): no trait, the baseline for the IPMs that give the trait.
## Stage 2 ('stage <- 2'): trait explains gamma, beta_x and beta_climate (one run per trait draw).
## Plant (individual) random effects in the non-centred form (setting 'plant_effects').
## Reads 'CLR_PartA_input.rds' (CLR_InputPartA.R); uses 'site_tag' and 'run_chains' (CLR_fns_nimble.R).


# Settings ---------------------------------------------------------------------
stage <- 1
  # Stage 2: matrix of trait draws (rows = draws, columns = species in the order of the species
  # table), standardised; one MCMC run per draw ('trait_draw')
trait_file <- NULL
trait_draw <- 1
  # Plant (individual) random effects: "both" = survival and growth, correlated (main run);
  # "survival" = survival only, "none" = no plant effects (sensitivity runs)
plant_effects <- "both"
  # MCMC
n_chains <- 4
n_iter   <- 20000
n_burn   <- 10000
n_thin   <- 10
  # Cores for the chains: 1 = one chain after another (progress bar in the console); 2-4 = chains
  # in parallel (each core builds and compiles its own copy of the model, so memory grows with
  # it; the progress bars go to a log file in 'data/data.derived/Idaho')
n_cores  <- 2
  # Test run: a few species and a short chain to check the code (FALSE for the real run)
test_run <- FALSE
test_species <- c("PSSP", "HECO", "ARTR")

der_dir <- here("data", "data.derived", site_tag)
stopifnot(stage %in% 1:2, plant_effects %in% c("survival", "both", "none"))


# ═══════════════════════════════════════════════════════════════════════════
# 模型
# ═══════════════════════════════════════════════════════════════════════════
# 生命率编号：1 存活、2 生长、3 补充（v.surv / v.grow / v.recr）
# 擁擠度编号 x：1 种内、2 禾草、3 灌木、4 forb（与 W 的列同序）
# ★ 擁擠度的四个效应合成一个「有效交互强度」向量（图里的 β_x,eff）：
#   βeff = β_x[k,v] + θ_x[v]·z + ω_x[v]·E + δ_x[v]·climate + ε_t,x[v]
#   线性预测子里只需一次 inprod(βeff, W)，与示意图逐项对应。
#   （δ_x·climate 写在 ε_t,x 的平均里：yr.slope.W[t,v,x] ~ N(δ_x·clim[t], σ)，见年份层）
model <- nimbleCode({

  # ═══════════════════ 群落层：气候（固定）+ 年份偏离（随机） ═══════════════════
  # ★ clim 是降水距平：与测站长期平均的差，以长期 SD 为单位（保留正负；0 = 长期平均）
  #   研究期比长期平均偏干，所以 clim 在转换年份上的平均不是 0；截距对应长期平均的降水
  for (v in 1:2){
    com.size.x.clim[v] ~ dnorm(0, sd = 2)          # δ_size：降水 × 自身大小
    yr.slope.size.sd[v] ~ dexp(2)                  # ε_t,size
  }
  for (v in 1:nVR){
    for (x in 1:nX){
      com.W.x.clim[v,x] ~ dnorm(0, sd = 2)         # δ_x：降水 × 擁擠度（Q2 的核心）
      yr.slope.W.sd[v,x] ~ dexp(2)                 # ε_t,x：该年擁擠斜率偏离多少
    }
    yr.int.sd[v]   ~ dexp(2)                       # ε_t：该年整体好坏
    yrsp.int.sd[v] ~ dexp(2)                       # ε_t,k：该年该物种的残余
  }
  # ★ yr.*.sd 的先验决定 δ 的后验宽度：紧的 sd 先验把变异推给 com.*.x.clim
  #   必做敏感度：dexp(2) vs dexp(1)

  # ★ 年份效应「中心化」在降水上（hierarchical centring; Gelfand et al. 1995）：
  #   年份效应的平均 = 该年降水的预测，所以 yr.int、yrsp.int、yr.slope.* 是「该年的总效应」，
  #   降水的部分（β_climate、δ）在它们的平均里，不再单独写进线性预测子。模型与原写法相同，
  #   只是抽样快得多：降水只有 22 个年份值，原写法下 β_climate、δ 与年份效应互相交换，链走得慢。
  #   yrsp.int[t,k,v] 也包含物种截距 pop.int[k,v]（线性预测子里只出现 yrsp.int）。
  # ★ Part B：某年降水改为 clim* 时，年份效应 = 平均（用 clim*）+ 该年的残差（后验值 − 原平均）
  for (t in 1:nYr){
    for (v in 1:nVR){
      yr.int[t,v] ~ dnorm(pop.int.mu[v] + pop.slope.clim.mu[v] * clim[t], sd = yr.int.sd[v])
      for (x in 1:nX){
        yr.slope.W[t,v,x] ~ dnorm(com.W.x.clim[v,x] * clim[t], sd = yr.slope.W.sd[v,x])
      }
      for (k in 1:nSp){
        yrsp.int[t,k,v] ~ dnorm(yr.int[t,v] + pop.int[k,v] - pop.int.mu[v] +
                                  (pop.slope.clim[k,v] - pop.slope.clim.mu[v]) * clim[t], sd = yrsp.int.sd[v])
      }
    }
    for (v in 1:2){
      yr.slope.size[t,v] ~ dnorm(com.size.x.clim[v] * clim[t], sd = yr.slope.size.sd[v])
    }
  }

  # ═══════════════════ 物种层 · 先验 ═══════════════════
  # Y = α0 / α1 / α2：只有物种随机效应，不放 trait（trait 由它们导出）
  # ★ size² (α2) 按 AIC 决定（补充方法 S3）：存活 ΔAIC −6.6、生长 −123，两个都放
  #   （第一次正式跑的 α2：result/tables/Idaho/CLR_PartA_stage1_size2.csv）
  # X = γ / β_x / β_climate：stage 2 才放 trait
  for (v in 1:nVR){
    pop.int.mu[v] ~ dnorm(0, sd = 5)
    pop.int.sd[v] ~ dexp(1)

    pop.slope.graze.mu[v] ~ dnorm(0, sd = 2)       # γ 的三件套
    pop.slope.graze.sd[v] ~ dexp(1)
    pop.slope.clim.mu[v]  ~ dnorm(0, sd = 2)       # β_climate 的三件套
    pop.slope.clim.sd[v]  ~ dexp(1)
    for (x in 1:nX){
      pop.slope.W.mu[v,x] ~ dnorm(0, sd = 2)       # β_x 的三件套
      pop.slope.W.sd[v,x] ~ dexp(1)
    }
    if (useTrait == 1){
      pop.slope.graze.trait[v] ~ dnorm(0, sd = 2)  # ★ 什么性状的物种更受放牧影响
      pop.slope.clim.trait[v]  ~ dnorm(0, sd = 2)  # ★ 什么性状的物种对降水更敏感
      for (x in 1:nX){
        pop.slope.W.trait[v,x] ~ dnorm(0, sd = 2)  # ★ 什么性状的物种更怕擁擠
      }
    }
  }
  for (v in 1:2){
    pop.slope.size.mu[v]  ~ dnorm(0, sd = 2)
    pop.slope.size.sd[v]  ~ dexp(1)
  }
  for (v in 1:2){
    pop.slope.size2.mu[v] ~ dnorm(0, sd = 2)       # α2
    pop.slope.size2.sd[v] ~ dexp(1)
  }

  # ═══════════ 物种层 · 每个物种的系数 ═══════════
  for (k in 1:nSp){
    for (v in 1:nVR){
      pop.int[k,v] ~ dnorm(pop.int.mu[v], sd = pop.int.sd[v])
      if (useTrait == 1){
        pop.slope.graze[k,v] ~ dnorm(pop.slope.graze.mu[v] + pop.slope.graze.trait[v] * trait[k],
                                     sd = pop.slope.graze.sd[v])
        pop.slope.clim[k,v]  ~ dnorm(pop.slope.clim.mu[v] + pop.slope.clim.trait[v] * trait[k],
                                     sd = pop.slope.clim.sd[v])
        for (x in 1:nX){
          pop.slope.W[k,v,x] ~ dnorm(pop.slope.W.mu[v,x] + pop.slope.W.trait[v,x] * trait[k],
                                     sd = pop.slope.W.sd[v,x])
        }
      } else {
        pop.slope.graze[k,v] ~ dnorm(pop.slope.graze.mu[v], sd = pop.slope.graze.sd[v])
        pop.slope.clim[k,v]  ~ dnorm(pop.slope.clim.mu[v],  sd = pop.slope.clim.sd[v])
        for (x in 1:nX){
          pop.slope.W[k,v,x] ~ dnorm(pop.slope.W.mu[v,x], sd = pop.slope.W.sd[v,x])
        }
      }
    }
    for (v in 1:2){
      pop.slope.size[k,v]  ~ dnorm(pop.slope.size.mu[v],  sd = pop.slope.size.sd[v])
    }
    for (v in 1:2){
      pop.slope.size2[k,v] ~ dnorm(pop.slope.size2.mu[v], sd = pop.slope.size2.sd[v])
    }
  }

  # ═══════════════════ 个体层的交互（全物种共享） ═══════════════════
  for (x in 1:nX){
    for (v in 1:2){
      com.W.x.size[v,x] ~ dnorm(0, sd = 2)         # θ_x：擁擠度 × 自身大小
    }
    for (v in 1:nVR){
      com.W.x.graze[v,x] ~ dnorm(0, sd = 2)        # ω_x：擁擠度 × 放牧
    }
  }

  # ═══════════════════ 样方层：组 + 组内样方（随机） ═══════════════════
  for (v in 1:nVR){
    grp.sd[v]  ~ dexp(1)
    quad.sd[v] ~ dexp(1)
    for (g in 1:nGrp){  grp.int[g,v]  ~ dnorm(0, sd = grp.sd[v]) }
    for (q in 1:nQuad){ quad.int[q,v] ~ dnorm(0, sd = quad.sd[v]) }
  }

  # ═══════════════════ 个体层（植株随机效应，非中心化） ═══════════════════
  # ★ 非中心化：ind.raw ~ N(0, 1)，ind.eff = ind.sd × ind.raw。六成植株只有一笔存活纪录；
  #   中心化写法（ind.eff ~ N(0, ind.sd)）下 ind.sd 与上万个效应互相卡住，第一次正式跑
  #   ind.sd[1] 四条链停在 0.08、0.10、0.26、1.43（Rhat 4.1）
  # ★ 有存活的个体效应时，存活系数是「条件」系数（给定植株）：IPM 的存活要对 ind.eff 的
  #   常态分布取平均（logit 尺度积分），代 0 会低估族群平均存活
  # ★ 随机效应按抽样设计放（同一植株的重复纪录）：存活与生长的个体效应都放，两者相关 ρ
  #   （Bolker et al. 2009；glmer 也支持：存活 SD 1.14，AIC 低 361；生长 SD 0.06，AIC 低 4）
  #   敏感度：只放存活（"survival"）、都不放（"none"）
  if (useIndS == 1){
    ind.sd[1] ~ dexp(1)
    for (i in 1:nInd){
      ind.raw[i,1] ~ dnorm(0, sd = 1)
      ind.eff[i,1] <- ind.sd[1] * ind.raw[i,1]
    }
  } else {
    for (i in 1:nInd){ ind.eff[i,1] <- 0 }
  }
  if (useIndG == 1){
    ind.sd[2] ~ dexp(1)
    ind.rho01 ~ dbeta(2, 2)                        # ρ = 2·ρ01 − 1，即 2×2 的 LKJ(η = 2)
    ind.rho <- 2 * ind.rho01 - 1
    for (i in 1:nInd){
      ind.raw[i,2] ~ dnorm(0, sd = 1)
      ind.eff[i,2] <- ind.sd[2] * (ind.rho * ind.raw[i,1] + sqrt(1 - ind.rho^2) * ind.raw[i,2])
    }
  } else {
    for (i in 1:nInd){ ind.eff[i,2] <- 0 }
  }

  # ═══════════════════ 观测层 ① 存活 ═══════════════════
  for (r in 1:nSurv){
    # （截距、β_climate、δ 都在年份效应里：yrsp.int = 物种截距 + 降水效应 + 年份残差）
    betaS[r,1:nX] <- pop.slope.W[idSp[r],1,1:nX] + com.W.x.size[1,1:nX] * size[r] +
      com.W.x.graze[1,1:nX] * graze[r] + yr.slope.W[idYr[r],1,1:nX]  # βeff，存活
    lpS[r] <- yrsp.int[idYr[r],idSp[r],1] +
      pop.slope.size[idSp[r],1]  * size[r] +
      pop.slope.size2[idSp[r],1] * size2[r] +
      pop.slope.graze[idSp[r],1] * graze[r] +
      inprod(betaS[r,1:nX], W[r,1:nX]) +
      yr.slope.size[idYr[r],1] * size[r] +
      grp.int[idGrp[r],1] + quad.int[idQuad[r],1] +
      ind.eff[idInd[r],1]
    surv[r] ~ dbern(ilogit(lpS[r]))
  }

  # ═══════════════════ 观测层 ② 生长（存活且两年都有面积） ═══════════════════
  # ★ size 与 sizeNext 用同一组 mean / sd（输入文件 scaling$z），IPM 的核才能迭代
  for (s in 1:nGrow){
    betaG[s,1:nX] <- pop.slope.W[idSp[idGrow[s]],2,1:nX] + com.W.x.size[2,1:nX] * size[idGrow[s]] +
      com.W.x.graze[2,1:nX] * graze[idGrow[s]] + yr.slope.W[idYr[idGrow[s]],2,1:nX]  # βeff，生长
    lpG[s] <- yrsp.int[idYr[idGrow[s]],idSp[idGrow[s]],2] +
      pop.slope.size[idSp[idGrow[s]],2]  * size[idGrow[s]] +
      pop.slope.size2[idSp[idGrow[s]],2] * size2[idGrow[s]] +
      pop.slope.graze[idSp[idGrow[s]],2] * graze[idGrow[s]] +
      inprod(betaG[s,1:nX], W[idGrow[s],1:nX]) +
      yr.slope.size[idYr[idGrow[s]],2] * size[idGrow[s]] +
      grp.int[idGrp[idGrow[s]],2] + quad.int[idQuad[idGrow[s]],2] +
      ind.eff[idInd[idGrow[s]],2]
    sizeNext[s] ~ dnorm(lpG[s], sd = grow.sd[idSp[idGrow[s]]])
  }
  for (k in 1:nSp){ grow.sd[k] ~ dexp(1) }
  # ★ Adler 系的生长残差随预测大小变化（Var = a·exp(b·û)）；若残差图显示漏斗形再改

  # ═══════════════════ 观测层 ③ 补充（物种 × 样方 × 转换） ═══════════════════
  # 有效盖度 C′ = p·样方盖度 + (1 − p)·组平均盖度（Adler et al. 2010; Chu & Adler 2015）
  # 密度依赖用盖度平方根（Chu & Adler 2015），种内项用 √C′（模型内算，再按输入文件的
  # scaling$C_intra 在物种内标准化，同种内擁擠度），禾草 / 灌木 / forb 用样方的值（输入文件已标准化，CRo 三列）
  # ★ 背景补充 exp(rec.bg[k])：组内无亲本时仍有新株（15 格、48 株），没有它这些格子的似然为 0
  for (k in 1:nSp){
    p.eff[k] ~ dbeta(1, 1)
    rec.bg[k] ~ dnorm(rec.bg.mu, sd = rec.bg.sd)
    nb.size[k] ~ dexp(0.1)
  }
  rec.bg.mu ~ dnorm(-2, sd = 2)
  rec.bg.sd ~ dexp(1)
  for (j in 1:nRec){
    Ceff[j] <- p.eff[idSpR[j]] * covQ[j] + (1 - p.eff[idSpR[j]]) * covG[j]
    CRi[j]  <- (sqrt(Ceff[j]) - cIntraMean[idSpR[j]]) / cIntraSd[idSpR[j]]
    betaR[j,1:nX] <- pop.slope.W[idSpR[j],3,1:nX] + com.W.x.graze[3,1:nX] * grazeR[j] +
      yr.slope.W[idYrR[j],3,1:nX]                                   # βeff，补充
    lpR[j] <- yrsp.int[idYrR[j],idSpR[j],3] +
      pop.slope.graze[idSpR[j],3] * grazeR[j] +
      betaR[j,1] * CRi[j] + inprod(betaR[j,2:nX], CRo[j,1:3]) +
      grp.int[idGrpR[j],3] + quad.int[idQuadR[j],3]
    lambdaR[j] <- Ceff[j] * exp(lpR[j]) + exp(rec.bg[idSpR[j]])
    recruits[j] ~ dnegbin(prob = nb.size[idSpR[j]] / (nb.size[idSpR[j]] + lambdaR[j]),
                          size = nb.size[idSpR[j]])
  }

  # ═══════════════════ 新个体大小分布 c(z') ═══════════════════
  # ★ 只用有面积的新株；灌木新株九成以上记为点（0.25 cm2），IPM 里点的大小固定
  for (k in 1:nSp){
    pop.sizeRec[k] ~ dnorm(0, sd = 5)
    sizeRec.sd[k]  ~ dexp(1)
  }
  for (j in 1:nRecSize){
    sizeRec[j] ~ dnorm(pop.sizeRec[idSpRS[j]], sd = sizeRec.sd[idSpRS[j]])
  }
})


# ═══════════════════════════════════════════════════════════════════════════
# data / constants
# ═══════════════════════════════════════════════════════════════════════════
partA_input <- readRDS(file.path(der_dir, "CLR_PartA_input.rds"))
sp_tab <- partA_input$species
rec    <- partA_input$main$rec
recr   <- partA_input$main$recruit
recsz  <- filter(partA_input$main$recsize, !point)
if (test_run) {
    # Test: a subset of species, renumbered 1..n
  sp_tab <- filter(sp_tab, code %in% test_species) %>% mutate(k_new = row_number())
  rec  <- inner_join(rec,  select(sp_tab, k, k_new), by = "k") %>% mutate(k = k_new, i = match(trackID, unique(trackID)))
  recr <- inner_join(recr, select(sp_tab, k, k_new), by = "k") %>% mutate(k = k_new)
  recsz <- inner_join(recsz, select(sp_tab, k, k_new), by = "k") %>% mutate(k = k_new)
  sp_tab$k <- sp_tab$k_new
  n_chains <- 2; n_iter <- 200; n_burn <- 100; n_thin <- 1
}
W_z <- c("W_intra_z", "W_grass_z", "W_shrub_z", "W_forb_z")
grow <- which(rec$grow_ok)

  # Trait (stage 2): one draw, standardised over species
trait <- rep(0, nrow(sp_tab))
if (stage == 2) {
  stopifnot(!is.null(trait_file))
  trait <- readRDS(trait_file)[trait_draw, sp_tab$code]
  trait <- (trait - mean(trait)) / sd(trait)
}

constList <- list(
  nSp = nrow(sp_tab), nYr = nrow(partA_input$years), nQuad = nrow(partA_input$quadrats),
  nGrp = nrow(partA_input$groups), nInd = max(rec$i), nVR = 3, nX = 4,
  nSurv = nrow(rec), nGrow = length(grow), nRec = nrow(recr), nRecSize = nrow(recsz),
  useTrait = as.integer(stage == 2), trait = trait,
  useIndS = as.integer(plant_effects != "none"), useIndG = as.integer(plant_effects == "both"),
    # 存活纪录的下标（生长借用 idGrow 指回存活纪录）
  idSp = rec$k, idYr = rec$t, idQuad = rec$q, idGrp = rec$g, idInd = rec$i, idGrow = grow,
  size = rec$z, size2 = rec$z^2, graze = rec$E, W = as.matrix(rec[W_z]),
  clim = partA_input$years$clim_anom,
    # 补充
  idSpR = recr$k, idYrR = recr$t, idQuadR = recr$q, idGrpR = recr$g, grazeR = recr$E,
  covQ = recr$cover_q, covG = recr$cover_g,
  CRo = as.matrix(recr[c("C_grass_z", "C_shrub_z", "C_forb_z")]),
    # 种内盖度的标准化常数（每个物种一组，按 sp_tab 的顺序）
  cIntraMean = partA_input$scaling$C_intra$mean[match(sp_tab$code, partA_input$scaling$C_intra$code)],
  cIntraSd   = partA_input$scaling$C_intra$sd[match(sp_tab$code, partA_input$scaling$C_intra$code)],
  idSpRS = recsz$k
)
dataList <- list(surv = rec$survives_tplus1, sizeNext = rec$z_next[grow],
                 recruits = recr$recruits, sizeRec = recsz$z_rec)

# ── 建模前的检查 ──────────────────────────────────────────────────────────
with(constList, stopifnot(
  all(sort(unique(idSp)) == 1:nSp), all(idYr %in% 1:nYr), all(idQuad %in% 1:nQuad), all(idGrp %in% 1:nGrp),
  all(sort(unique(idInd)) == 1:nInd), length(idGrow) == nGrow, max(idGrow) <= nSurv,
  all(is.finite(size)), all(is.finite(W)), all(is.finite(clim)), all(graze %in% 0:1),
  all(covG > 0 | covQ == 0), all(is.finite(CRo)), all(idSpRS %in% 1:nSp),
  all(is.finite(cIntraMean)), all(cIntraSd > 0)))
stopifnot(all(dataList$surv %in% 0:1), all(is.finite(dataList$sizeNext)), all(dataList$recruits >= 0),
          all(rec$survives_tplus1[grow] == 1))
cat(sprintf("Part A stage %d: %d species, %d survival, %d growth, %d recruitment rows, %d recruit sizes, %d plants, %d years.\n",
            stage, constList$nSp, constList$nSurv, constList$nGrow, constList$nRec, constList$nRecSize,
            constList$nInd, constList$nYr))


# ═══════════════════════════════════════════════════════════════════════════
# inits
# ═══════════════════════════════════════════════════════════════════════════
# ★ logit / log 尺度上的起始值一律压在 sd ≤ 0.5；函数形式，各链从不同地方出发
inits <- function() {
  nSp <- constList$nSp; nYr <- constList$nYr
  l <- list(
    pop.int.mu = c(rnorm(2, 0, 0.5), rnorm(1, -1, 0.3)), pop.int.sd = rexp(3, 2),
    pop.slope.size.mu = rnorm(2, 0.5, 0.2), pop.slope.size.sd = rexp(2, 2),
    pop.slope.size2.mu = rnorm(2, 0, 0.1), pop.slope.size2.sd = rexp(2, 2),
    pop.slope.graze.mu = rnorm(3, 0, 0.2), pop.slope.graze.sd = rexp(3, 2),
    pop.slope.clim.mu = rnorm(3, 0, 0.2), pop.slope.clim.sd = rexp(3, 2),
    pop.slope.W.mu = matrix(rnorm(12, 0, 0.1), 3, 4), pop.slope.W.sd = matrix(rexp(12, 4), 3, 4),
    pop.int = cbind(rnorm(nSp, 0, 0.3), rnorm(nSp, 0, 0.3), rnorm(nSp, -1, 0.3)),
    pop.slope.size = matrix(rnorm(nSp * 2, 0.5, 0.2), nSp, 2),
    pop.slope.size2 = matrix(rnorm(nSp * 2, 0, 0.1), nSp, 2),
    pop.slope.graze = matrix(rnorm(nSp * 3, 0, 0.1), nSp, 3),
    pop.slope.clim = matrix(rnorm(nSp * 3, 0, 0.1), nSp, 3),
    pop.slope.W = array(rnorm(nSp * 12, 0, 0.1), c(nSp, 3, 4)),
    com.size.x.clim = rnorm(2, 0, 0.1),
    com.W.x.clim = matrix(rnorm(12, 0, 0.1), 3, 4), com.W.x.graze = matrix(rnorm(12, 0, 0.1), 3, 4),
    com.W.x.size = matrix(rnorm(8, 0, 0.1), 2, 4),
      # 年份层（sd 起始值给小的，免得一开始就把年份效应吃掉降水效应）
    yr.int.sd = rexp(3, 4), yrsp.int.sd = rexp(3, 4), yr.slope.size.sd = rexp(2, 4),
    yr.slope.W.sd = matrix(rexp(12, 8), 3, 4),
    yr.int = matrix(rnorm(nYr * 3, 0, 0.1), nYr, 3),
    yrsp.int = array(rnorm(nYr * nSp * 3, 0, 0.05), c(nYr, nSp, 3)),
    yr.slope.size = matrix(rnorm(nYr * 2, 0, 0.05), nYr, 2),
    yr.slope.W = array(rnorm(nYr * 12, 0, 0.05), c(nYr, 3, 4)),
    grp.sd = rexp(3, 4), quad.sd = rexp(3, 4),
    grp.int = matrix(rnorm(constList$nGrp * 3, 0, 0.1), constList$nGrp, 3),
    quad.int = matrix(rnorm(constList$nQuad * 3, 0, 0.1), constList$nQuad, 3),
    grow.sd = rexp(nSp, 2) + 0.3,
    p.eff = runif(nSp, 0.3, 0.7), rec.bg.mu = rnorm(1, -2, 0.3), rec.bg.sd = rexp(1, 2),
    rec.bg = rnorm(nSp, -2, 0.3), nb.size = rgamma(nSp, 4, 2),
    pop.sizeRec = rnorm(nSp, -1, 0.3), sizeRec.sd = rexp(nSp, 2) + 0.3
  )
    # 年份效应从它们的平均附近出发（年份效应中心化在截距与降水上）
  for (v in 1:3) {
    l$yr.int[, v] <- l$yr.int[, v] + l$pop.int.mu[v]
    l$yrsp.int[, , v] <- sweep(l$yrsp.int[, , v], 2, l$pop.int[, v] - l$pop.int.mu[v], "+") + l$yr.int[, v]
  }
    # 植株效应：各链的 SD 从不同大小出发（0.1 到 1.5），可以看出链是否走到同一处
  n_ind_v <- constList$useIndS + constList$useIndG
  if (n_ind_v > 0) {
    l$ind.sd  <- runif(n_ind_v, 0.1, 1.5)
    l$ind.raw <- matrix(rnorm(constList$nInd * n_ind_v), constList$nInd, n_ind_v)
  }
  if (n_ind_v == 2) l$ind.rho01 <- runif(1, 0.3, 0.7)
  if (stage == 2) {
    l$pop.slope.graze.trait <- rnorm(3, 0, 0.1); l$pop.slope.clim.trait <- rnorm(3, 0, 0.1)
    l$pop.slope.W.trait <- matrix(rnorm(12, 0, 0.1), 3, 4)
  }
  l
}

params <- c(
    # 物种层超参数
  "pop.int.mu", "pop.int.sd", "pop.slope.size.mu", "pop.slope.size.sd", "pop.slope.size2.mu", "pop.slope.size2.sd",
  "pop.slope.graze.mu", "pop.slope.graze.sd", "pop.slope.clim.mu", "pop.slope.clim.sd",
  "pop.slope.W.mu", "pop.slope.W.sd",
    # 交互（个体层、群落层）
  "com.size.x.clim", "com.W.x.clim", "com.W.x.graze", "com.W.x.size",
    # 年份层、样方层、个体层
  "yr.int.sd", "yrsp.int.sd", "yr.slope.size.sd", "yr.slope.W.sd",
  "yr.int", "yrsp.int", "yr.slope.size", "yr.slope.W",
  "grp.sd", "quad.sd", "grp.int", "quad.int",
    # 补充与新株大小
  "p.eff", "rec.bg", "rec.bg.mu", "rec.bg.sd", "nb.size", "pop.sizeRec", "sizeRec.sd",
    # IPM 后处理需要
  "pop.int", "pop.slope.size", "pop.slope.size2", "pop.slope.graze", "pop.slope.clim", "pop.slope.W", "grow.sd")
if (plant_effects != "none") params <- c(params, "ind.sd")
if (plant_effects == "both")  params <- c(params, "ind.rho")
if (stage == 2) params <- c(params, "pop.slope.graze.trait", "pop.slope.clim.trait", "pop.slope.W.trait")


# ═══════════════════════════════════════════════════════════════════════════
# 建模、编译、MCMC
# ═══════════════════════════════════════════════════════════════════════════
t0 <- Sys.time()
m_partA <- nimbleModel(model, constants = constList, data = dataList, inits = inits(), calculate = FALSE)
  # 起始值下的对数似然必须是有限值（-Inf 多半是起始值太大或资料有误）
ll0 <- m_partA$calculate()
cat(sprintf("Model built (%.1f min); log probability at the inits: %.1f\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), ll0))
stopifnot(is.finite(ll0))
rm(m_partA); invisible(gc())   # the check model leaves memory; the chains build their own
  # Build, compile and run the chains ('run_chains', CLR_fns_nimble.R): one seed per chain, so a
  # rerun with the same n_cores gives the same samples
samples <- run_chains(model, constList, dataList, inits, monitors = params, thin = n_thin, n_iter = n_iter,
                      n_burn = n_burn, n_chains = n_chains, n_cores = n_cores, seed = 2026,
                      log_file = file.path(der_dir, sprintf("CLR_PartA_stage%d_chains.log", stage)))
cat(sprintf("MCMC done (%.1f min).\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Save -------------------------------------------------------------------------
  # (sensitivity runs of the plant effects get their own file: '_plantsurvival', '_plantnone')
out_file <- file.path(der_dir, sprintf("CLR_PartA_stage%d%s%s%s.rds", stage,
                                       if (stage == 2) sprintf("_draw%03d", trait_draw) else "",
                                       if (plant_effects != "both") paste0("_plant", plant_effects) else "",
                                       if (test_run) "_test" else ""))
saveRDS(list(samples = samples, species = sp_tab, constants = constList[c("nSp", "nYr", "nQuad", "nGrp", "nInd")],
             scaling = partA_input$scaling, settings = list(stage = stage, plant_effects = plant_effects, n_iter = n_iter,
                                                            n_burn = n_burn, n_thin = n_thin, n_chains = n_chains,
                                                            test_run = test_run)),
        out_file)
cat("Saved:\n  data/data.derived/", site_tag, ": ", basename(out_file), "\n", sep = "")
rm(list = intersect(c("rec", "recr", "recsz", "grow", "W_z", "t0", "ll0"), ls()))


# ═══════════════════════════════════════════════════════════════════════════
# 要跑的版本（写进 Methods）
# ═══════════════════════════════════════════════════════════════════════════
# A. 完整版（上面这个）
# B. 去掉 yr.slope.W，只留 yr.int / yrsp.int
#    → 比较 com.W.x.clim（δ_x）的后验缩水多少 = 这个交互有多少是靠先验撑着的
# C. yr.*.sd 的先验换成 dexp(1)
#    → 比较 δ_x 的 CI 移动多少
# D. 放牧换成连续的放牧量（标准化的羊日数）；或去掉 E，只留组的随机效应
# E. 只用有幼苗标记的物种；合并（coalescence）当作死亡
# F. 植株随机效应：plant_effects = "survival"（只放存活的个体效应）、"none"（不放个体效应）
