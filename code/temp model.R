# ═══════════════════════════════════════════════════════════════════════════
# Alpine 4 vital-rate IPM
#
# 待确认
#   个体随机效应叫 ind.eff；对比量前缀用 con.
#   v.surv 等作为下标是 constants 里的整数，NIMBLE 建模时直接替换。
#      若报错改回字面量 1–4。nVR 作矩阵维度同理。
#   零截断 NB 目前用 T(...,1,)。封闭形式是 1 - p^size，
#      profile 之后可换成 nimbleFunction 自订分布。
#
# trait[k]（生活史性状）的来源与用法（2026-10 定案，方案 A）
#   Stage 1：只估 baseline demography（不含 trait），对每个 posterior draw s
#            建 IPM_k^(s) → 算 trait_k^(s)（GenTime、Iteroparity 等）。
#   Stage 2：就是本模型。trait 作为 constants 传入；抽 M 组完整的
#            trait^(s)（保留偏斜与物种间相关），每组跑一次，最后合并后验
#            （posterior-draw propagation / multiple imputation，Plummer 2015 的 cut）。
#   trait 只解释 γ（mic）、β_intra、β_inter、β_precip；
#   α0 / α1 / α2（pop.int、pop.slope.size、pop.slope.size2）不放 trait，
#   因为 trait 正是由它们算出来的（放了就是同义反复）。
# ═══════════════════════════════════════════════════════════════════════════

library(nimble)

model <- nimbleCode({
  
  # ═══════════════════ 年份层：降水（固定）+ 年份偏离（随机） ═══════════════════
  # ★ precip 必须先置中，否则与 yr.int 和两个交互项绞在一起
  for (v in 1:nVR){
    com.size.x.precip[v]  ~ dnorm(0, sd = 2)  # 降水 × 自身大小（个体层调节降水响应）
    com.INTRA.x.precip[v] ~ dnorm(0, sd = 2)  # 降水 × 种内
    com.INTER.x.precip[v] ~ dnorm(0, sd = 2)  # 降水 × 种间（时间版 SGH）
    
    yr.int.sd[v]         ~ dexp(2)   # ★ 这五个先验直接决定四个降水调节项的后验宽度
    yr.slope.size.sd[v]  ~ dexp(2)   #   紧的 sd 先验会把变异推给 com.*.x.precip
    yr.slope.INTRA.sd[v] ~ dexp(2)   #   必做敏感度：dexp(2) vs dexp(1)
    yr.slope.INTER.sd[v] ~ dexp(2)
    yrsp.int.sd[v]       ~ dexp(2)   #   对应 pop.slope.precip：年×物种的残余
  }
  
  for (t in 1:nYr){
    for (v in 1:nVR){
      yr.int[t,v]         ~ dnorm(0, sd = yr.int.sd[v])          # 该年整体好坏
      yr.slope.size[t,v]  ~ dnorm(0, sd = yr.slope.size.sd[v])   # 该年大小斜率偏离直线多少
      yr.slope.INTRA[t,v] ~ dnorm(0, sd = yr.slope.INTRA.sd[v])  # 该年种内斜率偏离直线多少
      yr.slope.INTER[t,v] ~ dnorm(0, sd = yr.slope.INTER.sd[v])  # 该年种间斜率偏离直线多少
    }
  }
  
  # ═══════════════════ 年 × 物种层：该年该物种的残余好坏 ═══════════════════
  # 四个降水调节者各自对应一个年份层的竞争对手，缺一个那一项就会被高估。
  # size / INTRA / INTER 对应上面三条随机斜率，pop.slope.precip 对应这一条。
  for (t in 1:nYr){
    for (k in 1:nSp){
      for (v in 1:nVR){
        yrsp.int[t,k,v] ~ dnorm(0, sd = yrsp.int.sd[v])
      }
    }
  }
  
  # ═══════════════════ 物种层 · 先验 ═══════════════════
  for (v in 1:nVR){
    # α0 / α1 / α2：只有物种随机效应，不放 trait（trait 由它们导出）
    pop.int.mu[v]    ~ dnorm(0, sd = 5)   # 物种截距的总均值
    pop.int.sd[v]    ~ dexp(1)            # 物种间差异
    
    pop.slope.size.mu[v]    ~ dnorm(0, sd = 2)   # 物种大小斜率的均值
    pop.slope.size.sd[v]    ~ dexp(1)
    
    pop.slope.size2.mu[v]    ~ dnorm(0, sd = 2)  # 物种大小平方项的均值
    pop.slope.size2.sd[v]    ~ dexp(1)           # ★ 若压到 0 → 曲率一致，可收回 com
    
    # γ / β_intra / β_inter / β_precip：三件套，trait 只在这里出现
    
    pop.slope.INTRA.mu[v]    ~ dnorm(0, sd = 2)  # 物种种内斜率的三件套
    pop.slope.INTRA.trait[v] ~ dnorm(0, sd = 2)  # ★ 什么性状的物种更怕种内竞争
    pop.slope.INTRA.sd[v]    ~ dexp(1)
    
    pop.slope.INTER.mu[v]    ~ dnorm(0, sd = 2)  # 物种种间斜率的三件套
    pop.slope.INTER.trait[v] ~ dnorm(0, sd = 2)  # ★ 什么性状的物种更怕种间竞争
    pop.slope.INTER.sd[v]    ~ dexp(1)
    
    pop.slope.precip.mu[v]    ~ dnorm(0, sd = 2) # 物种降水斜率的三件套
    pop.slope.precip.trait[v] ~ dnorm(0, sd = 2) # ★ 什么性状的物种对降水更敏感
    pop.slope.precip.sd[v]    ~ dexp(1)          # ★ 压到 0 → 各物种降水响应一致
    
    for (m in 1:nMic){                           # 微环境三轴，每轴一套
      pop.slope.mic.mu[v,m]    ~ dnorm(0, sd = 2)
      pop.slope.mic.trait[v,m] ~ dnorm(0, sd = 2)
      pop.slope.mic.sd[v,m]    ~ dexp(1)         # ★ 压到 0 → 各物种环境响应一致
    }
  }
  
  # ═══════════ 物种层 · 每个物种的系数 ═══════════
  for (k in 1:nSp){
    for (v in 1:nVR){
      pop.int[k,v]         ~ dnorm(pop.int.mu[v],         sd = pop.int.sd[v])
      pop.slope.size[k,v]  ~ dnorm(pop.slope.size.mu[v],  sd = pop.slope.size.sd[v])
      pop.slope.size2[k,v] ~ dnorm(pop.slope.size2.mu[v], sd = pop.slope.size2.sd[v])
      pop.slope.INTRA[k,v] ~
        dnorm(pop.slope.INTRA.mu[v] + pop.slope.INTRA.trait[v] * trait[k],
              sd = pop.slope.INTRA.sd[v])
      pop.slope.INTER[k,v] ~
        dnorm(pop.slope.INTER.mu[v] + pop.slope.INTER.trait[v] * trait[k],
              sd = pop.slope.INTER.sd[v])
      pop.slope.precip[k,v] ~
        dnorm(pop.slope.precip.mu[v] + pop.slope.precip.trait[v] * trait[k],
              sd = pop.slope.precip.sd[v])
      
      for (m in 1:nMic){
        pop.slope.mic[k,v,m] ~
          dnorm(pop.slope.mic.mu[v,m] + pop.slope.mic.trait[v,m] * trait[k],
                sd = pop.slope.mic.sd[v,m])
      }
    }
  }
  # ★ pop.slope.mic 是 nSp × nVR × nMic = 17×4×3 = 204 个节点，
  #   现在是模型里最大的物种层区块，混合问题会最先从这里显现。
  
  # ═══════════════════ 全局层 ═══════════════════
  # 注意这里的不对称是有意的：微环境【主效应】按物种估（控制变量，要estimate准），
  # 微环境【× 竞争】仍然共享（核心结果，要集中全部样本量）。
  # 降水的部分相反：主效应按物种估（pop.slope.precip），三个调节者都在这一层。
  for (v in 1:nVR){
    com.INTRA.x.size[v] ~ dnorm(0, sd = 2)   # 种内 × 自身大小
    com.INTER.x.size[v] ~ dnorm(0, sd = 2)   # 种间 × 自身大小
    com.INTRA.x.mic[v]  ~ dnorm(0, sd = 2)   # 种内 × 微环境轴1
    com.INTER.x.mic[v]  ~ dnorm(0, sd = 2)   # ★ 种间 × 微环境轴1（空间版 SGH）
  }
  
  # ═══════════════════ 个体层（相关版：四个生命率共用一个 4×4 相关矩阵） ═══════════════════
  # 同一个体在四个生命率上的偏离可正可负相关；LKJ(eta = 2) 以 0 为中心，方向由资料决定。
  # ★ 37% 个体只有一笔纪录、69% 从未开花：涉及 v.flow / v.infl 的相关基本是先验在说话，
  #   surv–grow 较有资料支撑。报告前先用模拟确认哪些相关可恢复。
  # ★ 单笔纪录的个体，ind.eff[i,v.grow] 与 grow.sd 仍然混杂（相关版不解决这点）。
  for (v in 1:nVR){ ind.sd[v] ~ dexp(1) }
  ind.chCor[1:nVR,1:nVR] ~ dlkj_corr_cholesky(eta = 2, p = nVR)
  ind.chCov[1:nVR,1:nVR] <- uppertri_mult_diag(ind.chCor[1:nVR,1:nVR], ind.sd[1:nVR])
  for (i in 1:nInd){
    ind.eff[i,1:nVR] ~ dmnorm(zeros[1:nVR],
                              cholesky = ind.chCov[1:nVR,1:nVR], prec_param = 0)
  }
  ind.cor[1:nVR,1:nVR] <- t(ind.chCor[1:nVR,1:nVR]) %*% ind.chCor[1:nVR,1:nVR]   # 相关矩阵本身
  
  # ── 对角版（退路：相关版混合不良时换回）──────────────────────────────
  # for (v in 1:nVR){ ind.sd[v] ~ dexp(1) }
  # for (i in 1:nInd){
  #   for (v in 1:nVR){ ind.eff[i,v] ~ dnorm(0, sd = ind.sd[v]) }
  # }
  
  # ═══════════════════ Linear Predictor ═══════════════════
  # lp 对【所有活体纪录】都算，观测层再各自挑子集
  for (r in 1:nRecord){
    for (v in 1:nVR){
      lp[r,v] <-
        pop.int[idSp[r],v] +
        pop.slope.size[idSp[r],v]  * size[r] +
        pop.slope.size2[idSp[r],v] * size2[r] +
        pop.slope.mic[idSp[r],v,1] * mic[idInd[r],1] +
        pop.slope.mic[idSp[r],v,2] * mic[idInd[r],2] +
        pop.slope.mic[idSp[r],v,3] * mic[idInd[r],3] +
        pop.slope.INTRA[idSp[r],v] * INTRA[r] +
        pop.slope.INTER[idSp[r],v] * INTER[r] +
        com.INTRA.x.size[v]   * INTRA[r] * size[r] +
        com.INTER.x.size[v]   * INTER[r] * size[r] +
        com.INTRA.x.mic[v]    * INTRA[r] * mic[idInd[r],1] +
        com.INTER.x.mic[v]    * INTER[r] * mic[idInd[r],1] +
        pop.slope.precip[idSp[r],v] * precip[idYr[r]] +
        com.size.x.precip[v]  * size[r]  * precip[idYr[r]] +
        com.INTRA.x.precip[v] * INTRA[r] * precip[idYr[r]] +
        com.INTER.x.precip[v] * INTER[r] * precip[idYr[r]] +
        yr.int[idYr[r],v] +
        yr.slope.size[idYr[r],v]  * size[r] +
        yr.slope.INTRA[idYr[r],v] * INTRA[r] +
        yr.slope.INTER[idYr[r],v] * INTER[r] +
        yrsp.int[idYr[r],idSp[r],v] +
        ind.eff[idInd[r],v]
    }
  }
  
  # ═══════════════════ 对比量：种间 − 种内 ═══════════════════
  for (v in 1:nVR){
    con.slope[v]  <- pop.slope.INTER.mu[v] - pop.slope.INTRA.mu[v]  # 平均竞争强度
    con.mic[v]    <- com.INTER.x.mic[v]    - com.INTRA.x.mic[v]     # ★ 空间：环境调节
    con.precip[v] <- com.INTER.x.precip[v] - com.INTRA.x.precip[v]  # 时间：降水调节
  }
  
  # ═══════════════════ 观测层 · 四个子集各有自己的索引 ═══════════════════
  grow.sd ~ dexp(1)     # 生长的残差 sd，标准化 log 尺度上的量级
  nb.size ~ dexp(0.1)   # 负二项的 size 参数
  
  for (s in 1:nSurv){                              # ① 存活：t→t+1，最后一年没有
    surv[s] ~ dbern(ilogit(lp[idSurv[s], v.surv]))
  }
  
  for (g in 1:nGrow){                              # ② 生长：存活且量到 t+1 大小
    sizeNext[g] ~ dnorm(lp[idGrow[g], v.grow], sd = grow.sd)
  }
  # ★ size 与 sizeNext 必须走同一个变换（见下方 data 区）。两边尺度不同时，
  #   IPM 的 g(z'|z) 进去的 z 和出来的 z' 对不上，核没法迭代。
  
  for (r in 1:nRecord){                            # ③ 是否开花：t 当年，全部活体
    flow[r] ~ dbern(ilogit(lp[r, v.flow]))
  }
  
  for (f in 1:nFlow){                              # ④ 花序数：仅开花者
    nb.mu[f] <- exp(lp[idFlow[f], v.infl])         # 未截断的均值
    nb.p[f]  <- nb.size / (nb.size + nb.mu[f])
    infl[f] ~ T(dnegbin(nb.p[f], nb.size), 1, )    # 零截断
  }
  # ★ 后处理：E[infl | 开花] = nb.mu / (1 - nb.p^nb.size)，不是 exp(lp[,v.infl])
  
  # ═══════════════════ 补充子模型（物种 × 年） ═══════════════════
  pop.seed.mu ~ dnorm(0, sd = 2);   pop.seed.sd ~ dexp(1)
  pop.veg.mu  ~ dnorm(0, sd = 2);   pop.veg.sd  ~ dexp(1)
  pop.seed0.mu ~ dnorm(0, sd = 2);  pop.seed0.sd ~ dexp(1)  # 背景补充：与当年花序数无关的那一份
  
  for (k in 1:nSp){
    pop.seed.log[k] ~ dnorm(pop.seed.mu, sd = pop.seed.sd)
    pop.seed[k] <- exp(pop.seed.log[k])   # 每花序产生多少幼苗
    pop.veg.log[k]  ~ dnorm(pop.veg.mu,  sd = pop.veg.sd)
    pop.veg[k]  <- exp(pop.veg.log[k])    # 每单位同种盖度产生多少无性繁殖体
    pop.seed0.log[k] ~ dnorm(pop.seed0.mu, sd = pop.seed0.sd)
    pop.seed0[k] <- exp(pop.seed0.log[k]) # 种子库与样方外传入，每种每年一个常数速率
  }
  # ★ 连接函数不能放在 ~ 左边（NIMBLE），所以写成 log 尺度的随机节点 + exp 的确定性节点
  
  # ★ obsInfl[k,t] 是 t 年的花序总数，recSeed[k,t] 是 t+1 年出现的种子实生苗。
  #   差一格不会报错，只会让 pop.seed 整体偏掉，建资料时用 stopifnot 钉死。
  # ★ 零花年份是常态：Pseudocymopterus 5/8 年、Penstemon 3/8 年、Viola 2/8 年无花，
  #   119 个物种×transition 格子里有 10 格 obsInfl = 0。其中 Penstemon 2015 年 0 朵花
  #   而 2016 年有 3 株实生苗，所以补充不能只由花序数驱动；pop.seed0 承接这一份。
  #   没有它时这些格子的 rate 不依赖 pop.seed，对参数完全无资讯。
  for (k in 1:nSp){
    for (t in 1:nTrans){
      recSeed[k,t] ~ dpois(pop.seed0[k] + pop.seed[k] * obsInfl[k,t])
      recVeg[k,t]  ~ dpois(pop.veg[k]   * obsCov[k,t] + 1e-8)
    }
  }
  
  # ═══════════════════ 新个体大小分布 c(z') ═══════════════════
  sizeRec.sd ~ dexp(1)
  for (k in 1:nSp){ pop.sizeRec[k] ~ dnorm(0, sd = 5) }
  for (j in 1:nRecruit){
    sizeRec[j] ~ dnorm(pop.sizeRec[idSpRec[j]], sd = sizeRec.sd)
  }
  
})


# ═══════════════════════════════════════════════════════════════════════════
# data / constants
# ═══════════════════════════════════════════════════════════════════════════

# dataList <- list(
#   surv     = surv,       # 0/1,    长度 nSurv     = 8196
#   sizeNext = sizeNext,   # 连续,   长度 nGrow     = 5720
#   flow     = flow,       # 0/1,    长度 nRecord   = 8842
#   infl     = infl,       # 正整数, 长度 nFlow     = 2879
#   recSeed  = recSeed,    # nSp × nTrans
#   recVeg   = recVeg,     # nSp × nTrans
#   sizeRec  = sizeRec     # 连续,   长度 nRecruit  = 2115
# )

# constList <- list(
#   # 维度（数字是你 review 里算出来的，跑之前 stopifnot 再核一次）
#   nSp = 17,              # ★ 不是 18：NA NA 和 Indet indet 已排除
#   nYr = nYr, nTrans = nYr - 1, nVR = 4, nMic = 3,
#   nInd = 3126, nRecord = 8842,
#   nSurv = 8196, nGrow = 5720, nFlow = 2879, nRecruit = 2115,
#   # ★ nRecruit 不是 nRecordruit —— 别被 nRec → nRecord 的全局替换改到
#   # 生命率编号
#   v.surv = 1, v.grow = 2, v.flow = 3, v.infl = 4,
#   zeros = rep(0, 4),     # 个体随机效应 dmnorm 的均值向量
#   # 下标向量（都指回 1:nRecord 里的某一笔；必须是连续正整数）
#   idSp = idSp, idInd = idInd, idYr = idYr,
#   idSurv = idSurv, idGrow = idGrow, idFlow = idFlow, idSpRec = idSpRec,
#   # 协变量
#   size = size, size2 = size2, mic = mic,
#   trait = trait_draw,    # ★ Stage 1 的一组 posterior draw（长度 nSp，已置中/标准化），每个 imputation 换一组
#   precip = precip,       # ★ 必须置中： precip <- precip - mean(precip)
#   INTRA = INTRA, INTER = INTER,
#   obsInfl = obsInfl, obsCov = obsCov
# )

# ── 建模前的检查 ──────────────────────────────────────────────────────────
# stopifnot(all(sort(unique(idSp))  == 1:nSp))
# stopifnot(all(sort(unique(idYr))  == 1:nYr))
# stopifnot(length(surv)     == nSurv,  max(idSurv) <= nRecord)
# stopifnot(length(sizeNext) == nGrow,  max(idGrow) <= nRecord)
# stopifnot(length(infl)     == nFlow,  max(idFlow) <= nRecord)
# stopifnot(length(flow)     == nRecord)
# stopifnot(all(infl >= 1))                      # 零截断：不能有 0
# stopifnot(abs(mean(precip)) < 1e-8)            # precip 已置中
# # 存活子集应该正好排除最后一年
# stopifnot(!any(idYr[idSurv] == nYr))
# # size 与 sizeNext 同尺度、无 -Inf（log(0) 会从这里漏进来）
# stopifnot(all(is.finite(size)), all(is.finite(sizeNext)), all(is.finite(size2)))
# stopifnot(abs(mean(size)) < 1e-8)
# # obsInfl 取 t 年、recSeed 取 t+1 年，列必须错开一格
# stopifnot(colnames(obsInfl) == yrs[1:nTrans])
# stopifnot(colnames(recSeed) == yrs[2:nYr])

# ── size 与 sizeNext 的变换 ──────────────────────────────────────────────
# Length 是正值右偏，原始 cm 尺度上固定变异数的常态会对小个体预测出负的大小，
# IPM 的 g(z'|z) 在负值域有质量就没法正规化。所以在 log 尺度上建模。
# ★ 两边套同一组 mean / sd，并把 z.mean 与 z.sd 存下来跟着模型走：
#   建 IPM 核时 z' 网格要用同一个变换才能反推回 cm。
#   z.mean   <- mean(log(Length)); z.sd <- sd(log(Length))
#   size     <- (log(Length)     - z.mean) / z.sd
#   sizeNext <- (log(LengthNext) - z.mean) / z.sd
#   size2    <- size^2
# 物种内正交化（让 size 与 size2 不相关）能改善混合，但建核时 z' 网格要套同一个
# poly()，否则系数对不上。嫌麻烦就先用 size^2，看 traceplot 再决定。

# params <- c(
#   # 物种层超参数
#   "pop.int.mu","pop.int.sd",
#   "pop.slope.size.mu","pop.slope.size.sd",
#   "pop.slope.size2.mu","pop.slope.size2.sd",
#   "pop.slope.mic.mu","pop.slope.mic.trait","pop.slope.mic.sd",
#   "pop.slope.INTRA.mu","pop.slope.INTRA.trait","pop.slope.INTRA.sd",
#   "pop.slope.INTER.mu","pop.slope.INTER.trait","pop.slope.INTER.sd",
#   "pop.slope.precip.mu","pop.slope.precip.trait","pop.slope.precip.sd",
#   # 全局层
#   "com.size.x.precip","com.INTRA.x.size","com.INTER.x.size",
#   "com.INTRA.x.mic","com.INTER.x.mic",
#   "com.INTRA.x.precip","com.INTER.x.precip",
#   # 年份层
#   "yr.int.sd","yr.slope.size.sd","yr.slope.INTRA.sd","yr.slope.INTER.sd","yrsp.int.sd",
#   "yr.int","yr.slope.size","yr.slope.INTRA","yr.slope.INTER",
#   # 个体层
#   "ind.sd", "ind.cor",   # ★ ind.cor 就是四个生命率之间的个体层相关
#   # 对比量（★ 不写就不保存）
#   "con.slope","con.mic","con.precip",
#   # 观测层与补充
#   "grow.sd","nb.size",
#   "pop.seed.mu","pop.seed.sd","pop.veg.mu","pop.veg.sd",
#   "pop.seed0.mu","pop.seed0.sd",
#   "pop.sizeRec","sizeRec.sd",
#   # IPM 后处理需要
#   "pop.int","pop.slope.size","pop.slope.size2","pop.slope.mic",
#   "pop.slope.INTRA","pop.slope.INTER","pop.slope.precip","pop.seed","pop.veg","pop.seed0"
# )


# ═══════════════════════════════════════════════════════════════════════════
# inits
# ═══════════════════════════════════════════════════════════════════════════
# NIMBLE 缺 inits 时从先验随机抽。pop.int.mu 的先验是 dnorm(0, sd = 5)，
# 抽到 ±10 很正常，ilogit(±10) 就是 0 或 1，配上 Bernoulli 资料直接 -Inf。
# ★ 重点不是给值，是不能给大值：logit 尺度上的起始值一律压在 sd ≤ 0.5。
# ★ 必须是函数不是固定值，四条链才会从不同地方出发，Rhat 才有意义。

# inits <- function() list(
#   # 物种层超参数
#   pop.int.mu    = rnorm(4, 0, 0.5),   pop.int.sd    = rexp(4, 2),
#   pop.slope.size.mu   = rnorm(4, 0, 0.3), pop.slope.size.sd   = rexp(4, 2),
#   pop.slope.size2.mu  = rnorm(4, 0, 0.2), pop.slope.size2.sd  = rexp(4, 2),
#   pop.slope.INTRA.mu  = rnorm(4, 0, 0.3), pop.slope.INTRA.trait = rnorm(4, 0, 0.2),
#   pop.slope.INTRA.sd  = rexp(4, 2),
#   pop.slope.INTER.mu  = rnorm(4, 0, 0.3), pop.slope.INTER.trait = rnorm(4, 0, 0.2),
#   pop.slope.INTER.sd  = rexp(4, 2),
#   pop.slope.precip.mu = rnorm(4, 0, 0.2), pop.slope.precip.trait = rnorm(4, 0, 0.2),
#   pop.slope.precip.sd = rexp(4, 2),
#   pop.slope.mic.mu    = matrix(rnorm(4*3, 0, 0.2), 4, 3),
#   pop.slope.mic.trait = matrix(rnorm(4*3, 0, 0.2), 4, 3),
#   pop.slope.mic.sd    = matrix(rexp(4*3, 2), 4, 3),
#   # 物种层系数（给了超参数还得给这一层，否则 NIMBLE 仍旧从先验抽）
#   pop.int         = matrix(rnorm(nSp*4, 0, 0.3), nSp, 4),
#   pop.slope.size  = matrix(rnorm(nSp*4, 0, 0.3), nSp, 4),
#   pop.slope.size2 = matrix(rnorm(nSp*4, 0, 0.2), nSp, 4),
#   pop.slope.INTRA = matrix(rnorm(nSp*4, 0, 0.2), nSp, 4),
#   pop.slope.INTER = matrix(rnorm(nSp*4, 0, 0.2), nSp, 4),
#   pop.slope.precip= matrix(rnorm(nSp*4, 0, 0.2), nSp, 4),
#   pop.slope.mic   = array(rnorm(nSp*4*3, 0, 0.2), c(nSp, 4, 3)),
#   # 全局层
#   com.size.x.precip  = rnorm(4, 0, 0.2),
#   com.INTRA.x.precip = rnorm(4, 0, 0.2), com.INTER.x.precip = rnorm(4, 0, 0.2),
#   com.INTRA.x.size   = rnorm(4, 0, 0.2), com.INTER.x.size   = rnorm(4, 0, 0.2),
#   com.INTRA.x.mic    = rnorm(4, 0, 0.2), com.INTER.x.mic    = rnorm(4, 0, 0.2),
#   # 年份层（sd 起始值给小的，免得一开始就把年份效应吃掉降水效应）
#   yr.int.sd = rexp(4, 4), yr.slope.size.sd = rexp(4, 4),
#   yr.slope.INTRA.sd = rexp(4, 4), yr.slope.INTER.sd = rexp(4, 4),
#   yrsp.int.sd = rexp(4, 4),
#   yr.int         = matrix(rnorm(nYr*4, 0, 0.2), nYr, 4),
#   yr.slope.size  = matrix(rnorm(nYr*4, 0, 0.2), nYr, 4),
#   yr.slope.INTRA = matrix(rnorm(nYr*4, 0, 0.2), nYr, 4),
#   yr.slope.INTER = matrix(rnorm(nYr*4, 0, 0.2), nYr, 4),
#   yrsp.int       = array(rnorm(nYr*nSp*4, 0, 0.2), c(nYr, nSp, 4)),
#   # 个体层
#   ind.sd  = rexp(4, 2),
#   ind.chCor = diag(4),              # 相关矩阵从单位矩阵（零相关）出发
#   ind.eff = matrix(rnorm(nInd*4, 0, 0.2), nInd, 4),
#   # 观测层
#   grow.sd = rexp(1, 1),
#   nb.size = rgamma(1, 4, 1),        # 避开 dexp(0.1) 抽到接近 0 的情形
#   # 补充子模型
#   pop.seed.mu  = rnorm(1, 0, 0.3), pop.seed.sd  = rexp(1, 2),
#   pop.veg.mu   = rnorm(1, 0, 0.3), pop.veg.sd   = rexp(1, 2),
#   pop.seed0.mu = rnorm(1, 0, 0.3), pop.seed0.sd = rexp(1, 2),
#   pop.seed.log  = rnorm(nSp, 0, 0.3),
#   pop.veg.log   = rnorm(nSp, 0, 0.3),
#   pop.seed0.log = rnorm(nSp, 0, 0.3),
#   pop.sizeRec = rnorm(nSp, 0, 0.3), sizeRec.sd = rexp(1, 1)
# )


# ═══════════════════════════════════════════════════════════════════════════
# 要跑的三个版本（review 建议，写进 Methods）
# ═══════════════════════════════════════════════════════════════════════════
# A. 完整版（上面这个）
# B. 去掉 yr.slope.INTRA / yr.slope.INTER，只留 yr.int
#    → 比较 com.*.x.precip 的后验缩水多少 = 这个交互有多少是靠先验撑着的
# C. yr.*.sd 的先验换成 dexp(1)
#    → 比较 com.*.x.precip 和 con.precip 的 CI 移动多少
# 三个结果一起报。审稿人一定会问 7 个年份凭什么估三阶交互。