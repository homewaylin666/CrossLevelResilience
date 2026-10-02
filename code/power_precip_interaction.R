# ═══════════════════════════════════════════════════════════════════════════
# 种内/种间交互 × macro 气候：能不能估得出来
#
# 问题的结构：precip 每年只有一个值，所以 com.INTRA.x.precip 实际上是
#   「每年一个 INTRA 斜率，对该年降水回归」——一条穿过 nYr 个点的线。
#   每年多少个体决定每个点有多准；年份数决定有几个点。
#
# 两个引擎：
#   两阶段（快）——逐年估斜率，再对 precip 回归。power loop 用这个。
#   glmer（慢）——联合模型，确认两阶段没低估，顺便看收敛状况。
# ═══════════════════════════════════════════════════════════════════════════

set.seed(1)

# ── 场景参数（照 Baldy 的实际规模）──────────────────────────────────────
PAR <- list(
  nYr      = 8,        # 2014-2021 = 7 个 transition，这里用 8 个 census 年
  nPerYr   = 1100,     # 8842 / 8
  nSp      = 17,
  b_size   = 0.50,     # 大小主效应
  b_intra  = -0.30,    # 种内主效应（★ gamma 要跟这个比才有意义）
  b_inter  = -0.20,    # 种间主效应
  sd_sp    = 0.40,     # 物种截距的 sd
  sd_yrInt = 0.25,     # yr.int 的 sd
  gamma    = 0.15,     # ★ 真值：每 1 SD 降水，INTRA 斜率移动多少
  sd_yrSlp = 0.15      # ★ 年份自身的斜率擺盪（气候解释不掉的那部分）
)

# ── 造一份资料 ──────────────────────────────────────────────────────────
sim_data <- function(P, family = c("binomial","gaussian")){
  family <- match.arg(family)
  precip <- scale(rnorm(P$nYr))[,1]                  # 每年一个值，已置中
  sp.int <- rnorm(P$nSp, 0, P$sd_sp)
  yr.int <- rnorm(P$nYr, 0, P$sd_yrInt)
  yr.bA  <- rnorm(P$nYr, 0, P$sd_yrSlp)              # 该年 INTRA 斜率的偏离
  yr.bE  <- rnorm(P$nYr, 0, P$sd_yrSlp)

  d <- do.call(rbind, lapply(seq_len(P$nYr), function(t){
    n  <- P$nPerYr
    data.frame(
      yr    = t,
      sp    = sample.int(P$nSp, n, replace = TRUE),
      z     = rnorm(n),
      INTRA = rnorm(n),
      INTER = rnorm(n),
      prec  = precip[t]
    )
  }))
  bA <- P$b_intra + P$gamma * d$prec + yr.bA[d$yr]   # 实现的 INTRA 斜率
  bE <- P$b_inter + P$gamma * d$prec + yr.bE[d$yr]
  eta <- sp.int[d$sp] + P$b_size * d$z + bA * d$INTRA + bE * d$INTER + yr.int[d$yr]
  d$y <- if (family == "binomial") rbinom(nrow(d), 1, plogis(eta)) else eta + rnorm(nrow(d))
  attr(d, "precip") <- precip
  d
}

# ── 引擎 1：两阶段 ──────────────────────────────────────────────────────
# 阶段一：每年单独估 INTRA 斜率。阶段二：那 nYr 个斜率对 precip 回归。
fit_twostage <- function(d, family = "binomial"){
  yrs <- sort(unique(d$yr))
  bhat <- vapply(yrs, function(t){
    s <- d[d$yr == t, ]
    m <- if (family == "binomial")
           glm(y ~ z + INTRA + INTER, data = s, family = binomial())
         else lm(y ~ z + INTRA + INTER, data = s)
    coef(m)[["INTRA"]]
  }, numeric(1))
  pr <- attr(d, "precip")[yrs]
  s2 <- summary(lm(bhat ~ pr))$coefficients
  c(est = s2["pr","Estimate"], se = s2["pr","Std. Error"],
    p   = s2["pr","Pr(>|t|)"])
}

# ── 引擎 2：联合模型（需要 lme4）────────────────────────────────────────
# ★ 8 个年份上放随机斜率，singular fit 很常见 —— 这本身就是结果之一。
fit_joint <- function(d, family = "binomial"){
  if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)
  f <- y ~ z + I(z^2) + INTRA * prec + INTER * prec +
           (1 | sp) + (1 + INTRA + INTER | yr)
  m <- suppressMessages(suppressWarnings(
    if (family == "binomial")
      lme4::glmer(f, data = d, family = binomial(),
                  control = lme4::glmerControl(optimizer = "bobyqa"))
    else lme4::lmer(f, data = d, REML = FALSE)
  ))
  cf <- summary(m)$coefficients
  list(est = cf["INTRA:prec","Estimate"], se = cf["INTRA:prec","Std. Error"],
       singular = lme4::isSingular(m))
}

# ── power loop ──────────────────────────────────────────────────────────
power_twostage <- function(P, R = 300, family = "binomial"){
  out <- replicate(R, { r <- fit_twostage(sim_data(P, family), family)
                        c(r[["est"]], r[["se"]]) })
  est <- out[1,]; se <- out[2,]
  c(power = mean(abs(est) > 1.96 * se), mean_est = mean(est), sd_est = sd(est))
}

# ═══════════════════════════════════════════════════════════════════════════
# A) 个体数 vs 年份数 —— 这是核心结论
# ═══════════════════════════════════════════════════════════════════════════
cat("\nA) 加个体 vs 加年份   (gamma =", PAR$gamma, ", 年擺盪 =", PAR$sd_yrSlp, ")\n")
cat(sprintf("%4s %8s %10s | %6s\n", "nYr", "n/yr", "total", "power"))
grid <- list(c(8,1100), c(8,11000), c(12,1100), c(16,1100), c(20,1100), c(30,1100))
for (g in grid){
  P <- modifyList(PAR, list(nYr = g[1], nPerYr = g[2]))
  r <- power_twostage(P, R = 300)
  cat(sprintf("%4d %8d %10d | %6.2f\n", g[1], g[2], g[1]*g[2], r[["power"]]))
}

# ═══════════════════════════════════════════════════════════════════════════
# B) 8 年的情况下，效应要多大才看得见
# ═══════════════════════════════════════════════════════════════════════════
cat("\nB) nYr = 8, n/yr = 1100。行 = gamma，列 = 年擺盪\n")
gammas <- c(0.05, 0.10, 0.20, 0.40); sds <- c(0.05, 0.15, 0.30)
cat(sprintf("%7s |", "gamma")); cat(sprintf(" sd=%4.2f", sds)); cat("\n")
for (gm in gammas){
  cat(sprintf("%7.2f |", gm))
  for (sd in sds){
    P <- modifyList(PAR, list(gamma = gm, sd_yrSlp = sd))
    cat(sprintf(" %7.2f", power_twostage(P, R = 300)[["power"]]))
  }
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# C) 联合模型对照：两阶段有没有低估，以及 8 个年份撑不撑得住随机斜率
# ═══════════════════════════════════════════════════════════════════════════
cat("\nC) 联合模型 (lme4) 对照，20 次重复\n")
if (requireNamespace("lme4", quietly = TRUE)){
  res <- t(replicate(20, {
    d  <- sim_data(PAR, "binomial")
    j  <- fit_joint(d, "binomial"); t2 <- fit_twostage(d, "binomial")
    c(joint = j$est, joint_se = j$se, sing = as.numeric(j$singular),
      two = t2[["est"]], two_se = t2[["se"]])
  }))
  cat(sprintf("  真值 gamma          = %.3f\n", PAR$gamma))
  cat(sprintf("  联合模型  平均估计  = %+.3f   显著比例 = %.2f\n",
      mean(res[,"joint"]), mean(abs(res[,"joint"]) > 1.96*res[,"joint_se"])))
  cat(sprintf("  两阶段    平均估计  = %+.3f   显著比例 = %.2f\n",
      mean(res[,"two"]),   mean(abs(res[,"two"])   > 1.96*res[,"two_se"])))
  cat(sprintf("  ★ singular fit 比例 = %.2f  （8 个年份放随机斜率的代价）\n",
      mean(res[,"sing"])))
} else cat("  lme4 没装，跳过。install.packages('lme4')\n")

# ═══════════════════════════════════════════════════════════════════════════
# 怎么用这份脚本做真正的决策
# ═══════════════════════════════════════════════════════════════════════════
# 1. 先跑真资料的 Part A，拿到 yr.slope.INTRA.sd 的后验中位数 —— 那就是 sd_yrSlp。
# 2. 拿 com.INTRA.x.precip 的后验中位数当 gamma。
# 3. 把这两个数填回 PAR，重跑 A 和 B。
#    → 得到的 power 就是「在我真实的效应量下，8 年到底够不够」。
# ★ 这是事后 power，不能拿来判断已有结果的显著性，只能拿来判断
#   「不显著」到底是真的没效应，还是年份不够。这两件事在写作上完全不同。
