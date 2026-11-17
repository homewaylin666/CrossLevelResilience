################################################################################
## Helper functions for the NIMBLE models of Part A
##
## Written by Hung-wei Lin Oct 2026


# MCMC chains, one after another or on several cores ---------------------------
  # One group of chains: builds and compiles the model once, then runs the chains of the group.
  # Self-contained (it also runs inside a parallel worker): loads nimble and this file itself.
run_chain_group <- function(chain_ids, code, constants, data, chain_inits, monitors, thin, n_iter, n_burn,
                            seed, fns_file, shifts = NULL) {
  suppressPackageStartupMessages(library(nimble))
  source(fns_file)
  m  <- nimbleModel(code, constants = constants, data = data, inits = chain_inits[[chain_ids[1]]],
                    calculate = FALSE)
  cm <- compileNimble(m)
  conf <- configureMCMC(m, monitors = monitors, thin = thin, print = FALSE)
    # Extra shift samplers (sampler_shift below), added to the default samplers
  for (s in shifts) conf$addSampler(target = s$target, type = "sampler_shift", control = list(weights = s$weights))
  mc <- buildMCMC(conf)
  cmc <- compileNimble(mc, project = m)
  out <- runMCMC(cmc, niter = n_iter, nburnin = n_burn, nchains = length(chain_ids),
                 inits = chain_inits[chain_ids], setSeed = seed + chain_ids, progressBar = TRUE,
                 samplesAsCodaMCMC = TRUE)
  if (length(chain_ids) == 1) coda::mcmc.list(out) else out
}

  # All chains of a model. n_cores = 1: in this R session, one chain after another, with NIMBLE's
  # progress bar in the console. n_cores > 1: the chains are split over that many parallel
  # workers (each builds and compiles its own copy of the model, so memory use grows with
  # n_cores); the progress bars of the workers go to 'log_file'. The inits of every chain are
  # drawn here, and each chain has its own seed (seed + chain number), so a rerun gives the same
  # samples for the same n_cores. 'shifts': node sets that also get a shift sampler (list(target, weights)).
run_chains <- function(code, constants, data, inits, monitors, thin, n_iter, n_burn, n_chains,
                       n_cores = 1, seed = 2026, log_file = NULL, shifts = NULL,
                       fns_file = here::here("code", "CLR_fns_nimble.R")) {
  set.seed(seed)
  chain_inits <- lapply(seq_len(n_chains), function(ch) inits())
  args <- list(code = code, constants = constants, data = data, chain_inits = chain_inits,
               monitors = monitors, thin = thin, n_iter = n_iter, n_burn = n_burn, seed = seed,
               fns_file = fns_file, shifts = shifts)
  n_cores <- min(n_cores, n_chains)
  if (n_cores <= 1) return(do.call(run_chain_group, c(list(chain_ids = seq_len(n_chains)), args), quote = TRUE))

  groups <- split(seq_len(n_chains), rep_len(seq_len(n_cores), n_chains))
  if (is.null(log_file)) log_file <- tempfile(fileext = ".log")
  if (file.exists(log_file)) file.remove(log_file)
  cat(sprintf("Running %d chains on %d cores; progress of the chains in '%s'.\n", n_chains, n_cores, log_file))
  cl <- parallel::makeCluster(n_cores, outfile = log_file)
  on.exit(parallel::stopCluster(cl), add = TRUE)
    # (argument names chosen to avoid partial matching with the arguments of parLapply; quote = TRUE
    # passes the model code as code, without evaluating it)
  worker <- function(ids, args_list, chain_fun) do.call(chain_fun, c(list(chain_ids = ids), args_list), quote = TRUE)
  environment(worker) <- globalenv()   # keeps the data of this call out of the copy sent to the workers
  res <- parallel::parLapply(cl, groups, worker, args_list = args, chain_fun = run_chain_group)
    # Chains back in their order (1, 2, ..., n_chains)
  chains <- unlist(lapply(res, as.list), recursive = FALSE)
  coda::mcmc.list(chains[order(unlist(groups))])
}


# Shift sampler ----------------------------------------------------------------
  # Moves a community mean and the yearly effects built on it together: target + c * weights. With
  # the year effects centred on the community mean (CLR_ModelPartA.R), this move leaves the
  # species-by-year effects unchanged, so the data stay out of it and the step is cheap. Single-node
  # samplers move along this direction very slowly (first full model: Rhat 1.25, growth intercept).
  # One-dimensional random-walk Metropolis step (a translation, so the proposal is symmetric); the step
  # size adapts towards 44% acceptance with shrinking adjustments, as NIMBLE's own RW sampler.
sampler_shift <- nimbleFunction(
  name = "sampler_shift",
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    weights <- control$weights
    calc_nodes <- model$getDependencies(target)
    scale <- 0.1
    n_tried <- 0
    n_accepted <- 0
    times_adapted <- 0
  },
  run = function() {
    step <- rnorm(1, 0, scale)
    lp_old <- model$getLogProb(calc_nodes)
    values(model, target) <<- values(model, target) + step * weights
    lp_new <- model$calculate(calc_nodes)
    if (decide(lp_new - lp_old)) {
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calc_nodes, logProb = TRUE)
      n_accepted <<- n_accepted + 1
    } else {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calc_nodes, logProb = TRUE)
    }
    n_tried <<- n_tried + 1
    if (n_tried == 200) {
      times_adapted <<- times_adapted + 1
      scale <<- scale * exp(10 / (times_adapted + 3)^0.8 * (n_accepted / 200 - 0.44))
      n_tried <<- 0
      n_accepted <<- 0
    }
  },
  methods = list(reset = function() {})
)
