################################################################################
## Helper functions for the NIMBLE models of Part A
##
## Written by Hung-wei Lin Oct 2026


# MCMC chains, one after another or on several cores ---------------------------
  # One group of chains: builds and compiles the model once, then runs the chains of the group.
  # Self-contained (it also runs inside a parallel worker): loads nimble and this file itself.
run_chain_group <- function(chain_ids, code, constants, data, chain_inits, monitors, thin, n_iter, n_burn,
                            seed, fns_file) {
  suppressPackageStartupMessages(library(nimble))
  source(fns_file)
  m  <- nimbleModel(code, constants = constants, data = data, inits = chain_inits[[chain_ids[1]]],
                    calculate = FALSE)
  cm <- compileNimble(m)
  mc <- buildMCMC(configureMCMC(m, monitors = monitors, thin = thin, print = FALSE))
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
  # samples for the same n_cores.
run_chains <- function(code, constants, data, inits, monitors, thin, n_iter, n_burn, n_chains,
                       n_cores = 1, seed = 2026, log_file = NULL,
                       fns_file = here::here("code", "CLR_fns_nimble.R")) {
  set.seed(seed)
  chain_inits <- lapply(seq_len(n_chains), function(ch) inits())
  args <- list(code = code, constants = constants, data = data, chain_inits = chain_inits,
               monitors = monitors, thin = thin, n_iter = n_iter, n_burn = n_burn, seed = seed,
               fns_file = fns_file)
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
