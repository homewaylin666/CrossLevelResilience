################################################################################
## Model selection of Part A (Idaho): fit the glmmTMB models of one phase
##
## Written by Hung-wei Lin Oct 2026
##
## Phase 'select_phase' (set in the runscript): 1 climate variable (276 fits), 2 interaction terms
## (66048 fits, run on ARC: code/arc/CLR_select_array.sh), 3 model_after (192 fits).
## One small summary per fit in 'data/data.derived/Idaho/selection'; fits already saved are skipped.
## Functions in CLR_fns_select.R.


# Settings ---------------------------------------------------------------------
n_cores <- 4   # fits run in parallel (1-2 GB of memory each)
trait_file <- here("data", "data.derived", site_tag, "CLR_PartA_traits.csv")   # phase 3 (columns code, trait)

# Set-up: p of C' and the random effects left out -----------------------------
dir.create(sel_dir(), showWarnings = FALSE)
if (!file.exists(file.path(sel_dir(), "CLR_select_setup.rds"))) select_setup()
setup <- readRDS(file.path(sel_dir(), "CLR_select_setup.rds"))
cat("Random effects at zero in the full glmmTMB models (left out):\n"); print(setup$re_drop)

# Fit --------------------------------------------------------------------------
tasks <- make_tasks(select_phase)
todo <- which(!file.exists(file.path(sel_dir(), tasks$file)))
cat(sprintf("Phase %d: %d fits, %d still to run.\n", select_phase, nrow(tasks), length(todo)))
  # Each core is a fresh R session (forked copies of an RStudio session can fail on a Mac)
t0 <- Sys.time()
cl <- parallel::makeCluster(n_cores)
parallel::clusterExport(cl, c("site_tag", "select_phase", "trait_file"))
parallel::clusterCall(cl, setwd, here())
invisible(parallel::clusterEvalQ(cl, {library(here); library(glmmTMB); source(here("code", "CLR_fns_select.R"))}))
msg <- unlist(parallel::parLapplyLB(cl, todo, function(i)
  tryCatch({run_task(select_phase, i, trait_file); "ok"}, error = function(e) conditionMessage(e))))
parallel::stopCluster(cl)
if (any(msg != "ok")) { cat("Errors of the failed fits:\n"); print(table(msg[msg != "ok"])) }
cat(sprintf("Phase %d: %.0f min, %d fits still missing.\n", select_phase,
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            sum(!file.exists(file.path(sel_dir(), tasks$file)))))
rm(setup, tasks, todo, t0, cl, msg)
