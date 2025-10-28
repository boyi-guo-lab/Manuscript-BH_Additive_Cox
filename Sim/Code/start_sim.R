library(dplyr)


# Project root  -----------------------------------------------
PROJ_ROOT <- normalizePath("~/Documents/SSL", mustWork = TRUE)
JOB_SCRIPT <- file.path(PROJ_ROOT, "Manuscript-BH_Additive_Cox","Sim", "Code", "bcam_sim.job")

# Simulation Parameters ---------------------------------------------------
sim_prmt <- expand.grid(
  n_train = c(200),
  p = c(4,10), # c(4, 10, 50, 100, 200),        # Number of Predictors
  rho = c(0.5), # c(0, 0.5),                   # X Cov Structure AR(rho)
  pi_cns = c(0.15) #c(0.15, 0.3, 0.45)        # Proportional of Censoring
)

# ---- Cluster & job defaults for smoke test ----
ACCOUNT   <- "guo"
PARTITION <- "kingspeak-shared"
ARRAY     <- "1-10"         # 2 reps; bump later for production
TIME      <- "00:10:00"
MEM       <- "4G"
CPUS      <- "1"



start.sim <- function(n_train, p, rho, pi_cns) {
  # Unique job name per scenario
  job.name <- paste0("bcam_sim_p=", p, ",",
                     "rho=", rho, ",",
                     "pi_cns=", pi_cns)

  # R output repos
  res_dir   <- file.path(PROJ_ROOT,  "Manuscript-BH_Additive_Cox", "Sim", "Res",   job.name)
  scale_dir <- file.path(PROJ_ROOT,  "Manuscript-BH_Additive_Cox", "Sim", "scale", job.name)
  log_dir   <- file.path(PROJ_ROOT,  "Manuscript-BH_Additive_Cox", "Sim", "Log",   job.name)

  dir.create(res_dir,   recursive = TRUE, showWarnings = FALSE)
  dir.create(scale_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir,   recursive = TRUE, showWarnings = FALSE)

  # Slurm output repos
  slurm_out_dir <- file.path(PROJ_ROOT, "Manuscript-BH_Additive_Cox", "Sim", "Log", "slurm_out")
  dir.create(slurm_out_dir, recursive = TRUE, showWarnings = FALSE)
  slurm_out_path <- file.path(slurm_out_dir, "%x-%A_%a.out")

  # sbatch command ( pass parameters + paths via --export)
  cmd <- paste(
    "sbatch",
    paste0("--account=", ACCOUNT),
    paste0("--partition=", PARTITION),
    paste0("--array=",     ARRAY),
    paste0("--time=",      TIME),
    paste0("--mem=",       MEM),
    paste0("--cpus-per-task=", CPUS),
    paste0("--job-name=",  job.name),
    paste0("--output=",  slurm_out_path),
    paste0("--error=",  slurm_out_path),
    paste0("--export=ALL,",
           "n_train=", n_train, ",",
           "p=",       p,       ",",
           "rho=",     rho,     ",",
           "pi_cns=",  pi_cns,  ",",
           "resPath=",  res_dir,   ",",
           "scalePath=",scale_dir, ",",
           "logPath=",  log_dir),
    shQuote(JOB_SCRIPT)
  )

  message("Submitting: ", cmd)
  system(cmd)
}

# Submit all combinations (2 combos × array 1-2 = 4 tiny tasks)
apply(sim_prmt, 1, function(row) do.call(start.sim, as.list(row)))


