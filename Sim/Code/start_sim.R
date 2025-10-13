library(dplyr)


# Project root with symlink -----------------------------------------------
PROJ_ROOT <- normalizePath("~/Manuscript-BH_Additive_Cox", mustWork = TRUE)
JOB_SCRIPT <- file.path(PROJ_ROOT, "Sim", "Code", "bcam_sim.job")

# Simulation Parameters ---------------------------------------------------
sim_prmt <- expand.grid(
  n_train = c(200),
  p = c(4,10), # c(4, 10, 50, 100, 200),        # Number of Predictors
  rho = c(0.5), # c(0, 0.5),                   # X Cov Structure AR(rho)
  pi_cns = c(0.15) #c(0.15, 0.3, 0.45)        # Proportional of Censoring
)

# ---- Cluster & job defaults for smoke test ----
ACCOUNT   <- "guo"
PARTITION <- "lonepeak-shared"
ARRAY     <- "1-2"         # 2 reps; bump later for production
TIME      <- "00:15:00"
MEM       <- "4G"
CPUS      <- "1"
MAIL_TYPE <- "END,FAIL"
MAIL_USER <- "sophhuebler@gmail.com"

start.sim <- function(n_train, p, rho, pi_cns) {
  # Unique job name per scenario
  job.name <- paste0("bcam_sim_p=", p, ",",
                     "rho=", rho, ",",
                     "pi_cns=", pi_cns)

  # Repo-based output directories (visible in your project)
  res_dir   <- file.path(PROJ_ROOT, "Sim", "Res",   job.name)
  scale_dir <- file.path(PROJ_ROOT, "Sim", "scale", job.name)
  log_dir   <- file.path(PROJ_ROOT, "Sim", "Log",   job.name)

  # Create folders (job script will also mkdir -p)
  dir.create(res_dir,   recursive = TRUE, showWarnings = FALSE)
  dir.create(scale_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir,   recursive = TRUE, showWarnings = FALSE)

  # sbatch command (we pass parameters + paths via --export)
  cmd <- paste(
    "sbatch",
    paste0("--account=",   ACCOUNT),
    paste0("--partition=", PARTITION),
    paste0("--array=",     ARRAY),
    paste0("--time=",      TIME),
    paste0("--mem=",       MEM),
    paste0("--cpus-per-task=", CPUS),
    paste0("--mail-type=", MAIL_TYPE),
    paste0("--mail-user=", MAIL_USER),
    paste0("--job-name=",  job.name),
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

# Helper Function for Setting Up Job --------------------------------------
#start.sim <- function(
#  n_train,
#  p, rho,
#  pi_cns
#) {
#  #Compose job name
#  job.name <- paste0("bcam_sim_p=", p, ",",
#                     "rho=", rho, ",",
#                     "pi_cns=", pi_cns)
#
#  # NOTE:
#  ## Job name has to be unique for each of your simulation settings
#  ## DO NOT USE GENERIC JOB NAME FOR CONVENIENCE
#  job.flag <- paste0("--job-name=",job.name)
#
#  err.flag <- paste0("--error=",job.name,".err")
#
#  out.flag <- paste0("--output=",job.name,".out")
#
#  # Pass simulation parameters to jobs using export flag
#  arg.flag <- paste0("--export=n_train=", n_train, ",",
#                     "p=", p, ",",
#                     "rho=", rho, ",",
#                     "pi_cns=", pi_cns)
#
#  # Create Jobs
#  system(
#    paste("sbatch", job.flag, err.flag, out.flag, arg.flag,
#          "~/Manuscript-BH_Additive_Cox/Sim/Code/bcam_sim.job")
#  )
#}





# Set up job for all simulation settings ----------------------------------
for(i in 1:NROW(sim_prmt)){
  do.call(start.sim, sim_prmt[i,])
}
