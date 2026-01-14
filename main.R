### INITIAL CONFIGURATION ###

# Set up personal library path
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- path.expand("~/R/x86_64-pc-linux-gnu-library/4.5")
}

# Create personal library if it doesn't exist
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  cat("Created personal R library at:", user_lib, "\n")
}

# Add to library paths
.libPaths(c(user_lib, .libPaths()))

# Install and load packages
required_pkgs <- c("yaml", "parallelly", "this.path", "mvtnorm", 
                   "matrixcalc", "Matrix", "rstan", "loo")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

rstan_options(auto_write = TRUE)

# Cluster
options(mc.cores = parallelly::availableCores())

# Working environment
setwd(this.path::this.dir())
getwd()

config <- yaml.load_file("config.yaml") # Production
# config <- yaml.load_file("config_test.yaml") # Local testing

source("src/dgp_sem.R")
source("src/fitindices_utils.R")

### GET SETTINGS FROM CONFIG FILE ###

seed <- config$seed
N.all <- config$settings$N
Nt.all <- config$settings$N_t
size_miss <- config$settings$M
T_event <- config$settings$T_event
R <- config$settings$R  # replications per condition

# MCMC
iter_n <- if(!is.null(config$mcmc$iter)) config$mcmc$iter else 2000
chains_n <- if(!is.null(config$mcmc$chains)) config$mcmc$chains else 4
warmup_n <- if(!is.null(config$mcmc$warmup)) config$mcmc$warmup else 1000
cores_n <- if(!is.null(config$mcmc$cores)) config$mcmc$cores else 1

##########################
# Model parameters
##########################

phi0 <- diag(2)*.7+.3 # cov(eta)
mu0 <- c(0,0) # mean(eta)
ar0 <- c(.3,.3) # ar(1) structure
ly00 <- .6

ly0 <- matrix(c(ly00,ly00,ly00,0,0,0,
                0,0,0,ly00,ly00,ly00),6,2,byrow=F) # factor loadings
td <- diag(6)*(1-ly00^2) # cond. var(y) -> res var

##########################
# Load Stan models
##########################

rt1 <- stanc("models/DSEM_VAR1_RI.stan")
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)
params2 <- c("ly","sigmaeps","sigmaeta1","sigmaeta2","ty","ka","beta1")

rt1.null <- stanc("models/DSEM_NULL.stan")
sm1.null <- stan_model(stanc_ret = rt1.null, verbose=FALSE)

##########################
# Helper function
##########################

calculate_T_event <- function(N_t, proportion) {
  t0 <- 1 + floor((N_t - 1) * proportion)
  t0 <- max(2, min(N_t, t0))
  t0
}

##########################
# Generate conditions grid
##########################

conditions_grid <- expand.grid(
  N = N.all,
  Nt = Nt.all,
  M = size_miss,
  T_event_prop = T_event,
  stringsAsFactors = FALSE
)

##########################
# SLURM task ID
##########################

rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(rslurm_id)) rslurm_id <- 1 # local testing

# Validate task ID
if (rslurm_id > nrow(conditions_grid)) {
  stop(paste("SLURM_ARRAY_TASK_ID", rslurm_id, "exceeds number of conditions", nrow(conditions_grid)))
}

# Get condition for this task
cond <- conditions_grid[rslurm_id, ]
N <- cond$N
Nt <- cond$Nt
M <- cond$M
T_event_prop <- cond$T_event_prop
t0 <- calculate_T_event(Nt, T_event_prop)

# Create output directory
dir.create("results", showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("Condition ID:", rslurm_id, "\n")
cat("N =", N, "| Nt =", Nt, "| M =", M, "| T_event =", T_event_prop, "| t0 =", t0, "\n")
cat("Running", R, "replications...\n")
cat("========================================\n")

##########################
# Main loop: R replications for this condition
##########################

for (rep in 1:R) {
  
  # Set unique seed for this replication
  set.seed(seed + rslurm_id * 100000 + rep)
  
  # Generate data
  ly1 <- matrix(0, 6, 2)
  ly1[4, 1] <- M
  
  dat1 <- gendata02(N, Nt, phi0, mu0, ar0, ly0, ly1, td, t0 = t0)
  y <- dat1[["y"]]
  y.2 <- dat1[["y0"]]
  
  data.stan1 <- list(N = N, Nt = Nt, y = y)
  
  # Fit models
  fit1 <- sampling(sm1, data = data.stan1,
                   iter = iter_n,
                   chains = chains_n,
                   warmup = warmup_n,
                   cores = cores_n,
                   refresh = 0)
  
  fit1.null <- sampling(sm1.null, data = data.stan1,
                        iter = iter_n,
                        chains = chains_n,
                        warmup = warmup_n,
                        cores = cores_n,
                        refresh = 0)
  
  # Extract convergence
  fit_summary <- summary(fit1)
  param_names <- rownames(fit_summary$summary)
  beta_indices <- c()
  for (t1 in 1:length(params2)) {
    beta_indices <- c(beta_indices, grep(params2[t1], param_names))
  }
  conv <- max(fit_summary$summary[beta_indices, "Rhat"])
  
  # Compute fit indices
  modelfit <- getfit(fit1, y, y.2, nvar = 6, Nt = Nt, cfi = TRUE, fit1.null = fit1.null)
  modelfit[[3]] <- conv
  
  # Save results
  nomfile <- paste0(
    "results/N", N,
    "_Nt", Nt,
    "_M", M,
    "_t0", t0,
    "_prop", T_event_prop,
    "_rep", rep,
    ".RDS"
  )
  saveRDS(modelfit, file = nomfile)
  
  cat("Replication", rep, "/", R, "completed.\n")
  
} # end replications loop

cat("========================================\n")
cat("Condition", rslurm_id, "finished successfully.\n")
cat("========================================\n")