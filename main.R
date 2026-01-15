### INITIAL CONFIGURATION ###

# 1. Setup Library Path (Mantengo esto porque es crítico para tu cluster)
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") user_lib <- path.expand("~/R/x86_64-pc-linux-gnu-library/4.5")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

# 2. Load Packages (Loop simplificado)
required_pkgs <- c("yaml","parallelly","this.path","mvtnorm","matrixcalc",
                   "Matrix","rstan","loo","DBI","RSQLite")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# 3. Environment & Config
rstan_options(auto_write = TRUE)
options(mc.cores = parallelly::availableCores())
setwd(this.path::this.dir())

config <- yaml.load_file("config.yaml") 
source("src/dgp_sem.R")
source("src/fitindices_utils.R")

# 4. Global Constants
SEED_BASE <- config$seed
R_REPS    <- config$settings$R

# 5. Define Grid Directly (Robustez aquí: unlist + as.numeric previene errores de yaml)
conditions_grid <- expand.grid(
  N            = as.numeric(unlist(config$settings$N)),
  Nt           = as.numeric(unlist(config$settings$N_t)), # Unificamos nombre a 'Nt'
  M            = as.numeric(unlist(config$settings$M)),
  T_event_prop = as.numeric(unlist(config$settings$T_event)),
  stringsAsFactors = FALSE
)

# 6. SLURM Task Identification
rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(rslurm_id)) rslurm_id <- 1

# Validation
if (rslurm_id > nrow(conditions_grid)) {
  stop(paste("Task ID", rslurm_id, "exceeds grid size", nrow(conditions_grid)))
}

# Extract Condition Parameters (Directamente del grid)
current_cond <- conditions_grid[rslurm_id, ]
N            <- current_cond$N
Nt           <- current_cond$Nt
M            <- current_cond$M
T_event_prop <- current_cond$T_event_prop

# Calculate Derived Parameter
t0 <- max(2, min(Nt, 1 + floor((Nt - 1) * T_event_prop)))

# 7. Prepare Output Directory & DB
run_id <- Sys.getenv("RUN_ID")
if (run_id == "") run_id <- format(Sys.time(), "%Y%m%d_%H%M%S") # Fallback local
job_id <- as.integer(Sys.getenv("SLURM_JOB_ID")); if(is.na(job_id)) job_id <- -1

run_dir <- file.path(config$output$results_dir, paste0("run_", run_id))
if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

# Database Connection (Con retry pause para evitar colisiones en cluster)
Sys.sleep(runif(1, 0.1, 2.0)) 
db_file <- file.path(run_dir, "fit_indices.sqlite")
con <- DBI::dbConnect(RSQLite::SQLite(), db_file)

# Optimize SQLite for Cluster
DBI::dbExecute(con, "PRAGMA journal_mode = WAL;")
DBI::dbExecute(con, "PRAGMA synchronous = NORMAL;")
DBI::dbExecute(con, "PRAGMA busy_timeout = 60000;")

# Create Tables (Solo si no existen)
DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS fit_overall (
    run_id TEXT, job_id INTEGER, array_id INTEGER, rep INTEGER,
    N INTEGER, Nt INTEGER, M REAL, T_event_prop REAL, t0 INTEGER,
    conv_max_rhat REAL,
    BRMSEA REAL, BGamma REAL, BadjGamma REAL, BMc REAL, BCFI REAL, BTLI REAL,
    created_at TEXT
);")

DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS fit_timepoint (
    run_id TEXT, job_id INTEGER, array_id INTEGER, rep INTEGER, time INTEGER,
    N INTEGER, Nt INTEGER, M REAL, T_event_prop REAL, t0 INTEGER,
    BRMSEA REAL, BGamma REAL, BadjGamma REAL, BMc REAL, BCFI REAL, BTLI REAL
);")

# 8. Compile Models (Fuera del loop)
sm1      <- stan_model(file = config$stan$model_file, verbose = FALSE)
sm1.null <- stan_model(file = config$stan$null_model_file, verbose = FALSE)

# Model Constants
phi0 <- diag(2)*.7+.3
mu0  <- c(0,0)
ar0  <- c(.3,.3)
ly00 <- .6
ly0  <- matrix(c(ly00,ly00,ly00,0,0,0, 0,0,0,ly00,ly00,ly00), 6, 2)
td   <- diag(6)*(1-ly00^2)

cat(sprintf("=== STARTING: ID %d | N=%d, Nt=%d, M=%.2f, T_event_prop=%.2f ===\n", rslurm_id, N, Nt, M, T_event_prop))

# 9. Main Replication Loop
db_append_retry <- function(con, table, df) {
  for (k in 1:20) {
    if (!inherits(try(DBI::dbAppendTable(con, table, df), silent=TRUE), "try-error")) return(TRUE)
    Sys.sleep(0.5 * k)
  }
  stop("DB locked")
}

for (rep in 1:R_REPS) {
  set.seed(SEED_BASE + rslurm_id * 100000 + rep)
  
  # Data Generation
  ly1 <- matrix(0, 6, 2); ly1[4, 1] <- M
  dat1 <- gendata02(N, Nt, phi0, mu0, ar0, ly0, ly1, td, t0 = t0)
  
  # Stan Data List (Aseguramos nombres exactos para Stan)
  stan_data <- list(N = N, Nt = Nt, y = dat1[["y"]])
  
  # Sampling
  mcmc_args <- list(data = stan_data, iter = config$mcmc$iter, 
                    chains = config$mcmc$chains, warmup = config$mcmc$warmup, 
                    cores = config$mcmc$cores, refresh = 0)
  
  fit1      <- do.call(sampling, c(list(object = sm1), mcmc_args))
  fit1.null <- do.call(sampling, c(list(object = sm1.null), mcmc_args))
  
  # Convergence Check
  sry <- summary(fit1)$summary
  # Buscamos Rhats solo de parámetros relevantes (evitando greps complejos si posible)
  # Usamos una regex simple que capture los params principales
  params_regex <- "ly|sigmaeps|sigmaeta|ty|ka|beta1"
  rhats <- sry[grep(params_regex, rownames(sry)), "Rhat"]
  conv_stat <- if(length(rhats) > 0) max(rhats, na.rm=TRUE) else NA
  
  # Fit Indices
  fits <- getfit(fit1, dat1[["y"]], dat1[["y0"]], nvar=6, Nt=Nt, cfi=TRUE, fit1.null=fit1.null)
  
  # Prepare Data for DB
  # Overall
  overall_df <- data.frame(
    run_id = run_id, job_id = job_id, array_id = rslurm_id, rep = rep,
    N = N, Nt = Nt, M = M, T_event_prop = T_event_prop, t0 = t0,
    conv_max_rhat = conv_stat,
    BRMSEA = fits[[1]][1], BGamma = fits[[1]][2], BadjGamma = fits[[1]][3],
    BMc = fits[[1]][4], BCFI = fits[[1]][5], BTLI = fits[[1]][6],
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  # Timepoint
  tp_df <- as.data.frame(fits[[2]])
  colnames(tp_df) <- c("BRMSEA","BGamma","BadjGamma","BMc","BCFI","BTLI")
  tp_df$time <- 1:Nt
  tp_df$run_id <- run_id; tp_df$job_id <- job_id; tp_df$array_id <- rslurm_id; 
  tp_df$rep <- rep; tp_df$N <- N; tp_df$Nt <- Nt; tp_df$M <- M; 
  tp_df$T_event_prop <- T_event_prop; tp_df$t0 <- t0
  
  # Write to DB
  DBI::dbWithTransaction(con, {
    db_append_retry(con, "fit_overall", overall_df)
    db_append_retry(con, "fit_timepoint", tp_df)
  })
  
  cat(sprintf("Rep %d/%d OK.\n", rep, R_REPS))
}

DBI::dbDisconnect(con)
cat("=== COMPLETED ===\n")
