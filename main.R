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
required_pkgs <- c(
  "yaml","parallelly","this.path","mvtnorm","matrixcalc","Matrix",
  "rstan","loo",
  "DBI","RSQLite"
)

for (pkg in required_pkgs) {
  library(pkg, character.only = TRUE)
}

rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1)

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

seed_base <- if(!is.null(config$seed)) config$seed else 42
overall_fits <- if(!is.null(config$overall_fits)) config$overall_fits else FALSE
N.all <- unlist(config$settings$N)
Nt.all <- unlist(config$settings$N_t)
M_values <- unlist(config$settings$M)
T_event_proportions <- unlist(config$settings$T_event)
R <- config$settings$R  # Total replications per condition

# --- CHUNKING CONFIGURATION ---
# Default to 5 reps per chunk if not specified in config
reps_per_chunk <- if(!is.null(config$mcmc$reps_per_chunk)) config$mcmc$reps_per_chunk else 5

# MCMC
iter_n <- if(!is.null(config$mcmc$iter)) config$mcmc$iter else 2000
chains_n <- if(!is.null(config$mcmc$chains)) config$mcmc$chains else 4
warmup_n <- if(!is.null(config$mcmc$warmup)) config$mcmc$warmup else 1000
cores_n <- if(!is.null(config$mcmc$cores)) config$mcmc$cores else 1

##############################
# Helpers
##############################

# Helper for appending to DB with retries
db_append_retry <- function(con, table, df, retries = 40, sleep_sec = 0.2) {
  for (k in seq_len(retries)) {
    ok <- try({
      DBI::dbAppendTable(con, table, df)
      TRUE
    }, silent = TRUE)
    if (!inherits(ok, "try-error")) return(invisible(TRUE))
    Sys.sleep(sleep_sec * k)
  }
  stop("DB write failed repeatedly (database locked?)")
}

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
params2 <- c("ly","sigmaeps","sigmaeta1","sigmaeta2","ty","ka","beta1","phi1")

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
# Generate conditions grid (WITH CHUNKING)
##########################

# 1. Base grid (Conditions only)
base_grid <- expand.grid(
  N = N.all,
  Nt = Nt.all,
  M = M_values,
  T_event_prop = T_event_proportions,
  stringsAsFactors = FALSE
)

# Avoid bug with list and SQLite
for (col in names(base_grid)) {
  if (is.list(base_grid[[col]])) {
    base_grid[[col]] <- unlist(base_grid[[col]])
  }
}

# 2. Generate Chunk Map
num_chunks <- ceiling(R / reps_per_chunk)
chunks_df <- data.frame(
  chunk_id = 1:num_chunks,
  start_rep = seq(1, R, by = reps_per_chunk)
)
chunks_df$end_rep <- pmin(chunks_df$start_rep + reps_per_chunk - 1, R)

# 3. Expand Grid (Cartesian Join: Conditions x Chunks)
# Each condition is now split into multiple rows
conditions_grid <- merge(base_grid, chunks_df, by = NULL) 

cat(sprintf("Total Conditions: %d\n", nrow(base_grid)))
cat(sprintf("Replications per Chunk: %d\n", reps_per_chunk))
cat(sprintf("Total Atomic Tasks: %d\n", nrow(conditions_grid)))


##########################
# Cluster Detection & Workload Manager
##########################

is_cluster <- (!is.na(Sys.getenv("SLURM_JOB_ID", NA)) || 
               !is.na(Sys.getenv("SLURM_ARRAY_TASK_ID", NA)))

# Recuperar RUN_ID (vendrá de la variable de entorno, sea nueva o manual)
run_id <- Sys.getenv("RUN_ID", "local_run")

if (is_cluster) {
  cat("Running on HPC cluster (SLURM detected)\n")
  rslurm_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
  job_id <- as.integer(Sys.getenv("SLURM_ARRAY_JOB_ID", "1"))
} else {
  cat("Running locally (no SLURM environment)\n")
  rslurm_id <- 1
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  job_id <- 1
}

cat(sprintf("Job ID: %d | Run ID: %s\n", rslurm_id, run_id))

##########################
# Configuration
##########################

# Routes
base_dir <- config$output$results_dir
run_dir  <- file.path(base_dir, paste0("run_", run_id))

# LÓGICA DE REANUDACIÓN ROBUSTA
db_file <- file.path(run_dir, "fit_indices.sqlite")
config$db_path <- db_file

# Si la carpeta no existe, la creamos (Nuevo Run)
if (!dir.exists(run_dir)) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created NEW run directory:", run_dir, "\n")
  init_db_needed <- TRUE
} else {
  cat(">>> RESUMING from EXISTING directory:", run_dir, "\n")
  # Si reanudamos, NO queremos reinicializar la tabla tasks si ya tiene datos
  init_db_needed <- FALSE 
  
  if (!file.exists(db_file)) {
      cat("Warning: Directory exists but DB missing. Initializing as fresh.\n")
      init_db_needed <- TRUE
  }
}

Sys.sleep(runif(1, min=0.01, max=0.2))

con <- DBI::dbConnect(RSQLite::SQLite(), db_file)

DBI::dbExecute(con, "PRAGMA journal_mode = WAL;")
DBI::dbExecute(con, "PRAGMA synchronous = NORMAL;")
DBI::dbExecute(con, "PRAGMA busy_timeout = 60000;")

#Overall results table
DBI::dbExecute(con, "
  CREATE TABLE IF NOT EXISTS fit_overall (
    run_id TEXT,
    job_id INTEGER,
    array_id INTEGER,
    rep INTEGER,
    chunk_id INTEGER, -- Added for traceability
    N INTEGER,
    Nt INTEGER,
    M REAL,
    T_event_prop REAL,
    t0 INTEGER,
    id0 INTEGER,
    -- DIAGNOSTICS --
    min_n_eff REAL,       
    n_divergent INTEGER,  
    exec_time_sec REAL,   
    conv_max_rhat REAL,
    -- PARAMETER RECOVERY --
    est_ar1 REAL,      
    est_ar2 REAL,      
    est_ri_var1 REAL,  
    est_ri_var2 REAL,  
    
    -- PARAMETERS FOR BIAS CHECK: Factor loadings and correlation --
    est_ly1 REAL,      
    est_ly2 REAL,      
    est_ly3 REAL,      
    est_ly4 REAL,      
    est_corr_factor REAL, 

    BRMSEA REAL,
    BGamma REAL,
    BadjGamma REAL,
    BMc REAL,
    BCFI REAL,
    BTLI REAL,
    created_at TEXT
  );
")

DBI::dbExecute(con, "
  CREATE TABLE IF NOT EXISTS fit_timepoint (
    run_id TEXT,
    job_id INTEGER,
    array_id INTEGER,
    rep INTEGER,
    time INTEGER,
    N INTEGER,
    Nt INTEGER,
    M REAL,
    T_event_prop REAL,
    t0 INTEGER,
    id0 INTEGER,
    conv_max_rhat REAL,
    BRMSEA REAL,
    BGamma REAL,
    BadjGamma REAL,
    BMc REAL,
    BCFI REAL,
    BTLI REAL
  );
")

##########################
# Task farm
##########################

# Only load workload manager on cluster to avoid local dependencies
if (is_cluster) {
  source("src/workload_manager.R")
  cat("Loaded workload_manager.R\n")
  init_tasks_db(config$db_path, conditions_grid)
}

############################
# Helper functions
############################

calculate_seed <- function(cond_row, seed_base, rep) {
  N <- as.numeric(cond_row$N)
  Nt <- as.numeric(cond_row$Nt)
  M <- as.numeric(cond_row$M)
  T_prop <- as.numeric(cond_row$T_event_prop)
  
  condition_hash <- (
    (N * 1) + 
    (Nt * 1000) + 
    (round(M * 100) * 100000) + 
    (round(T_prop * 100) * 10000000)
  ) %% 100000
  
  # Crucial: rep identifies unique stream. Rep must be global index.
  val <- (seed_base + condition_hash * 100000 + rep)
  return(as.integer(val %% 2147483647))
}

clean_val <- function(x) {
  if (is.null(x) || length(x) == 0 || is.nan(x)) return(NA_real_)
  val <- as.numeric(unlist(x))[1]
  return(val)
}

############################
# While Loop
############################
keep_working <- TRUE
while(keep_working) {

  ############################
  # Get bundle of conditions
  ############################
  if (is_cluster) {
    bundle <- get_bundle(
      db_path = config$db_path,
      job_id = rslurm_id,
      batch_size = 1, # Fetch 1 chunk at a time (Atomic) to prevent timeout
      n_jobs = 900
    )
    
    if (bundle$status == "empty") {
      cat("✓ No pending tasks (chunks). Worker finished.\n")
      keep_working <- FALSE
      break
    }
    
    conditions_to_run <- bundle$conditions
    cat(sprintf(">>> Starting bundle: %s (%d chunks)\n", bundle$bundle_id, nrow(conditions_to_run)))
    
  } else {
    # LOCAL MODE: Run 1 chunk for testing
    if (rslurm_id > nrow(conditions_grid)) stop("rslurm_id out of bounds")
    conditions_to_run <- conditions_grid[rslurm_id, , drop = FALSE]
    keep_working <- FALSE
  }

  ############################
  # Process conditions in the current bundle
  ############################
  
  for (cond_idx in seq_len(nrow(conditions_to_run))) {
    
    cond <- conditions_to_run[cond_idx, ]
    
    # Extract params
    N <- as.numeric(cond$N)
    Nt <- as.numeric(cond$Nt)
    M <- as.numeric(cond$M)
    T_event_prop <- as.numeric(cond$T_event_prop)
    
    # --- GLOBAL REP RANGE ---
    start_rep <- as.integer(cond$start_rep)
    end_rep   <- as.integer(cond$end_rep)
    chunk_id  <- as.integer(cond$chunk_id)
    
    # Task ID from DB
    condition_id <- if("condition_id" %in% names(cond)) cond$condition_id else NA
    if (!is.na(condition_id) && length(condition_id) > 0) {
      condition_id <- as.integer(unlist(condition_id))[1]
    }

    t0 <- calculate_T_event(Nt, T_event_prop)
    
    cat(sprintf("\n>>> Processing Task (Chunk) %d: N=%d, Nt=%d, Reps %d-%d\n", 
                condition_id, N, Nt, start_rep, end_rep))
    
    # Inner loop: iterate over the SPECIFIC RANGE for this chunk
    # This replaces the old 1:R loop
    for (rep in start_rep:end_rep) {
        
        tryCatch({
            # Ensure we use GLOBAL rep for seed
            current_seed <- calculate_seed(cond, seed_base, rep)
            set.seed(current_seed)
          
            cat("  -> Rep", rep, "of global", R, "...\n")
            
            # Generate data
            ly1 <- matrix(0, 6, 2)
            ly1[4, 1] <- M
            
            dat1 <- gendata02(N, Nt, phi0, mu0, ar0, ly0, ly1, td, t0 = t0)
            y <- dat1[["y"]]
            y.2 <- dat1[["y0"]]
            
            data.stan1 <- list(N = N, Nt = Nt, y = y)
            
            stan_control_list <- list(adapt_delta = 0.95, max_treedepth = 12)
            
            # Fit models
            fit1 <- sampling(sm1, data = data.stan1,
                            iter = iter_n,
                            chains = chains_n,
                            warmup = warmup_n,
                            cores = cores_n,
                            refresh = 0,
                            control = stan_control_list,
                            pars = c("ly","ty","sigmaeps","ka","beta1","sigmaeta1","sigmaeta2","phi1"),
                            include = TRUE)
            
            fit1.null <- sampling(sm1.null, data = data.stan1,
                                  iter = iter_n,
                                  chains = chains_n,
                                  warmup = warmup_n,
                                  cores = cores_n,
                                  refresh = 0,
                                  control = stan_control_list,
                                  pars = c("ty","sigmaeps"),
                                  include = TRUE)
            
            # Extract convergence
            fit_summary <- summary(fit1)

            # Diagnostics
            sampler_params <- get_sampler_params(fit1, inc_warmup = FALSE)
            n_div_val <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
            time_mat <- get_elapsed_time(fit1)
            exec_time_val <- sum(time_mat)
            min_n_eff_val <- min(fit_summary$summary[, "n_eff"], na.rm = TRUE)

            # Extract params
            est_ar1_val <- fit_summary$summary["beta1[1,1]", "mean"]
            est_ar2_val <- fit_summary$summary["beta1[2,2]", "mean"]
            est_ri_var1_val <- fit_summary$summary["sigmaeta2[1,1]", "mean"]
            est_ri_var2_val <- fit_summary$summary["sigmaeta2[2,2]", "mean"]
            est_ly1_val <- fit_summary$summary["ly[1]", "mean"]
            est_ly2_val <- fit_summary$summary["ly[2]", "mean"]
            est_ly3_val <- fit_summary$summary["ly[3]", "mean"]
            est_ly4_val <- fit_summary$summary["ly[4]", "mean"]
            est_phi_val <- fit_summary$summary["phi1[1,2]", "mean"]

            param_names <- rownames(fit_summary$summary)
            beta_indices <- c()
            for (t1 in 1:length(params2)) {
              beta_indices <- c(beta_indices, grep(params2[t1], param_names))
            }
            
            if(length(beta_indices) > 0) {
              conv <- max(fit_summary$summary[beta_indices, "Rhat"], na.rm = TRUE)
            } else {
              conv <- NA
            }
            
            # Compute fit indices
            modelfit <- getfit(fit1, y, y.2, nvar = 6, Nt = Nt, cfi = TRUE, 
                              fit1.null = fit1.null, calc_overall = overall_fits)

            val_id0 <- as.integer(modelfit$id0)
            idx_names <- c("BRMSEA","BGamma","BadjGamma","BMc","BCFI","BTLI")
            
            raw_overall <- modelfit$overall
            if(is.null(names(raw_overall))) names(raw_overall) <- idx_names
            
            overall_row <- data.frame(
              run_id = as.character(run_id),
              job_id = as.integer(job_id),
              array_id = as.integer(rslurm_id),
              rep = as.integer(rep),
              chunk_id = as.integer(chunk_id),
              N = as.integer(N),
              Nt = as.integer(Nt),
              M = as.numeric(M),
              T_event_prop = as.numeric(T_event_prop),
              t0 = as.integer(t0),
              id0 = as.integer(val_id0),
              conv_max_rhat = clean_val(conv),
              min_n_eff = clean_val(min_n_eff_val),
              n_divergent = clean_val(n_div_val),
              exec_time_sec = clean_val(exec_time_val),
              est_ar1 = clean_val(est_ar1_val),
              est_ar2 = clean_val(est_ar2_val),
              est_ri_var1 = clean_val(est_ri_var1_val),
              est_ri_var2 = clean_val(est_ri_var2_val),
              est_ly1 = clean_val(est_ly1_val),
              est_ly2 = clean_val(est_ly2_val),
              est_ly3 = clean_val(est_ly3_val),
              est_ly4 = clean_val(est_ly4_val),
              est_corr_factor = clean_val(est_phi_val),
              BRMSEA = clean_val(raw_overall["BRMSEA"]),
              BGamma = clean_val(raw_overall["BGamma"]),
              BadjGamma = clean_val(raw_overall["BadjGamma"]),
              BMc = clean_val(raw_overall["BMc"]),
              BCFI = clean_val(raw_overall["BCFI"]),
              BTLI = clean_val(raw_overall["BTLI"]),
              created_at = as.character(format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
              stringsAsFactors = FALSE
            )
            
            tp_mat <- modelfit$timepoint
            if (is.list(tp_mat)) {
                tp_mat <- do.call(cbind, lapply(seq_len(ncol(tp_mat)), function(i) unlist(tp_mat[,i])))
            }
            tp_df <- as.data.frame(tp_mat, stringsAsFactors = FALSE)
            colnames(tp_df) <- idx_names
            for(col in names(tp_df)) tp_df[[col]] <- as.numeric(unlist(tp_df[[col]]))
            
            tp_df$time <- seq_len(Nt)
            tp_df$run_id <- as.character(run_id)
            tp_df$job_id <- as.integer(job_id)
            tp_df$array_id <- as.integer(rslurm_id)
            tp_df$rep <- as.integer(rep)
            tp_df$N <- as.integer(N)
            tp_df$Nt <- as.integer(Nt)
            tp_df$M <- as.numeric(M)
            tp_df$T_event_prop <- as.numeric(T_event_prop)
            tp_df$t0 <- as.integer(t0)
            tp_df$id0 <- as.integer(val_id0)
            tp_df$conv_max_rhat <- clean_val(conv)
            
            tp_df <- tp_df[, c("run_id","job_id","array_id","rep","time","N","Nt","M","T_event_prop","t0","id0","conv_max_rhat", idx_names)]
            
            DBI::dbWithTransaction(con, {
              db_append_retry(con, "fit_overall", overall_row)
              db_append_retry(con, "fit_timepoint", tp_df)
            })
            
            cat("  -> Rep", rep, "saved OK.\n")

        }, error = function(e) {
            cat(sprintf("!!! ERROR in Rep %d (skipping): %s\n", rep, e$message))
        })
      } # End Inner Loop
      
      # Mark Chunk Completed
      if (is_cluster && !is.na(condition_id)) {
        mark_condition_completed(config$db_path, condition_id)
        cat(sprintf("✓ Task/Chunk %d marked as completed.\n", condition_id))
      }
    
  } # End Bundle Loop
  
  if(is_cluster) Sys.sleep(1)
  
} # End While

# Close DB connection
DBI::dbDisconnect(con)

cat("(Ignore if running locally) worker finished all assigned conditions.\n")