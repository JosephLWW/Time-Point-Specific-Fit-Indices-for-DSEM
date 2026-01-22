# ============================================================================
# workload_manager.R
# ============================================================================
# Task Farming System for DSEM Simulation Study
# Distributes 324 conditions × 250 replications across 50 SLURM array jobs
# 
# Role: Assign bundles (packages) of conditions to individual workers
# Each worker gets a batch of conditions to execute, reducing DB contention
# ============================================================================

# ============================================================================
# 1. INITIALIZATION: Create Tasks Database Table
# ============================================================================

init_tasks_db <- function(db_path, conditions_grid, batch_size = 3) {
  require(RSQLite)
  
  db <- dbConnect(RSQLite::SQLite(), db_path)
  
  tryCatch({
    dbExecute(db, "PRAGMA journal_mode = WAL")
    dbExecute(db, "PRAGMA busy_timeout = 60000")  # ← Aumentar a 60s
    
    dbExecute(db, "
      CREATE TABLE IF NOT EXISTS tasks (
        condition_id INTEGER PRIMARY KEY,
        N INTEGER,
        Nt INTEGER,
        M REAL,
        T_event_prop REAL,
        status INTEGER DEFAULT 0,
        job_assigned_to INTEGER,
        started_at TEXT,
        completed_at TEXT,
        error_msg TEXT
      )
    ")
    
    n_existing <- dbGetQuery(db, "SELECT COUNT(*) FROM tasks")[1, 1]
    
    if (n_existing == 0) {
      conditions_with_id <- cbind(
        condition_id = seq_len(nrow(conditions_grid)),
        conditions_grid
      )
      
      dbWriteTable(db, "tasks", conditions_with_id, append = TRUE, row.names = FALSE)
      cat(sprintf("Created tasks table with %d conditions\n", nrow(conditions_grid)))
    } else {
      cat(sprintf("Tasks table already exists with %d conditions\n", n_existing))
    }
    
    dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_status ON tasks(status)")
    
  }, finally = {
    dbDisconnect(db)  # ← SIEMPRE cierra
  })
  
  invisible(TRUE)  # ← Devuelve TRUE en lugar de conexión
}


# ============================================================================
# 2. WORKER FUNCTION: Request and Execute Bundle
# ============================================================================

get_bundle <- function(db_path, job_id, batch_size = 3, 
                       n_jobs = 50, conditions_grid = NULL) {
  require(RSQLite)
  
  db <- dbConnect(RSQLite::SQLite(), db_path)
  
  # ← AÑADIR tryCatch AQUÍ para asegurar cierre
  tryCatch({
    dbExecute(db, "PRAGMA busy_timeout = 60000")
    
    max_retries <- 20
    retry_count <- 0
    
    while (retry_count < max_retries) {
      tryCatch({
        dbBegin(db)
        
        pending <- dbGetQuery(db, "
          SELECT condition_id, N, Nt, M, T_event_prop
          FROM tasks 
          WHERE status = 0
          ORDER BY condition_id
          LIMIT ?
        ", params = list(batch_size))
        
        if (nrow(pending) == 0) {
          dbRollback(db)
          return(list(
            conditions = NULL,
            bundle_id = NA,
            job_id = job_id,
            status = "empty",
            n_tasks = 0
          ))
        }
        
        condition_ids <- pending$condition_id
        
        dbExecute(db, "
          UPDATE tasks 
          SET status = 1, 
              job_assigned_to = ?, 
              started_at = datetime('now')
          WHERE condition_id IN (
            SELECT condition_id FROM (
              SELECT condition_id FROM tasks 
              WHERE status = 0
              ORDER BY condition_id
              LIMIT ?
            )
          )
        ", params = list(job_id, batch_size))
        
        dbCommit(db)
        
        bundle_id <- sprintf("job%03d_batch%04d", job_id, min(condition_ids))
        
        cat(sprintf(
          " Job %2d: Assigned %d conditions (IDs %d-%d) [bundle: %s]\n",
          job_id, nrow(pending), min(condition_ids), max(condition_ids), bundle_id
        ))
        
        return(list(
          conditions = pending,
          bundle_id = bundle_id,
          job_id = job_id,
          status = "success",
          n_tasks = nrow(pending)
        ))
        
      }, error = function(e) {
        retry_count <<- retry_count + 1
        
        if (retry_count < max_retries) {
          wait_time <- 0.5 * (2 ^ (retry_count - 1))
          cat(sprintf("DB lock (retry %d/%d, waiting %.1fs): %s\n", 
                      retry_count, max_retries, wait_time, conditionMessage(e)))
          Sys.sleep(wait_time)
        } else {
          stop(sprintf("Failed to get bundle after %d retries: %s", 
                       max_retries, conditionMessage(e)))
        }
      })
    }
  }, finally = {
    # ← SIEMPRE cierra la conexión
    if (exists("db") && dbIsValid(db)) dbDisconnect(db)
  })
}



# ============================================================================
# 3. MARK COMPLETION: Update Task Status
# ============================================================================

mark_condition_completed <- function(db_path, condition_id, error_msg = NULL) {
  require(RSQLite)
  
  max_retries <- 10
  
  for (attempt in 1:max_retries) {
    tryCatch({
      db <- dbConnect(RSQLite::SQLite(), db_path)
      dbExecute(db, "PRAGMA busy_timeout = 60000")
      
      status_val <- 2  # Always mark as completed (even if error)
      
      dbExecute(db, "
        UPDATE tasks 
        SET status = ?, 
            completed_at = datetime('now'),
            error_msg = ?
        WHERE condition_id = ?
      ", params = list(status_val, error_msg, condition_id))
      
      dbDisconnect(db)
      return(invisible(TRUE))
      
    }, error = function(e) {
      if (exists("db") && dbIsValid(db)) dbDisconnect(db)
      
      if (attempt < max_retries) {
        wait_time <- 0.05 * (2 ^ (attempt - 1))  # Exponential backoff
        Sys.sleep(wait_time)
      } else {
        cat(sprintf("Failed to mark condition %d after %d retries: %s\n", 
                    condition_id, max_retries, e$message))
        return(invisible(FALSE))
      }
    })
  }
}


# ============================================================================
# 4. MONITORING: Get Status Summary
# ============================================================================

get_task_status <- function(db_path) {
  # """
  # Get summary of task completion status.
  
  # Returns:
  #   Data frame with columns:
  #     - status: "pending", "in_progress", "completed"
  #     - count: number of tasks in that status
  #     - pct_total: percentage of total tasks
  # """
  
  require(RSQLite)
  
  db <- dbConnect(RSQLite::SQLite(), db_path)
  
  status_map <- c("0" = "pending", "1" = "in_progress", "2" = "completed")
  
  summary <- dbGetQuery(db, "
    SELECT status, COUNT(*) as count
    FROM tasks
    GROUP BY status
  ")
  
  dbDisconnect(db)
  
  # Map numeric status to labels
  summary$status_label <- status_map[as.character(summary$status)]
  total <- sum(summary$count)
  summary$pct_total <- round(100 * summary$count / total, 1)
  
  return(summary[, c("status_label", "count", "pct_total")])
}


# ============================================================================
# 5. CONDITION ID MAPPING: Create Global Unique ID
# ============================================================================

# calculate_seed in main.R

# ============================================================================
# 6. BUNDLE EXECUTION WRAPPER
# ============================================================================

# In main.R

# ============================================================================
# 7. MAIN ENTRY POINT: Task Farming Worker
# ============================================================================

# in main.R

# ============================================================================
# END OF workload_manager.R
# ============================================================================
