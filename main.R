### INITIAL CONFIGURATION ###
# Install packages if needed and load
sapply(c("yaml", "parallelly", "this.path" ,"mvtnorm", "matrixcalc", "Matrix", "rstan", "loo"), 
       function(pkg) if (!require(pkg, character.only = TRUE)) {
         install.packages(pkg)
         library(pkg, character.only = TRUE)
       })

rstan_options(auto_write = TRUE)

# Cluster
options(mc.cores = parallelly::availableCores())
rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(rslurm_id)) rslurm_id <- 1  # pruebas locales

# Working environment
setwd(this.path::this.dir())
getwd()
config <- yaml.load_file("config.yaml") # Production
# config <- yaml.load_file("config_test.yaml") # Local testing
source("src/dgp_sem.R")
source("src/fitindices_utils.R")

#setwd("C:\\holger\\SEM\\dsem_modelfit_timepoint\\sim")
#setwd("C:\\holger\\SEM\\modelfit\\stanversion")
#source("C:/holger/SEM/modelfit/stanversion/funs/lavaan_dsem_models_randomintercept.R")
#source("C:/holger/SEM/modelfit/stanversion/funs/lavaan_dsem_models_randomintercept_ar.R")
#source("C:/holger/SEM/modelfit/stanversion/funs/lavaan_dsem_models_randomintercept_tt1.R")
#source("C:/holger/SEM/modelfit/stanversion/funs/lavaan_dsem_models_comp1.R")
#source("C:/holger/SEM/modelfit/stanversion/funs/lavaan_dsem_models_comp2.R")
#source("C:/holger/SEM/modelfit/stanversion/funs/gen_data_version04.R")

### GET SETTINGS FROM CONFIG FILE ###
seed <- config$seed
N.all <- config$settings$N
Nt.all <- config$settings$N_t
size_miss <- config$settings$M  # misspecification magnitude
T_event <- config$settings$T_event

#MCMC
iter_n <- if(!is.null(config$mcmc$iter)) config$mcmc$iter else 2000
chains_n <- if(!is.null(config$mcmc$chains)) config$mcmc$chains else 4
warmup_n <- if(!is.null(config$mcmc$warmup)) config$mcmc$warmup else 1000
cores_n <- if(!is.null(config$mcmc$cores)) config$mcmc$cores else 1


#type_TRUE_MISS_SIMULATE <- c("tt","tt1") # within time points, between time points #"none",
#run_Samples_SIMULATE <- 50
##########################


##########################
phi0 <- diag(2)*.7+.3 # cov(eta)
mu0  <- c(0,0)        #mean(eta)
ar0  <- c(.3,.3)      # ar(1) structure
ly00 <- .6
ly0  <- matrix(c(ly00,ly00,ly00,0,0,0,
                 0,0,0,ly00,ly00,ly00),6,2,byrow=F) # factor loadings
td   <- diag(6)*(1-ly00^2) # cond. var(y) -> res var
##########################


##########################
# load stan file
rt1 <- stanc("models/DSEM_VAR1_RI.stan") 
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)
params2 <- c("ly","sigmaeps","sigmaeta1","sigmaeta2","ty","ka","beta1")

rt1.null <- stanc("models/DSEM_NULL.stan") 
sm1.null <- stan_model(stanc_ret = rt1.null, verbose=FALSE)
##########################

# Helper function to calculate T_event based on N_t and proportion
calculate_T_event <- function(N_t, proportion) {
  t0 <- 1 + floor((N_t - 1) * proportion)
  t0 <- max(2, min(N_t, t0))
  t0
}


#diag((ly0)%*%phi0%*%t(ly0)+td)

# Create output directory
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# Main loop
for (i in 1:length(N.all)) {
  for (j in 1:length(Nt.all)) {
    
    for (l in 1:length(size_miss)) {
      for (k in 1:length(T_event)) {
        
        # Set different seed for each iteration -< avoids producing same numbers in each iteration
        set.seed(seed + i*1 + j*10 + l*100 + k*1000 + rslurm_id*10000)
        
        N <- N.all[i]
        Nt <- Nt.all[j]
        t0 <- calculate_T_event(Nt, T_event[k])
        
        ly1 <- matrix(0, 6, 2)
        ly1[4, 1] <- size_miss[l]
        
        dat1 <- gendata02(N, Nt, phi0, mu0, ar0, ly0, ly1, td, t0 = t0)
        
        y <- dat1[["y"]]
        y.2 <- dat1[["y0"]]
        
        data.stan1 <- list(N = N, Nt = Nt, y = y)
        fit1 <- sampling(sm1, data = data.stan1, 
                         iter = iter_n, 
                         chains = chains_n, 
                         warmup = warmup_n,
                         cores = cores_n,
                         refresh = 0) # refresh 0 evita que llene la consola de logs
        
        fit1.null <- sampling(sm1.null, data = data.stan1, 
                              iter = iter_n, 
                              chains = chains_n, 
                              warmup = warmup_n,
                              cores = cores_n,
                              refresh = 0)
        
        fit_summary <- summary(fit1)
        param_names <- rownames(fit_summary$summary)
        
        beta_indices <- c()
        for (t1 in 1:length(params2)) {
          beta_indices <- c(beta_indices, grep(params2[t1], param_names))
        }
        conv <- max(fit_summary$summary[beta_indices, "Rhat"])
        
        modelfit <- getfit(fit1, y, y.2, nvar = 6, Nt = Nt, cfi = TRUE, fit1.null = fit1.null)
        modelfit[[3]] <- conv
        
        nomfile <- paste0(
          "results/N", N,
          "_Nt", Nt,
          "_M", size_miss[l],
          "_t0", t0,
          "_prop", T_event[k], 
          "_", rslurm_id,
          ".RDS"
        )
        saveRDS(modelfit, file = nomfile)
        
      } # end T_event loop
    } # end M loop
  }
}


#plot(modelfit[[2]][,6],type="l")


#################
# compare to blavaan as check

#res1 <- bsem(dsem[[2]], data=y.2,std.lv = TRUE,
#             burnin = 1000, n.chains = 4, sample = 1000)
#
#res0 <- bsem(null.dsem[[2]], data=y.2,
#             burnin = 1000, n.chains = 4, sample = 1000)
#summary(res0)
#
#a1 <- blavFitIndices(res1, baseline.model = res0)
#fitm1 <- matrix(unlist(lapply(a1@indices,mean)),nrow=1)


#########################################################
#null.dsem[[2]] <- '
#y1t1 ~~ y1t1
#y1t2 ~~ y1t2
##
#y2t1 ~~ y2t1
#y2t2 ~~ y2t2
##
#y3t1 ~~ y3t1
#y3t2 ~~ y3t2
##
#y4t1 ~~ y4t1
#y4t2 ~~ y4t2
##
#y5t1 ~~ y5t1
#y5t2 ~~ y5t2
##
#y6t1 ~~ y6t1
#y6t2 ~~ y6t2
##
#y1t1 ~ 1
#y2t1 ~ 1
#y3t1 ~ 1
#y4t1 ~ 1
#y5t1 ~ 1
#y6t1 ~ 1
#y1t2 ~ 1
#y2t2 ~ 1
#y3t2 ~ 1
#y4t2 ~ 1
#y5t2 ~ 1
#y6t2 ~ 1
#'


resultado <- readRDS("results/N10_Nt10_M0.3_t05_prop0.5_1.RDS")
print(resultado)


