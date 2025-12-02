# DO NOT CHANGE THIS SECTION
## Receive the simulation parameters from job file
## Evaluate the simulation parameters in the R global environment
## For the Toy Example
## It is equivalent to run
 # n_train <- 200
 # p <- c(4, 10, 50, 100, 200)[2]
 # rho <- c(0, 0.5)[2]
 # pi_cns <- c(0.15, 0.3, 0.4)[2]

library(argparse)
parser <- ArgumentParser(description="Run a simulation study")

parser$add_argument("-n", "--n_train", type="integer", default=200,
                    help="Number of training samples")

parser$add_argument("-p", "--p", type="integer", default=10,
                    help="Number of predictors")

parser$add_argument("-r", "--rho", type="double", default=0.5,
                    help="Correlation parameter")

parser$add_argument("-c", "--pi_cns", type="double", default=0.1,
                    help="Censoring proportion")

parser$add_argument("--resPath", type="character", required=TRUE,
                    help="Path to save results files")

parser$add_argument("--scalePath", type="character", required=TRUE,
                    help="Path to save scaling info")


args <- parser$parse_args()
n_train <- args$n_train
p <- args$p
rho <- args$rho
pi_cns <- args$pi_cns
resPath <- args$resPath
scalePath <- args$scalePath

#args=(commandArgs(TRUE))
#
#if(length(args)==0){
#  print("No arguments supplied.")
#}else{
#  for(i in 1:length(args)){
#    eval(parse(text=args[[i]]))
#  }
#}



# Library & Helper Functions ----------------------------------------------
## Required Libraries
library(tidyverse)
library(mgcv)
library(cosso)
#library(BhGLM)
library(BHAM)
library(survival)
library(simsurv)
library(glmnet)

## Helper Functions
 source("Sim/Code/find_censor_parameter.R")
 source("Sim/Code/create_HD_formula.R")
 source("Sim/Code/make_null_res.R")



# Data Generating Process -------------------------------------------------
# * Simulation Parameters -------------------------------------------------
source("Sim/Code/sim_pars_funs.R")

## Job Name
job_name <- Sys.getenv('SLURM_JOB_NAME')
print(job_name)

## Use Array ID as random seed ID
it <- Sys.getenv('SLURM_ARRAY_TASK_ID') |> as.numeric()
# it <- 1
set.seed(it)

# * Generate Data -------------------------------------------------
x_all <- MASS::mvrnorm(n_train+n_test, rep(0, p), AR(p, rho)) |>
  data.frame()
eta_all <- with(x_all, f_1(X1) + f_2(X2) + f_3(X3) + f_4(X4))
# Adding adjustable signal to noise ratio


dat_all <- simsurv::simsurv(dist = "weibull",
                            lambdas = scale.t,
                            gammas = shape.t,
                            x = data.frame(eta = eta_all) ,
                            beta = c(eta = 1)) %>%
  data.frame( x_all, eta = eta_all, .)

train_dat <- dat_all[1:n_train, ]
test_dat <- dat_all[(n_train+1):n_total, ]


## Censoring Distribution, Weibull(shape.c, scale.c)
scale.c <- tryCatch({
  find_censor_parameter(lambda = exp(-1*train_dat$eta/shape.t),
                        pi.cen = pi_cns,
                        shape_hazard = shape.t, shape_censor = shape.c)
},
error = function(err) {
  if(!file.exists("Sim/Code/scale_vec.RDS"))
    stop("Please Generate scale_vec, and use 'R/calculate_scales' to generates scale_vec.RDS")
  scale_vec <- readRDS("Sim/Code/scale_vec.RDS")
  temp <- gsub("-", ",", job_name)
  scale.c <- scale_vec[[temp]]
  #scale.c <- scale_vec[["bcam_sim_p=10,rho=0.5,pi_cns=0.3"]]
  if(is.null(scale.c)) stop("No scale for this scenario")
  return(scale.c)
})



# Save Scale Parameter ----------------------------------------------------
# TODO: if you need to generate scale vector, please uncomment the following two line
# saveRDS(scale.c,
#         paste0("/data/user/boyiguo1/bcam/scale/", job_name,"/it_",it,".rds"))


train_dat <-  train_dat |>
  data.frame(
    c_time = rweibull(n = n_train, shape = shape.c, scale = scale.c)
  ) |>
  mutate(
    cen_ind = (c_time < eventtime),
    status = (!cen_ind)*1
  ) |>
  rowwise() |>
  mutate(time = min(c_time, eventtime)) |>
  ungroup()


# * Spline Specification --------------------------------------------------

mgcv_df <- data.frame(
  Var = grep("X", names(train_dat), value = TRUE),
  Func = "s",
  Args = paste0("bs='cr', k=", k)
)

# Set Up Parallelization -------------------------------------------------
library(parallel)
library(future)
library(furrr)
library(doParallel)

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
plan(multisession, workers = n_cores)
registerDoParallel(cores = n_cores)

# Fit Models------------------------------------------------------------------

#### Oracle ####
start_time_oracle <- Sys.time()
oracle_mdl <- coxph(Surv(time, status) ~ eta,
                    data = train_dat)



#Prediction
oracle_train <- BHAM::measure_cox(Surv(train_dat$time, train_dat$status),
                                  predict(oracle_mdl,
                                          newdata = train_dat |> select(eta),
                                          type = "lp")
)

oracle_test <- BHAM::measure_cox(Surv(test_dat$eventtime, test_dat$status),
                                 predict(oracle_mdl,
                                         newdata = test_dat |> select(eta),
                                         type = "lp")
)
end_time_oracle <- Sys.time()
oracle_time = end_time_oracle-start_time_oracle

#### Linear Lasso ####

start_time_lasso <- Sys.time()
lasso_mdl <- cv.glmnet(x = data.matrix(train_dat |> select(starts_with("X"))),
                       y = train_dat %>% select(time, status) |> data.matrix(),
                       nfolds = 5, family = "cox")

lasso_fnl_mdl <- glmnet(x = data.matrix(train_dat %>% select(starts_with("X"))),
                        y = train_dat %>% select(time, status) |> data.matrix(),
                        family = "cox", lambda = lasso_mdl$lambda.min)

# Prediction
lasso_train <- BHAM::measure_cox(Surv(train_dat$time, train_dat$status),
                           predict(lasso_fnl_mdl,
                                   newx = data.matrix(train_dat |> select(starts_with("X"))),
                                   type = "link")
)
lasso_test <- BHAM::measure_cox(Surv(test_dat$eventtime, test_dat$status),
                          predict(lasso_fnl_mdl,
                                  newx = data.matrix(test_dat |> select(starts_with("X"))),
                                  type = "link")
)

# Variable Selection
lasso_var <- ((lasso_fnl_mdl$beta %>% as.vector())!=0) |>
  `names<-`(names(test_dat %>% select(starts_with("X"))))


end_time_lasso <- Sys.time()
lasso_time <- end_time_lasso - start_time_lasso

# * mgcv --------------------------------------------------------------------

start_time_mgcv<- Sys.time()
mgcv_converge <- TRUE
mgcv_mdl <- tryCatch({
  gam(create_HD_formula(time~1, spl_df = mgcv_df), data = train_dat,
      family = cox.ph(), weight = status)
},
error = function(err) {
  mgcv_mdl <- NULL
  return(NULL)
})

mgcv_train <- make_null_res("cox")
mgcv_test <- make_null_res("cox")
mgcv_var <- rep(NA, p) |> `names<-`(names(test_dat %>% select(starts_with("X"))))
mgcv_plot <- NULL

if(!is.null(mgcv_mdl)){
  # Prediction Results
  mgcv_train <- BHAM::measure_cox(Surv(train_dat$time, train_dat$status) , mgcv_mdl$linear.predictors)
  mgcv_test <- BHAM::measure_cox(Surv(test_dat$eventtime, test_dat$status),
                           predict(mgcv_mdl, newdata=test_dat, type = "link"))

  # Variable Selection Results
  mgcv_var <- (summary(mgcv_mdl)$s.table[,"p-value"] < 0.05) |>
    `names<-`(names(test_dat |> select(starts_with("X"))))

  # Plotting Results
}

end_time_mgcv<- Sys.time()
mgcv_time <- end_time_mgcv - start_time_mgcv

if(all(is.na(mgcv_var))){
  mgcv_time <- NA
}


# * COSSO -------------------------------------------------------------------
start_time_cosso <- Sys.time()
cosso_mdl <-  tryCatch({cosso(x = train_dat |> select(starts_with("X")) |> data.matrix(),
                              y = train_dat |> select(time, status) |> data.matrix(), family = "Cox",
                              nbasis = k, scale = F)
},
error = function(err) {
  cosso_mdl <- NULL
  return(NULL)
}
)

if(!is.null(cosso_mdl)){
  cosso_tn_mdl <- tryCatch({
    tune.cosso(cosso_mdl, plot.it = FALSE)
  },
  error = function(err) {
    return(NULL)
  }
  )
}

# If tuning failed, increase Kfold
if(is.null(cosso_tn_mdl)){
  cosso_tn_mdl <- tryCatch({
    tune.cosso(cosso_mdl, plot.it = FALSE, folds = 10)
  },
  error = function(err) {
    return(NULL)
  }
  )

}

cosso_var <- rep(NA, p) |> `names<-`(names(test_dat |> select(starts_with("X"))))
if(!is.null(cosso_mdl) && !is.null(cosso_tn_mdl)){
  # Prediction
  cosso_train_lp <- predict.cosso(cosso_mdl,
                                  xnew=train_dat |> select(starts_with("X")) |> data.matrix(),
                                  M=ifelse(!is.null(cosso_tn_mdl), cosso_tn_mdl$OptM, 2), type = "fit")
  cosso_train <- BHAM::measure_cox(Surv(train_dat$time, train_dat$status), cosso_train_lp)

  cosso_test_lp <- predict.cosso(cosso_mdl,
                                 xnew=test_dat |> select(starts_with("X")) |> data.matrix(),
                                 M=ifelse(!is.null(cosso_tn_mdl), cosso_tn_mdl$OptM, 2), type = "fit")
  cosso_test <- BHAM::measure_cox(Surv(test_dat$eventtime, test_dat$status), cosso_test_lp)

  # Variable Selection
  cosso_var <- rep(FALSE, p) %>% `names<-`(names(test_dat |> select(starts_with("X"))))
  cosso_var[predict.cosso(cosso_mdl, M=ifelse(!is.null(cosso_tn_mdl), cosso_tn_mdl$OptM, 2), type = "nonzero")] <- TRUE

  # Effect Plotting
  cosso_plot <- NULL


} else{
  cosso_train <- make_null_res("cox")
  cosso_test <- make_null_res("cox")
}

end_time_cosso <- Sys.time()
cosso_time <- end_time_cosso - start_time_cosso

if(all(is.na(cosso_var))){
  cosso_time <- NA
}

#### Fit ACOSSO Models ####

start_time_acosso <- Sys.time()
acosso_mdl <-  tryCatch({
  wt_mdl <- SSANOVAwt(x = train_dat |> select(starts_with("X")) |> data.matrix(),
                      y = train_dat |> select(time, status) |> data.matrix(), family = "Cox", nbasis=k)
  cosso(x = train_dat |> select(starts_with("X")) |> data.matrix(),
        y = train_dat |> select(time, status) |> data.matrix(), family = "Cox",
        wt= wt_mdl, scale = F, nbasis=k)
},
error = function(err) {
  acosso_mdl <- NULL
  return(NULL)
})

if(!is.null(acosso_mdl)){
  acosso_tn_mdl <- tryCatch({

    if(acosso_mdl$tune$Mgrid[1]<0.1){
      acosso_mdl$tune$Mgrid <- acosso_mdl$tune$Mgrid[2:length(acosso_mdl$tune$Mgrid)]
      acosso_mdl$tune$ACV <- acosso_mdl$tune$ACV[2:length(acosso_mdl$tune$ACV)]
      acosso_mdl$tune$L2norm <- acosso_mdl$tune$L2norm[2:nrow(acosso_mdl$tune$L2norm),
                                                       2:ncol(acosso_mdl$tune$L2norm)]
    }

    tune.cosso(acosso_mdl, plot.it = FALSE)
  },
  error = function(err) {
    return(NULL)
  }
  )
}


# If tuning failed, increase Kfold
if(is.null(acosso_tn_mdl)){
  acosso_tn_mdl <- tryCatch({
    tune.cosso(acosso_mdl, plot.it = FALSE, folds = 10)
  },
  error = function(err) {
    return(NULL)
  }
  )

}


acosso_var <- rep(NA, p) |> `names<-`(names(test_dat |> select(starts_with("X"))))
if(!is.null(acosso_mdl) && !is.null(acosso_tn_mdl)){

  acosso_train_lp <- predict.cosso(acosso_mdl,
                                   xnew = train_dat %>% select(starts_with("X")) %>% data.matrix(),
                                   M = ifelse(!is.null(acosso_tn_mdl), acosso_tn_mdl$OptM, 2), type = "fit")
  acosso_train <- BHAM::measure_cox(Surv(train_dat$time, train_dat$status), acosso_train_lp)

  acosso_test_lp <- predict.cosso(acosso_mdl,
                                  xnew=test_dat %>% select(starts_with("X")) %>% data.matrix,
                                  M=ifelse(!is.null(acosso_tn_mdl), acosso_tn_mdl$OptM, 2), type = "fit")
  acosso_test <- BHAM::measure_cox(Surv(test_dat$eventtime, test_dat$status), acosso_test_lp)

  acosso_var <- rep(FALSE, p) %>% `names<-`(names(test_dat %>% select(starts_with("X"))))
  acosso_var[predict.cosso(acosso_mdl, M=ifelse(!is.null(acosso_tn_mdl), acosso_tn_mdl$OptM, 2), type = "nonzero")] <- TRUE

} else {
  acosso_train <- acosso_test <- make_null_res("cox")
  acosso_var <- acosso_plot <- NULL
}

end_time_acosso <- Sys.time()
acossso_time <- end_time_acosso - start_time_acosso


if(all(is.na(acosso_var))){
  acosso_time <- NA
}

# * BHAM ----------------------------------------------------------

start_time_bham <- Sys.time()
train_sm_dat <- construct_smooth_data(mgcv_df, train_dat)
train_smooth_data <- train_sm_dat$data

test_sm_dat <- BHAM::make_predict_dat(train_sm_dat$Smooth, dat = test_dat)

bam_group <- BHAM::make_group(names(train_smooth_data))


#** bmlasso -----------------------------------------------------------------
bamlasso_raw_mdl <- BHAM::bamlasso( x = train_smooth_data, y = Surv(train_dat$time, event = train_dat$status),
                              family = "cox", group = BHAM::make_group(names(train_smooth_data)),
                              ss = c(0.04, 0.5))

blasso_s0_seq <- seq(0.005, 0.1, length.out = 20)    # TODO: need to be optimized
blasso_cv_res <- BHAM::tune.bgam(bamlasso_raw_mdl, nfolds = 5, s0= blasso_s0_seq, verbose = FALSE)

blasso_s0_min <- blasso_cv_res$s0[which.min(blasso_cv_res$deviance)]
bamlasso_mdl <- BHAM::bamlasso( x = train_smooth_data, y = Surv(train_dat$time, event = train_dat$status),
                          family = "cox", group = BHAM::make_group(names(train_smooth_data)),
                          ss = c(blasso_s0_min, 0.5))

# Prediction
bamlasso_train <- BHAM::measure_cox(Surv(train_dat$time, train_dat$status) , bamlasso_mdl$linear.predictors)

bamlasso_mdl$offset <- 0
bamlasso_test_lp <- predict(bamlasso_mdl,
                            newx = as.matrix(test_sm_dat)[,colnames(bamlasso_mdl$x), drop = F],
                            type = "link")
bamlasso_test <- BHAM::measure_cox(Surv(test_dat$eventtime, test_dat$status) , bamlasso_test_lp)
#bamlasso_test <- BHAM::measure_cox(Surv(test_dat$eventtime, test_dat$status) , bamlasso_mdl)
#bamlasso_test <- measure.bh(bamlasso_mdl, test_sm_dat, Surv(test_dat$eventtime, test_dat$status))


# Variable Selection
bamlasso_vs_part <- bamlasso_var_selection(bamlasso_mdl)
bamlasso_vs_part$`Non-parametric` <- bamlasso_vs_part$`Non-parametric` %>%
  rowwise()%>%
  mutate(Selected = ifelse(Linear|Nonlinear, TRUE, FALSE))
bamlasso_var <- rep(FALSE, p) %>% `names<-`(names(test_dat %>% select(starts_with("X"))))
bamlasso_var[bamlasso_vs_part$`Non-parametric`$Variable] <- bamlasso_vs_part$`Non-parametric`$Selected

end_time_bham <- Sys.time()
bham_time <- end_time_bham - start_time_bham

if(all(is.na(bamlasso_var))){
  bham_time <- NA
}

# Effect Plotting

# Save Simulation Results -------------------------------------------------

# Overall
ret <- list(
  train_res = list(
    oracle = oracle_train,
    lasso = lasso_train,
    mgcv = mgcv_train,
    cosso = cosso_train,
    acosso = acosso_train,
    # bacox = bacox_train,
    bamlasso = bamlasso_train
  ),
  test_res = list(
    oracle = oracle_test,
    lasso = lasso_test,
    mgcv = mgcv_test,
    cosso = cosso_test,
    acosso = acosso_test,
    # bacox = bacox_test,
    bamlasso = bamlasso_test
  ),

  censoring = list(
   scale.c = scale.c,                             # The scale parameter
  p.cen_train = mean(train_dat$status==0)              # Censoring proportion in training data
 ),

  var_slct = list(
    lasso = lasso_var,
    mgcv = mgcv_var,
    cosso = cosso_var,
    acosso = acosso_var,
    bamlasso = bamlasso_var
  ),

  bam_select = bamlasso_vs_part,

 timing = list(
   oracle = oracle_time,
   lasso = lasso_time,
   mgcv = mgcv_time,
   cosso = cosso_time,
   acosso = acosso_time,
   bamlasso = bham_time
 )
)

out_dir <- if (exists("resPath")) resPath else "Sim/Res"
saveRDS(ret, file.path(out_dir, sprintf("it_%s.rds", it)))

# Recommendation: to save the results in individual rds files
#saveRDS(ret,
#        paste0("/data/user/boyiguo1/bcam/Res/", job_name,"/it_",it,".rds"))
