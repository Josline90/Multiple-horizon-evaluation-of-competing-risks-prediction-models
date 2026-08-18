# Nested repeated cross-validation strategy: nested resampling design
# - Hyperparameters tuned on the training set of each outer split (inner CV)
# - Performance is computed only on held-out outer test folds

# tuned RSF hyperparameters by max mean C-index_t across the eval times

# load libraries
pacman::p_load(here, pec, mlr3, mlr3cmprsk, mlr3proba, mlr3tuning, mlr3extralearners, 
               bbotk, paradox, randomForestSRC, survival, riskRegression, prodlim, splines, ggplot2) #lgr, cmprsk

# ICI + calibration plot functions
source(here::here('integrated_calibration_index.R'))

# Data set; PBC
source(here::here('pbc_data.R'))
cause <- 2 #death as event of interest
time_horizon <- c(22, 39, 75) #evaluation time points

# Basic checks
stopifnot(all(c("time", "status") %in% names(dat)))
stopifnot(all(time_horizon > 0))
stopifnot(cause %in% unique(dat$status))
stopifnot(!any(sapply(dat, is.ordered)))

# Formula for the final model fit inside @outer fold
rsf_formula <- Surv(time, status) ~ 
  bili + albumin + copper + alk.phos + protime + ast +
  chol + age + trig + platelet + trt + sex +
  hepato + spiders + edema + stage

set.seed(1357) 
# Outer folds
outer_folds   <- 5
outer_repeats <- 50 # more iterations better for small data

# Inner folds/ tuning CV; inside each outer training set
inner_folds <- 5
n_evals_tuning <- 200 #number of hyper. combinations tried during inner CV

# Creating mlr3 competing risk task 
task <- mlr3cmprsk::as_task_cmprsk(dat, time = "time", event = "status", id = "PBC_nested")

# and stratification for the outer resampling
# - stratified across event types and time points; 
# - making strata from event status and quantile-based time bins
strata_dt <- task$data()
strata_dt[, time_bin := as.integer(cut(
  time,
  breaks = unique(quantile(time, probs = seq(0, 1, length.out = 5), na.rm = TRUE)),
  labels = FALSE,
  include.lowest = TRUE
))]
strata_dt[, strata := interaction(status, time_bin, drop = TRUE)]

# Adding strata column to the mlr3 task
task$cbind(strata_dt[, .(time_bin, strata)])
task$set_col_roles(c("time_bin", "strata"), roles = character(0))
task$set_col_roles("strata", add_to = "stratum")


# Inner tuning measure; pec::cindex evaluated at multiple horizons
# For @ hyper. config. in the inner CV;
# - train RSF on inner-training folds
# - predict CIF for the event of interest on inner-validation folds
# - compute pec::cindex at each time_horizon
# - average C-index across horizons
# - choose the hyperparameter set with the max mean C-index
# all tunable models optimized for discrimination; 
# other metrics on the held-out outer test folds

MeasurePEC_cindex_multi <- R6::R6Class(
  "MeasurePEC_cindex_multi",
  inherit = mlr3::Measure,
  public = list(
    time_col = NULL,
    status_col = NULL,
    cause = NULL,
    t0_vec = NULL,
    
    # measure properties
    initialize = function(time_col = "time", status_col = "status", cause, t0_vec) {
      super$initialize(
        id = "pec_cindex_multi",
        minimize = FALSE,
        predict_type = "cif",
        properties = c("requires_task")
      )
      
      # store parameters
      self$time_col <- time_col
      self$status_col <- status_col
      self$cause <- cause
      self$t0_vec <- t0_vec
    }
  ),
  
  private = list(
    .score = function(prediction, task, ...) { #called once @fold
      
      stopifnot(!is.null(task))
      dat_fold <- task$data(rows = prediction$row_ids)
      cif_list <- prediction$cif # CIF predictions: a list of matrices
      cause_id <- as.character(self$cause)
      
      if (!cause_id %in% names(cif_list)) {
        stop("The requested cause is not present in prediction$cif.")
      }
      
      pred_times <- as.numeric(colnames(cif_list[[cause_id]]))
      
      # Safe evaluation times (trim to max follow-up)
      eval_times <- pmin(self$t0_vec, max(dat_fold[[self$time_col]], na.rm = TRUE))
      
      f <- as.formula(paste0("Hist(", self$time_col, ",", self$status_col, ") ~ 1"))
      cindex_values <- rep(NA_real_, length(eval_times))
      
      for (j in seq_along(eval_times)) {
        t0 <- eval_times[j]
        idx_time <- which.min(abs(pred_times - t0))
        risk_vec <- cif_list[[cause_id]][, idx_time]
        risk_mat <- matrix(risk_vec, ncol = 1)
        
        cindex_values[j] <- tryCatch({
          as.numeric(pec::cindex(
            object = risk_mat,
            formula = f,
            data = dat_fold,
            cause = self$cause,
            eval.times = t0,
            pred.times = t0,
            cens.model = "marginal",
            tiedPredictionsIn = TRUE,
            tiedOutcomeIn = FALSE,
            tiedMatchIn = FALSE
          )$AppCindex)[1]
        }, error = function(e) NA_real_)
      }
      
      out = mean(cindex_values, na.rm = TRUE)
      cat("[Measure] Fold C-index (mean over times):", round(out, 4), "\n")
      return(out)
    }
  )
)

# The C-index
inner_measure <- MeasurePEC_cindex_multi$new(
  time_col = "time",
  status_col = "status",
  cause = cause,
  t0_vec = time_horizon
)

# RSF learner and search space for inner loop
rsf_learner <- lrn(
  "cmprsk.rfsrc",
  splitrule = "logrankCR",
  samptype = "swr",
  cause = cause,
  predict_type = "cif"
)

p <- length(task$feature_names)
search_space <- ps(
  mtry      = p_int(lower = max(1L, floor(0.5 * sqrt(p))), upper = min(p, ceiling(1.75 * sqrt(p)))),
  nodesize  = p_int(lower = 15L, upper = 40L),
  nsplit    = p_int(lower = 0L,  upper = 5L),
  nodedepth = p_int(lower = 3L,  upper = 5L)
)

# tuning setup
inner_resampling <- rsmp("cv", folds = inner_folds)
tuner <- tnr("random_search")
terminator <- trm("evals", n_evals = n_evals_tuning)

# Create outer repeated CV splits/ resampling: should be same for all models
outer_resampling <- rsmp("repeated_cv", folds = outer_folds, repeats = outer_repeats)
outer_resampling$instantiate(task)

n_outer <- outer_resampling$iters
message("Number of outer evaluation iterations: ", n_outer)

# Save the same splits for other models
n_outer2 <- data.table(
  outer_iteration = seq_len(outer_resampling$iters),
  train_ids = lapply(
    seq_len(outer_resampling$iters),
    function(i) outer_resampling$train_set(i)
  ),
  test_ids = lapply(
    seq_len(outer_resampling$iters),
    function(i) outer_resampling$test_set(i)
  )
)
#saveRDS(n_outer2, 'shared_outer_splits.rds')


# Storage objects.
outer_results <- list()
outer_best_params <- list()
outer_predictions <- list()

# Main nested repeated CV loop

for (i in seq_len(n_outer)) {
  
  message("\n================ OUTER ITERATION ", i, " / ", n_outer, " ================")
  
  # - Define outer training and outer test row ids
  outer_train_ids <- outer_resampling$train_set(i)
  outer_test_ids  <- outer_resampling$test_set(i)
  
  train_df <- as.data.table(task$data(rows = outer_train_ids))
  test_df  <- as.data.table(task$data(rows = outer_test_ids))
  
  message("Outer train n = ", nrow(train_df), "; outer test n = ", nrow(test_df))
  message("Outer test event counts: ",
          paste(names(table(test_df$status)), table(test_df$status), collapse = "; "))
  
  # - Inner hyperparameter tuning RSF on outer_training set
  learner_i <- rsf_learner$clone(deep = TRUE)
  
  autotuner_i <- AutoTuner$new(
    learner = learner_i,
    resampling = inner_resampling,
    measure = inner_measure,
    search_space = search_space,
    terminator = terminator,
    tuner = tuner,
    store_tuning_instance = TRUE
  )
  
  # mlr3 performs the inner CV and selects hyperparameters
  task_outer_train <- task$clone(deep = TRUE)
  task_outer_train$filter(outer_train_ids)
  
  autotuner_i$train(task_outer_train) # Train AutoTuner on outer_train rows only
  
  # Extract the best hyperparameter values selected by inner CV
  best_params_i <- autotuner_i$archive$best()$x_domain[[1]]
  
  outer_best_params[[i]] <- data.table::as.data.table(best_params_i)
  outer_best_params[[i]][, outer_iteration := i]
  
  message("Best inner-CV parameters: ", paste(names(best_params_i), best_params_i, sep = "=", collapse = ", "))
  
  
  # Refit final RSF on the full outer_train data with selected parameters
  rsf_fit_i <- randomForestSRC::rfsrc(
    formula = rsf_formula,
    data = train_df,
    cause = cause,
    splitrule = "logrankCR",
    samptype = "swr",
    save.memory = TRUE,
    mtry = as.integer(best_params_i$mtry),
    nodesize = as.integer(best_params_i$nodesize),
    nsplit = as.integer(best_params_i$nsplit),
    nodedepth = as.integer(best_params_i$nodedepth),
    ntree = 1000,
    seed = 1357 + i
  )
  
  # Predict risk/CIF on outer_test data at all requested time horizons
  # - pred_rsf_i is an n_test x length(time_horizon) matrix
  pred_rsf_i <- riskRegression::predictRisk(
    object = rsf_fit_i,
    newdata = test_df,
    times = time_horizon,
    cause = cause
  )
  
  #pred_rsf_i <- as.matrix(pred_rsf_i)
  colnames(pred_rsf_i) <- paste0("risk_t", time_horizon)
  
  # patient-level outer-test predictions: for other analysis
  pred_dt_i <- data.table(
    outer_iteration = i,
    row_id = outer_test_ids,
    time = test_df$time,
    status = test_df$status
  )
  pred_dt_i <- cbind(pred_dt_i, as.data.table(pred_rsf_i))
  outer_predictions[[i]] <- pred_dt_i
  
  
  
  # Compute time-dependent C-index on outer_test
  cindex_i <- rep(NA_real_, length(time_horizon))
  
  for (j in seq_along(time_horizon)) {
    
    t0 <- time_horizon[j]
    
    cindex_i[j] <- tryCatch({
      
      ci_obj <- pec::cindex(
        object = matrix(pred_rsf_i[, j], ncol = 1),
        formula = Hist(time, status) ~ 1,
        data = test_df,
        cause = cause,
        eval.times = t0,
        pred.times = t0,
        cens.model = "marginal",
        tiedPredictionsIn = TRUE,
        tiedOutcomeIn = FALSE,
        tiedMatchIn = FALSE
      )
      
      as.numeric(ci_obj$AppCindex)[1]
      
    }, error = function(e) {
      message("C-index failed: outer iteration ", i,
              ", time = ", t0,
              ". Error: ", conditionMessage(e))
      NA_real_
    })
  }
  
  # Compute Brier score on outer_test
  brier_i <- tryCatch({
    sc <- riskRegression::Score(
      object = list(RSF = rsf_fit_i),
      formula = Hist(time, status) ~ 1,
      data = test_df,
      times = time_horizon,
      metrics = "brier",
      cause = cause,
      null.model = FALSE,
      cens.model = "marginal",
      se.fit = FALSE,
      conf.int = FALSE
    )
    
    bdt <- as.data.table(sc$Brier$score)
    # Typical columns include model, times, Brier.
    if ("Brier" %in% names(bdt)) {
      bdt[model == "RSF"][match(time_horizon, times), Brier]
    } else if ("brier" %in% names(bdt)) {
      bdt[model == "RSF"][match(time_horizon, times), brier]
    } else {
      rep(NA_real_, length(time_horizon))
    }
  }, error = function(e) rep(NA_real_, length(time_horizon)))
  
  
  # Compute ICI on outer_test
  # - computed separately at each time horizon using the corresponding
  # - predicted risk column
  
  ici_i <- rep(NA_real_, length(time_horizon))
  for (j in seq_along(time_horizon)) {
    ici_i[j] <- compute_ici(
      test_data = test_df,
      pred_risk = pred_rsf_i[, j],
      t0 = time_horizon[j],
      cause = cause,
      k = 3
    )
  }
  
  # Save one row per outer iteration and time horizon
  
  outer_results[[i]] <- data.table(
    model = "RSF",
    tuning_metric = "mean time-dependent C-index across horizons",
    outer_iteration = i,
    outer_train_n = nrow(train_df),
    outer_test_n = nrow(test_df),
    time = time_horizon,
    cindex = cindex_i,
    brier = as.numeric(brier_i),
    ici = as.numeric(ici_i)
  )
  
  print(outer_results[[i]])
  
}

# Combine outer-fold results

results_long <- rbindlist(outer_results, fill = TRUE)
best_params_all <- rbindlist(outer_best_params, fill = TRUE)
predictions_all <- rbindlist(outer_predictions, fill = TRUE)## Combine all outer-fold patient predictions

# Mean estimates and CI from outer-resampling distribution
source(here::here('summary_metric.R'))
print(summary_table)

# Calibration plot
calibration_outputs <- list()

for (t0 in time_horizon) {
  
  pred_col <- paste0("risk_t", t0)
  
  calib_data_t <- data.frame(
    time = predictions_all$time,
    status = predictions_all$status
  )
  
  calibration_outputs[[as.character(t0)]] <- calib_plot(
    test_data = calib_data_t,
    pred_risk = predictions_all[[pred_col]],
    t0 = t0,
    cause = cause,
    k = 3,
    eps = 1e-8,
    grid_points = 200
  )
  
}

# Marginal histogram of pred
ggplot(
  predictions_all,
  aes(x = risk_t22)
) +
  geom_histogram(
    bins = 30,
    color = "black",
    fill = "grey80"
  ) +
  labs(
    title = "RSF predicted risk distribution at time 22",
    x = "Predicted risk",
    y = "Count"
  ) +
  theme_minimal(base_size = 13)
