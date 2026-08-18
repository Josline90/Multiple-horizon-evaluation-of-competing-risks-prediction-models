#  Calibration / ICI function; Austine et. al 2022
#------------------------------------------------------------------------------
# - Computes integrated calibration index and Builds smoothed calibration curves 
#   by fitting a Fine–Gray subdistribution hazard model 
#   (covariate is a spline of the complementary log–log transform of the predicted risk).

# - For nested repeated cross-validation strategy, calibration model is fitted 
#   within the outer test fold to estimate observed risk as a smooth function of predicted risk.

.harrell_probs <- function(k) {
  switch(
    as.character(k),
    "3" = c(0.10, 0.50, 0.90),
    "4" = c(0.05, 0.35, 0.65, 0.95),
    "5" = c(0.05, 0.275, 0.50, 0.725, 0.95),
    stop("Supported k are 3, 4, or 5.")
  )
}

compute_ici <- function(test_data, pred_risk, t0, cause, k = 3, eps = 1e-8) {
  stopifnot(nrow(test_data) == length(pred_risk))
  
  # Avoid 0 or 1 before complementary log-log transformation.
  p_trans <- pmin(pmax(as.numeric(pred_risk), eps), 1 - eps)
  x <- log(-log(1 - p_trans))
  
  probs <- .harrell_probs(k)
  qs <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
  
  # If predictions are almost constant, calibration spline cannot be fit reliably.
  if (length(unique(round(qs, 10))) < length(qs)) {
    return(NA_real_)
  }
  
  boundary <- c(qs[1], qs[length(qs)])
  internal <- if (length(qs) > 2) qs[2:(length(qs) - 1)] else NULL
  
  spline_basis <- splines::ns(x, knots = internal, Boundary.knots = boundary)
  colnames(spline_basis) <- paste0("spl", seq_len(ncol(spline_basis)))
  
  df_fit <- cbind(as.data.frame(test_data[, .(time, status)]), as.data.frame(spline_basis))
  rhs <- paste(colnames(spline_basis), collapse = " + ")
  fml <- as.formula(paste0("Hist(time, status) ~ ", rhs))
  
  fg_cal <- tryCatch({
    riskRegression::FGR(formula = fml, data = df_fit, cause = cause)
  }, error = function(e) NULL)
  
  if (is.null(fg_cal)) return(NA_real_)
  
  obs_hat <- tryCatch({
    as.numeric(riskRegression::predictRisk(fg_cal, newdata = df_fit, times = t0))
  }, error = function(e) rep(NA_real_, length(p_trans)))
  
  if (any(!is.finite(obs_hat))) return(NA_real_)
  
  mean(abs(obs_hat - p_trans), na.rm = TRUE)
  
}

# For calibration plot

calib_plot <- function(test_data, pred_risk, t0, cause, k = 3, eps = 1e-8, grid_points = 200) {
  stopifnot(nrow(test_data) == length(pred_risk))
  
  # Avoid 0 or 1 before complementary log-log transformation.
  p_trans <- pmin(pmax(as.numeric(pred_risk), eps), 1 - eps)
  x <- log(-log(1 - p_trans))
  
  probs <- .harrell_probs(k)
  qs <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
  
  # If predictions are almost constant, calibration spline cannot be fit reliably.
  if (length(unique(round(qs, 10))) < length(qs)) {
    return(NA_real_)
  }
  
  boundary <- c(qs[1], qs[length(qs)])
  internal <- if (length(qs) > 2) qs[2:(length(qs) - 1)] else NULL
  
  spline_basis <- splines::ns(x, knots = internal, Boundary.knots = boundary)
  colnames(spline_basis) <- paste0("spl", seq_len(ncol(spline_basis)))
  
  df_fit <- cbind(test_data[, c("time", "status")], as.data.frame(spline_basis))
  rhs <- paste(colnames(spline_basis), collapse = " + ")
  fml <- as.formula(paste0("Hist(time, status) ~ ", rhs))
  
  fg_cal <- tryCatch({
    riskRegression::FGR(formula = fml, data = df_fit, cause = cause)
  }, error = function(e) NULL)
  
  if (is.null(fg_cal)) return(NA_real_)
  
  #  Calibration plot on grid
  p_grid <- seq(max(min(p_trans), eps), min(max(p_trans), 1 - eps),
                length.out = grid_points) 
  x_grid <- log(-log(1 - p_grid))
  spline_grid <- ns(x_grid, knots = internal, Boundary.knots = boundary)
  colnames(spline_grid) <- colnames(spline_basis)
  new_grid <- as.data.frame(spline_grid)
  y_grid <- as.numeric(riskRegression::predictRisk(
    fg_cal, 
    newdata = new_grid, 
    times = t0
  ))
  curve_df <- data.frame(pred = p_grid, obs = y_grid)
   
  # # Plot
  plt <- ggplot2::ggplot(curve_df, aes(x = pred, y = obs)) + 
    geom_line(linewidth = 1) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + 
    labs(
      x = paste0("Predicted risk at time = ", t0),
      y = paste0("Observed risk at time = ", t0)
    ) + 
    theme_minimal(base_size = 13)
}





