# Mean estimates and CI from outer-resampling distribution

# For each metric and time horizon:
#   mean = average over outer iterations
#   SE   = SD across outer iterations / sqrt(number of valid outer iterations)
#   CI   = percentile-based 95% CI for the mean outer-fold performance


summarise_metric <- function(x, probs = c(0.025, 0.975)) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n == 0L) {
    return(list(mean = NA_real_, sd = NA_real_,
                lower = NA_real_, upper = NA_real_, n = 0L))
  }
  
  if (n == 1L) {
    return(list(mean = mean(x), sd = NA_real_,
                lower = NA_real_, upper = NA_real_, n = 1L))
  }
  
  qs <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 6)
  
  list(
    mean = mean(x),
    sd = stats::sd(x),
    lower = qs[1],
    upper = qs[2],
    n = n
  )
}

summary_table <- results_long[, {
  cind <- summarise_metric(cindex)
  br   <- summarise_metric(brier)
  ic   <- summarise_metric(ici)
  
  .(
    cindex_mean = cind$mean,
    cindex_ci_lower = cind$lower,
    cindex_ci_upper = cind$upper,
    cindex_n_outer = cind$n,
    
    brier_mean = br$mean,
    brier_ci_lower = br$lower,
    brier_ci_upper = br$upper,
    brier_n_outer = br$n,
    
    ici_mean = ic$mean,
    ici_ci_lower = ic$lower,
    ici_ci_upper = ic$upper,
    ici_n_outer = ic$n
  )
}, by = .(model, time)]