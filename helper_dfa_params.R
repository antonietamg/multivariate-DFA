# Helper utilities for choosing DFA parameters (pol and n_scales).
#
# Usage:
#   source("get_scale_function.R")
#   source("Multi_DFA.R")
#   source("helper_dfa_params.R")
#   x <- read.csv("todasvocesRESC.csv")
#   out <- choose_dfa_params(x)
#   print(out$results)
#   print(out$best)

choose_dfa_params <- function(
  x,
  pol_values = c(1, 2, 3),
  n_scales_values = c(12, 20, 30),
  signal_col = 2
) {
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (ncol(x) < 2) stop("x must include time column + at least one signal column")
  if (!signal_col %in% seq_len(ncol(x))) stop("signal_col is out of bounds")

  rows <- list()
  idx <- 1

  for (pol in pol_values) {
    for (n_scales in n_scales_values) {
      sec <- get_scale(x[[signal_col]], n_scales)
      result <- dfa_function(x, pol = pol, sec = sec)

      dplot <- result$dplot
      alpha <- result$alpha

      fit <- lm(fx ~ sx, data = dplot)
      r2 <- summary(fit)$r.squared

      rows[[idx]] <- data.frame(
        pol = pol,
        n_scales = n_scales,
        unique_scales = length(sec),
        alpha = alpha,
        r2 = r2
      )
      idx <- idx + 1
    }
  }

  results <- do.call(rbind, rows)
  results <- results[order(-results$r2), ]
  rownames(results) <- NULL

  best <- results[1, ]
  return(list(results = results, best = best))
}
