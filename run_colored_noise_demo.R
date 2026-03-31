# End-to-end demo: generate colored noises, estimate DFA alpha, and
# create diagnostics with:
#   - LEFT  panel: log-log power spectrum
#   - RIGHT panel: DFA log-log fluctuation plot

source("colored_noise_utils.R")
source("get_scale_function.R")
source("Multi_DFA.R")

noise_df <- build_colored_noise_dataframe(n = 4096, seed = 123)
write.csv(noise_df, "colored_noise_test.csv", row.names = FALSE)

# Use one shared set of scales for fair comparison across signals.
pol <- 1
n_scales <- 20
sec <- get_scale(noise_df$white, n_scales)

signals <- c("white", "pink", "brown", "blue")
expected_alpha <- c(white = 0.5, pink = 1.0, brown = 1.5, blue = 0.0)

power_spectrum_df <- function(signal) {
  n <- length(signal)
  centered <- signal - mean(signal)
  fft_vals <- stats::fft(centered)

  # Keep positive frequencies, excluding DC.
  idx <- 2:floor(n / 2)
  freq <- (idx - 1) / n
  power <- Mod(fft_vals[idx])^2 / n

  data.frame(freq = freq, power = power)
}

plot_diagnostics <- function(signal_name, signal, dplot, expected_alpha, estimated_alpha) {
  spec <- power_spectrum_df(signal)
  out_file <- sprintf("colored_noise_%s_diagnostic.png", signal_name)
  beta_expected <- c(white = 0, pink = 1, brown = 2, blue = -1)[[signal_name]]

  grDevices::png(out_file, width = 1400, height = 600, res = 140)
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.2, 1.5))

  # LEFT: log-log power spectrum
  graphics::plot(
    log10(spec$freq), log10(spec$power),
    pch = 16, cex = 0.5, col = "steelblue4",
    xlab = "log10(frequency)", ylab = "log10(power)",
    main = sprintf("%s noise - Power spectrum", signal_name)
  )
  # Fit on central spectral band to reduce edge effects.
  lo <- stats::quantile(spec$freq, probs = 0.05)
  hi <- stats::quantile(spec$freq, probs = 0.95)
  spec_mid <- spec[spec$freq >= lo & spec$freq <= hi, , drop = FALSE]
  spec_fit <- stats::lm(log10(power) ~ log10(freq), data = spec_mid)
  graphics::abline(spec_fit, col = "firebrick3", lwd = 2)
  beta_est <- -as.numeric(stats::coef(spec_fit)[2])
  graphics::legend(
    "topright",
    legend = c(
      sprintf("beta_est = %.3f", beta_est),
      sprintf("beta_exp = %.3f", beta_expected)
    ),
    bty = "n", text.col = "firebrick3"
  )

  # RIGHT: DFA log-log fluctuation
  graphics::plot(
    dplot$sx, dplot$fx,
    pch = 16, cex = 0.8, col = "black",
    xlab = "log10(s)", ylab = "log10(F(s))",
    main = sprintf("%s noise - DFA", signal_name)
  )
  dfa_fit <- stats::lm(fx ~ sx, data = dplot)
  graphics::abline(dfa_fit, col = "firebrick3", lwd = 2)
  graphics::legend(
    "topleft",
    legend = c(
      sprintf("alpha_est = %.3f", estimated_alpha),
      sprintf("alpha_exp = %.3f", expected_alpha)
    ),
    bty = "n",
    text.col = c("firebrick3", "gray30")
  )

  grDevices::dev.off()
  out_file
}

rows <- lapply(signals, function(sig) {
  x_sig <- data.frame(time = noise_df$time, value = noise_df[[sig]])
  result <- dfa_function(x_sig, pol = pol, sec = sec)
  alpha <- result$alpha
  dplot <- result$dplot
  fig_file <- plot_diagnostics(
    signal_name = sig,
    signal = noise_df[[sig]],
    dplot = dplot,
    expected_alpha = expected_alpha[[sig]],
    estimated_alpha = alpha
  )

  data.frame(
    signal = sig,
    expected_alpha = expected_alpha[[sig]],
    estimated_alpha = alpha,
    abs_error = abs(alpha - expected_alpha[[sig]]),
    diagnostic_figure = fig_file
  )
})

alpha_results <- do.call(rbind, rows)
write.csv(alpha_results, "colored_noise_alpha_results.csv", row.names = FALSE)

print(alpha_results)
message("Saved: colored_noise_test.csv")
message("Saved: colored_noise_alpha_results.csv")
message("Saved: colored_noise_<color>_diagnostic.png (one per color)")
