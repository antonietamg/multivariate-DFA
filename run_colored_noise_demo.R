# End-to-end demo: generate colored noises and estimate DFA alpha for each.

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

rows <- lapply(signals, function(sig) {
  x_sig <- data.frame(time = noise_df$time, value = noise_df[[sig]])
  result <- dfa_function(x_sig, pol = pol, sec = sec)
  alpha <- result$alpha

  data.frame(
    signal = sig,
    expected_alpha = expected_alpha[[sig]],
    estimated_alpha = alpha,
    abs_error = abs(alpha - expected_alpha[[sig]])
  )
})

alpha_results <- do.call(rbind, rows)
write.csv(alpha_results, "colored_noise_alpha_results.csv", row.names = FALSE)

print(alpha_results)
message("Saved: colored_noise_test.csv")
message("Saved: colored_noise_alpha_results.csv")
