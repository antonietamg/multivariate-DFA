# Utilities to generate colored-noise test signals for DFA validation.
#
# Convention: power spectral density S(f) ~ 1 / f^beta
#   beta =  0 -> white
#   beta =  1 -> pink
#   beta =  2 -> brown
#   beta = -1 -> blue

normalize_series <- function(x) {
  x <- as.numeric(x)
  x <- x - mean(x)
  s <- sd(x)
  if (is.na(s) || s == 0) return(x)
  x / s
}

colored_noise_fft <- function(n, beta) {
  if (n < 8) stop("n must be >= 8")

  # Frequency-domain shaping. Amplitude scales as f^(-beta/2).
  freqs <- 1:floor(n / 2)
  amp <- freqs ^ (-beta / 2)

  # Random phases in [0, 2*pi).
  phases <- runif(length(freqs), min = 0, max = 2 * pi)

  # Build full Hermitian spectrum to obtain real-valued signal.
  spectrum <- complex(length = n)
  spectrum[1] <- 0 + 0i

  pos <- amp * exp(1i * phases)
  spectrum[freqs + 1] <- pos

  if (n %% 2 == 0) {
    nyquist_idx <- n / 2 + 1
    spectrum[nyquist_idx] <- Re(spectrum[nyquist_idx]) + 0i
    mirror_src <- 2:(nyquist_idx - 1)
  } else {
    mirror_src <- 2:(floor(n / 2) + 1)
  }

  mirror_tgt <- n: (n - length(mirror_src) + 1)
  spectrum[mirror_tgt] <- Conj(spectrum[mirror_src])

  x <- Re(stats::fft(spectrum, inverse = TRUE) / n)
  normalize_series(x)
}

generate_colored_noise <- function(n, color = c("white", "pink", "brown", "blue")) {
  color <- match.arg(color)
  beta <- switch(
    color,
    white = 0,
    pink = 1,
    brown = 2,
    blue = -1
  )

  if (color == "white") {
    return(normalize_series(stats::rnorm(n)))
  }

  colored_noise_fft(n, beta)
}

build_colored_noise_dataframe <- function(n = 4096, seed = 123) {
  set.seed(seed)
  data.frame(
    time = seq_len(n),
    white = generate_colored_noise(n, "white"),
    pink = generate_colored_noise(n, "pink"),
    brown = generate_colored_noise(n, "brown"),
    blue = generate_colored_noise(n, "blue")
  )
}
