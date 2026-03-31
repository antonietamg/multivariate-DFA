# multivariate-DFA in R

Multivariate DFA (detrended fluctuation analysis) in R, adapted from:
https://github.com/spiralizing/InfoSeries.jl/blob/master/src/func_series.jl

## Repository files

- `Multi_DFA.R`: core `dfa_function(x, pol, sec)` implementation.
- `get_scale_function.R`: helper that builds log-spaced window scales.
- `plot_dfa.R`: helper function to plot DFA output (`log(s)` vs `log(F(s))`).
- `run_dfa.R`: complete runnable example using the sample CSV.
- `helper_dfa_params.R`: helper to compare `(pol, n_scales)` combinations and rank them.
- `colored_noise_utils.R`: utilities to generate white/pink/brown/blue noise signals.
- `run_colored_noise_demo.R`: demo that generates colored-noise test data and estimates alpha for each noise type.
- `todasvocesRESC.csv`: sample multivariate time-series input.

## Input format

Your input must be a `data.frame` where:

1. Column 1 is time/index.
2. Columns 2..N are signals (one signal per column).

Example header from sample data:

```text
time,V1,V2,...
```

## Quick start

```r
source("get_scale_function.R")
source("Multi_DFA.R")
source("plot_dfa.R")

x <- read.csv("todasvocesRESC.csv")
pol <- 1
sec <- get_scale(x[[2]], 20)

result <- dfa_function(x, pol, sec)
dplot <- result$dplot
alpha <- result$alpha

library(ggplot2)
p <- plot_dfa(dplot, "Multivariate DFA")
print(p)
```

Or just run:

```bash
Rscript run_dfa.R
```

This will create:

- `alphadfa.txt` (estimated alpha)
- `dfa_plot.png` (scatter plot)

## Notes

- `sec` should contain increasing window sizes; `get_scale()` handles this.
- `pol` is polynomial detrending order (common first try: `pol = 1`).

## Optional helper: choose `pol` and `n_scales`

```r
source("get_scale_function.R")
source("Multi_DFA.R")
source("helper_dfa_params.R")

x <- read.csv("todasvocesRESC.csv")
out <- choose_dfa_params(
  x,
  pol_values = c(1, 2, 3),
  n_scales_values = c(12, 20, 30),
  signal_col = 2
)

out$results  # all combinations sorted by R^2
out$best     # top-ranked row
```

If one combination fails (for example, a singular fit inside DFA), the helper now marks that row with `status = "error"` and stores the message in `error` instead of stopping the full search.

## Colored-noise validation demo

You can generate test signals (white, pink, brown, blue noise) and run DFA on each:

```bash
Rscript run_colored_noise_demo.R
```

This script:

- creates `colored_noise_test.csv` with columns `time, white, pink, brown, blue`
- runs DFA on each signal with shared `pol = 1` and shared scales (`n_scales = 20`)
- writes `colored_noise_alpha_results.csv` with expected vs estimated alpha values
- saves one diagnostic figure per color:
  - `colored_noise_white_diagnostic.png`
  - `colored_noise_pink_diagnostic.png`
  - `colored_noise_brown_diagnostic.png`
  - `colored_noise_blue_diagnostic.png`
    (left panel: log-log power spectrum, right panel: DFA log-log plot)
  - spectrum fit now uses the central frequency band (5% to 95%) for a more stable `beta_est`

Expected alpha reference (approximate):

- white noise: ~0.5
- pink noise: ~1.0
- brown noise: ~1.5
- blue noise: ~0.0

Note: blue noise can show larger alpha deviations in finite-length samples because anti-correlated signals are more sensitive to scale range and detrending choices.
Also, colored-noise generation uses random complex Fourier coefficients (not fixed magnitudes), so pink and brown spectra should have distinct slopes but realistic variability between runs.
