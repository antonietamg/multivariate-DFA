# multivariate-DFA in R

Multivariate DFA (detrended fluctuation analysis) in R, adapted from:
https://github.com/spiralizing/InfoSeries.jl/blob/master/src/func_series.jl

## Repository files

- `Multi_DFA.R`: core `dfa_function(x, pol, sec)` implementation.
- `get_scale_function.R`: helper that builds log-spaced window scales.
- `plot_dfa.R`: helper function to plot DFA output (`log(s)` vs `log(F(s))`).
- `run_dfa.R`: complete runnable example using the sample CSV.
- `helper_dfa_params.R`: helper to compare `(pol, n_scales)` combinations and rank them.
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
