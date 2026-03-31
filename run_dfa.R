################## RUN #################################
# Example end-to-end runner for multivariate DFA.
#
# Expected input format:
#   - First column: time index
#   - Remaining columns: one signal per column
#
# This script writes:
#   - alphadfa.txt (single alpha estimate)
#   - dfa_plot.png (log(s) vs log(F(s)) scatter)

library(ggplot2)

source("get_scale_function.R")
source("Multi_DFA.R")
source("plot_dfa.R")

# 1) Load sample data (replace with your own CSV/data.frame if needed)
x <- read.csv("todasvocesRESC.csv")

# 2) Set DFA parameters
pol <- 1               # polynomial detrending order
n_scales <- 20         # number of log-spaced scales
sec <- get_scale(x[[2]], n_scales)

# 3) Run DFA
result <- dfa_function(x, pol, sec)
dplot <- result$dplot
alpha <- result$alpha

# 4) Save alpha
write(alpha, file = "alphadfa.txt")

# 5) Plot and save
plotdfa <- plot_dfa(dplot, title = "Multivariate DFA")
ggsave("dfa_plot.png", plot = plotdfa, width = 7, height = 5)

message(sprintf("Done. Alpha = %.6f", alpha))
