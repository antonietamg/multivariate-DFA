# Plot helper for DFA output data.
# Input dplot must contain:
#   - sx: log10(scale)
#   - fx: log10(fluctuation)

plot_dfa <- function(dplot, title = "DFA plot") {
  ggplot(dplot, aes(x = sx, y = fx)) +
    geom_point(col = "firebrick2") +
    geom_point(shape = 1, colour = "black") +
    labs(title = title, x = "log(s)", y = "log(F(s))") +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
}
