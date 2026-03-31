# Build log-spaced DFA window sizes.
# x: numeric vector from a representative signal column (or any vector with desired length)
# y: number of scale points desired
get_scale <- function(x, y){
  n <- 4
  tam <- length(x)
  tmax <- trunc((tam / n), digits = 0)
  tmin <- 10

  window.size.range <- c(tmin, tmax)
  npoints <- y

  log.window.sizes <- seq(
    log10(window.size.range[[1]]),
    log10(window.size.range[[2]]),
    len = npoints
  )

  scale <- unique(round(10 ^ log.window.sizes))
  return(scale)
}
