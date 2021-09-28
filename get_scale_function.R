get_scale <- function(x,y){   # x = data,  y = number of points I want
  n = 4                    
  tam = length(x)
  tmax = trunc((tam/n), digits=0)
  tmin = 10
  window.size.range = c(tmin, tmax); npoints = y
  log.window.sizes = seq(log10(window.size.range[[1]]), 
                         log10(window.size.range[[2]]), len = npoints)
  scale = unique(round(10^log.window.sizes))
  return(scale)
}

sec <- get_scale(x, y)

#### Now, you can use "sec" for DFA and MagDFA Functions 
## dplot <- dfa_function(x, pol, sec) 
