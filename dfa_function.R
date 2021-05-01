##################################
##################################
##################################
##################################
# MULTIVARIATE DFA FUNCTION
##################################
##################################
##################################

library(compiler)
library(pracma)

dfa_function <- function(x, pol, sec){ # x = dataframe; pol = polynomial; sec = scale values (for scale use get_scale function)
  
  #profile series
  
  x_series <- data.frame(x[ ,2 : ncol(x)])
  x_time <- as.numeric(x[,1])
  
  integrate.serie <- function(x) {      
    p <- cumsum((x) - mean(x))
    return(p)
  }
  
  integrate.serie <- cmpfun(integrate.serie)
  
  serie_int <- lapply(x_series, function(x) integrate.serie(x))
  
  
  #define variables
  
  variables.dfa <- function(x, y){  # x = x_series; y = x_time
    size = nrow(x)  
      if(is.null(size) == TRUE){
        print("nrow to length")
        size = length(x)
    }
    
  n = 4                     # min number of scales
  tmax = trunc((size / n), digits = 0)
  tmin = 10
  si = (tmax - tmin + 1)
  fu = (tmax - tmin + 1)
  N = size
  nu = 0
  dats <- sec
  ndats <- length(sec)
  dats <- as.numeric(sec)
    
  for (s in 1 : ndats){
    nu[s] = (trunc(size / dats[s]) - 1) 
    Fn = size
  }
    
  my_list <- list("tmax"=tmax, "tmin"=tmin, "N"=N, "dats"=dats, "ndats"=ndats, "Fn"=Fn, "si"=si, "fu"=fu, "nu"=nu, "s"=s )
  return(my_list)
  
  }
  
  variables.dfa <- cmpfun(variables.dfa)
  variables <- variables.dfa(x=x_series, y=x_time)
  
  #split series in chunks of n size and do polyfit
  
  blockf  <- function(w, x, s = s, k = k){
    block = x[(k * w[s] + 1): (k * w[s] + w[s])]
  }
  
  blockf <- cmpfun(blockf)
  
  blockft  <- function(w, y, s = s, k = k){
    blockt = y[(k * w[s] + 1): (k * w[s] + w[s])]
  }
  
  blockft <- cmpfun(blockft)
  
  dfa_fluctuation <- function(x){  
    Fs <- list()
    Fss = numeric()
      for (j in 1 : length(variables[["dats"]])){
        q <- c(1 : length(variables[["dats"]]))
        nu = variables[["nu"]]
        knu = nu[j]
        y = x_time
        z = knu
        w = variables[["dats"]]
        s = q[j]
        dfb <- list()
        dfbt <- list()
        for (k in 0 : z){
          bl <- blockf(w, x, s = s, k = k)
          nam = sprintf("B%s_%s_%s", 1, w[s], k)
          dfb[[nam]] <- bl
        }
        for (k in 0:z){
          blt <- blockf(w, y, s = s, k = k)
          nam = sprintf("B%s_%s_%s", 1, w[s], k)
          dfbt[[nam]] <- blt
        }
    coef1 <- mapply(x = dfbt, y = dfb, FUN = function(x, y) polyfit(x, y, pol)) 
    coef1 <- data.frame(coef1)
    coef1 <- as.list(coef1)
    evalp <- mapply(x = coef1, y = dfbt, FUN = function(x, y) polyval(x, y))
    evalp <- data.frame(evalp)
    evalp <- as.list(evalp)
    k_var <- mapply(x = dfb, y = evalp, FUN = function(x, y) mean(( x - y) ^ 2))
    yevalp <- as.vector(k_var)
    Fnx <- yevalp
    Fss[[j]] <- sum(Fnx)
    Fs <- Fss
    }
    return(Fs)
  }
  
  dfa_fluctuation <- cmpfun(dfa_fluctuation)
  funcfn <- lapply(serie_int, function(x) dfa_fluctuation(x))

  #get fluctuation function value
  
  sifu.alpha.dfa <- function(x,v,vari,l) {
    
    lon = length(x)
    funcntot = 0
    
    for(i in 1:lon){
      serie_f <- x[[i]]
      funcntot <- funcntot + sapply(serie_f, as.numeric)
    }
    
    v=v-1
    funcntot <- funcntot/(vari[["nu"]])
    Fnn = sqrt(funcntot)
    t1 <- vari[["dats"]][1]
    tt <- tail(vari[["dats"]], n=1)
    si = vari[["si"]]
    fu = vari[["fu"]]
    
    for (s in 1:length(sec)){
      si[s] = log10(sec[s])
      fu = log10(Fnn)
    }
    
    A = polyfit(si, fu) # alpha value
    my_list <- list("A"=A, "si"=si, "fu"=fu)
    return(my_list)
   }
  
  l = nrow(x)
  sifu.alpha.dfa <- cmpfun(sifu.alpha.dfa)
  sifu <- sifu.alpha.dfa(x=funcfn, v=ncol(x), vari=variables, l)

  # get the final values
  
  bin_func <- function(x){
     alpha_val <<- x[[1]]
     sx = x[[2]]
     fx = x[[3]]
     dplot <- data.frame()
     dplot <- data.frame(cbind(sx, fx))
     return(dplot)
  }
  
   bin_func <- cmpfun(bin_func)
   dplot_serie <- bin_func(sifu)
   return(dplot_serie)
  }

dfa_function <- cmpfun(dfa_function)

##################################
##################################
##################################
# END DFA FUNCTION
##################################
##################################
##################################
