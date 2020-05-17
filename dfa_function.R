ibrary(OneR)
library(phonTools)
library(pracma)
library(dplyr)
library(ggplot2)
library(qdap)
library(lessR)
library(compiler)
library(batchtools)


# FUNCTION


dfa_function <- function(x){
  
  #profile series
  
  notas <- x[,2:ncol(x)]
  temp <- x[,1]
  
  integrate.serie <- function(x) {       #integrate series
    p <- cumsum((x)-mean(x))
    return(p)
  }
  integrate.serie <- cmpfun(integrate.serie)
  
  serie_int <- lapply(notas, function(x) integrate.serie(x))
  
  
  #define variables
  
  variables.dfa <- function(x, y){
    
    #x = notas
    #y = temp
    
    tam = nrow(x)            
    
    n = 4                     #min value boxes - how many scales i'll get
    
    tmax = trunc((tam/n), digits=0)
    tmin = 6
    
    si = (tmax-tmin+1)
    fu = (tmax-tmin+1)
    N = tam
    nu=0
    
    dats <- data.frame(tmin:tmax)
    ndats <- length(dats$tmin.tmax)
    dats <- as.numeric(dats$tmin.tmax)
    
    for (s in 1:ndats){
      nu[s] = (trunc(tam/dats[s])-1) # computes how many chunks of "s" elements I can get. 
      Fn = tam
    }
    my_list <- list("tmax"=tmax, "tmin"=tmin, "N"=N, "dats"=dats, "ndats"=ndats, "Fn"=Fn, "si"=si, "fu"=fu, "nu"=nu, "s"=s )
    return(my_list)
  }
  variables.dfa <- cmpfun(variables.dfa)
  variables <- variables.dfa(x=notas, y=temp)
  
  #split series in chunks of n size and do polyfit
  
  blockf  <- function(w, x, s = s, k = k){
    
    block = x[(k * w[s] + 1): (k *  w[s] +  w[s])]
  }
  blockf <- cmpfun(blockf)
  
  blockft  <- function(w, x, s = s, k = k){
    
    blockt = y[(k * w[s] + 1): (k *  w[s] +  w[s])]
  }
  blockft <- cmpfun(blockft)
  
  dfa_fluctuation <- function(x){  
    
    Fs <- list()
    Fss =numeric()
    
    for (j in 1: length(variables[["dats"]])){
      q <- c(1:length(variables[["dats"]]))
      nu = variables[["nu"]]
      knu=nu[j]
      y = temp
      z = knu
      w = variables[["dats"]]
      s = q[j]
      
      
      dfb <- list()
      dfbt <- list()
      
      for (k in 0:z){
        bl <- blockf(w, x, s = s, k = k)
        nam = sprintf("B%s_%s_%s", 1, w[s], k)
        dfb[[nam]] <- bl
        
      }
      
      for (k in 0:z){
        blt <- blockf(w, y, s = s, k = k)
        nam = sprintf("B%s_%s_%s", 1, w[s], k)
        dfbt[[nam]] <- blt
        
      }
      
      #POLYNOMIAL ADJUSTMENT & SUBSTRACTION
      
      coef1 <- mapply(x = dfbt, y = dfb, FUN = function(x, y) polyfit(x, y, 1)) # el numero es el grado del polinomio
      coef1 <- data.frame(coef1)
      coef1 <- as.list(coef1)
      evalp <- mapply(x = coef1, y = dfbt, FUN = function(x, y) polyval(x, y))
      evalp <- data.frame(evalp)
      evalp <- as.list(evalp)
      yevalp <- mapply(x = dfb, y = evalp, FUN = function(x, y) x-y)
      yevalp <- as.vector(yevalp)
      Fnx <- yevalp
      Fss[[j]] <- dot(Fnx, Fnx)
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
      funcntot <- funcntot+ sapply(serie_f, as.numeric)
    }
    funcntot <- funcntot/(l*v)
    Fnn = sqrt(funcntot)
    
    t1 <- vari[["dats"]][1]
    tt <- tail(vari[["dats"]], n=1)
    
    si = vari[["si"]]
    fu = vari[["fu"]]
    
    for (s in t1:tt){
      si[s-vari[["tmin"]]+1] = log10(s)
      fu = log10(Fnn)
    }
    
    #ALPHA
    A = polyfit(si, fu)
    A
    
    my_list <- list("A"=A, "si"=si, "fu"=fu)
    return(my_list)
  }
  
  l = nrow(x)
  
  sifu.alpha.dfa <- cmpfun(sifu.alpha.dfa)
  sifu <- sifu.alpha.dfa(x=funcfn, v=ncol(x), vari=variables, l)
  
  
  
  ################# get the final values
  
  bin_func <- function(x){
    
    destring <- function(x,keep="0-9.-") {
      return( as.numeric(gsub(paste("[^",keep,"]+",sep=""),"",x)) )
    }
    
    si <- x[[2]]
    fu <- x[[3]]
    alpha_val <<- sifu[["A"]]
    tmin = 6
    
    tbinsx = si - log10(tmin)
    tbinsx = tail(tbinsx, n=1)
    
    binsx = round(tbinsx/0.05)
    
    losk <- bin(si, nbins= binsx, method = "clusters") #binning
    
    levels_split <- strsplit(levels(losk), ",")
    destring(levels_split[[1]])
    
    dlev=0
    for (ds in 1:binsx){
      dlev[ds] = destring(levels_split[[ds]])
    }
    
    tryCatch(
      expr = {
        sx = 0
        fx = 0
        for(bn in 1:length(si)){
          sibn = si[bn]
          
          for(nbn in 1:binsx){
            if(sibn > dlev[nbn] & sibn < dlev[nbn+1]){
              sx[nbn] = si[bn]
              fx[nbn] = fu[bn]
            }
            
          }
        }
        
      },
      error = function(e){
        message('Error with NA')
        print(e)
      },
      warning = function(w){
        message('Warning with NA')
        print(w)
      },
      finally = {
        message('All done succesfully... quitting.')
      }
    )
    dplot <- data.frame()
    dplot <- data.frame(cbind(sx, fx))
    return(dplot)
  }
  
  bin_func <- cmpfun(bin_func)
  dplot_serie <- bin_func(sifu)
  return(dplot_serie)
  
}
dfa_function <- cmpfun(dfa_function)
