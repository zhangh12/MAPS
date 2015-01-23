
maps <- function(data, formula = NULL, subset = NULL, nperm = 1e5, 
                   rho=seq(-1, 1, length.out=21), kappa=seq(0, 1, length.out=21), 
                   na.rm = FALSE, seed = 0, nthread=NULL, plot.pval = FALSE){
  
  if(length(rho) == 1 && rho == 0 && length(kappa) == 1 && kappa == .5){ #svchom
    rho <- NULL
    kappa <- NULL
    nperm <- NULL
    seed <- NULL
    nthread <- NULL
    plot.pval <- NULL
    maps.obj <- svchom(data, formula, subset, na.rm)
  }else if(length(rho) == 1 && rho == 0 && !is.null(kappa)){ #svchet
    rho <- NULL
    nperm <- NULL
    nthread <- NULL
    plot.pval <- NULL
    maps.obj <- svchet(data, formula, subset, kappa, na.rm, seed)
  }else if(length(kappa) == 1 && kappa == .5 && !is.null(rho)){ #svccor
    kappa <- NULL
    plot.pval <- NULL
    maps.obj <- svccor(data, formula, subset, nperm, rho, na.rm, seed, nthread)
  }else if(length(rho) >= 1 && length(kappa) >= 1){ #svcopt
    maps.obj <- svcopt(data, formula, subset, nperm, rho, kappa, na.rm, seed, nthread, plot.pval)
  }else{
    stop("kappa or rho is not configured correctly")
  }
  
  class(maps.obj) <- c("maps", class(maps.obj))
  maps.obj
  
}





