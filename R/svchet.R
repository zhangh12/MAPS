svchet <-
function(data, formula = NULL, subset = NULL, kappa=seq(0, 1, length.out=21), na.rm = FALSE, seed = 0){
  
  set.seed(seed)
  
  if(!("ini" %in% class(data))){
    if(is.null(formula)){
      stop("formula cannot be NULL if data is not a object of class ini")
    }
    ini <- data.parse(formula, data, subset, na.rm)
  }else{
    ini <- data
  }
  
  S.B <- ini$S.B
  S.G <- ini$S.G
  V.B <- ini$V.B
  V.G <- ini$V.G
  
  ei.B <- eigen(V.B)
  lam.B<-ei.B$values
  ei.G <-eigen(V.G)
  lam.G<-ei.G$values
  
  stat.B<-sum(S.B^2)
  stat.G<-sum(S.G^2)
  
  #this function is modified from liu() in package "CompQuadForm"
  liu.mod <- function(q, lambda, h = rep(1, length(lambda)), delta = rep(0,length(lambda))){
    
    r <- length(lambda)
    if (length(h) != r) 
      stop("lambda and h should have the same length!")
    if (length(delta) != r) 
      stop("lambda and delta should have the same length!")
    c1 <- sum(lambda * h) + sum(lambda * delta)
    c2 <- sum(lambda^2 * h) + 2 * sum(lambda^2 * delta)
    c3 <- sum(lambda^3 * h) + 3 * sum(lambda^3 * delta)
    c4 <- sum(lambda^4 * h) + 4 * sum(lambda^4 * delta)
    s1 <- c3/(c2^(3/2))
    s2 <- c4/c2^2
    muQ <- c1
    sigmaQ <- sqrt(2 * c2)
    tstar <- (q - muQ)/sigmaQ
    if (s1^2 > s2) {
      a <- 1/(s1 - sqrt(s1^2 - s2))
      delta <- s1 * a^3 - a^2
    }
    else {
      a <- 1/sqrt(s2)
      delta <- 0
    }
    
    l <- a^2 - 2*delta
    
    muX <- l + delta
    sigmaX <- sqrt(2) * a
    Qq <- pchisq(tstar * sigmaX + muX, df = l, ncp = delta, lower.tail = FALSE)
    list(Qq = Qq, mu.Q = muQ, sigma.Q = sigmaQ, mu.X=muX, sigma.X=sigmaX, df=l, ncp=delta)
    
  }
  
  
  np.B <- nrow(V.B)
  np.G <- nrow(V.G)
  p.kappa <- NULL
  lam <- list()
  var.comp <- NULL
  #   if(0 %in% kappa){
  #     kappa[kappa==0] <- 1e-12
  #   }
  for(i in 1:length(kappa)){
    r <- kappa[i]
    var.comp <- c(var.comp, r * stat.B + (1 - r) * stat.G)
    lam[[i]] <- c(r*lam.B, (1-r)*lam.G)
    p.kappa <- c(p.kappa, pchisqsum(var.comp[i], rep(1, length(lam[[i]])), lam[[i]], lower.tail=FALSE, method="saddlepoint"))
  }
  
  
  stat <- min(p.kappa)
  
  if(stat == 1.0){
    pval.VC.Het <- 1.0
  }else{
    
    find.perc <- function(stat, lam, q0){
      
      df <- rep(1, length(lam))
      pr <- pchisqsum(q0, df, lam, lower.tail=FALSE, method="saddlepoint")
      if(pr == stat){
        return(q0)
      }
      
      v <- sqrt(2 * sum(lam^2))
      q1 <- NULL
      for(i in 1:100){
        q1 <- q0+i*v
        pr <- pchisqsum(q1, df, lam, lower.tail=FALSE, method="saddlepoint")
        if(pr < stat){
          break
        }
      }
      
      q <- NULL
      while(abs(q1-q0)>1e-6){
        q <- (q0+q1)/2
        pr <- pchisqsum(q, df, lam, lower.tail=FALSE, method="saddlepoint")
        if(pr < stat){
          q1 <- q
        }else if(pr > stat){
          q0 <- q
        }else{
          return(q)
        }
      }
      
      return(q)
      
    }
    perc <- NULL
    for(i in 1:length(lam)){
      perc <- c(perc, find.perc(stat, lam[[i]], var.comp[i]))
    }
    
    par.G <- liu.mod(stat.G, lam.G)
    par.B <- liu.mod(stat.B, lam.B)
    
    #t(S.B)%*%S.B <= min((perc-(1-kappa) * t(S.G) %*% S.G)/kappa)
    
    ### try the build-in integration routine first
    
    fn1 <- function(x, perc, kappa, par.G, par.B){
      
      perc <- perc[kappa>0]
      kappa <- kappa[kappa>0]
      
      x.star <- (x-par.G$mu.Q)/par.G$sigma.Q * par.G$sigma.X + par.G$mu.X
      x.star[x.star == 0] <- -1e-12
      pdf.G <- par.G$sigma.X / par.G$sigma.Q * dchisq(x.star, df = par.G$df, ncp = par.G$ncp)
      #q <- min((perc-(1-kappa) * x)/kappa)
      q <- apply(rep(1,length(x))%*%t(perc/kappa) - x %*% t(1/kappa-1), 1, min)
      #cdf.B <- pchisq((q-par.B$mu.Q)/par.B$sigma.Q * par.B$sigma.X + par.B$mu.X, df = par.B$df, ncp = par.B$ncp, lower.tail=TRUE)
      cdf.B <- pchisqsum(q, rep(1,length(lam.B)), lam.B, lower.tail=TRUE, method="saddlepoint")
      fval <- pdf.G * cdf.B
      
      fval
      
    }
    
    fn <- function(x, perc, kappa, par.G, par.B){
      
      perc <- perc[kappa>0]
      kappa <- kappa[kappa>0]
      
      x.star <- (x-par.G$mu.Q)/par.G$sigma.Q * par.G$sigma.X + par.G$mu.X
      x.star[x.star == 0] <- -1e-12
      pdf.G <- log(par.G$sigma.X / par.G$sigma.Q) + dchisq(x.star, df = par.G$df, ncp = par.G$ncp, log=TRUE)
      #q <- min((perc-(1-kappa) * x)/kappa)
      q <- apply(rep(1,length(x))%*%t(perc/kappa) - x %*% t(1/kappa-1), 1, min)
      cdf.B <- pchisq((q-par.B$mu.Q)/par.B$sigma.Q * par.B$sigma.X + par.B$mu.X, df = par.B$df, ncp = par.B$ncp, lower.tail=TRUE, log.p=TRUE)
      fval <- exp(pdf.G + cdf.B)
      
      fval
      
    }
    
    #this function decides the interval of integration
    search.interval <- function(perc, kappa, par.G, par.B){
      
      x <- seq(0,perc[kappa==0][1],length.out=1e5+1)
      f <- fn(x, perc,kappa,par.G,par.B)
      cf <- cumsum(f)
      x.lower <- x[head(which(cf>1e-6),1)]
      x <- rev(x)
      f <- rev(f)
      cf <- cumsum(f)
      x.upper <- x[head(which(cf>1e-6),1)]
      
      if(is.na(x.lower) || is.nan(x.lower)){
        x.lower <- 0
      }
      if(is.na(x.upper) || is.nan(x.upper)){
        x.upper <- perc[kappa==0][1]
      }
      
      c(x.lower, x.upper)
      
    }
    
    interval <- search.interval(perc, kappa, par.G, par.B)
    interval[1] <- 0
    
    stat.thr <- 1e-4
    if(stat < stat.thr){
      try.int <- try(pval.VC.Het <- 1-
                     integrate(fn1, lower=interval[1], upper=interval[2], perc=perc, kappa=kappa, par.G=par.G, par.B=par.B, subdivisions=2000, rel.tol=1e-8)$value, silent=TRUE)
    }else{
      try.int <- try(pval.VC.Het <- 1-
                     integrate(fn, lower=interval[1], upper=interval[2], perc=perc, kappa=kappa, par.G=par.G, par.B=par.B, subdivisions=2000, rel.tol=1e-8)$value, silent=TRUE)
    }
    
    #class(try.int) <- "try-error"
    
    if(class(try.int) == "try-error"){
      
      mcmc <- TRUE
      c1 <- seq(interval[1], interval[2], length.out=7)[-c(1,7)]
      ec1 <- NULL
      int1 <- NULL
      for(ii in length(c1):1){
        if(stat < stat.thr){
          t1 <- try(int1 <- integrate(fn1, lower=interval[1], upper=c1[ii], perc=perc, kappa=kappa, par.G=par.G, par.B=par.B, subdivisions=2000, rel.tol=1e-8)$value, silent=TRUE)
        }else{
          t1 <- try(int1 <- integrate(fn, lower=interval[1], upper=c1[ii], perc=perc, kappa=kappa, par.G=par.G, par.B=par.B, subdivisions=2000, rel.tol=1e-8)$value, silent=TRUE)
        }
        if(class(t1) != "try-error"){
          ec1 <- c1[ii]
          break
        }
      }
      
      if(!is.null(int1)){
        c2 <- seq(ec1, interval[2], length.out=6)[-1]
        ec2 <- NULL
        int2 <- NULL
        for(ii in length(c2):1){
          if(stat < stat.thr){
            t2 <- try(int2 <- integrate(fn1, lower=ec1, upper=c2[ii], perc=perc, kappa=kappa, par.G=par.G, par.B=par.B, subdivisions=2000, rel.tol=1e-8)$value, silent=TRUE)
          }else{
            t2 <- try(int2 <- integrate(fn, lower=ec1, upper=c2[ii], perc=perc, kappa=kappa, par.G=par.G, par.B=par.B, subdivisions=2000, rel.tol=1e-8)$value, silent=TRUE)
          }
          if(class(t2) != "try-error"){
            ec2 <- c2[ii]
            break
          }
        }
        
        if(!is.null(int2)){
          if(ec2<interval[2]){
            if(stat < stat.thr){
              t3 <- try(int3 <- integrate(fn1, lower=ec2, upper=interval[2], perc=perc, kappa=kappa, par.G=par.G, par.B=par.B, subdivisions=2000, rel.tol=1e-8)$value, silent=TRUE)
            }else{
              t3 <- try(int3 <- integrate(fn, lower=ec2, upper=interval[2], perc=perc, kappa=kappa, par.G=par.G, par.B=par.B, subdivisions=2000, rel.tol=1e-8)$value, silent=TRUE)
            }
            if(class(t3) != "try-error"){
              pval.VC.Het <- 1-int1-int2-int3
              mcmc <- FALSE
              method="int3"
            }else{
              mcmc <- TRUE
              method <- "MCMC"
            }
          }else{
            pval.VC.Het <- 1-int1-int2
            mcmc <- FALSE
            method="int2"
          }
        }
      }else{
        mcmc <- TRUE
        method <- "MCMC"
      }
    }else{
      mcmc <- FALSE
      method <- "int1"
    }
    
    ##### this MCMC method is not accurated as desired
    ##### if integrate() fails, try MCMC
    
    if(mcmc){
      
      print("Fail to compute the p-value with build-in integration function. Trying MCMC")
      kappa[kappa==0] <- 1e-12
      pval.VC.Het <- NULL
      NP <- ifelse(stat >= 1e-5, 1, ifelse(stat >=1e-6, 10, 100))
      for(i in 1:NP){
        nperm <- 1e5
        nsnp <- nrow(V.G)
        u <- matrix(rnorm(nperm * nsnp), nrow = nperm)
        stat.G.null <- apply(u %*% V.G * u, 1, sum)
        rm(u)
        gc()
        
        upper <- apply(perc/kappa - (1/kappa-1) %*% t(stat.G.null), 2, min)
        
        #variance component test for heterogeneity (svchet)
        p <- liu.mod(upper, lam.B)$Qq
        pval.VC.Het <- c(pval.VC.Het, mean(p))
        rm(p)
        rm(upper)
        gc()
      }
      
      pval.VC.Het <- mean(pval.VC.Het)
      method <- "MCMC"
    }
    
  }
  
  kappa.opt <- mean(kappa[p.kappa == stat])
  if(kappa.opt<=1e-12){
    kappa.opt <- 0
  }
  
  
  ############
  
  pval <- c(VC.Het=pval.VC.Het)
  
  svchet.obj <- list()
  svchet.obj$pval <- pval
  svchet.obj$method <- method
  svchet.obj$kappa.opt <- kappa.opt
  class(svchet.obj) <- "svchet"
  svchet.obj
  
}
