
svcopt <-
function(data, formula = NULL, subset = NULL, nperm = 1e5, rho=seq(-1, 1, length.out=21), kappa=seq(0, 1, length.out=21), na.rm = FALSE, seed = 0, nthread=NULL, plot.pval = FALSE){
  
  if(!("ini" %in% class(data))){
    if(is.null(formula)){
      stop("formula cannot be NULL if data is not a object of class ini")
    }
    ini <- data.parse(formula, data, subset, na.rm)
  }else{
    ini <- data
  }
  
  #data <- ini$data
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
  stat.BG<-sum(S.B * S.G)
  
  sqrt.V.B <- ei.B$vectors %*% diag(sqrt(abs(ei.B$values))) %*% t(ei.B$vectors)
  sqrt.V.G <- ei.G$vectors %*% diag(sqrt(abs(ei.G$values))) %*% t(ei.G$vectors)
  
  
  if(is.null(nthread)){
    nthread <- 0
  }
  
  if(is.null(seed)){
    seed <- 0
  }
  
  pval <- -1.0
  obs.rank <- rep(-1, length(rho) * length(kappa))
  stat <- rep(-1, nperm+1)
  refine <- -1
  #dyn.load("/home/zhangh12/vc/code/evalp.so")
  tmp <- .C("eval_pval_opt", as.double(S.B), as.double(S.G), 
            as.double(as.vector(sqrt.V.B)), as.double(as.vector(sqrt.V.G)), 
            as.double(rho), as.double(kappa), 
            as.integer(ncol(V.G)), as.integer(nperm), 
            as.integer(length(rho)), as.integer(length(kappa)), 
            as.integer(seed), as.integer(nthread), 
            pval = as.double(pval), obs.rank = as.integer(obs.rank), 
            stat = as.integer(stat), refine = as.integer(refine))
  #dyn.unload("/home/zhangh12/vc/code/evalp.so")
  
  pval<-c(VC.Opt = tmp$pval)
  obs.rank <- tmp$obs.rank
  stat <- tmp$stat
  refine <- tmp$refine
  if(refine == -1){
    stop("Invalid tag of $refine, please debug")
  }
  
  min.id<-which(obs.rank==min(obs.rank))
  k<-0
  kappa.opt <- NULL
  rho.opt <- NULL
  for(kappa0 in kappa){
    for(rho0 in rho){
      k<-k+1
      if(k %in% min.id){
        kappa.opt<-c(kappa.opt, kappa0)
        rho.opt<-c(rho.opt, rho0)
      }
    }
  }
  
  svcopt.obj <- list()
  svcopt.obj$pval <- pval
  svcopt.obj$nperm <- nperm
  svcopt.obj$rho.opt <- mean(rho.opt)
  svcopt.obj$kappa.opt <- mean(kappa.opt)
  svcopt.obj$obs.rank <- matrix(obs.rank, nrow = length(kappa), byrow = TRUE)
  rownames(svcopt.obj$obs.rank) <- paste0("kappa_", 1:length(kappa))
  colnames(svcopt.obj$obs.rank) <- paste0("rho_", 1:length(rho))
  svcopt.obj$rho <- rho
  svcopt.obj$kappa <- kappa
  svcopt.obj$stat <- stat
  svcopt.obj$refine <- ifelse(refine == 1, TRUE, FALSE)
  class(svcopt.obj) <- "svcopt"
  
  if(plot.pval){
    plot(svcopt.obj)
  }
  
  svcopt.obj
  
}


