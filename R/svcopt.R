svcopt <-
function(formula, data, subset = NULL, nperm = 1e5, rho=seq(0, 1, length.out=21), kappa=seq(0, 1, length.out=21), seed = 0, nthread=NULL, plot.pval = FALSE){
  
  
  formula<-Formula(formula)
  
  if(is.null(subset)){
    subset <- 1:nrow(data)
  }else{
    subset <- sort(intersect(subset, 1:nrow(data)))
  }
  
  data <- data[subset, ]
  
  mf<-model.frame(formula, na.action = na.pass, data = data, lhs=1:2)
  null.B<-model.part(formula, mf, lhs=1, rhs=1, drop=F)
  null.G<-model.part(formula, mf, lhs=2, rhs=1, drop=F)
  
  alt.B<-model.part(formula, mf, lhs=1, rhs=1:2, drop=F)
  alt.G<-model.part(formula, mf, lhs=2, rhs=1:2, drop=F)
  
  resp.B<-colnames(null.B)[1]
  covar.B<-colnames(null.B)[-1]
  snp.B<-setdiff(colnames(alt.B), c(resp.B, covar.B))
  
  resp.G<-colnames(null.G)[1]
  covar.G<-colnames(null.G)[-1]
  snp.G<-setdiff(colnames(alt.G), c(resp.G, covar.G))
  
  null.G<-null.G[!is.na(null.G[, resp.G]), ,drop=F]
  alt.G<-alt.G[!is.na(alt.G[, resp.G]), ,drop=F]
  
  n.B<-nrow(null.B)
  n.G<-nrow(null.G)
  n<-n.B
  
  X.B<-as.matrix(cbind(Intercept=rep(1,n.B), null.B[, covar.B, drop=F]))
  X.G<-as.matrix(cbind(Intercept=rep(1,n.G), null.G[, covar.G, drop=F]))
  
  G.B<-as.matrix(alt.B[, snp.B, drop=F])
  G.G<-as.matrix(alt.G[, snp.G, drop=F])
  
  mdl0.B<-glm(as.formula(paste(resp.B, "~ .")), data=null.B, family="binomial")
  mdl1.B<-glm(as.formula(paste(resp.B, "~ .")), data=alt.B, family="binomial")
  
  mdl0.G<-glm(as.formula(paste(resp.G, "~ .")), data=null.G, family="gaussian")
  mdl1.G<-glm(as.formula(paste(resp.G, "~ .")), data=alt.G, family="gaussian")
  
  
  y.hat<-mdl0.B$fitted.values
  res.B<-null.B[, resp.B] - y.hat
  #A<-diag(y.hat * (1-y.hat))
  #V.B<-t(G.B)%*%A%*%G.B - t(G.B)%*%A %*% X.B %*% solve(t(X.B) %*% A %*% X.B) %*% t(X.B) %*% A%*%G.B #much faster
  A<-y.hat * (1-y.hat)
  V.B<-t(G.B)%*% (A * G.B) - t(G.B)%*% (A * X.B) %*% solve(t(X.B) %*% (A * X.B)) %*% t(X.B) %*% (A * G.B) #much faster
  V.B<-V.B/n.B #standarize
  ei.B <- eigen(V.B)
  lam.B<-ei.B$values
  
  res.G<-null.G[, resp.G]-mdl0.G$fitted.values
  s2<-sum(res.G^2)/(n.G-ncol(X.G))
  V.G<-(t(G.G) %*% G.G - t(G.G) %*% X.G %*% solve(t(X.G) %*% X.G) %*% t(X.G) %*% G.G)/s2
  V.G<-V.G/n.G #standarize
  ei.G <-eigen(V.G)
  lam.G<-ei.G$values
  
  S.B<-t(G.B)%*%res.B/sqrt(n.B)
  S.G<-t(G.G)%*%res.G/s2/sqrt(n.G)
  
  stat.B<-t(S.B)%*%S.B
  stat.G<-t(S.G)%*%S.G
  stat.BG<-t(S.B)%*%S.G
  
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
  tmp <- .C("eval_pval_opt", as.double(as.vector(S.B)), as.double(as.vector(S.G)), 
            as.double(as.vector(sqrt.V.B)), as.double(as.vector(sqrt.V.G)), 
            as.double(rho), as.double(kappa), 
            as.integer(ncol(G.B)), as.integer(nperm), 
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
    stop("Invalid tag of refine, please debug")
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
  #svcopt.obj$stat <- stat
  svcopt.obj$refine <- ifelse(refine == 1, TRUE, FALSE)
  class(svcopt.obj) <- "svcopt"
  
  if(plot.pval){
    plot(svcopt.obj)
  }
  
  svcopt.obj
  
}


