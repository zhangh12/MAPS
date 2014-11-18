
# load("./data/data.rda")
# subset <- NULL
# nperm <- 1e5
# rho <- seq(0, 1, length.out=11)
# kappa <- seq(0, 1, length.out=11)
# seed <- 0
# nthread <- NULL
# plot.pval <- FALSE


svcopt <-
function(formula, data, subset = NULL, nperm = 1e5, rho=seq(0, 1, length.out=11), kappa=seq(0, 1, length.out=11), seed = 0, nthread=NULL, plot.pval = FALSE){
  
  if(!("Formula" %in% class(formula))){
    if("formula" %in% class(formula)){
      formula<-Formula(formula)
    }else{
      stop("formula should be of class \"formula\"")
    }
  }
  
  if(is.null(subset)){
    subset <- 1:nrow(data)
  }else{
    subset <- sort(intersect(subset, 1:nrow(data)))
  }
  
  data <- data[subset, ]
  
  mf <- model.frame(formula, na.action = na.pass, data = data, rhs=1:2, lhs=1:2, drop=FALSE)
  null.data <- model.matrix(formula, mf, rhs=1, drop=F)
  alt.data <- model.matrix(formula, mf, rhs=1:2, drop=F)
  
  outcome.B <- model.part(formula, mf, lhs=1, drop=F)
  outcome.G <- model.part(formula, mf, lhs=2, drop=F)
  
  null.B <- data.frame(outcome.B, null.data)
  null.G <- data.frame(outcome.G, null.data)
  
  alt.B <- data.frame(outcome.B, alt.data)
  alt.G <- data.frame(outcome.G, alt.data)
  
  
  resp.B <- colnames(null.B)[1]
  resp.G <- colnames(null.G)[1]
  covar <- colnames(null.B)[-1]
  
  snp.B <- setdiff(colnames(alt.data), colnames(null.data))
  snp.G <- snp.B
  
  null.G<-null.G[!is.na(null.G[, resp.G]), ,drop=F]
  alt.G<-alt.G[!is.na(alt.G[, resp.G]), ,drop=F]
  
  n.B<-nrow(null.B)
  n.G<-nrow(null.G)
  n<-n.B
  
  X.B<-as.matrix(null.B[, covar, drop=F])
  X.G<-as.matrix(null.G[, covar, drop=F])
  
  G.B<-as.matrix(alt.B[, snp.B, drop=F])
  G.G<-as.matrix(alt.G[, snp.G, drop=F])
  
  
  ########
  
  s2.marg <- rep(NA, ncol(G.B))
  
  if(sum(is.na(G.B)) + sum(is.na(G.G)) != 0){
    nsnp <- ncol(G.B)
    nx <- ncol(X.B)
    obs.id.B <- !is.na(G.B)
    suff.n.B <- t(obs.id.B) %*% obs.id.B
    
    mdl0.B <- glm(paste0(resp.B, " ~ .-1"), data = null.B, family="binomial")
    y.hat<-mdl0.B$fitted.values
    A<-y.hat * (1-y.hat)
    inv.I.aa.B <- solve(t(X.B[,,drop=FALSE]) %*% (A * X.B[,,drop=FALSE])/length(y.hat))
    
    I.ab.B <- matrix(NA, nx, nsnp)
    V.B <- matrix(NA, nsnp, nsnp)
    S.B <- rep(NA, nsnp)
    for(k in 1:nsnp){
      id.B <- which(obs.id.B[,k])
      mdl0.B <- glm(paste0(resp.B, " ~ .-1"), data = null.B, subset = id.B, family="binomial")
      y.hat<-mdl0.B$fitted.values
      A<-y.hat * (1-y.hat)
      I.ab.B[, k] <- t(X.B[id.B,,drop=FALSE]) %*% (A * G.B[id.B, k])/length(y.hat)
      V.B[k, k] <- t(G.B[id.B, k])%*% (A * G.B[id.B, k])/length(y.hat)  - t(I.ab.B[, k]) %*% inv.I.aa.B %*% I.ab.B[, k]
      S.B[k] <- sum(G.B[id.B,k] * mdl0.B$residuals)
    }
    
    
    for(k in 1:(nsnp-1)){
      for(l in (k+1):nsnp){
        id.B <- which(obs.id.B[, k] & obs.id.B[, l])
        mdl0.B <- glm(paste0(resp.B, " ~ .-1"), data = null.B, subset = id.B, family="binomial")
        y.hat<-mdl0.B$fitted.values
        A<-y.hat * (1-y.hat)
        V.B[k, l] <- t(G.B[id.B, k])%*% (A * G.B[id.B, l])/length(y.hat) - t(I.ab.B[, k]) %*% inv.I.aa.B %*% I.ab.B[, l] #much faster
        V.B[l, k] <- V.B[k, l]
      }
    }
    
    V.B <- V.B * suff.n.B
    
    ########
    
    nsnp <- ncol(G.G)
    nx <- ncol(X.G)
    obs.id.G <- !is.na(G.G)
    suff.n.G <- t(obs.id.G) %*% obs.id.G
    
    mdl0.G <- lm(paste0(resp.G, " ~ .-1"), data = null.G)
    s2<-sum(mdl0.G$residuals^2)/(n.G-nx)
    inv.I.aa.G <- solve(t(X.G[,,drop=FALSE]) %*% X.G[,,drop=FALSE] / s2/n.G)
    
    I.ab.G <- matrix(NA, nx, nsnp)
    V.G <- matrix(NA, nsnp, nsnp)
    S.G <- rep(NA, nsnp)
    for(k in 1:nsnp){
      id.G <- which(obs.id.G[,k])
      mdl0.G <- lm(paste0(resp.G, " ~ .-1"), data = null.G, subset = id.G)
      s2<-sum(mdl0.G$residuals^2)/(n.G-nx)
      s2.marg[k] <- s2
      print(s2)
      I.ab.G[, k] <- t(X.G[id.G,,drop=FALSE]) %*% G.G[id.G, k] / s2/length(id.G)
      V.G[k, k] <- t(G.G[id.G, k])%*% G.G[id.G, k] / s2/length(id.G)  - t(I.ab.G[, k]) %*% inv.I.aa.G %*% I.ab.G[, k]
      S.G[k] <- sum(G.G[id.G,k] * mdl0.G$residuals)/s2
    }
    
    
    for(k in 1:(nsnp-1)){
      for(l in (k+1):nsnp){
        id.G <- which(obs.id.G[, k] & obs.id.G[, l])
        mdl0.G <- lm(paste0(resp.G, " ~ .-1"), data = null.G, subset = id.G)
        s2<-sum(mdl0.G$residuals^2)/(n.G-nx)
        print(s2)
        V.G[k, l] <- t(G.G[id.G, k])%*% G.G[id.G, l] / s2/length(id.G) - t(I.ab.G[, k]) %*% inv.I.aa.G %*% I.ab.G[, l] #much faster
        V.G[l, k] <- V.G[k, l]
      }
    }
    
    V.G <- V.G * suff.n.G
  }else{
    
    mdl0.B<-glm(as.formula(paste(resp.B, "~ .-1")), data=null.B, family="binomial")
    mdl1.B<-glm(as.formula(paste(resp.B, "~ .-1")), data=alt.B, family="binomial")
    
    mdl0.G<-glm(as.formula(paste(resp.G, "~ .-1")), data=null.G, family="gaussian")
    mdl1.G<-glm(as.formula(paste(resp.G, "~ .-1")), data=alt.G, family="gaussian")
    
    
    y.hat<-mdl0.B$fitted.values
    res.B<-mdl0.B$residuals
    A<-y.hat * (1-y.hat)
    V.B<-t(G.B)%*% (A * G.B) - t(G.B)%*% (A * X.B) %*% solve(t(X.B) %*% (A * X.B)) %*% t(X.B) %*% (A * G.B) #much faster
    #V.B<-V.B/n.B #standarize
    
    res.G<-mdl0.G$residuals
    s2<-sum(res.G^2)/(n.G-ncol(X.G))
    V.G<-(t(G.G) %*% G.G - t(G.G) %*% X.G %*% solve(t(X.G) %*% X.G) %*% t(X.G) %*% G.G)/s2
    #V.G<-V.G/n.G #standarize
    
    S.B<-as.vector(t(G.B)%*%res.B)
    S.G<-as.vector(t(G.G)%*%res.G/s2)
    
  }
  
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


