
data.parse <- function(formula, data, subset = NULL, na.rm = FALSE){
  
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
  
  snp <- setdiff(colnames(alt.data), colnames(null.data))
  
  null.G<-null.G[!is.na(null.G[, resp.G]), ,drop=F]
  alt.G<-alt.G[!is.na(alt.G[, resp.G]), ,drop=F]
  
  n.B<-nrow(null.B)
  n.G<-nrow(null.G)
  n<-n.B
  
  X.B<-as.matrix(null.B[, covar, drop=F])
  X.G<-as.matrix(null.G[, covar, drop=F])
  
  G.B<-as.matrix(alt.B[, snp, drop=F])
  G.G<-as.matrix(alt.G[, snp, drop=F])
  
  
  ########
  
  if(na.rm){
    message("We suggest and set na.rm as FALSE")
    na.rm <- FALSE
  }
  
  if(sum(is.na(G.B)) + sum(is.na(G.G)) != 0){
    message("Missing genotypes detected")
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
      id.B <- as.vector(which(obs.id.B[,k]))
      mdl0.B <- glm(paste0(resp.B, " ~ .-1"), data = null.B[id.B,,drop=FALSE], family="binomial")
      y.hat<-as.vector(mdl0.B$fitted.values)
      res.B <- null.B[id.B, resp.B] - y.hat
      A<-y.hat * (1-y.hat)
      I.ab.B[, k] <- t(X.B[id.B,,drop=FALSE]) %*% (A * G.B[id.B, k])/length(y.hat)
      V.B[k, k] <- t(G.B[id.B, k])%*% (A * G.B[id.B, k])/length(y.hat)  - t(I.ab.B[, k]) %*% inv.I.aa.B %*% I.ab.B[, k]
      S.B[k] <- sum(G.B[id.B,k] * res.B)
    }
    
    for(k in 1:(nsnp-1)){
      for(l in (k+1):nsnp){
        id.B <- as.vector(which(obs.id.B[, k] & obs.id.B[, l]))
        mdl0.B <- glm(paste0(resp.B, " ~ .-1"), data = null.B[id.B,,drop=FALSE], family="binomial")
        y.hat<-mdl0.B$fitted.values
        A<-y.hat * (1-y.hat)
        V.B[k, l] <- t(G.B[id.B, k]) %*% (A * G.B[id.B, l])/length(y.hat) - t(I.ab.B[, k]) %*% inv.I.aa.B %*% I.ab.B[, l] #much faster
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
    res <- mdl0.G$residuals
    s2<-sum(res^2)/(length(res)-nx)
    inv.I.aa.G <- solve(t(X.G[,,drop=FALSE]) %*% X.G[,,drop=FALSE] / s2/n.G)
    
    I.ab.G <- matrix(NA, nx, nsnp)
    V.G <- matrix(NA, nsnp, nsnp)
    S.G <- rep(NA, nsnp)
    for(k in 1:nsnp){
      id.G <- as.vector(which(obs.id.G[,k]))
      mdl0.G <- lm(paste0(resp.G, " ~ .-1"), data = null.G[id.G,,drop=FALSE])
      res <- mdl0.G$residuals
      s2<-sum(res^2)/(length(res)-nx)
      I.ab.G[, k] <- t(X.G[id.G,,drop=FALSE]) %*% G.G[id.G, k] / s2/length(res)
      V.G[k, k] <- t(G.G[id.G, k])%*% G.G[id.G, k] / s2/length(id.G)  - t(I.ab.G[, k]) %*% inv.I.aa.G %*% I.ab.G[, k]
      S.G[k] <- sum(G.G[id.G,k] * res)/s2
    }
    
    
    for(k in 1:(nsnp-1)){
      for(l in (k+1):nsnp){
        id.G <- as.vector(which(obs.id.G[, k] & obs.id.G[, l]))
        mdl0.G <- lm(paste0(resp.G, " ~ .-1"), data = null.G[id.G,, drop=FALSE ])
        res <- mdl0.G$residuals
        s2<-sum(res^2)/(length(res)-nx)
        V.G[k, l] <- t(G.G[id.G, k])%*% G.G[id.G, l] / s2/length(res) - t(I.ab.G[, k]) %*% inv.I.aa.G %*% I.ab.G[, l] #much faster
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
    res.B<-null.B[, resp.B] - y.hat
    A<-y.hat * (1-y.hat)
    V.B<-t(G.B)%*% (A * G.B) - t(G.B)%*% (A * X.B) %*% solve(t(X.B) %*% (A * X.B)) %*% t(X.B) %*% (A * G.B) #much faster
    
    res.G<-mdl0.G$residuals
    s2<-sum(res.G^2)/(length(res.G)-ncol(X.G))
    V.G<-(t(G.G) %*% G.G - t(G.G) %*% X.G %*% solve(t(X.G) %*% X.G) %*% t(X.G) %*% G.G)/s2
    
    S.B<-as.vector(t(G.B)%*%res.B)
    S.G<-as.vector(t(G.G)%*%res.G/s2)
    
  }
  
  S.B <- S.B/sqrt(n)
  S.G <- S.G/sqrt(n)
  V.B <- V.B/n
  V.G <- V.G/n
  
  ini <- list(data = data, S.B = S.B, S.G = S.G, V.B = V.B, V.G = V.G)
  class(ini) <- "ini"
  ini
  
}
