
vchom <- function(formula, data){
  
  #library("Formula")
  formula <- Formula(formula)
  mf <- model.frame(formula, na.action = na.pass, data = data, lhs=1:2)
  null.B <- model.part(formula, mf, lhs=1, rhs=1, drop=F)
  null.G <- model.part(formula, mf, lhs=2, rhs=1, drop=F)
  
  alt.B <- model.part(formula, mf, lhs=1, rhs=1:2, drop=F)
  alt.G <- model.part(formula, mf, lhs=2, rhs=1:2, drop=F)
  
  resp.B <- colnames(null.B)[1]
  covar.B <- colnames(null.B)[-1]
  snp.B <- setdiff(colnames(alt.B), c(resp.B, covar.B))
  
  resp.G <- colnames(null.G)[1]
  covar.G <- colnames(null.G)[-1]
  snp.G <- setdiff(colnames(alt.G), c(resp.G, covar.G))
  
  null.G <- null.G[!is.na(null.G[, resp.G]), ,drop=F]
  alt.G <- alt.G[!is.na(alt.G[, resp.G]), ,drop=F]
  
  n.B <- nrow(null.B)
  n.G <- nrow(null.G)
  n <- n.B
  
  X.B <- as.matrix(cbind(Intercept=rep(1,n.B), null.B[, covar.B, drop=F]))
  X.G <- as.matrix(cbind(Intercept=rep(1,n.G), null.G[, covar.G, drop=F]))
  
  G.B <- as.matrix(alt.B[, snp.B, drop=F])
  G.G <- as.matrix(alt.G[, snp.G, drop=F])
  
  mdl0.B <- glm(as.formula(paste(resp.B, "~ .")), data=null.B, family="binomial")
  mdl1.B <- glm(as.formula(paste(resp.B, "~ .")), data=alt.B, family="binomial")
  
  mdl0.G <- glm(as.formula(paste(resp.G, "~ .")), data=null.G, family="gaussian")
  mdl1.G <- glm(as.formula(paste(resp.G, "~ .")), data=alt.G, family="gaussian")
  
  
  y.hat <- mdl0.B$fitted.values
  res.B <- null.B[, resp.B] - y.hat
  #A <- diag(y.hat * (1-y.hat))
  #V.B <- t(G.B)%*%A%*%G.B - t(G.B)%*%A %*% X.B %*% solve(t(X.B) %*% A %*% X.B) %*% t(X.B) %*% A%*%G.B #much faster
  A <- y.hat * (1-y.hat)
  V.B <- t(G.B)%*% (A * G.B) - t(G.B)%*% (A * X.B) %*% solve(t(X.B) %*% (A * X.B)) %*% t(X.B) %*% (A * G.B) #much faster
  V.B <- V.B/n #standarize
  lam.B <- eigen(V.B)$values
  
  res.G <- null.G[, resp.G]-mdl0.G$fitted.values
  s2 <- sum(res.G^2)/(n.G-ncol(X.G))
  V.G <- (t(G.G) %*% G.G - t(G.G) %*% X.G %*% solve(t(X.G) %*% X.G) %*% t(X.G) %*% G.G)/s2
  V.G <- V.G/n #standarize
  lam.G <- eigen(V.G)$values
  
  S.B <- t(G.B)%*%res.B/sqrt(n)
  S.G <- t(G.G)%*%res.G/s2/sqrt(n)
  
  stat.B <- t(S.B)%*%S.B
  stat.G <- t(S.G)%*%S.G
  
  #variance component test assuming equal variance components in two traits
  stat.BG <- stat.B + stat.G
  lam.BG <- c(lam.B, lam.G)
  
  pval.VC.Hom <- pchisqsum(stat.BG, rep(1, length(lam.BG)), lam.BG, lower.tail=FALSE, method = "saddlepoint")
  pval.SKAT.B <- pchisqsum(stat.B, rep(1, length(lam.B)), lam.B, lower.tail=FALSE, method = "saddlepoint")
  pval.SKAT.G <- pchisqsum(stat.G, rep(1, length(lam.G)), lam.G, lower.tail=FALSE, method = "saddlepoint")
  
  stat.Fisher <- -2*(log(pval.SKAT.B) + log(pval.SKAT.G))
  pval.Fisher <- pchisq(stat.Fisher, df=4, lower.tail=FALSE)
  
  stat.minp <- min(pval.SKAT.B, pval.SKAT.G)
  pval.minp <- 1-(1-stat.minp)^2
  
  pval <- c(VC.Hom=pval.VC.Hom, 
            SKAT.B=pval.SKAT.B, SKAT.G=pval.SKAT.G, 
            Fisher=pval.Fisher, minp=pval.minp)
  
  vchom.obj <- list()
  vchom.obj$pval <- pval
  class(vchom.obj) <- "vchom"
  vchom.obj
  
}

