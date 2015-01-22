svchom <-
function(data, formula = NULL, subset = NULL, na.rm = FALSE){
  
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
  
  svchom.obj <- list()
  svchom.obj$pval <- pval
  class(svchom.obj) <- "svchom"
  svchom.obj
  
}
