wu <-
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
  
  ##################
  
  np.B <- nrow(V.B)
  np.G <- nrow(V.G)
  
  #score tests using information from single trait
  stat.sc.B <- t(S.B) %*% solve(V.B) %*% S.B
  pval.Score.B <- pchisq(stat.sc.B, df = np.B, lower.tail= FALSE)
  
  stat.sc.G <- t(S.G) %*% solve(V.G) %*% S.G
  pval.Score.G <- pchisq(stat.sc.G, df = np.G, lower.tail= FALSE)
  
  #Wu's method, a score test for joint analysis
  V.sc <- matrix(0, nrow=np.B+np.G, ncol=np.B+np.G)
  V.sc[1:np.B, 1:np.B] <- V.B
  V.sc[(np.B+1):(np.B+np.G), (np.B+1):(np.B+np.G)] <- V.G
  score <- c(S.B, S.G)
  stat.sc <- t(score) %*% solve(V.sc) %*% score
  pval.Wu <- pchisq(stat.sc, df=np.B+np.G, lower.tail = FALSE)
  
  
  ############
  
  pval <- c(Wu=pval.Wu, Score.B=pval.Score.B, Score.G=pval.Score.G)
  
  wu.obj <- list()
  wu.obj$pval <- pval
  class(wu.obj) <- "wu"
  wu.obj
  
}
