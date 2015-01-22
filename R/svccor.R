svccor <-
function(data, formula = NULL, subset = NULL, nperm = 1e5, rho=seq(-1, 1, length.out=21), na.rm = FALSE, seed = 0, nthread = NULL){
  
  obj <- svcopt(data, formula, subset, nperm, rho, kappa = c(.5), na.rm, seed, nthread, plot.pval = FALSE)
  
  svccor.obj <- list()
  svccor.obj$pval <- obj$pval
  names(svccor.obj$pval) <- "VC.Cor"
  svccor.obj$nperm <- nperm
  svccor.obj$rho.opt <- obj$rho.opt
  svccor.obj$obs.rank <- as.vector(obj$obs.rank)
  names(svccor.obj$obs.rank) <- paste0("rho_", 1:length(rho))
  svccor.obj$stat <- obj$stat
  svccor.obj$rho <- rho
  svccor.obj$refine <- ifelse(obj$refine == 1, TRUE, FALSE)
  class(svccor.obj) <- "svccor"
  svccor.obj
  
}
