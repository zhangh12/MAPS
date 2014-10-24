svccor <-
function(formula, data, subset = NULL, nperm = 1e5, rho=seq(0, 1, length.out=101), seed = 0, nthread = NULL){
  
  obj <- svcopt(formula, data, subset, nperm, rho, kappa = c(.5), seed, nthread, plot.pval = FALSE)
  
  svccor.obj <- list()
  svccor.obj$pval <- obj$pval
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
