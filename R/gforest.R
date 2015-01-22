gforest <- 
function(formula, data, subset = NULL, pdf.file = NULL, main = NULL, width = NULL, height = NULL, ...){
  
  formula<-Formula(formula)
  
  if(is.null(subset)){
    subset <- 1:nrow(data)
  }else{
    subset <- sort(intersect(subset, 1:nrow(data)))
  }
  
  data <- data[subset, ]
  
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
  
  n.B <- nrow(null.B)
  n.G <- nrow(null.G)
  n <- n.B
  
  or.B <- NULL
  coef.G <- NULL
  ci.B <- NULL
  ci.G <- NULL
  pv.B <- NULL
  pv.G <- NULL
  for(rs in snp.G){
    subdata.B <- data[, c(resp.B, covar.B, rs)]
    subdata.G <- data[, c(resp.G, covar.G, rs)]
    mdl.B <- glm(as.formula(paste(resp.B, "~ .")), data=subdata.B, family="binomial")
    mdl.G <- lm(as.formula(paste(resp.G, "~ .")), data=subdata.G)
    
    ci.B <- rbind(ci.B, exp(confint.default(mdl.B, rs)))
    ci.G <- rbind(ci.G, confint.default(mdl.G, rs))
    
    coef.G <- c(coef.G, coef(mdl.G)[rs])
    or.B <- c(or.B, exp(coef(mdl.B)[rs]))
    
    pv.B <- c(pv.B, summary(mdl.B)$coefficients[rs, "Pr(>|z|)"])
    pv.G <- c(pv.G, summary(mdl.G)$coefficients[rs, "Pr(>|t|)"])
    
  }
  
  names(pv.B)<-snp.B
  names(pv.G)<-snp.G
  lci.B <- ci.B[,1]
  rci.B <- ci.B[,2]
  lci.G <- ci.G[,1]
  rci.G <- ci.G[,2]
  
  mg <- coef.G
  mb <- or.B
  
  lg <- lci.G
  ug <- rci.G
  lb <- lci.B
  ub <- rci.B
  
  min.G <- min(lg) - (max(ug)-min(lg)) * 4
  max.G <- max(ug) + (max(ug)-min(lg)) * 3
  
  if(!is.null(pdf.file)){
    rr <- length(snp.G)/20
    if(is.null(height) && is.null(width)){
      rr <- length(snp.G)/20
      height <- max(as.integer(6 * rr), 6)
      width <- 2 * height
    }else if(is.null(height) && !is.null(width)){
      height <- as.integer(width/2)
    }else if(!is.null(height) && is.null(width)){
      width <- as.integer(height*2)
    }else{
      # do nothing
    }
    
    height <- max(height, 6)
    width <- 2 * height
    
    pdf(pdf.file, width = width, height = height, ...)
  }
  
  par(mfrow=c(1,2), cex = .9, cex.lab=.9)
  
  forest.default(mg, ci.lb=lg, ci.ub=ug, refline =0,slab=snp.G, 
         xlab=expression(theta), xlim = c(min.G, max.G))
  
  text(min.G, 1.5+length(snp.G), labels="SNP", pos=4)
  text(max.G, 1.5+length(snp.G), labels=expression(theta * "  [ 95% CI ]"), pos=2)
  text(min(lg) - (max(ug)-min(lg)) * 2, 1.5+length(snp.G), labels="P-Value", pos=4)
  text(min(lg) - (max(ug)-min(lg)) * 2, length(snp.G):1, labels = signif(pv.G,2), pos=4)
  
  min.B <- min(lb) - (max(ub)-min(lb)) * 4
  max.B <- max(ub) + (max(ub)-min(lb)) * 3
  forest.default(mb, ci.lb=lb, ci.ub =ub, refline=1,slab=snp.B, 
         xlab=expression(e^beta), xlim = c(min.B, max.B))
  text(min.B, 1.5+length(snp.B), labels="SNP", pos=4)
  text(max.B, 1.5+length(snp.B), labels=expression(e^beta * "  [ 95% CI ]"), pos=2)
  text(min(lb) - (max(ub)-min(lb)) * 2, 1.5+length(snp.B), labels="P-Value", pos=4)
  text(min(lb) - (max(ub)-min(lb)) * 2, length(snp.B):1, labels = signif(pv.B,2), pos=4)
  
  if(!is.null(main)){
    mtext(main, 3, -4, TRUE)
    #text((min.G + max.G)/2, 3 + length(snp.G), labels = main, pos=4)
  }
  
  if(!is.null(pdf.file)){
    dev.off()
  }
  
}

