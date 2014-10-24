

plot.svcopt <- function(obj){
  
  mat <- -log10(obj$obs.rank/(obj$nperm + 1))
  mat0 <- cut(mat, breaks=c(-Inf, seq(0, 7, by = .1), Inf), labels = FALSE)
  mat1 <- matrix(1:length(mat0), nrow(mat), ncol(mat))
  
  col <- heat.colors(71)[mat0]
  
  xlim <- nrow(mat1) * c(-.2, 1.4)
  ylim <- ncol(mat1) * c(-.2, 1.1)
  image(x = 1:nrow(mat1), y = 1:ncol(mat1), mat1, col = col, axes = FALSE, ylim = ylim, xlim = xlim, xlab=NA, ylab=NA)
  text(1:nrow(mat1), 0, rownames(mat), srt = -45, pos = 4)
  text(0.5, 1:ncol(mat1), colnames(mat), pos = 2)
  tempY <- seq(0.5, ncol(mat1) + 0.5, length = 71)
  rect(xleft = nrow(mat1) * 1.1, ybottom = tempY[-length(tempY)], xright = nrow(mat1) * 1.15, ytop = tempY[-1], col = heat.colors(71), border = NA)
  legendIdx <- cut(0:7, breaks=c(-Inf, seq(0, 7, by = .1), Inf), labels = FALSE)
  
  expr <- c(1, expression(10^-1), expression(10^-2), expression(10^-3), expression(10^-4), expression(10^-5), expression(10^-6), expression(""<=10^-7))
  text(nrow(mat1) * 1.15, tempY[legendIdx], expr, pos = 4)
  
  
}


