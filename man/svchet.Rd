\name{svchet}
\alias{svchet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
svchet(formula, data, subset = NULL, kappa = seq(0, 1, length.out = 101))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
  \item{kappa}{
%%     ~~Describe \code{kappa} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data, subset = NULL, kappa = seq(0, 1, length.out = 101)) 
{
    formula <- Formula(formula)
    if (is.null(subset)) {
        subset <- 1:nrow(data)
    }
    else {
        subset <- sort(intersect(subset, 1:nrow(data)))
    }
    data <- data[subset, ]
    mf <- model.frame(formula, na.action = na.pass, data = data, 
        lhs = 1:2)
    null.B <- model.part(formula, mf, lhs = 1, rhs = 1, drop = F)
    null.G <- model.part(formula, mf, lhs = 2, rhs = 1, drop = F)
    alt.B <- model.part(formula, mf, lhs = 1, rhs = 1:2, drop = F)
    alt.G <- model.part(formula, mf, lhs = 2, rhs = 1:2, drop = F)
    resp.B <- colnames(null.B)[1]
    covar.B <- colnames(null.B)[-1]
    snp.B <- setdiff(colnames(alt.B), c(resp.B, covar.B))
    resp.G <- colnames(null.G)[1]
    covar.G <- colnames(null.G)[-1]
    snp.G <- setdiff(colnames(alt.G), c(resp.G, covar.G))
    null.G <- null.G[!is.na(null.G[, resp.G]), , drop = F]
    alt.G <- alt.G[!is.na(alt.G[, resp.G]), , drop = F]
    n.B <- nrow(null.B)
    n.G <- nrow(null.G)
    n <- n.B
    X.B <- as.matrix(cbind(Intercept = rep(1, n.B), null.B[, 
        covar.B, drop = F]))
    X.G <- as.matrix(cbind(Intercept = rep(1, n.G), null.G[, 
        covar.G, drop = F]))
    G.B <- as.matrix(alt.B[, snp.B, drop = F])
    G.G <- as.matrix(alt.G[, snp.G, drop = F])
    mdl0.B <- glm(as.formula(paste(resp.B, "~ .")), data = null.B, 
        family = "binomial")
    mdl1.B <- glm(as.formula(paste(resp.B, "~ .")), data = alt.B, 
        family = "binomial")
    mdl0.G <- glm(as.formula(paste(resp.G, "~ .")), data = null.G, 
        family = "gaussian")
    mdl1.G <- glm(as.formula(paste(resp.G, "~ .")), data = alt.G, 
        family = "gaussian")
    y.hat <- mdl0.B$fitted.values
    res.B <- null.B[, resp.B] - y.hat
    A <- y.hat * (1 - y.hat)
    V.B <- t(G.B) \%*\% (A * G.B) - t(G.B) \%*\% (A * X.B) \%*\% solve(t(X.B) \%*\% 
        (A * X.B)) \%*\% t(X.B) \%*\% (A * G.B)
    V.B <- V.B/n
    lam.B <- eigen(V.B)$values
    res.G <- null.G[, resp.G] - mdl0.G$fitted.values
    s2 <- sum(res.G^2)/(n.G - ncol(X.G))
    V.G <- (t(G.G) \%*\% G.G - t(G.G) \%*\% X.G \%*\% solve(t(X.G) \%*\% 
        X.G) \%*\% t(X.G) \%*\% G.G)/s2
    V.G <- V.G/n
    lam.G <- eigen(V.G)$values
    S.B <- t(G.B) \%*\% res.B/sqrt(n)
    S.G <- t(G.G) \%*\% res.G/s2/sqrt(n)
    stat.B <- t(S.B) \%*\% S.B
    stat.G <- t(S.G) \%*\% S.G
    liu.mod <- function(q, lambda, h = rep(1, length(lambda)), 
        delta = rep(0, length(lambda))) {
        r <- length(lambda)
        if (length(h) != r) 
            stop("lambda and h should have the same length!")
        if (length(delta) != r) 
            stop("lambda and delta should have the same length!")
        c1 <- sum(lambda * h) + sum(lambda * delta)
        c2 <- sum(lambda^2 * h) + 2 * sum(lambda^2 * delta)
        c3 <- sum(lambda^3 * h) + 3 * sum(lambda^3 * delta)
        c4 <- sum(lambda^4 * h) + 4 * sum(lambda^4 * delta)
        s1 <- c3/(c2^(3/2))
        s2 <- c4/c2^2
        muQ <- c1
        sigmaQ <- sqrt(2 * c2)
        tstar <- (q - muQ)/sigmaQ
        if (s1^2 > s2) {
            a <- 1/(s1 - sqrt(s1^2 - s2))
            delta <- s1 * a^3 - a^2
        }
        else {
            a <- 1/sqrt(s2)
            delta <- 0
        }
        l <- a^2 - 2 * delta
        muX <- l + delta
        sigmaX <- sqrt(2) * a
        Qq <- pchisq(tstar * sigmaX + muX, df = l, ncp = delta, 
            lower.tail = FALSE)
        list(Qq = Qq, mu.Q = muQ, sigma.Q = sigmaQ, mu.X = muX, 
            sigma.X = sigmaX, df = l, ncp = delta)
    }
    np.B <- nrow(V.B)
    np.G <- nrow(V.G)
    p.kappa <- NULL
    lam <- list()
    var.comp <- NULL
    for (i in 1:length(kappa)) {
        r <- kappa[i]
        var.comp <- c(var.comp, r * stat.B + (1 - r) * stat.G)
        lam[[i]] <- c(r * lam.B, (1 - r) * lam.G)
        p.kappa <- c(p.kappa, pchisqsum(var.comp[i], rep(1, length(lam[[i]])), 
            lam[[i]], lower.tail = FALSE, method = "saddlepoint"))
    }
    stat <- min(p.kappa)
    if (stat == 1) {
        pval.VC.Het <- 1
    }
    else {
        find.perc <- function(stat, lam, q0) {
            df <- rep(1, length(lam))
            pr <- pchisqsum(q0, df, lam, lower.tail = FALSE, 
                method = "saddlepoint")
            if (pr == stat) {
                return(q0)
            }
            v <- sqrt(2 * sum(lam^2))
            q1 <- NULL
            for (i in 1:100) {
                q1 <- q0 + i * v
                pr <- pchisqsum(q1, df, lam, lower.tail = FALSE, 
                  method = "saddlepoint")
                if (pr < stat) {
                  break
                }
            }
            q <- NULL
            while (abs(q1 - q0) > 1e-06) {
                q <- (q0 + q1)/2
                pr <- pchisqsum(q, df, lam, lower.tail = FALSE, 
                  method = "saddlepoint")
                if (pr < stat) {
                  q1 <- q
                }
                else if (pr > stat) {
                  q0 <- q
                }
                else {
                  return(q)
                }
            }
            return(q)
        }
        perc <- NULL
        for (i in 1:length(lam)) {
            perc <- c(perc, find.perc(stat, lam[[i]], var.comp[i]))
        }
        par.G <- liu.mod(stat.G, lam.G)
        par.B <- liu.mod(stat.B, lam.B)
        fn1 <- function(x, perc, kappa, par.G, par.B) {
            perc <- perc[kappa > 0]
            kappa <- kappa[kappa > 0]
            x.star <- (x - par.G$mu.Q)/par.G$sigma.Q * par.G$sigma.X + 
                par.G$mu.X
            x.star[x.star == 0] <- -1e-12
            pdf.G <- par.G$sigma.X/par.G$sigma.Q * dchisq(x.star, 
                df = par.G$df, ncp = par.G$ncp)
            q <- apply(rep(1, length(x)) \%*\% t(perc/kappa) - 
                x \%*\% t(1/kappa - 1), 1, min)
            cdf.B <- pchisqsum(q, rep(1, length(lam.B)), lam.B, 
                lower.tail = TRUE, method = "saddlepoint")
            fval <- pdf.G * cdf.B
            fval
        }
        fn <- function(x, perc, kappa, par.G, par.B) {
            perc <- perc[kappa > 0]
            kappa <- kappa[kappa > 0]
            x.star <- (x - par.G$mu.Q)/par.G$sigma.Q * par.G$sigma.X + 
                par.G$mu.X
            x.star[x.star == 0] <- -1e-12
            pdf.G <- log(par.G$sigma.X/par.G$sigma.Q) + dchisq(x.star, 
                df = par.G$df, ncp = par.G$ncp, log = TRUE)
            q <- apply(rep(1, length(x)) \%*\% t(perc/kappa) - 
                x \%*\% t(1/kappa - 1), 1, min)
            cdf.B <- pchisq((q - par.B$mu.Q)/par.B$sigma.Q * 
                par.B$sigma.X + par.B$mu.X, df = par.B$df, ncp = par.B$ncp, 
                lower.tail = TRUE, log.p = TRUE)
            fval <- exp(pdf.G + cdf.B)
            fval
        }
        search.interval <- function(perc, kappa, par.G, par.B) {
            x <- seq(0, perc[kappa == 0][1], length.out = 1e+05 + 
                1)
            f <- fn(x, perc, kappa, par.G, par.B)
            cf <- cumsum(f)
            x.lower <- x[head(which(cf > 1e-06), 1)]
            x <- rev(x)
            f <- rev(f)
            cf <- cumsum(f)
            x.upper <- x[head(which(cf > 1e-06), 1)]
            if (is.na(x.lower) || is.nan(x.lower)) {
                x.lower <- 0
            }
            if (is.na(x.upper) || is.nan(x.upper)) {
                x.upper <- perc[kappa == 0][1]
            }
            c(x.lower, x.upper)
        }
        interval <- search.interval(perc, kappa, par.G, par.B)
        interval[1] <- 0
        stat.thr <- 1e-04
        if (stat < stat.thr) {
            try.int <- try(pval.VC.Het <- 1 - integrate(fn1, 
                lower = interval[1], upper = interval[2], perc = perc, 
                kappa = kappa, par.G = par.G, par.B = par.B, 
                subdivisions = 2000, rel.tol = 1e-08)$value, 
                silent = TRUE)
        }
        else {
            try.int <- try(pval.VC.Het <- 1 - integrate(fn, lower = interval[1], 
                upper = interval[2], perc = perc, kappa = kappa, 
                par.G = par.G, par.B = par.B, subdivisions = 2000, 
                rel.tol = 1e-08)$value, silent = TRUE)
        }
        if (class(try.int) == "try-error") {
            mcmc <- TRUE
            c1 <- seq(interval[1], interval[2], length.out = 7)[-c(1, 
                7)]
            ec1 <- NULL
            int1 <- NULL
            for (ii in length(c1):1) {
                if (stat < stat.thr) {
                  t1 <- try(int1 <- integrate(fn1, lower = interval[1], 
                    upper = c1[ii], perc = perc, kappa = kappa, 
                    par.G = par.G, par.B = par.B, subdivisions = 2000, 
                    rel.tol = 1e-08)$value, silent = TRUE)
                }
                else {
                  t1 <- try(int1 <- integrate(fn, lower = interval[1], 
                    upper = c1[ii], perc = perc, kappa = kappa, 
                    par.G = par.G, par.B = par.B, subdivisions = 2000, 
                    rel.tol = 1e-08)$value, silent = TRUE)
                }
                if (class(t1) != "try-error") {
                  ec1 <- c1[ii]
                  break
                }
            }
            if (!is.null(int1)) {
                c2 <- seq(ec1, interval[2], length.out = 6)[-1]
                ec2 <- NULL
                int2 <- NULL
                for (ii in length(c2):1) {
                  if (stat < stat.thr) {
                    t2 <- try(int2 <- integrate(fn1, lower = ec1, 
                      upper = c2[ii], perc = perc, kappa = kappa, 
                      par.G = par.G, par.B = par.B, subdivisions = 2000, 
                      rel.tol = 1e-08)$value, silent = TRUE)
                  }
                  else {
                    t2 <- try(int2 <- integrate(fn, lower = ec1, 
                      upper = c2[ii], perc = perc, kappa = kappa, 
                      par.G = par.G, par.B = par.B, subdivisions = 2000, 
                      rel.tol = 1e-08)$value, silent = TRUE)
                  }
                  if (class(t2) != "try-error") {
                    ec2 <- c2[ii]
                    break
                  }
                }
                if (!is.null(int2)) {
                  if (ec2 < interval[2]) {
                    if (stat < stat.thr) {
                      t3 <- try(int3 <- integrate(fn1, lower = ec2, 
                        upper = interval[2], perc = perc, kappa = kappa, 
                        par.G = par.G, par.B = par.B, subdivisions = 2000, 
                        rel.tol = 1e-08)$value, silent = TRUE)
                    }
                    else {
                      t3 <- try(int3 <- integrate(fn, lower = ec2, 
                        upper = interval[2], perc = perc, kappa = kappa, 
                        par.G = par.G, par.B = par.B, subdivisions = 2000, 
                        rel.tol = 1e-08)$value, silent = TRUE)
                    }
                    if (class(t3) != "try-error") {
                      pval.VC.Het <- 1 - int1 - int2 - int3
                      mcmc <- FALSE
                      method = "int3"
                    }
                    else {
                      mcmc <- TRUE
                      method <- "MCMC"
                    }
                  }
                  else {
                    pval.VC.Het <- 1 - int1 - int2
                    mcmc <- FALSE
                    method = "int2"
                  }
                }
            }
            else {
                mcmc <- TRUE
                method <- "MCMC"
            }
        }
        else {
            mcmc <- FALSE
            method <- "int1"
        }
        if (mcmc) {
            print("Fail to compute the p-value with build-in integration function. Trying MCMC")
            kappa[kappa == 0] <- 1e-12
            pval.VC.Het <- NULL
            NP <- ifelse(stat >= 1e-05, 1, ifelse(stat >= 1e-06, 
                10, 100))
            for (i in 1:NP) {
                nperm <- 1e+05
                nsnp <- nrow(V.G)
                u <- matrix(rnorm(nperm * nsnp), nrow = nperm)
                stat.G.null <- apply(u \%*\% V.G * u, 1, sum)
                rm(u)
                gc()
                upper <- apply(perc/kappa - (1/kappa - 1) \%*\% 
                  t(stat.G.null), 2, min)
                p <- liu.mod(upper, lam.B)$Qq
                pval.VC.Het <- c(pval.VC.Het, mean(p))
                rm(p)
                rm(upper)
                gc()
            }
            pval.VC.Het <- mean(pval.VC.Het)
            method <- "MCMC"
        }
    }
    kappa.opt <- mean(kappa[p.kappa == stat])
    if (kappa.opt <= 1e-12) {
        kappa.opt <- 0
    }
    pval <- c(VC.Het = pval.VC.Het)
    svchet.obj <- list()
    svchet.obj$pval <- pval
    svchet.obj$method <- method
    svchet.obj$kappa.opt <- kappa.opt
    class(svchet.obj) <- "svchet"
    svchet.obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
