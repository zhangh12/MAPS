\name{svcopt}
\alias{svcopt}
\title{
Optimal Stratified Variance Component Test
}
\description{
\code{svcopt} is used to perform the stratified variance component test with optimal tuning on parameters of variance component proportion and effect correlation.
}
\usage{
svcopt(formula, data, subset = NULL, nperm = 1e+05, 
       rho = seq(0, 1, length.out = 21), 
       kappa = seq(0, 1, len = 21), na.rm = FALSE, 
       seed = 0, nthread = NULL, plot.pval = FALSE)
}
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
an object of clas "\code{Formula}". The user can specify two outcomes on the left-hand side. The covariates and the variables of interest can be specified on the right-hand side. See 'Details' below.
}
  \item{data}{
%%     ~~Describe \code{data} here~~
a data frame containing all variables specified in \code{formula}.
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
an optional vector specifying a subset of observations to be used in the testing process.
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
an integer specifying the number of replicates used in Monte Carlo test.
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
a vector of effect correlation parameter used in tuning process.
}
  \item{kappa}{
%%     ~~Describe \code{kappa} here~~
a vector of variance component proportion parameter used in tuning process.
}
  \item{na.rm}{
a logical indicating whether to allow individuals with missing values in variables of interest. The current version always set it as \code{FALSE} as using more samples can increase power. See 'Details' below.
}
  \item{seed}{
an integer used as the random seed. 
}
  \item{nthread}{
an interger specifying the number of threads used in parallelizing the Monte Carlo test. Use all available threads by default. 
}
  \item{plot.pval}{
a logical indicating whether to draw a plot of the unadjusted p-values for each pair of \code{rho} and \code{kappa}. See 'Details' below.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The \code{formula} is parsed by the package \code{\link[Formula]{Formula}}. On the left-hand side, the binary outcome must be specified first, then the continuous outcome is specified, separated by |. One the right-hand side, the covariates must be specified first, then the variables of interest are specified, separated by |. A valid formula, e.g.,

\code{SMOKE | CIG_PER_DAY ~ AGE + SEX | SNP1 + SNP2 + SNP3}

means that the variable \code{SMOKE} is the binary outcome, \code{CIG_PER_DAY} is the continuous outcome. Both outcomes should be adjusted by \code{AGE} and \code{SEX}. The function will test the joint association effect of three SNPs, i.e., \code{SNP1}, \code{SNP2}, \code{SNP3}.

Ge's minp algorithm is used in evaluating the final p-value accounting for multiple-comparison in tuning parameters of \code{rho} and \code{kappa}. It will produce unadjusted statistics which are returned as \code{$obs.rank}. See 'Value' below. A generic function \code{plot} can be used to visualize these statistics (standardized as p-values), which gives intuition of the optimal chosen \code{$rho.opt} and \code{$kappa.opt}.

One of the major problems for gene-level test is that it doesn't allow missing genotypes. In practice, the user has to exclude individuals even with one missing entry. Although generally the SNPs included in a gene pass the quality control, e.g., missing rate < 2\%, however, a substantial proportion of individuals can be excluded in testing the association, especially for large gene. This can reduce the statistical power or more seriously, bias the inference. We propose to use the modified scores defined on all observed genotypes to generalize the score tests, which provides more flexibility to in real application. Please refer to our paper for more details.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{svcopt} returns an object of class "\code{svcopt}" containing p-value and other useful information.

An object of class "\code{svcopt}" is a list containing the following components:

\item{pval}{the final p-value for the test after accounting for multiple-comparison.}
\item{rho.opt}{the optimal chosen \code{rho}.}
\item{kappa.opt}{the optimal chosen \code{kappa}.}
\item{nperm}{the number of replicates used in Monte Carlo test.}
\item{rho}{the vector of effect correlation parameter used in tuning process.}
\item{kappa}{the vector of variance component proportion parameter used in tuning process.}
\item{refine}{if \code{TRUE}, larger \code{nperm} should be try to obtain a more stable estimate of p-value.}
\item{obs.rank}{a vector containing unadjusted statistics produced by Ge's minp algorithm.}
\item{stat}{a vector of statistics used to estimate the final p-value \code{$pval}.}
}
\references{
Zhang, H., Wu, C.O., Yang, Y., and Yu K. (2014) Optimal variance component test for stratified sampling design with its application to the genome-wide association study of aggressive prostate cancer. Submitted.

Wu, C.O., Zheng, G., and Kwak, M. (2013) A joint regression analysis for genetic association studies with outcome stratified samples. \emph{Biometrics} \bold{69}, 417--426.

Ge, Y., Dudoit, S., and Speed, T. (2003) Resampling-based multiple testing for microarray data analysis, \emph{Test} \bold{12}, 1--77. 
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
data(test)
obj <- svcopt(formula, data, nperm = 1e5, seed = 1)
plot(obj)
}
