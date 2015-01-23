\name{wu}
\alias{wu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A multi-locus score test directly generalized from Wu's single-marker test.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This test is generalized from Wu's test. It is modified so that it allows missing genotypes (MCAR) as well. 
}
\usage{
wu(data, formula = NULL, subset = NULL, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
a data frame containing all variables specified in \code{formula}.
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
an object of clas "\code{Formula}". The user can specify two outcomes on the left-hand side. The covariates and the variables of interest can be specified on the right-hand side. See 'Details' in page of function \code{maps}. 
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
an optional vector specifying a subset of observations to be used in the testing process.
}

  \item{na.rm}{
a logical indicating whether to allow individuals with missing values on variables of interest. The current version always set it as \code{FALSE} as using more samples can increase power. See 'Details' in page of \code{maps}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
If only one SNP is included in \code{data}, \code{wu} is actually the test proposed in Wu et al. (2013).
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
Zhang, H., Wu, C.O., Yang, Y., Berndt, S., and Yu K. (2015) A multi-locus genetic association test for a dichotomous trait and its secondary phenotype. Submitted.

Wu, C.O., Zheng, G., and Kwak, M. (2013) A joint regression analysis for genetic association studies with outcome stratified samples. \emph{Biometrics} \bold{69}, 417--426.
}
\author{
%%  ~~who you are~~
Han Zhang <han.zhang2@nih.gov>
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
obj <- wu(data, formula)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
