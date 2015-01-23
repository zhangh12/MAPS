\name{MAPS-package}
\alias{MAPS-package}
\alias{MAPS}
\docType{package}
\title{
A multi-locus genetic association test for a dichotomous trait and its secondary phenotype
}
\description{
The dichotomous trait could be the disease status from a prospective cohort study or from a retrospective case-control study. The quantitative secondary phenotype is assumed to be only collected on the subjects with disease condition. This test is only allow one secondary phenotype. 
}
\details{
\tabular{ll}{
Package: \tab MAPS\cr
Type: \tab Package\cr
Version: \tab 0.1.10\cr
Date: \tab 2015-01-22\cr
License: \tab GPL-3\cr
}


The proposed test can be performed through the function \code{maps}. A multi-locus score test generalized from Wu's single-marker test is implemented in function \code{wu}. Notice that when there is only one SNP, \code{wu} is actually the test in Wu et al. (2013). The function \code{gforest} can generate the forest plot for each SNP specified in the \code{formula}
}
\author{
Han Zhang and Yifan Yang

Maintainer: Han Zhang <han.zhang2@nih.gov>
}
\references{
Zhang, H., Wu, C.O., Yang, Y., Berndt, S., and Yu K. (2015) A multi-locus genetic association test for a dichotomous trait and its secondary phenotype. Submitted.

Wu, C.O., Zheng, G., and Kwak, M. (2013) A joint regression analysis for genetic association studies with outcome stratified samples. \emph{Biometrics} \bold{69}, 417--426.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{

}
\examples{
%~~ simple examples of the most important functions ~~
}
