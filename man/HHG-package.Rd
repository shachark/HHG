\name{HHG-package}
\alias{HHG-package}
\alias{HHG}
\docType{package}
\title{Heller-Heller-Gorfine Tests of Independence, Equality of Distributions, and Goodness-of-fit}
\description{Heller-Heller-Gorfine (HHG) tests are a set of powerful statistical tests of independnece between two random vectors of arbitrary dimensions, and equality of two or more multivariate distributions given independent samples. For testing independence between two scalar random variables, k-sample testing, and goodness-of-fit testing, the package also contains implementations of the data-derived partitions (DDP) and all data partitions (ADP) tests, which are distribution-free and thus are much faster to apply.}
\details{
\tabular{ll}{
Package: \tab HHG\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2014-03-25\cr
License: \tab GPL-2\cr
}
The package contains four major functions: hhg.test, hhg.test.k.sample, hhg.test.2.sample, xdp.test, xdp.test.k.sample, xdp.test.gof.
}
\author{
Shachar Kaufman, based in part on an earlier implementation of the original HHG test by Ruth Heller <ruheller@post.tau.ac.il> and Yair Heller <heller.yair@gmail.com>.
Maintainer: Shachar Kaufman <shachark@post.tau.ac.il>
}
\references{
Heller, R., Heller, Y., & Gorfine, M. (2013). A consistent multivariate test of association based on ranks of distances. \emph{Biometrika}, 100(2), 503-510.

Heller, R., Heller, Y., Kaufman S., & Gorfine, M. (2014). Consistent distribution-free tests of association between univariate random variables. \emph{arXiv:1308.1559}.
}
%\seealso{
%}
\examples{
# See examples in the documentation for hhg.test and xdp.test.
}
