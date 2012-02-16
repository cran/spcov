\name{ProxADMM}
\alias{ProxADMM}
\title{Solving penalized Frobenius problem.}
\description{
  This function solves the optimization problem

  Minimize_X { (1/2)||X - A||_F^2 + lam||P*X||_1} s.t. X >= del * I.

  This is the prox function for the generalized gradient descent of Bien
  & Tibshirani 2011 (see full reference below).}
\usage{
ProxADMM(A, del, lam, P, rho=.1, tol=1e-6, maxiters=100, verb=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{A symmetric matrix.}
  \item{del}{A non-negative scalar. Lower bound on eigenvalues.}
  \item{lam}{A non-negative scalar. L1 penalty parameter.}
  \item{P}{Matrix with non-negative elements and dimension of A. Allows
    for differing L1 penalty parameters.}
  \item{rho}{ADMM parameter.  Can affect rate of convergence a lot.}
  \item{tol}{Convergence threshold.}
  \item{maxiters}{Maximum number of iterations.}
  \item{verb}{Controls whether to be verbose.}
}
\details{
  This is the R implementation of the algorithm in Appendix 3 of
  Bien, J., and Tibshirani, R. (2011), "Sparse Estimation of a Covariance Matrix," Biometrika. 98(4). 807--820.  It uses an ADMM approach to solve the problem

  Minimize_X { (1/2)||X - A||_F^2 + lam||P*X||_1} s.t. X >= del * I.

  Here, the multiplication between P and X is elementwise.  The inequality in the constraint is a lower bound on the minimum eigenvalue of the matrix X.

  Note that there are two variables X and Z that are outputted.  Both are estimates of the optimal X.  However, Z has exact zeros whereas X has eigenvalues at least del.  Running the ADMM algorithm long enough, these two are guaranteed to converge.
}
\value{
  \item{X}{Estimate of optimal X.}
  \item{Z}{Estimate of optimal X.}
  \item{obj}{Objective values.}
}
\references{
  Bien, J., and Tibshirani, R. (2011), "Sparse Estimation of a Covariance Matrix," Biometrika. 98(4). 807--820.
  }
\author{
Jacob Bien and Rob Tibshirani
}

\seealso{spcov}

\examples{
set.seed(1)
n <- 100
p <- 200
# generate a covariance matrix:
model <- GenerateCliquesCovariance(ncliques=4, cliquesize=p / 4, 1)

# generate data matrix with x[i, ] ~ N(0, model$Sigma):
x <- matrix(rnorm(n * p), ncol=p) \%*\% model$A
S <- var(x)

# compute sparse, positive covariance estimator:
P <- matrix(1, p, p)
diag(P) <- 0
lam <- 0.1
aa <- ProxADMM(S, 0.01, lam, P)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

