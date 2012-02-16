\name{spcov}
\alias{spcov}
\title{Sparse Covariance Estimation}
\description{
Provides a sparse and positive definite estimate of a covariance
matrix.  This function performs the majorize-minimize algorithm described in Bien & Tibshirani
2011 (see full reference below).}
\usage{
spcov(Sigma, S, lambda, step.size, nesterov = TRUE, n.outer.steps =
10000, n.inner.steps = 10000, tol.outer = 1e-04, thr.inner = 0.01,
backtracking = 0.2, trace = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{an initial guess for Sigma (suggestions: \code{S} or \code{diag(diag(S)))}.}
  \item{S}{the empirical covariance matrix of the data.  Must be
    positive definite (if it is not, add a small constant to the diagonal).}
  \item{lambda}{penalty parameter.  Either a scalar or a matrix of the
    same dimension as \code{Sigma}.  This latter choice should be used
    to penalize only off-diagonal elements.  All elements of \code{lambda} must be non-negative.}
  \item{step.size}{the step size to use in generalized gradient
    descent.  Affects speed of algorithm.}
  \item{nesterov}{indicates whether to use Nesterov's modification of
    generalized gradient descent. Default: \code{TRUE}.}
  \item{n.outer.steps}{maximum number of majorize-minimize steps to take
  (recall that MM is the outer loop).}
  \item{n.inner.steps}{maximum number of generalized gradient steps to
    take (recall that generalized gradient descent is the inner loop).}
  \item{tol.outer}{convergence threshold for outer (MM) loop.  Stops
    when drop in objective between steps is less than \code{tol.outer}.}
  \item{thr.inner}{convergence threshold for inner (i.e. generalized
    gradient) loop.  Stops when mean absolute change in \code{Sigma} is less than \code{thr.inner * mean(abs(S))}.}
  \item{backtracking}{if \code{FALSE}, then fixed step size used.  If
    numeric and in (0,1), this is the parameter of backtracking that
    multiplies \code{step.size} on each step.  Usually, in range of
    (0.1, 0.8). Default: \code{0.2}.}
  \item{trace}{controls how verbose output should be.}
}
\details{
  This is the R implementation of Algorithm 1 in
  Bien, J., and Tibshirani, R. (2011), "Sparse Estimation of a
  Covariance Matrix," Biometrika. 98(4). 807--820.  The goal is to approximately
  minimize (over Sigma) the following non-convex optimization problem:
  
  minimize logdet(Sigma) + trace(S Sigma^{-1}) + || lambda*Sigma ||_1
  subject to Sigma positive definite.

  Here, the L1 norm and matrix multiplication between lambda and Sigma
  are elementwise.  The empirical covariance matrix must be positive
  definite for the optimization problem to have bounded objective (see
  Section 3.3 of paper). We suggest adding a small constant to the diagonal of S
  if it is not.  Since the above problem is not convex, the returned
  matrix is not guaranteed to be a global minimum of the problem.
  
  In Section 3.2 of the paper, we mention a
  simple modification of gradient descent due to Nesterov.  The argument
  \code{nesterov} controls whether to use this modification (we suggest
  that it be used).  We also strongly recommend using backtracking.
  This allows the algorithm to begin by taking large steps (the initial
  size is determined by the argument \code{step.size}) and then to
  gradually reduce the size of steps.

  At the start of the algorithm, a lower bound (\code{delta} in the paper) on the
  eigenvalues of the solution is calculated.  As shown in Equation (3) of the paper, the prox function for our generalized gradient
  descent amounts to minimizing (over a matrix X) a problem of the form

  minimize (1/2)|| X-A ||_F^2 + || lambda*X ||_1 subject to X >= delta I

  This is implemented using an alternating direction method of
  multipliers approach given in Appendix 3.
}
\value{
  \item{Sigma}{the sparse covariance estimate}
  \item{n.iter}{a vector giving the number of generalized gradient steps
  taken on each step of the MM algorithm}
  \item{obj}{a vector giving the objective values after each step of the
    MM algorithm}
}
\references{
  Bien, J., and Tibshirani, R. (2011), "Sparse Estimation of a Covariance Matrix," Biometrika. 98(4). 807--820.
  }
\author{
Jacob Bien and Rob Tibshirani
}

\seealso{ProxADMM}

\examples{
set.seed(1)
n <- 100
p <- 20
# generate a covariance matrix:
model <- GenerateCliquesCovariance(ncliques=4, cliquesize=p / 4, 1)

# generate data matrix with x[i, ] ~ N(0, model$Sigma):
x <- matrix(rnorm(n * p), ncol=p) \%*\% model$A
S <- var(x)

# compute sparse, positive covariance estimator:
step.size <- 100
tol <- 1e-3
P <- matrix(1, p, p)
diag(P) <- 0
lam <- 0.06
mm <- spcov(Sigma=S, S=S, lambda=lam * P,
            step.size=step.size, n.inner.steps=200,
            thr.inner=0, tol.outer=tol, trace=1)
sqrt(mean((mm$Sigma - model$Sigma)^2))
sqrt(mean((S - model$Sigma)^2))
\dontrun{image(mm$Sigma!=0)}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

