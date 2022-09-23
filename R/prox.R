# Minimize_X { (1/2)||X - A||_F^2 + lam||P*X||_1} s.t. X >= del * I
# ...using ADMM



#' Solving penalized Frobenius problem.
#' 
#' This function solves the optimization problem
#' 
#' Minimize_X (1/2)||X - A||_F^2 + lam||P*X||_1 s.t. X >= del * I.
#' 
#' This is the prox function for the generalized gradient descent of Bien &
#' Tibshirani 2011 (see full reference below).
#' 
#' This is the R implementation of the algorithm in Appendix 3 of Bien, J., and
#' Tibshirani, R. (2011), "Sparse Estimation of a Covariance Matrix,"
#' Biometrika. 98(4). 807--820.  It uses an ADMM approach to solve the problem
#' 
#' Minimize_X (1/2)||X - A||_F^2 + lam||P*X||_1 s.t. X >= del * I.
#' 
#' Here, the multiplication between P and X is elementwise.  The inequality in
#' the constraint is a lower bound on the minimum eigenvalue of the matrix X.
#' 
#' Note that there are two variables X and Z that are outputted.  Both are
#' estimates of the optimal X.  However, Z has exact zeros whereas X has
#' eigenvalues at least del.  Running the ADMM algorithm long enough, these two
#' are guaranteed to converge.
#' 
#' @param A A symmetric matrix.
#' @param del A non-negative scalar. Lower bound on eigenvalues.
#' @param lam A non-negative scalar. L1 penalty parameter.
#' @param P Matrix with non-negative elements and dimension of A. Allows for
#' differing L1 penalty parameters.
#' @param rho ADMM parameter.  Can affect rate of convergence a lot.
#' @param tol Convergence threshold.
#' @param maxiters Maximum number of iterations.
#' @param verb Controls whether to be verbose.
#' @return \item{X}{Estimate of optimal X.} \item{Z}{Estimate of optimal X.}
#' \item{obj}{Objective values.}
#' @author Jacob Bien and Rob Tibshirani
#' @seealso spcov
#' @references Bien, J., and Tibshirani, R. (2011), "Sparse Estimation of a
#' Covariance Matrix," Biometrika. 98(4). 807--820.
#' @keywords multivariate
#' @export
#' @examples
#' 
#' set.seed(1)
#' n <- 100
#' p <- 200
#' # generate a covariance matrix:
#' model <- GenerateCliquesCovariance(ncliques=4, cliquesize=p / 4, 1)
#' 
#' # generate data matrix with x[i, ] ~ N(0, model$Sigma):
#' x <- matrix(rnorm(n * p), ncol=p) %*% model$A
#' S <- var(x)
#' 
#' # compute sparse, positive covariance estimator:
#' P <- matrix(1, p, p)
#' diag(P) <- 0
#' lam <- 0.1
#' aa <- ProxADMM(S, 0.01, lam, P)
#' 
ProxADMM <- function(A, del, lam, P, rho=.1, tol=1e-6, maxiters=100, verb=FALSE) {
  # Minimize_X { (1/2)||X - A||_F^2 + lam||P*X||_1} s.t. X >= del * I
  #
  # ADMM approach

  # first, check if simple soft-thesholding works... if so, skip the ADMM!
  soft <- SoftThreshold(A, lam * P)
  minev <- min(eigen(soft, symmetric=T, only.values=T)$val)
  if (minev >= del) {
    return(list(X=soft, Z=soft, obj=ComputeProxObjective(soft, A, lam, P)))
  }
  
  p <- nrow(A)
  obj <- NULL
  
  # initialize Z, Y
  Z <- soft
  Y <- matrix(0, p, p)

  # main loop
  for (i in seq(maxiters)) {    
    # update X:
    B <- (A + rho * Z - Y) / (1 + rho)
    if (min(eigen(B, symmetric=T, only.values=T)$val) < del) {
      # note: even though eigen is called twice, only.values=T is
      #       much faster, making this worthwhile.
      eig <- eigen(B, symmetric=T)
      X <- eig$vec %*% diag(pmax(eig$val, del)) %*% t(eig$vec)
    }
    else {
      X <- B
    }
    # check for convergence:
    obj <- c(obj, ComputeProxObjective(X, A, lam, P))
    if (verb)
      cat(" ", obj[i], fill=T)
    if (i > 1)
      if (obj[i] > obj[i - 1] - tol) {
        if (verb)
          cat(" ADMM converged after ", i, " steps.", fill=T)
        break
      }
      
    # update Z:
    Z <- SoftThreshold(X + Y / rho, lam * P / rho)

    # update Y:
    Y <- Y + rho * (X - Z)
  }

  list(X=X, Z=Z, obj=obj)
}

SoftThreshold <- function(x, lam) {
  # note: this works also if lam is a matrix of the same size as x.
  sign(x) * (abs(x) - lam) * (abs(x) > lam)
}

ComputeProxObjective <- function(X, A, lam, P) {
  sum((X-A)^2) / 2 + lam * sum(abs(P*X))
}
