# Implements MM algorithm with generalized gradient descent (and Nesterov's
# method) using ADMM to solve prox function.



#' Sparse Covariance Estimation
#' 
#' Provides a sparse and positive definite estimate of a covariance matrix.
#' This function performs the majorize-minimize algorithm described in Bien &
#' Tibshirani 2011 (see full reference below).
#' 
#' This is the R implementation of Algorithm 1 in Bien, J., and Tibshirani, R.
#' (2011), "Sparse Estimation of a Covariance Matrix," Biometrika. 98(4).
#' 807--820.  The goal is to approximately minimize (over Sigma) the following
#' non-convex optimization problem:
#' 
#' minimize logdet(Sigma) + trace(S Sigma^-1) + || lambda*Sigma ||_1 subject to
#' Sigma positive definite.
#' 
#' Here, the L1 norm and matrix multiplication between lambda and Sigma are
#' elementwise.  The empirical covariance matrix must be positive definite for
#' the optimization problem to have bounded objective (see Section 3.3 of
#' paper). We suggest adding a small constant to the diagonal of S if it is
#' not.  Since the above problem is not convex, the returned matrix is not
#' guaranteed to be a global minimum of the problem.
#' 
#' In Section 3.2 of the paper, we mention a simple modification of gradient
#' descent due to Nesterov.  The argument \code{nesterov} controls whether to
#' use this modification (we suggest that it be used).  We also strongly
#' recommend using backtracking.  This allows the algorithm to begin by taking
#' large steps (the initial size is determined by the argument
#' \code{step.size}) and then to gradually reduce the size of steps.
#' 
#' At the start of the algorithm, a lower bound (\code{delta} in the paper) on
#' the eigenvalues of the solution is calculated.  As shown in Equation (3) of
#' the paper, the prox function for our generalized gradient descent amounts to
#' minimizing (over a matrix X) a problem of the form
#' 
#' minimize (1/2)|| X-A ||_F^2 + || lambda*X ||_1 subject to X >= delta I
#' 
#' This is implemented using an alternating direction method of multipliers
#' approach given in Appendix 3.
#' 
#' @param Sigma an initial guess for Sigma (suggestions: \code{S} or
#' \code{diag(diag(S)))}.
#' @param S the empirical covariance matrix of the data.  Must be positive
#' definite (if it is not, add a small constant to the diagonal).
#' @param lambda penalty parameter.  Either a scalar or a matrix of the same
#' dimension as \code{Sigma}.  This latter choice should be used to penalize
#' only off-diagonal elements.  All elements of \code{lambda} must be
#' non-negative.
#' @param step.size the step size to use in generalized gradient descent.
#' Affects speed of algorithm.
#' @param nesterov indicates whether to use Nesterov's modification of
#' generalized gradient descent. Default: \code{TRUE}.
#' @param n.outer.steps maximum number of majorize-minimize steps to take
#' (recall that MM is the outer loop).
#' @param n.inner.steps maximum number of generalized gradient steps to take
#' (recall that generalized gradient descent is the inner loop).
#' @param tol.outer convergence threshold for outer (MM) loop.  Stops when drop
#' in objective between steps is less than \code{tol.outer}.
#' @param thr.inner convergence threshold for inner (i.e. generalized gradient)
#' loop.  Stops when mean absolute change in \code{Sigma} is less than
#' \code{thr.inner * mean(abs(S))}.
#' @param backtracking if \code{FALSE}, then fixed step size used.  If numeric
#' and in (0,1), this is the parameter of backtracking that multiplies
#' \code{step.size} on each step.  Usually, in range of (0.1, 0.8). Default:
#' \code{0.2}.
#' @param trace controls how verbose output should be.
#' @return \item{Sigma}{the sparse covariance estimate} \item{n.iter}{a vector
#' giving the number of generalized gradient steps taken on each step of the MM
#' algorithm} \item{obj}{a vector giving the objective values after each step
#' of the MM algorithm}
#' @author Jacob Bien and Rob Tibshirani
#' @seealso ProxADMM
#' @references Bien, J., and Tibshirani, R. (2011), "Sparse Estimation of a
#' Covariance Matrix," Biometrika. 98(4). 807--820.
#' @keywords multivariate
#' @export
#' @examples
#' 
#' set.seed(1)
#' n <- 100
#' p <- 20
#' # generate a covariance matrix:
#' model <- GenerateCliquesCovariance(ncliques=4, cliquesize=p / 4, 1)
#' 
#' # generate data matrix with x[i, ] ~ N(0, model$Sigma):
#' x <- matrix(rnorm(n * p), ncol=p) %*% model$A
#' S <- var(x)
#' 
#' # compute sparse, positive covariance estimator:
#' step.size <- 100
#' tol <- 1e-3
#' P <- matrix(1, p, p)
#' diag(P) <- 0
#' lam <- 0.06
#' mm <- spcov(Sigma=S, S=S, lambda=lam * P,
#'             step.size=step.size, n.inner.steps=200,
#'             thr.inner=0, tol.outer=tol, trace=1)
#' sqrt(mean((mm$Sigma - model$Sigma)^2))
#' sqrt(mean((S - model$Sigma)^2))
#' \dontrun{image(mm$Sigma!=0)}
#' 
spcov <- function(Sigma, S, lambda, step.size, nesterov=TRUE,
               n.outer.steps=1e4, n.inner.steps=1e4,
               tol.outer=1e-4, thr.inner=1e-2,
               backtracking=0.2, trace=0) {
  # Performs MM to optimize the non-convex objective function
  # Args:
  #  Sigma: an initial guess for Sigma
  #  lambda: either a scalar or a matrix of the same dimensions as Sigma
  #  nesterov: indicates whether to use Nesterov or standard generalized gradient
  #                    descent to perform inner loop optimization
  #  tol.outer: convergence threshold for outer (MM) loop.
  #  thr.inner: convergence threshold for inner (GG/Nesterov) loop.
  #             stops when mean absolute change in Sigma is
  #             less than thr.inner * mean(abs(S))
  #  backtracking: see "GGDescent"
  if (all(lambda == 0)) {
    cat("Skipping MM.  Solution is S!", fill=T)
    return(list(n.iter=0, Sigma=S, obj=ComputeObjective(S, S, lambda)))
  }
  stopifnot(lambda >= 0)
  if (trace > 0) {
    cat("---", fill=T)
    cat(ifelse(nesterov, "using Nesterov, ", ""))
    cat(ifelse(backtracking, "backtracking line search", ""), fill=T)
    cat("---", fill=T)
  }
  mean.abs.S <- mean(abs(S))
  if (min(eigen(Sigma, symmetric=T, only.values=T)$val) < 1e-5)
    warning("Starting value is nearly singular.")

  del <- ComputeDelta(S, lambda, trace=trace-1) # get a lower bound on minimum eval

  objective <- ComputeObjective(Sigma, S, lambda)
  if (trace > 0)
    cat("objective: ", objective, fill=T)
  n.iter <- NULL # number of inner iterations on each step
  for (i in seq(n.outer.steps)) {
    Sigma0 <- Sigma
    if (trace > 0)
      cat("step size given to GGDescent/Nesterov:", step.size, fill=T)
    gg <- GGDescent(Sigma=Sigma, Sigma0=Sigma0, S=S, lambda=lambda,
                    del=del, nsteps=n.inner.steps,
                    step.size=step.size,
                    nesterov=nesterov,
                    tol=thr.inner * mean.abs.S,
                    trace=trace - 1,
                    backtracking=backtracking)
    Sigma <- gg$Sigma
    objective <- c(objective, ComputeObjective(Sigma, S, lambda))
    if (trace > 0) {
      cat("objective: ", objective[length(objective)],
          " (", gg$niter, "iterations, max step size:",
          max(gg$step.sizes), ")",
          fill=T)
    }
    if (backtracking) {
      if (max(gg$step.sizes) < step.size * backtracking ^ 2) {
        step.size <- step.size * backtracking
        if (trace > 0)
          cat("Reducing step size to", step.size, fill=T)
      }
    }
    n.iter <- c(n.iter, gg$niter)
    if(objective[i + 1] > objective[i] - tol.outer) {
       cat("MM converged in", i, "steps!", fill=T)
       break
    }
  }
  list(n.iter=n.iter, Sigma=gg$Sigma, obj=objective)
}

GGDescent <- function(Sigma, Sigma0, S, lambda, del, nsteps,
                                     step.size, nesterov=FALSE, backtracking=FALSE,
                                     tol = 1e-3, trace=0) {
  # solves the problem
  # Min_{Sigma} tr(solve(Sigma0,Sigma)) + tr(solve(Sigma,S))
  #                                     + ||lambda * Sigma||_1
  # using Nesterov's method with backtracking.
  # Note: We wish to solve this with the constraint that Sigma pd.
  # However, this algorithm does not impose this constraint.
  # Args:
  #  Sigma: an initial guess for Sigma
  #  Sigma0, S, lambda, del: parameters of the optimization problem
  #  nsteps: number of generalized gradient steps to take
  #  nesterov: TRUE/FALSE, indicates whether to take Nesterov vs. standard gen grad steps.
  #  backtracking: if FALSE, then fixed step size used.  If numeric and in
  #                (0,1), this is the beta parameter of backtracking.
  #                Usually, beta is in (0.1, 0.8).
  #  tol:  convergence threshold.  Stops when mean(abs(Sigma-Sigma.last)) < tol
  if (backtracking) {
    beta <- backtracking
    if (beta <= 0 | beta >= 1)
      stop("Backtracking parameter beta must be in (0,1).")
  }
  tt <- step.size
  converged <- FALSE
  exit <- FALSE
  obj.starting <- ComputeObjective(Sigma, S, lambda)
  Sigma.starting <- Sigma
  Omega <- Sigma
  Sigma.last <- Sigma
  ttts <- NULL
  ttt <- tt # note: as in Beck & Teboulle, step size only gets smaller.
  for (i in seq(nsteps)) {
    inv.Sigma0 <- solve(Sigma0)
    log.det.Sigma0 <- LogDet(Sigma0)
    grad.g <- ComputeGradientOfg(Omega, S, Sigma0, inv.Sigma0=inv.Sigma0)
    grad.g <- (grad.g + t(grad.g)) / 2 # make sure this stays symmetric
    g.omega <- g(Omega, S, Sigma0,
                 inv.Sigma0=inv.Sigma0, log.det.Sigma0=log.det.Sigma0)
    # backtracking line search:
    while (backtracking) {
      #soft.thresh <- SoftThreshold(Omega - ttt * grad.g, lambda * ttt)
      soft.thresh <- ProxADMM(Omega - ttt * grad.g, del, 1, P=lambda*ttt, rho=.1)$X
      gen.grad.g <- (Omega - soft.thresh) / ttt
      left <- g(soft.thresh, S, Sigma0,
                inv.Sigma0=inv.Sigma0, log.det.Sigma0=log.det.Sigma0)
      right <- g.omega - ttt * sum(grad.g * gen.grad.g) + ttt * sum(gen.grad.g ^ 2) / 2
      if (is.na(left) || is.na(right)) {
        print("left or right is NA.")
        browser()
      }
      if (left <= right) {
        # accept this step size
        Sigma <- soft.thresh
        ttts <- c(ttts, ttt)
        # check for convergence
        if (mean(abs(Sigma - Sigma.last)) < tol) {
          converged <- TRUE
          break # note: this break only ends the backtracking loop
        }
        if (nesterov)
          Omega <- Sigma + (i - 1) / (i + 2) * (Sigma - Sigma.last)
        else
          Omega <- Sigma
        Sigma.last <- Sigma
        if (trace > 0)
          cat("--true objective:", ComputeObjective(Sigma, S, lambda), fill=T)
        if (trace > 0)
          cat(i, ttt, " ")
        break
      }
      ttt <- beta * ttt
      if (ttt < 1e-15) {
        cat("Step size too small: no step taken", fill=T)
        exit <- TRUE
        break
      }
    }
    if (!backtracking) {
      #Sigma <- SoftThreshold(Sigma - ttt * grad.g, lambda * ttt)
      Sigma <- ProxADMM(Sigma - ttt * grad.g, del, 1, P=lambda*ttt, rho=.1)$X
      # check for convergence:
      if (mean(abs(Sigma - Sigma.last)) < tol)
        converged <- TRUE
      if (nesterov)
        Omega <- Sigma + (i - 1)/(i + 2) * (Sigma - Sigma.last)
      else
        Omega <- Sigma
      Sigma.last <- Sigma
    }
    if (converged) {
      if(trace > 0) {
        cat("--GG converged in", i, "steps!")
        if (backtracking)
          cat(" (last step size:", ttt, ")", fill=T)
        else
          cat(fill=T)
      }
      break
    }
    if (exit) {
      break
    }
  }
  obj.end <- ComputeObjective(Sigma, S, lambda)
  if (obj.starting < obj.end) {
    if (nesterov) {
      cat("Objective rose with Nesterov.  Using generalized gradient instead.", fill=T)
      return(GGDescent(Sigma=Sigma.starting, Sigma0=Sigma0, S=S, lambda=lambda,
                       del=del, nsteps=nsteps, step.size=step.size,
                       nesterov=FALSE,
                       backtracking=backtracking,
                       tol = tol, trace=trace))
    }
    
    browser()
    cat("--Returning initial Sigma since GGDescent/Nesterov did not decrease objective", fill=T)
    Sigma <- Sigma.starting
  }

  list(Sigma=Sigma, niter=i, step.sizes=ttts)
}


ComputeObjective <- function(Sigma, S, lambda) {
  # the original non-convex problem's objective function
  -2 * ComputeLikelihood(Sigma, S) + ComputePenalty(Sigma, lambda)
}
