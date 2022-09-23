#' Generate a block diagonal covariance matrix
#' 
#' This function is included in the package so that it can be used in the
#' example code provided in \code{\link{spcov}}.
#' 
#' This function generates a block diagonal positive definite matrix with
#' randomly-signed, non-zero elements.  A shift is added to the diagonal of the
#' matrix so that its condition number equals \code{p}, the number of
#' variables.
#' 
#' @param ncliques number of blocks
#' @param cliquesize size of each block
#' @param theta magnitude of non-zeros
#' @return \item{Sigma}{the covariance matrix} \item{A}{symmetric square root
#' of \code{Sigma}} \item{shift}{how much the eigenvalues were shifted.  See
#' details.}
#' @author Jacob Bien and Rob Tibshirani
#' @references Bien, J., and Tibshirani, R. (2011), "Sparse Estimation of a
#' Covariance Matrix," accepted for publication in Biometrika.
#' @keywords multivariate
#' @export
GenerateCliquesCovariance <- function(ncliques, cliquesize, theta) {
  # Generates a block diagonal positive definite matrix with randomly signed
  # non-zero elements and condition number equal to p=ncliques*cliquesize
  #
  # Args:
  #   ncliques: number of cliques
  #   cliquesize: size of clique
  #   theta: magnitude of non-zeros
  p <- ncliques * cliquesize
  sizes <- rep(cliquesize, ncliques)
  Sigma <- matrix(0, p, p)
  lims <- c(0, cumsum(sizes))
  for (i in seq(ncliques)) {
    ii <- (lims[i] + 1):lims[i+1]
    signs <- 2 * (stats::rnorm(sizes[i] * (sizes[i] - 1) / 2) < 0) - 1
    Sigma[ii, ii][upper.tri(Sigma[ii, ii])] <- signs * theta
  }
  Sigma <- Sigma + t(Sigma)
  eig <- eigen(Sigma, symmetric=T)
  shift <- (max(eig$val) - p * min(eig$val)) / (p - 1)
  cat("Shifting eigenvalues by ", shift, fill=T)
  diag(Sigma) <- diag(Sigma) + shift

  # compute symmetric square root for generating data
  A <- eig$vect %*% diag(sqrt(eig$val + shift)) %*% t(eig$vect)
  list(Sigma=Sigma, A=A, shift=shift)
}
