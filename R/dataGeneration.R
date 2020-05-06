#' Generate a high-dimensional mediation dataset
#'
#' This function generates a dataset from an x -> M -> y model, where M may
#' be of any size with any correlation matrix.
#'
#' @param n Sample size
#' @param a Vector of a path coefficients within 0 and 1
#' @param b Vector of b path coefficients within 0 and 1
#' @param r2y Proportion of explained variance in y. Set to
#' \code{b \%*\% Sigma \%*\% b} for \code{var(y) == 1}.
#' @param dir Direct path from x to y
#' @param Sigma Desired true covariance matrix between the mediators M
#' @param residual Whether Sigma indicates residual or marginal covariance
#' @param empirical Ensure observed data matrix has exactly the requested covmat
#' (only if Sigma is specified)
#' @param scaley Whether to standardise y (changes b path coefficients)
#' @param forma Functional form of the a paths. Function that accepts a matrix
#' as input and transforms each column to the desired form.
#' @param formb Functional form of the b paths. Function that accepts a vector.
#'
#' @return A data frame with columns x, M.1 - M.p, y
#'
#' @examples
#' # Generate a suppression dataset where M.2 is suppressed
#' sup <- generateMed(n = 100,
#'                    a = c(-0.4, 0.4),
#'                    b = c(0.8, 0.48),
#'                    Sigma = matrix(c(1, -0.6, -0.6, 1), 2))
#'
#' @importFrom stats rnorm cov
#' @import Matrix
#'
#' @export
generateMed <- function(n = 1e2L,
                        a = 0.3, b = 0.3,
                        r2y = 0.5, dir = 0,
                        Sigma, residual = FALSE, empirical = FALSE,
                        scaley = FALSE,
                        forma = identity, formb = identity) {

  # Input checking
  if (!inherits(a, "dsparseVector")) {
    # user knows what they're doing if sparse vectors, skip input check
    if (r2y > 1 || r2y <= 0) {
      stop("Arg r2y should be in range (0,1]")
    }
    if (length(a) != length(b)) {
      stop("There should be as many a paths as b paths!")
    }
    if (!missing(Sigma)) {
      if (!is.matrix(Sigma) ||
          ncol(Sigma) != nrow(Sigma) ||
          !isSymmetric(Sigma)) {
        stop("Sigma needs to be a square symmetric covariance matrix.")
      }
      if (ncol(Sigma) != length(a)) {
        stop("Sigma should have as many cols & rows as a & b paths.")
      }
    }
  }


  if (missing(Sigma)) {
    if (any(c(a > 1, a < -1))) {
      stop("Elements in 'a' out of range, should all be between 0 and 1",
           "if Sigma is not specified. (Assume variance of 1 for M)")
    }

    # No residual covariance, fast to generate
    # Calculate residual variance of y
    vary <- b %*% b + dir^2 # variance of y = sum(b^2) because M std & no cov
    resy <- vary/r2y - vary # residual variance of y

    # Simulate the residuals of M
    resM <- sapply(1 - a^2, function(x) rnorm(n, sd = x))

  } else {

    # There is residual covariance

    # Switch between residual and marginal covariance of M.
    # This uses the Schur complement which simplifies to Sigma +/- tcrossprod(a)
    # because C = B' and A = 1. See Abadir & Magnus p. 102.
    if (residual) {
      psi <- Sigma
      Sigma <- psi + tcrossprod(a)
    } else {
      psi <- Sigma - tcrossprod(a)
    }

    # Marginal covariance matrix of x & M
    if (inherits(Sigma, "dgCMatrix")) {
      # use sparse matrices. empirical unavailable
      if (empirical) stop("Empirical not available with sparse matrices")


      # generate residuals of M
      prec <- Matrix::solve(psi, sparse = TRUE)
      resM <- sparseMVN::rmvn.sparse(n, numeric(length(a)), Cholesky(prec))


      # calculate residual variance of y
      Sigma2 <- cbind(as.matrix(a), Sigma)
      SigmaXM <- rbind(t(as.matrix(Matrix::c.sparseVector(1, a))), Sigma2)
      pathsToY <- Matrix::c.sparseVector(dir, b)
      vary <- as.numeric(t(pathsToY) %*% SigmaXM %*% pathsToY)
      resy <- vary/r2y - vary # residual variance of y

    } else {

      sigmaXM <- diag(ncol(Sigma) + 1)
      sigmaXM[-1,-1] <- Sigma
      sigmaXM[ 1,-1] <- a
      sigmaXM[-1, 1] <- a

      # Calculate residual variance of y
      pathsToY <- c(dir, b)
      vary <- t(pathsToY) %*% sigmaXM %*% pathsToY # propagation of error
      resy <- vary/r2y - vary # residual variance of y

      # Simulate the residuals of M
      if (!empirical) resM <- MASS::mvrnorm(n, numeric(ncol(psi)), psi)

    }
  }


  # simulate the dataset
  if (!empirical) {
    x <- rnorm(n)
    M <- forma(x %*% t(a)) + resM
    y <- formb(M %*% b) + formb(dir * x) + rnorm(n, sd = sqrt(resy))
    if (scaley) y <- scale(y)
    if (inherits(M, "dgeMatrix")) M <- matrix(M, n)
    if (inherits(y, "dgeMatrix")) y <- as.numeric(y)
    return(data.frame(x = x, M = M, y = y))
  } else {
    if (missing(Sigma)) stop("Empirical not allowed without Sigma")
    # calculate full sigma
    nvars <- length(a) + 2
    fullSigma <- diag(nvars)
    fullSigma[-nvars, -nvars] <- sigmaXM
    fullSigma[nvars,-nvars] <- fullSigma[-nvars, nvars] <- pathsToY %*% sigmaXM
    fullSigma[nvars, nvars] <- vary + resy

    # now it's simple to generate
    dmat <- matrix(rnorm(n*nvars), n)

    d <- as.data.frame(empiricalTransform(dmat, fullSigma))
    colnames(d) <- c("x", paste0("M.", 1:length(a)), "y")
    return(d)
  }
}

#' @keywords internal
empiricalTransform <- function(X, Sigma) {
  # Force covariance matrix Sigma on X as per Mair, Satorra et al.
  return(X %*% expm::sqrtm(solve(cov(X))) %*% expm::sqrtm(Sigma))
}
