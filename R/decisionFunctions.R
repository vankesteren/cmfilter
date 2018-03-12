# This file contains decision functions (not exported)

#' the difference in coefficients decision function (correlation - partial cor)
#'
#' @importFrom stats cor qnorm
#'
#' @keywords internal
corMinusPartCor <- function(x, m, y, p.value = 0.1) {
  n <- length(x)
  rxy <- cor(x, y)
  rym <- cor(y, m)
  rxm <- cor(x, m)
  rxy2 <- rxy^2
  rym2 <- rym^2
  rxm2 <- rxm^2
  rxy.m <- (rxy - rym*rxm) / sqrt((1 - rym2)*(1 - rxm2))
  rdif <- rxy - rxy.m

  if (rdif == 0) return(FALSE)

  partder <- c(
    (rym - rxm*rxy) / (sqrt(1 - rym2) * (1 - rxm2)^(3/2)),
    1 - 1/sqrt((1 - rym2)*(1 - rxm2)),
    (rxm - rxy*rym) / (sqrt(1 - rxm2) * (1 - rym2)^(3/2))
  )

  # Then the large-sample variances
  rvars <- c(
    (1 - rxm2)^2 / n, # var(rxm)
    (1 - rxy2)^2 / n, # var(rxy)
    (1 - rym2)^2 / n # var(rmy)
  )

  # Create the variance-covariance matrix
  vcov <- diag(rvars)

  c <- (1 - rym2 - rxm2 - rxy2) / 2 # constant term in all covs
  vcov[1,2] <- vcov[2,1] <- ((2*rym - rxm*rxy)*c + rym^3) / n
  vcov[1,3] <- vcov[3,1] <- ((2*rxy - rxm*rym)*c + rxy^3) / n
  vcov[2,3] <- vcov[3,2] <- ((2*rxm - rxy*rym)*c + rxm^3) / n

  seOlkinFinn <- sqrt(partder %*% vcov %*% partder)

  return(abs(rdif/seOlkinFinn) > qnorm(p.value/2, lower.tail = F))
}


#' the product of coefficients decision function
#'
#' @keywords internal
prodCoef <- function(x, m, y, p.value = 0.1, dir = TRUE) {
  n <- length(x)

  # first the alpha path
  cpx <-  crossprod(x)                                # cross product of x
  alpha <- solve(cpx, crossprod(x, m))                # alpha path
  res_m <- m - x * c(alpha)                           # residual of m~x+0
  var_m <- as.numeric(crossprod(res_m) / (n - 1))     # rss variance
  var_a <- var_m/cpx                                  # variance of alpha

  # then the beta path
  if (dir) {
    mm <- cbind(x, m)                                 # model matrix
  } else {
    mm <- cbind(m)
  }
  cpm <- crossprod(mm)                                # cross product of mm
  beta <- solve(cpm, crossprod(mm, y))                # beta
  res_y <- y - mm %*% c(beta)                         # residual of y~m+x+0
  var_y <- as.numeric(crossprod(res_y) / (n - 1))     # rss variance
  var_b <- diag(var_y * chol2inv(chol(cpm)))          # variance of beta

  stat <- alpha * beta[2] # product of coefficients
  se <- sqrt(alpha^2 * var_b[2] + beta[2]^2 * var_a) #- var_a * var_b[2])

  if (is.na(stat) || is.na(se) || !is.numeric(stat) || !is.numeric(se)) {
    return(FALSE)
  } else {
    return(abs(stat/se) > qnorm(p.value/2, lower.tail = F))
  }
}
