chiplot <- function(x, y, vlabs = c('X', 'Y'), matrix = 'NO') {
  n <- length(x)
  ind <- numeric(length = n)
  for (i in 1:n) {
    for (j in (1:n)[-i])  if (x[i] > x[j] & y[i] > y[j]) ind[i] <- ind[i] + 1
  }
  ind <- ind / (n-1)

  ind1 <- numeric(length=n)
  for (i in 1:n) {
    for (j in (1:n)[-i])  if (x[i] > x[j]) ind1[i] <- ind1[i] + 1
  }
  ind1 <- ind1 / (n-1)

  ind2 <- numeric(length=n)
  for (i in 1:n) {
    for (j in (1:n)[-i])  if (y[i] > y[j]) ind2[i] <- ind2[i] + 1
  }
  ind2 <- ind2/(n-1)

  s <- sign((ind1 - 0.5) * (ind2 - 0.5))
  
  chi <- (ind-ind1*ind2) / sqrt(ind1 * (1 - ind1) * ind2 * (1 - ind2))
  lambda <- 4 * s * pmax((ind1 - 0.5) ^ 2, (ind2 - 0.5) ^ 2)
  thresh <- 4 * (1 / (n - 1) - 0.5) ^ 2

  if (matrix == 'NO') {
    par(mfrow = c(1,2))
    plot(x, y, xlab = vlabs[1], ylab = vlabs[2])
    plot(lambda[abs(lambda) < thresh],
      chi[abs(lambda) < thresh],
      ylim = c(-1, 1),
      xlab = 'lambda',
      ylab='Chi')
    abline(h = 1.78 / sqrt(n))
    abline(h = -1.78 / sqrt(n))
  }

  if (matrix=='YES') {
    plot(lambda[abs(lambda) < thresh], chi[abs(lambda) < thresh], ylim = c(-1, 1))
    abline(h = 1.78 / sqrt(n))
    abline(h = -1.78 / sqrt(n))
  }
}
