biweight <- function(a, const1 = 9, const2 = 36, err = 0.0001) {
  #
  # a is data matrix with two cols.
  # const1 = common tuning constant
  # const2 = bivariate tuning constant
  # err = convergence criterion.
  #

  x <- a[,1]
  y <- a[,2]
  n <- length(x)
  mx <- median(x)
  my <- median(y)
  madx <- median(abs(x - mx))
  mady <- median(abs(y - my))
  if (madx != 0) {
    ux <- (x - mx) / (const1 * madx)
    ux1 <- ux[abs(ux) < 1]
    tx <- mx + (sum((x[abs(ux) < 1] - mx) * (1 - ux1 * ux1) ^ 2) / sum((1 - ux1 ^ 2) ^ 2))
    sx <- sqrt(n) * sqrt(sum((x[abs(ux) < 1] - mx) ^ 2 * (1 - ux1 * ux1) ^ 4)) / abs(sum((1 - ux1 * ux1) * (1 - 5 * ux1 * ux1)))
  } else {
    tx <- mx
    sx <- sum(abs(x - mx)) / n
  }
  if (mady != 0) {
    uy <- (y - my) / (const1 * mady)
    uy1 <- uy[abs(uy) < 1]
    ty <- my + (sum((y[abs(uy) < 1] - my) * (1 - uy1 * uy1) ^ 2) / sum((1 - uy1 ^ 2) ^ 2))
    sy <- sqrt(n) * sqrt(sum((y[abs(uy) < 1] - my) ^ 2 * (1 - uy1 * uy1) ^ 4)) / abs(sum((1 - uy1 * uy1) * (1 - 5 * uy1 * uy1)))
  } else {
    ty <- my
    sy <- sum(abs(y - my)) / n
  }

  z1 <- (y - ty) / sy + (x - tx) / sx
  z2 <- (y - ty) / sy - (x - tx) / sx
  mz1 <- median(z1)
  mz2 <- median(z2)
  madz1 <- median(abs(z1 - mz1))
  madz2 <- median(abs(z2 - mz2))
  
  if (madz1 != 0) {
    uz1 <- (z1 - mz1) / (const1 * madz1)
    uz11 <- uz1[abs(uz1) < 1]
    tz1 <- mz1 + (sum((z1[abs(uz1) < 1] - mz1) * (1 - uz11 * uz11) ^ 2) / sum((1 - uz11 ^ 2) ^ 2))
    sz1 <- sqrt(n) * sqrt(sum((z1[abs(uz1) < 1] - mz1) ^ 2 * (1 - uz11 * uz11) ^ 4)) / abs(sum((1 - uz11 * uz11) * (1 - 5 * uz11 * uz11)))
  } else {
    tz1 <- mz1
    sz1 <- sum(abs(z1 - mz1)) / n
  }
  
  if (mady != 0) {
    uz2 <- (z2 - mz2) / (const1 * madz2)
    uz21 <- uz2[abs(uz2) < 1]
    tz2 <- mz2 + (sum((z2[abs(uz2) < 1] - mz2) * (1 - uz21 * uz21) ^ 2) / sum((1 - uz21 ^ 2) ^ 2))
    sz2 <- sqrt(n) * sqrt(sum((z2[abs(uz2) < 1] - mz2) ^ 2 * (1 - uz21 * uz21) ^ 4)) / abs(sum((1 - uz21 * uz21) * (1 - 5 * uz21 * uz21)))
  } else {
    tz2 <- mz2
    sz2 <- sum(abs(z2 - mz2)) / n
  }

  esq <- ((z1 - tz1) / sz1) ^ 2 + ((z2 - tz2) / sz2) ^ 2
  w <- numeric(length = n)
  c2 <- const2

  for (i in 1:10) {
    w[esq < const2] <- (1 - esq[esq < const2] / const2) ^ 2
    w[esq >= const2] <- 0
    l <- length(w[w == 0])
    if (l < 0.5 * n) break
    else const2 <- 2 * const2
  }

  tx <- sum(w * x) / sum(w)
  sx <- sqrt(sum(w * (x - tx) ^ 2) / sum(w))
  ty <- sum(w * y) / sum(w)
  sy <- sqrt(sum(w * (y - ty) ^ 2) / sum(w))
  r <- sum(w * (x - tx) * (y - ty))/(sx * sy * sum(w))
  const2 <- c2
  wold <- w

  for (i in 1:100) {
    z1 <- ((y - ty) / sy + (x - tx) / sx) / sqrt(2 * (1 + r))
    z2 <- ((y - ty) / sy - (x - tx) / sx) / sqrt(2 * (1 + r))
    esq <- z1 * z1 + z2 * z2
    for (j in 1:10) {
      w[esq < const2] <- (1 - esq[esq < const2] / const2) ^ 2
      w[esq >= const2] <- 0
      l <- length(w[w == 0])
      if (l < 0.5 * n) break
      else const2 <- 2 * const2
    }
    tx <- sum(w * x) / sum(w)
    sx <- sqrt(sum(w * (x - tx) ^ 2) / sum(w))
    ty <- sum(w * y) / sum(w)
    sy <- sqrt(sum(w * (y - ty) ^ 2) / sum(w))
    r <- sum(w * (x - tx) * (y - ty)) / (sx * sy * sum(w))
    term <- sum((w - wold) ^ 2) / (sum(w) / n) ^ 2
    if (term <- err) break
    else {
      wold <- w
      const2 <- c2
    }
  }

  param <- c(tx, ty, sx, sy, r)
  param
}
