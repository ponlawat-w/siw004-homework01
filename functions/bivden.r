bivden <- function(x, y, ngridx = 30, ngridy = 30, constant.x = 1, constant.y = 1) {
  
  # x and y are vectors containing the bivariate data
  # ngridx and ngridy are the number of points in the grid
  mx <- mean(x)
  sdx <- sqrt(var(x))
  my <- mean(y)
  sdy <- sqrt(var(y))
  # scale x and y before estimation
  x <- scale(x)
  y <- scale(y)

  den <- matrix(0, ngridx, ngridy)
  
  # find possible value for bandwidth
  n <- length(x)
  
  hx <- constant.x * n^(-0.2)
  hy <- constant.y * n^(-0.2)
  h <- hx * hy
  hsqrt <- sqrt(h)
  
  seqx <- seq(range(x)[1], range(x)[2], length = ngridx)
  seqy <- seq(range(y)[1], range(y)[2], length = ngridy)
  
  for (i in 1:n) {
    X <- x[i]
    Y <- y[i]
    xx <- (seqx - X) / hsqrt
    yy <- (seqy - Y) / hsqrt
    den <- den + outer(xx, yy, function (x, y) {
      exp(-0.5 * (x ^ 2 + y ^ 2))
    })
  }

  den <- den / (n * 2 * pi * h)
  seqx <- sdx * seqx + mx
  seqy <- sdy * seqy + my
  result <- list(seqx = seqx, seqy = seqy, den = den)
  result
}
