bvbox <- function(a, d = 7, mtitle = 'Bivariate Boxplot',
  method = 'robust', xlab='X', ylab='Y') {
  #
  #a is data matrix
  #d is constant(usually 7)
  #
  p <- length(a[1,  ])
  if (method == 'robust') {
    param <- biweight(a[, 1:2])
    m1 <- param[1]
    m2 <- param[2]
    s1 <- param[3]
    s2 <- param[4]
    r <- param[5]
  } else {
    m1 <- mean(a[, 1])
    m2 <- mean(a[, 2])
    s1 <- sqrt(var(a[, 1]))
    s2 <- sqrt(var(a[, 2]))
    r <- cor(a[, 1:2])[1, 2]
  }
  x <- (a[, 1] - m1)/s1
  y <- (a[, 2] - m2)/s2
  e <- sqrt((x * x + y * y - 2 * r * x * y)/(1 - r * r))
  e2 <- e * e
  em <- median(e)
  emax <- max(e[e2 < d * em * em])
  r1 <- em * sqrt((1 + r)/2)
  r2 <- em * sqrt((1 - r)/2)
  theta <- ((2 * pi)/360) * seq(0, 360, 3)
  xp <- m1 + (r1 * cos(theta) + r2 * sin(theta)) * s1
  yp <- m2 + (r1 * cos(theta) - r2 * sin(theta)) * s2
  r1 <- emax * sqrt((1 + r)/2)
  r2 <- emax * sqrt((1 - r)/2)
  theta <- ((2 * pi)/360) * seq(0, 360, 3)
  xpp <- m1 + (r1 * cos(theta) + r2 * sin(theta)) * s1
  ypp <- m2 + (r1 * cos(theta) - r2 * sin(theta)) * s2
  maxxl <- max(xpp)
  minxl <- min(xpp)
  maxyl <- max(ypp)
  minyl <- min(ypp)
  b1 <- (r * s2)/s1
  a1 <- m2 - b1 * m1
  y1 <- a1 + b1 * minxl
  y2 <- a1 + b1 * maxxl
  b2 <- (r * s1)/s2
  a2 <- m1 - b2 * m2
  x1 <- a2 + b2 * minyl
  x2 <- a2 + b2 * maxyl
  maxx <- max(c(a[, 1], xp, xpp, x1, x2))
  minx <- min(c(a[, 1], xp, xpp, x1, x2))
  maxy <- max(c(a[, 2], yp, ypp, y1, y2))
  miny <- min(c(a[, 2], yp, ypp, y1, y2))
  plot(a[, 1], a[, 2],
    xlim = c(minx, maxx),
    ylim = c(miny, maxy),
    xlab = xlab,
    ylab = ylab,
    lwd = 2,
    pch = 1)
  lines(xp, yp, lwd = 2)
  lines(xpp, ypp, lty = 2, lwd = 2)
  segments(minxl, y1, maxxl, y2, lty = 3, lwd = 2)
  segments(x1, minyl, x2, maxyl, lty = 4, lwd = 2)
}
