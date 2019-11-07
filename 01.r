# Load Data
usair <- source('./data/chap3usair.dat')$value

# See column names of the data
names(usair)
# Result: [1] "SO2"      "Neg.Temp" "Manuf"    "Pop"      "Wind"     "Precip"   "Days"

attach(usair)

# We select SO2, Pop, Precip, Wind
summary(SO2)
summary(Pop)
summary(Precip)
summary(Wind)

sd(SO2)
sd(Pop)
sd(Precip)
sd(Wind)

# Scatter Plot

# SO2 and Population
par(fig = c(0, 0.7, 0, 0.7))
plot(SO2, Pop, lwd = 2)
abline(lm(SO2~Pop), lwd = 2)
lines(lowess(SO2, Pop), lwd = 2)

## Draw a histogram
par(fig = c(0, 0.7, 0.65, 1), new = TRUE)
hist(SO2, lwd = 2) 

## Draw a boxplot
par(fig = c(0.65, 1, 0, 0.7), new = TRUE)
boxplot(Pop, lwd = 2) 

dev.off()

# SO2 and Precipitation
par(fig = c(0, 0.7, 0, 0.7))
plot(SO2, Precip, lwd = 2)
abline(lm(SO2~Precip), lwd = 2)
lines(lowess(SO2, Precip), lwd = 2)

## Draw a histogram
par(fig = c(0, 0.7, 0.65, 1), new = TRUE)
hist(SO2, lwd = 2) 

## Draw a boxplot
par(fig = c(0.65, 1, 0, 0.7), new = TRUE)
boxplot(Precip, lwd = 2) 

dev.off()

# SO2 and Wind
par(fig = c(0, 0.7, 0, 0.7))
plot(SO2, Wind, lwd = 2)
abline(lm(SO2~Wind), lwd = 2)
lines(lowess(SO2, Wind), lwd = 2)

## Draw a histogram
par(fig = c(0, 0.7, 0.65, 1), new = TRUE)
hist(SO2, lwd = 2) 

## Draw a boxplot
par(fig = c(0.65, 1, 0, 0.7), new = TRUE)
boxplot(Wind, lwd = 2) 

dev.off()

# Convex Hulls

par(mfrow = c(1, 3))

hull <- chull(SO2, Pop)
plot(SO2, Pop, pch = 1)
polygon(SO2[hull], Pop[hull], density = 15, angle = 30)

hull <- chull(SO2, Precip)
plot(SO2, Precip, pch = 1)
polygon(SO2[hull], Precip[hull], density = 15, angle = 30)

hull <- chull(SO2, Wind)
plot(SO2, Wind, pch = 1)
polygon(SO2[hull], Wind[hull], density = 15, angle = 30)

dev.off()

# Chiplot

source('./functions/chiplot.r') # including source file of the function

chiplot(SO2, Pop, vlabs = c('SO2', 'Pop'))
chiplot(SO2, Precip, vlabs = c('SO2', 'Precip'))
chiplot(SO2, Wind, vlabs = c('SO2', 'Wind'))

# bvbox

source('./functions/bvbox.r') # including source file of the function
source('./functions/biweight.r')

par(mfrow = c(1, 3))
bvbox(cbind(SO2, Pop), xlab = 'SO2', ylab = 'Pop')
bvbox(cbind(SO2, Precip), xlab = 'SO2', ylab = 'Precip')
bvbox(cbind(SO2, Wind), xlab = 'SO2', ylab = 'Wind')

# bivden

source('./functions/bivden.r')

par(mfrow = c(1, 3))
den1 <- bivden(SO2, Pop)
persp(den1$seqx, den1$seqy, den1$den,
  xlab = 'SO2', ylab = 'Pop', zlab = 'Density',
  lwd = 2)
den2 <- bivden(SO2, Precip)
persp(den2$seqx, den2$seqy, den2$den,
  xlab = 'SO2', ylab = 'Precip', zlab = 'Density',
  lwd = 2)
den3 <- bivden(SO2, Wind)
persp(den3$seqx, den3$seqy, den3$den,
  xlab = 'SO2', ylab = 'Wind', zlab = 'Density',
  lwd = 2)

# pairs
usair2 <- data.frame(SO2, Pop, Precip, Wind)
usair2
pairs(usair2, panel = function(x, y) {
  abline(lsfit(x, y)$coef, lwd = 2)
  lines(lowess(x, y), lty = 2, lwd = 2)
  points(x, y)
})

# conditional plot

coplot(SO2~Precip | Wind, panel = function(x, y, col, pch) {
  panel.smooth(x, y, span = 1)
})

dev.off()

##----------##

# Euclidean distances
dist2full <- function(dis) {
  n <- attr(dis, 'Size')
  full <- matrix(0, n, n)
  full[lower.tri(full)] <- dis
  return(full + t(full))
}

dis.matrix <- dist2full(dist(usair))
dis.matrix

names(usair)

# Normalised Distances
## 1) find standard variation
std <- c(
  sd(usair[,1]),
  sd(usair[,2]),
  sd(usair[,3]),
  sd(usair[,4]),
  sd(usair[,5]),
  sd(usair[,6]),
  sd(usair[,7])
)
std

## 2) perform sweep() to normalise values in matrix
usair.std <- sweep(usair, 2, std, FUN = '/')
usair.std

## 3) Find the distance
dis <- dist(usair.std)
dis.matrix <- dist2full(dis)
round(dis.matrix, digits = 2)

# Mahalanobis Distance

mahalanobisdist <- function(data) {
  covmatrix <- cov(data)
  invcov <- solve(covmatrix)
  dataMahalanobis <- data

  ## Xi-mu
  for (i in 1:ncol(dataMahalanobis)) {
    dataMahalanobis[,i] <- dataMahalanobis[,i] - mean(dataMahalanobis[,i])
  }

  dataMahalanobismatrix <- data.matrix(dataMahalanobis)

  ## di^2 = x'S^-1x   , x = xi-mu
  mahalanobisA <- vector()
  for (i in 1:nrow(dataMahalanobismatrix)) {
    x <- t(dataMahalanobismatrix[i,]) %*% invcov %*% dataMahalanobismatrix[i,]
    mahalanobisA <- append(mahalanobisA, x)
  }

  names(mahalanobisA) <- rownames(data)
  mahalanobisA
}

usair.mahalanobis <- mahalanobisdist(usair)
usair.mahalanobis

# QQ Plot

par(mfrow = c(3, 3))
qqnorm(usair[,1], xlab = 'SO2', ylab = 'Ordered observations')
qqline(usair[,1])
qqnorm(usair[,2], xlab = 'Neg.Temp', ylab = 'Ordered observations')
qqline(usair[,2])
qqnorm(usair[,3], xlab = 'Manuf', ylab = 'Ordered observations')
qqline(usair[,3])
qqnorm(usair[,4], xlab = 'Pop', ylab = 'Ordered observations')
qqline(usair[,4])
qqnorm(usair[,5], xlab = 'Wind', ylab = 'Ordered observations')
qqline(usair[,5])
qqnorm(usair[,6], xlab = 'Precip', ylab = 'Ordered observations')
qqline(usair[,6])
plot(x = vector())
qqnorm(usair[,7], xlab = 'Days', ylab = 'Ordered observations')
qqline(usair[,7])

dev.off()

# Chiplot

chiplot(Precip, Neg.Temp, vlabs = c('Precip', 'Neg.Temp'))
chiplot(Precip, Wind, vlabs = c('Precip', 'Wind'))
chiplot(Manuf, Pop, vlabs = c('Manuf', 'Pop'))

dev.off()
