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

mahalanobisdist <- function(matrix) {
  matrix <- data.matrix(usair)
  covarmatrix <- var(matrix)
  inversematrix <- solve(covarmatrix)
  col <- ncol(matrix)
  row <- nrow(matrix)
  means <- array()
  for (i in 1:col) {
    means[i] <- mean(matrix[, i])
  }
  substractedmatrix <- data.matrix(matrix)
  for (i in 1:col) {
    for (j in 1:row) {
      substractedmatrix[j, i] <- matrix[j, i] - means[i]
    }
  }
  sqrt(t(substractedmatrix) %*% inversematrix %*% substractedmatrix)
}

m <- mahalanobisdist(usair)
nrow(m)


# Mahalanobis Distance

covmatrix <- cov(usair)
invcov <- solve(covmatrix)
usairMahalanobis <- usair


## Xi-mu
for (i in 1:ncol(usairMahalanobis)){
  usairMahalanobis[,i] <- usairMahalanobis[,i]-mean(usairMahalanobis[,i])
}

usairMahalanobismatrix<- data.matrix(usairMahalanobis)

## di^2 = x'S^-1x   , x = xi-mu
mahalanobisA <- vector()
for (i in 1:nrow(usairMahalanobismatrix)){
  x <- t(usairMahalanobismatrix[i,]) %*% invcov %*% usairMahalanobismatrix[i,]
  mahalanobisA <- append(mahalanobisA,x)
}

names(mahalanobisA) <- rownames(usair2)
mahalanobisA
