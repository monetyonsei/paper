source("./MEM.R")
source("./Functions_for_Monte_Carlo_Experiments.R")

###############################
# To prepare for Figure 2 & 3 #
###############################

set.seed(12345)

d <- 10
D <- 0.5*d*(d + 1)
N <- 20

# Hyperparameters
eta <- rnorm(D, sd = 0.1) 
alpha.h <- 5; alpha.J <- 25
alpha <- c(rep(alpha.h, d), rep(alpha.J, D - d)) 

# Model parameters
theta <- gen.theta(N, d, eta, alpha)


###########
# Tn = 25 #
###########

# Data
Tn <- 25
dat <- gen.data(N, d, Tn, theta)
x <- dat$x
subj <- dat$subj

# MLE
res.mle25 <- matrix(NA, N, D)
for ( n in 1:N ) 
  res.mle25[n,] <- MLE.Boltzmann(x[subj == n,])$theta

# VEM 
res.vem25 <- VEM.Boltzmann(x, N, subj, doplot = FALSE, dotrace = FALSE)


###########
# Tn = 50 #
###########

# Data
Tn <- 50
dat <- gen.data(N, d, Tn, theta)
x <- dat$x
subj <- dat$subj

# MLE
res.mle50 <- matrix(NA, N, D)
for ( n in 1:N ) 
  res.mle50[n,] <- MLE.Boltzmann(x[subj == n,])$theta

# VEM 
res.vem50 <- VEM.Boltzmann(x, N, subj, doplot = FALSE, dotrace = FALSE)


############
# Tn = 100 #
############

# Data
Tn <- 100
dat <- gen.data(N, d, Tn, theta)
x <- dat$x
subj <- dat$subj

# MLE
res.mle100 <- matrix(NA, N, D)
for ( n in 1:N ) 
  res.mle100[n,] <- MLE.Boltzmann(x[subj == n,])$theta

# VEM 
res.vem100 <- VEM.Boltzmann(x, N, subj, doplot = FALSE, dotrace = FALSE)



############
# Figure 2 #
############

with(res.vem25, plot(ELBO, type = "l", lwd = 3, col = "navy",
                     xlab = "Iteration", main = "Tn = 25",
                     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5))

with(res.vem50, plot(ELBO, type = "l", lwd = 3, col = "navy",
                     xlab = "Iteration", main = "Tn = 50",
                     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5))

with(res.vem100, plot(ELBO, type = "l", lwd = 3, col = "navy",
                     xlab = "Iteration", main = "Tn = 100",
                     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5))



############
# Figure 3 #
############

n <- 1

lbub <- range(c(res.mle25[n,], res.vem25$mu[n,], theta[n,]))
plot(theta[n,], res.mle25[n,], 
     col = 2, pch = "+", xlim = lbub, ylim = lbub,
     xlab = "Model parameters", ylab = "Estimated", 
     main = "Tn = 25", 
     sub = paste("Distance: MLE = ", round(sqrt(sum((res.mle25[n,] - theta[n,])^2)), 2), 
                 ", VEM = ", round(sqrt(sum((res.vem25$mu[n,] - theta[n,])^2)), 2), 
                 sep = ""),
     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.25)
points(theta[n,], res.vem25$mu[n,], pch = 19, col = "navy")
abline(0, 1, col = 2, lty = 2)

lbub <- range(c(res.mle50[n,], res.vem50$mu[n,], theta[n,]))
plot(theta[n,], res.mle50[n,], 
     col = 2, pch = "+", xlim = lbub, ylim = lbub,
     xlab = "Model parameters", ylab = "Estimated", 
     main = "Tn = 50",
     sub = paste("Distance: MLE = ", round(sqrt(sum((res.mle50[n,] - theta[n,])^2)), 2), 
                 ", VEM = ", round(sqrt(sum((res.vem50$mu[n,] - theta[n,])^2)), 2), 
                 sep = ""),
     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.25)
points(theta[n,], res.vem50$mu[n,], pch = 19, col = "navy")
abline(0, 1, col = 2, lty = 2)

lbub <- range(c(res.mle100[n,], res.vem100$mu[n,], theta[n,]))
plot(theta[n,], res.mle100[n,], 
     col = 2, pch = "+", xlim = lbub, ylim = lbub,
     xlab = "Model parameters", ylab = "Estimated", 
     main = "Tn = 100",
     sub = paste("Distance: MLE = ", round(sqrt(sum((res.mle100[n,] - theta[n,])^2)), 2), 
                 ", VEM = ", round(sqrt(sum((res.vem100$mu[n,] - theta[n,])^2)), 2), 
                 sep = ""),
     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.25)
points(theta[n,], res.vem100$mu[n,], pch = 19, col = "navy")
abline(0, 1, col = 2, lty = 2)


############
# Figure 4 #
############

load("sim_Tn25.RData")

VEM <- rowMeans(resN20d05$err.VEM)
MLE <- rowMeans(resN20d05$err.MLE)
lbub <- c(0.5, 500)
plot(MLE, VEM, xlim = lbub, ylim = lbub, pch = "x",
     cex.lab = 1.5, cex.axis = 1.5, log = "xy")
abline(0, 1, col = 2, lty = 2)

VEM <- rowMeans(resN40d05$err.VEM)
MLE <- rowMeans(resN40d05$err.MLE)
points(MLE, VEM, col = 2, pch = "x")

VEM <- rowMeans(resN60d05$err.VEM)
MLE <- rowMeans(resN60d05$err.MLE)
points(MLE, VEM, col = 4, pch = "x")

VEM <- rowMeans(resN20d10$err.VEM)
MLE <- rowMeans(resN20d10$err.MLE)
points(MLE, VEM, pch = "+")

VEM <- rowMeans(resN40d10$err.VEM)
MLE <- rowMeans(resN40d10$err.MLE)
points(MLE, VEM, col = 2, pch = "+")

VEM <- rowMeans(resN60d10$err.VEM)
MLE <- rowMeans(resN60d10$err.MLE)
points(MLE, VEM, col = 4, pch = "+")


load("sim_Tn50.RData")

VEM <- rowMeans(resN20d05$err.VEM)
MLE <- rowMeans(resN20d05$err.MLE)
points(MLE, VEM, pch = "x")

VEM <- rowMeans(resN40d05$err.VEM)
MLE <- rowMeans(resN40d05$err.MLE)
points(MLE, VEM, col = 2, pch = "x")

VEM <- rowMeans(resN60d05$err.VEM)
MLE <- rowMeans(resN60d05$err.MLE)
points(MLE, VEM, col = 4, pch = "x")

VEM <- rowMeans(resN20d10$err.VEM)
MLE <- rowMeans(resN20d10$err.MLE)
points(MLE, VEM, pch = "+")
VEM <- rowMeans(resN40d10$err.VEM)
MLE <- rowMeans(resN40d10$err.MLE)
points(MLE, VEM, col = 2, pch = "+")
VEM <- rowMeans(resN60d10$err.VEM)
MLE <- rowMeans(resN60d10$err.MLE)
points(MLE, VEM, col = 4, pch = "+")


load("sim_Tn100.RData")

VEM <- rowMeans(resN20d05$err.VEM)
MLE <- rowMeans(resN20d05$err.MLE)
points(MLE, VEM, pch = "x")

VEM <- rowMeans(resN40d05$err.VEM)
MLE <- rowMeans(resN40d05$err.MLE)
points(MLE, VEM, col = 2, pch = "x")

VEM <- rowMeans(resN60d05$err.VEM)
MLE <- rowMeans(resN60d05$err.MLE)
points(MLE, VEM, col = 4, pch = "x")

VEM <- rowMeans(resN20d10$err.VEM)
MLE <- rowMeans(resN20d10$err.MLE)
points(MLE, VEM, pch = "+")
VEM <- rowMeans(resN40d10$err.VEM)
MLE <- rowMeans(resN40d10$err.MLE)
points(MLE, VEM, col = 2, pch = "+")
VEM <- rowMeans(resN60d10$err.VEM)
MLE <- rowMeans(resN60d10$err.MLE)
points(MLE, VEM, col = 4, pch = "+")

legend("topleft", c("N = 20", "N = 40", "N = 60"), col = c(1, 2, 4), pch = rep(19, 3))


############
# Figure 5 #
############

# Tn = 25 ############################
load("sim_Tn25.RData")

# d = 5
VEM <- rowMeans(resN20d05$err.VEM)
MLE <- rowMeans(resN20d05$err.MLE)
plot(density(MLE/VEM), log = "x",
     lwd = 3, xlim = c(5, 100), ylim = c(0, 0.15),
     main = "Tn = 25, d = 5", xlab = "Ratio",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
VEM <- rowMeans(resN40d05$err.VEM)
MLE <- rowMeans(resN40d05$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 2)
VEM <- rowMeans(resN60d05$err.VEM)
MLE <- rowMeans(resN60d05$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 4)
legend("topleft", c("N = 20", "N = 40", "N = 60"), 
       col = c(1, 2, 4), lwd = c(3, 3))

# d = 10
VEM <- rowMeans(resN20d10$err.VEM)
MLE <- rowMeans(resN20d10$err.MLE)
plot(density(MLE/VEM), log = "x",
     lwd = 3, xlim = c(5, 100), ylim = c(0, 0.15),
     main = "Tn = 25, d = 10", xlab = "Ratio",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
VEM <- rowMeans(resN40d10$err.VEM)
MLE <- rowMeans(resN40d10$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 2)
VEM <- rowMeans(resN60d10$err.VEM)
MLE <- rowMeans(resN60d10$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 4)
legend("topleft", c("N = 20", "N = 40", "N = 60"), 
       col = c(1, 2, 4), lwd = c(3, 3))

# Tn = 50 ############################
load("sim_Tn50.RData")

# d = 5
VEM <- rowMeans(resN20d05$err.VEM)
MLE <- rowMeans(resN20d05$err.MLE)
plot(density(MLE/VEM), log = "x",
     lwd = 3, xlim = c(1, 30), ylim = c(0, 0.5),
     main = "Tn = 50, d = 5", xlab = "Ratio",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
VEM <- rowMeans(resN40d05$err.VEM)
MLE <- rowMeans(resN40d05$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 2)
VEM <- rowMeans(resN60d05$err.VEM)
MLE <- rowMeans(resN60d05$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 4)
legend("topleft", c("N = 20", "N = 40", "N = 60"), 
       col = c(1, 2, 4), lwd = c(3, 3))

# d = 10
VEM <- rowMeans(resN20d10$err.VEM)
MLE <- rowMeans(resN20d10$err.MLE)
plot(density(MLE/VEM), log = "x",
     lwd = 3, xlim = c(1, 30), ylim = c(0, 0.5),
     main = "Tn = 50, d = 10", xlab = "Ratio",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
VEM <- rowMeans(resN40d10$err.VEM)
MLE <- rowMeans(resN40d10$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 2)
VEM <- rowMeans(resN60d10$err.VEM)
MLE <- rowMeans(resN60d10$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 4)
legend("topleft", c("N = 20", "N = 40", "N = 60"), 
       col = c(1, 2, 4), lwd = c(3, 3))

# Tn = 100 ############################
load("sim_Tn100.RData")

# d = 5
VEM <- rowMeans(resN20d05$err.VEM)
MLE <- rowMeans(resN20d05$err.MLE)
plot(density(MLE/VEM), log = "x",
     lwd = 3, xlim = c(1, 10), ylim = c(0, 1.5),
     main = "Tn = 100, d = 5", xlab = "Ratio",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
VEM <- rowMeans(resN40d05$err.VEM)
MLE <- rowMeans(resN40d05$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 2)
VEM <- rowMeans(resN60d05$err.VEM)
MLE <- rowMeans(resN60d05$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 4)
legend("topleft", c("N = 20", "N = 40", "N = 60"), 
       col = c(1, 2, 4), lwd = c(3, 3))

# d = 10
VEM <- rowMeans(resN20d10$err.VEM)
MLE <- rowMeans(resN20d10$err.MLE)
plot(density(MLE/VEM), log = "x",
     lwd = 3, xlim = c(1, 10), ylim = c(0, 1.5),
     main = "Tn = 100, d = 10", xlab = "Ratio",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
VEM <- rowMeans(resN40d10$err.VEM)
MLE <- rowMeans(resN40d10$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 2)
VEM <- rowMeans(resN60d10$err.VEM)
MLE <- rowMeans(resN60d10$err.MLE)
den <- density(MLE/VEM)
lines(den$x, den$y, lwd = 3, col = 4)
legend("topleft", c("N = 20", "N = 40", "N = 60"), 
       col = c(1, 2, 4), lwd = c(3, 3))


############
# Figure 6 #
############

# res <- sim(M = 100, N = 20, d = 10, Tn = 50)
load("comparison.RData")


par(mfrow = c(1, 3))
lbub <- c(min(rowMeans(res$err.VBMEM), rowMeans(res$err.VEM)) - 0.5, 
          max(rowMeans(res$err.VBMEM), rowMeans(res$err.VEM)) + 0.5)
plot(rowMeans(res$err.VBMEM), rowMeans(res$err.VEM),
     xlim = lbub, ylim = lbub, 
     pch = "+", col = "navy", 
     cex = 2, cex.lab = 1.5, cex.axis = 1.5,
     xlab = "BMEM", ylab = "VEM", log = "xy")
abline(0, 1, col = 2, lty = 2)

lbub <- sqrt(c(min(res$err.VBMEM, res$err.VEM), 
               max(res$err.VBMEM, res$err.VEM)))
for ( r in 1:2 ) {
  plot(sqrt(res$err.VBMEM[r,]), sqrt(res$err.VEM[r,]), 
       xlim = lbub, ylim = lbub, 
       pch = 19, col = "gray",
       xlab = "BMEM", ylab = "VEM", log = "xy",
       cex = 2, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
       main = paste("Virtual Dataset #", r, sep = ""))
  abline(0, 1, col = 2, lty = 2)
}
