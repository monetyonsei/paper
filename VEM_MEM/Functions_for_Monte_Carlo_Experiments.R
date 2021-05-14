### R functions for Monte Carlo experiment

source("./MEM.R")

gen.theta <- function(N, d, eta, alpha)
{
  D <- d*(d+1)/2
  theta <- matrix(rnorm(N*D), N, D) * matrix(1/sqrt(alpha), 
                                             nrow = N, ncol = D,
                                             byrow = TRUE)
  theta <- matrix(eta, nrow = N, ncol = D, byrow = TRUE) + theta
  
  return(theta)
}

gen.data <- function(N, d, Tn, theta)
{
  x <- matrix(NA, N*Tn, d)
  subj <- rep(1:N, Tn)
  for ( n in 1:N ) {
    h <- theta[n, 1:d]
    J <- ivech(theta[n, -(1:d)], d)
    x[subj == n,] <- rBoltzmann(Tn, h, J)
  }
  
  return(list(x = x, subj = subj))
}

sim <- function(M, N, d, Tn)
{
  # M: the number of replication
  # N: the number of subjects
  # d: the number of nodes
  # Tn: the length of time series for each subject
  
  D <- 0.5*d*(d + 1)
  err.MLE <- err.BMEM <- err.VEM <- matrix(NA, M, N)
  
  for (r in 1:M) {
    
    cat("Iteration ", r, "\n")
    
    # Hyperparameters
    eta <- rnorm(D, sd = 0.1) 
    alpha.h <- 5; alpha.J <- 25
    alpha <- c(rep(alpha.h, d), rep(alpha.J, D - d)) 
    
    # Model parameters
    theta <- gen.theta(N, d, eta, alpha)
    
    # Data
    dat <- gen.data(N, d, Tn, theta)
    x <- dat$x
    subj <- dat$subj
    
    # MLE
    res.mle <- matrix(NA, N, D)
    for ( n in 1:N ) 
      res.mle[n,] <- MLE.Boltzmann(x[subj == n,])$theta
    se.mle <- (res.mle - theta)^2
    
    # BMEM
    res.BMEM <- BMEM.Boltzmann(x, N, subj)
    se.BMEM <- (theta - res.BMEM$mu)^2
    
    # VEM 
    res.VEM <- VEM.Boltzmann(x, N, subj, doplot = FALSE, dotrace = FALSE)
    se.VEM <- (theta - res.VEM$mu)^2
    
    err.MLE[r,] <- rowSums(se.mle)
    err.VEM[r,] <- rowSums(se.VEM)
    err.BMEM[r,] <- rowSums(se.BMEM)
  }
  
  return(list(err.MLE = err.MLE, err.BMEM = err.BMEM, err.VEM = err.VEM)) # list of M by N matrices
}
