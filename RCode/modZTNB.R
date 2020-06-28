
### initial settings of two parameters size and mu in a negative binomial
### distribution for the numeric optimal searching function optim in R
mod.SIZE.INIT <- 1
mod.MU.INIT <- 0.5

### termination conditions for EM algorithm
mod.TOLERANCE <- 1e-10
mod.ITER.TOLERANCE <- 1e5

### To quicken the algorithm - if it has run for more than 1000 iterations, increase tolerance to 1e-8
mod.TOLERANCE.INC <- 1e-8 #this higher value seems to only change output with one percent or so, doesn't matter much
mod.ITER.TOLERANCE.INC <- 1000


### check the input histogram in an appropriate format
mod.checking.hist <- function(n) 
{
  if (ncol(n)!=2 || is.numeric(n[,1])==FALSE || is.numeric(n[,2])==FALSE) {
    stop("Input must be a two-column matrix")
  }
  ## the first column is the frequency i
  ## the second column is the number of species represented i times 
  ## in the sample
  freq <- n[, 1]
  num <- n[, 2]
  
  ## check whether frequencies are at least one and the histogram is sorted
  ## based on frequencies
  for (i in 1:length(freq))
    if (freq[i] <= 0 || freq[i] != floor(freq[i])) {
      stop("The first column must be positive integers!")
    } else if (num[i] < 0) {
      stop("The second column must be non negative")
    }
  else {
    if (i > 1 && freq[i - 1] >= freq[i])
      stop("The first column is not sorted in the ascending order")
  }
  
  return(n)
}

### density function of a zero-truncated negative binomial distribution
### size and mu are two parameters for the negative binomial
mod.dztnb <- function(x, size, mu, log = FALSE)
{
  ## the density of x in negative binomial
  p <- dnbinom(x, size = size, mu = mu, log = log)
  
  ## set zeros in x with zero probability
  if (log == FALSE) {
    p[ which(x == 0) ] <- 0
  } else {
    p[ which(x == 0) ] <- -Inf
  }
  
  ## the density of non-zero in negative binomial
  q <- 1 - dnbinom(0, size = size, mu = mu)
  
  ## normalize all non-zero values in negrative binomial
  if (log == FALSE) {
    return( p/q )
  } else {
    return( p - log(q) )
  }
}


### zero-truncated negative loglikelihood
mod.ztnb.minus.loglikelihood <- function(n, size, mu)
{
  prob <- mod.dztnb(n[, 1], size, mu, log = TRUE)
  
  ## negative loglikelihood
  prob <- -prob
  return( prob %*% n[, 2] )
}


### calculate the negative binomial loglikelihood
### zero.count is the number of unobserved species
mod.nb.loglikelihood <- function(n, zero.count, size, mu)
{
  ## loglikelihood for nonzero counts
  log.prob <- dnbinom(n[, 1], size = size, mu = mu, log = TRUE)
  loglikelihood <- log.prob %*% n[, 2]
  
  ## add loglikelihood for zero count
  log.zero.prob <- dnbinom(0, size = size, mu = mu, log = TRUE)
  loglikelihood <- loglikelihood + zero.count * log.zero.prob
  
  return(loglikelihood)
}


### EM algorithm to fit the histogram with a negative binomial distribution
### n is the histogram for observed species
### the number of unobserved species is the missing data
mod.preseqR.ztnb.em <- function(n, size=mod.SIZE.INIT, mu=mod.MU.INIT, incTol = mod.TOLERANCE.INC, iterIncTol = mod.ITER.TOLERANCE.INC)
{
  mod.checking.hist(n)
  
  n[, 1] <- as.numeric(n[, 1])
  n[, 2] <- as.numeric(n[, 2])
  zero.prob <- exp(dnbinom(0, size = size, mu = mu, log = TRUE))
  
  S <- sum(n[, 2])
  ## estimate the total number of species
  L <- S / ( 1 - zero.prob )
  
  ## expected the number of zero counts
  zero.counts <- L*zero.prob
  
  ## estimated mean and variance
  m <- (n[, 1] %*% n[, 2]) / L
  v <- ( (n[, 1] - as.vector(m))^2 %*% n[, 2] + m^2 * zero.counts )/(L - 1)
  
  ## target function f
  f <- function(x) {
    return( -mod.nb.loglikelihood(n, zero.counts, size = x, mu = m)/L )
  }
  
  ## derivative of f
  ## zero.counts is an external variable that are updated by the EM algorithm
  ## CHECK IT!
  gr <- function(x)
  {
    first.term <- ( digamma(x) * zero.counts +
                      digamma(n[, 1] + x) %*% n[, 2] )/L
    second.term <- digamma(x)
    third.term <- log(x) - log(x + m)
    result <- first.term - second.term + third.term
    # f is negative loglikelihood
    return(-result)
  }
  
  ## estimate size and mu based on first and second moments
  if (v > m) {
    res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
                 lower = 0.0001, upper = 10000)
  } else {
    res <- optim(size, f, gr, method = "L-BFGS-B",
                 lower = 0.0001, upper = 10000)
  }
  
  ## count the times of iteration
  iter <- as.double(1)
  
  ## initialize the negative loglikelihood
  loglikelihood.pre <- Inf
  
  ## zero-truncated loglikelihood
  loglikelihood <- mod.ztnb.minus.loglikelihood(n, res$par, m)
  
  
  iterations = 0
  
  ## EM algorithm
  while (( loglikelihood.pre - loglikelihood ) / S > mod.TOLERANCE &&
         iter < mod.ITER.TOLERANCE && 
         !(iter >= iterIncTol && ( loglikelihood.pre - loglikelihood ) / S <= incTol))
  {
    iterations = iterations + 1
#    if (iterations %% 100 == 0) {
#      print(paste0("ITER: ", iterations, " Tolerance val: ", ( loglikelihood.pre - loglikelihood ) / S ))
#    }
    ## update negative loglikelihood
    loglikelihood.pre <- loglikelihood
    
    ## update parameters
    size <- res$par
    mu <- m
    
    ### E-step: estimate the number of unobserved species
    
    zero.prob <- exp(dnbinom(0, size = size, mu = mu, log = TRUE))
    L <- S / ( 1 - zero.prob )
    zero.counts <- L * zero.prob
    m <- (n[, 1] %*% n[, 2])/L
    v <- ( (n[, 1] - as.vector(m))^2 %*% n[, 2] + m^2 * zero.counts )/(L - 1)
    
    ### M step: estimate the parameters size and mu
    
    if (v > m) {
      res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
                   lower = 0.0001, upper = 10000)
    } else {
      res <- optim(size, f, gr, method = "L-BFGS-B",
                   lower = 0.0001, upper = 10000)
    }
    #print(res$par)
    iter <- iter + 1
    loglikelihood <- mod.ztnb.minus.loglikelihood(n, res$par, m)
    #print(( loglikelihood.pre - loglikelihood ) / S)
  }
  #print(paste0("iterations: ", iterations))
  return(list(size = size, mu = mu, loglik = -loglikelihood.pre))
}


## fitting the negative binoimal distribution to the data by EM algorithm
mod.ztnb.rSAC <- function(n, r=1, size=mod.SIZE.INIT, mu=mod.MU.INIT, incTol = mod.TOLERANCE.INC, iterIncTol = mod.ITER.TOLERANCE.INC)
{
  mod.checking.hist(n)
  
  n[, 2] <- as.numeric(n[, 2])
  S <- sum(n[, 2])
  
  ## estimate parameters
  opt <- mod.preseqR.ztnb.em(n, size, mu, incTol, iterIncTol)
  size <- opt$size
  mu <- opt$mu
  
  ## the probability of a species observed in the initial sample
  p <- 1 - dnbinom(0, size = size, mu = mu)
  
  ## L is the estimated number of species in total
  L <- S / p
  
  f.rSAC <- function(t) {
    L * pnbinom(r - 1, size=size, mu=mu*t, lower.tail=FALSE)
  }
  f.rSAC
}


## Best practice - changed this function to make it much faster
mod.preseqR.rSAC <- function(n, r=1, mt=20, size=mod.SIZE.INIT, mu=mod.MU.INIT)
{
  para <- mod.preseqR.ztnb.em(n, incTol = 1e-5, iterIncTol = 200)
  shape <- para$size
  mu <- para$mu
  ## the population is heterogeneous
  ## because the coefficient of variation is large $1 / sqrt(shape)$
  if (shape <= 1) {
    f.rSAC <- ds.rSAC(n=n, r=r, mt=mt)
  } else {
    ## the population is close to be homogeneous
    ## the ZTNB approach is applied
    
    ## the probability of a species observed in the initial sample
    p <- 1 - dnbinom(0, size = size, mu = mu)
    ## L is the estimated number of species in total
    L <- sum(as.numeric(n[, 2])) / p
    ## ZTNB estimator
    f.rSAC <- function(t) {
      L * pnbinom(r - 1, size=size, mu=mu*t, lower.tail=FALSE)
    }
  }
  return(f.rSAC)
}

## Best practice - changed this function to make it much faster
## Also changed it so the size param is now used in the ztnb estimation
mod.preseqR.rSAC.fixed <- function(n, r=1, mt=20, size=mod.SIZE.INIT, mu=mod.MU.INIT)
{
  para <- mod.preseqR.ztnb.em(n, incTol = 1e-5, iterIncTol = 200)
  size <- para$size
  mu <- para$mu
  ## the population is heterogeneous
  ## because the coefficient of variation is large $1 / sqrt(shape)$
  if (size <= 1) {
    f.rSAC <- ds.rSAC(n=n, r=r, mt=mt)
  } else {
    ## the population is close to be homogeneous
    ## the ZTNB approach is applied
    
    ## the probability of a species observed in the initial sample
    p <- 1 - dnbinom(0, size = size, mu = mu)
    ## L is the estimated number of species in total
    L <- sum(as.numeric(n[, 2])) / p
    ## ZTNB estimator
    f.rSAC <- function(t) {
      L * pnbinom(r - 1, size=size, mu=mu*t, lower.tail=FALSE)
    }
  }
  return(f.rSAC)
}

#Test code for the c++ implementation:
#freq = c(1,2,3,4)
#counts = c(28,16,5,1)
#n = data.frame(freq=freq, counts=counts)
#f = mod.ztnb.rSAC(n, incTol = 1e-5, iterIncTol = 200)
#f(10)#gives 80.1079
#counts = c(1,0,0,0,0,0,0,0,0,1)
#counts = c(1,0,0,0,0,0,0,0,1,0,0,1,0,1,1,1,3,0,1,0,0,0,0,0,0,1,0,0,0,1)#this one is problematic in cpp
#freq = 1:length(counts)
#n = data.frame(freq=freq, counts=counts)
#f = mod.ztnb.rSAC(n, incTol = 1e-5, iterIncTol = 50)
#f(10)

#freq = 1:16
#counts = c(11,2,1,1,0,0,0,0,0,0,0,0,0,0,0,1)
#n = data.frame(freq=freq, counts=counts)
#f = mod.ztnb.rSAC(n, incTol = 1e-5, iterIncTol = 200)
#f(10)
#counts = c(74,24,3,1)
#freq=1:4

#counts = c(3,3,2,3,3,3,2,0,0,3,0,2,3,5,3,3,3,4,0,1,0,2,3,1,1,0,1)
#freq = 1:length(counts)

