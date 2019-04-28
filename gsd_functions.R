library(mvtnorm)

# Function to generate the correlation matrix
## Input
### t: information time
## Output
### Correlation matrix
cr_function <- function(t){
  K <- length(t)
  cr <- diag(K)
  for (i in 1:K){
    for (j in i:K){
      cr[i, j] <- sqrt(t[i] / t[j])
    }
  }
  cr <- cr + t(cr) - diag(K)
  return(cr)
}

# Function to derive the stopping boundary using the cumulative error spent
## Input
### alpha: one-sided significance level
### t: information time
### cumulative: cumulative error spent with the same length as t
## Output
### Vector of boundary values
solver_boundary_esf_function <- function(alpha, t, cumulative){
  K <- length(t)
  cr <- cr_function(t)
  solver <- function(x, cumu, cr = cr, c_past){
    z <- c(c_past, x)
    return(1 - cumu - pmvnorm(upper = z, corr = cr[1:length(z), 1:length(z)],
                              algorithm = Miwa(steps = 128))
    )
  }
  c_boundary <- rep(0,K)
  c_boundary[1] <- min(qnorm(1 - cumulative[1]), 10)
  if (K>1){
    for (k in 2:K){
      a <- uniroot(solver, interval=c(0.001, 10), cumu = cumulative[k],
                   cr = cr, c_past = c_boundary[1:(k-1)])$root
      c_boundary[k] <- min(a, 10)
    }
  }
  return(c_boundary)
}

# Function to calculate the conditional error rate
## Input
### t: information time
### c_boundary: boundary values
## Output
### Vecotr of the conditional error rate
CE_function <- function(t, c_boundary) {
  K <- length(t)
  y <- 1 - pnorm((c_boundary[K] - c_boundary * sqrt(t)) / sqrt(1 - t))
  return(y[-K])
}
