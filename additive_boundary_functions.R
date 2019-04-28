# Boundary function for the first subfamily
## Input
### x: boundary value to be solved
### alpha: one-sided significance level
### t: information time
### gamma: parameter gamma
## Output
### Difference between alpha and the actual type I error
boundary_F1_function <- function(x, alpha, t, gamma){
  if (gamma < 0.5 | gamma >= 1) {
    stop("gamma should be between [0.5, 1)")
  } else {
    K <- length(t)
    cr <- cr_function(t)
    y <- 1 - alpha - pmvnorm(upper = x / sqrt(t) -
                               qnorm(1 - gamma) * sqrt(1 - t) / sqrt(t),
                             corr = cr, algorithm = Miwa(steps = 128))
    return(y)
  }
}

# Boundary function for the second subfamily
## Input
### x: boundary value to be solved
### alpha: one-sided significance level
### t: information time
### gamma: parameter gamma
## Output
### Difference between alpha and the actual type I error
boundary_F2_function <- function(x, alpha, t, gamma){
  if (gamma < 1 - pnorm(qnorm(1 - alpha / 2) / 2) | gamma >= 1) {
    stop(paste0("gamma should be between (",
                round(1 - pnorm(qnorm(1 - alpha / 2)/2), 3), ", 1)"))
  } else {
    K <- length(t)
    cr <- cr_function(t)
    y <- 1 - alpha - pmvnorm(upper = x / sqrt(t) -
                               qnorm(1 - gamma) * (1 - t) / sqrt(t),
                             corr = cr, algorithm = Miwa(steps = 128))
    return(y)
  }
}

# Boundary function for the third subfamily
## Input
### x: boundary value to be solved
### alpha: one-sided significance level
### t: information time
### gamma: parameter gamma
## Output
### Difference between alpha and the actual type I error
boundary_F3_function <- function(x, alpha, t, gamma){
  if (gamma <= alpha / 2 | gamma >= 1) {
    stop(paste0("gamma should be between (", alpha / 2, ", 1)"))
  } else {
    K <- length(t)
    cr <- cr_function(t)
    y <- 1 - alpha - pmvnorm(upper = x / sqrt(t) -
                               qnorm(1 - gamma) * (1 - sqrt(t)) / sqrt(t),
                             corr = cr, algorithm = Miwa(steps = 128))
    return(y)
  }
}

# Boundary function for the combined subfamily of the first and the second subfamilies
## Input
### x: boundary value to be solved
### alpha: one-sided significance level
### t: information time
### gamma: parameter gamma
## Output
### Difference between alpha and the actual type I error
boundary_combined_function <- function(x, alpha, t, gamma){
  if (gamma < 1 - pnorm(qnorm(1 - alpha / 2) / 2) | gamma >= 1) {
    stop(paste0("gamma should be between (",
                round(1 - pnorm(qnorm(1 - alpha / 2)/2), 3), ", 1)"))
  } else {
    K <- length(t)
    cr <- cr_function(t)
    if (gamma >= 0.5 & gamma < 1) {
      y <- 1 - alpha - pmvnorm(upper = x / sqrt(t) -
                                 qnorm(1 - gamma) * sqrt(1 - t) / sqrt(t),
                               corr = cr, algorithm = Miwa(steps = 128))
    } else {
      y <- 1 - alpha - pmvnorm(upper = x / sqrt(t) -
                                 qnorm(1 - gamma) * (1 - t) / sqrt(t),
                               corr = cr, algorithm = Miwa(steps = 128))
    }
    return(y)
  }
}

# Error spending function for the first subfamily
## Input
### alpha: one-sided significance level
### t: information time
### gamma: parameter gamma
## Output
### Vector of cumulative error spent
esf_F1_function <- function(alpha, t, gamma) {
  if (gamma < 0.5 | gamma >= 1) {
    stop("gamma should be between [0.5, 1)")
  } else {
    y <- 2 - 2 * pnorm(qnorm(1 - alpha / 2) / sqrt(t) -
                         qnorm(1 - gamma) * sqrt(1 - t) / sqrt(t))
    return(y)
  }
}

# Error spending function for the second subfamily
## Input
### alpha: one-sided significance level
### t: information time
### gamma: parameter gamma
## Output
### Vector of cumulative error spent
esf_F2_function <- function(alpha, t, gamma) {
  if (gamma < 1 - pnorm(qnorm(1 - alpha / 2) / 2) | gamma >= 1) {
    stop(paste0("gamma should be between (",
                round(1 - pnorm(qnorm(1 - alpha / 2)/2), 3), ", 1)"))
  } else {
    y <- 2 - 2 * pnorm(qnorm(1 - alpha / 2) / sqrt(t) -
                         qnorm(1 - gamma) * (1 - t) / sqrt(t))
    return(y)
  }
}

# Error spending function for the combined subfamily of the first and the second subfamilies
## Input
### alpha: one-sided significance level
### t: information time
### gamma: parameter gamma
## Output
### Vector of cumulative error spent
esf_combine_function <- function(alpha, t, gamma) {
  if (gamma < 1 - pnorm(qnorm(1 - alpha / 2) / 2) | gamma >= 1) {
    stop(paste0("gamma should be between (",
                round(1 - pnorm(qnorm(1 - alpha / 2)/2), 3), ", 1)"))
  } else {
    if (gamma < 0.5) {
      y <- 2 - 2 * pnorm(qnorm(1 - alpha / 2) / sqrt(t) -
                           qnorm(1 - gamma) * (1 - t) / sqrt(t))
    } else {
      y <- 2 - 2 * pnorm(qnorm(1 - alpha / 2) / sqrt(t) -
                           qnorm(1 - gamma) * sqrt(1 - t) / sqrt(t))
    }
    return(y)
  }
}

# Error spending function for the third subfamily
## Input
### alpha: one-sided significance level
### t: information time
### gamma: parameter gamma
## Output
### Vector of cumulative error spent
esf_F3_function <- function(alpha, t, gamma) {
  if (gamma <= alpha / 2 | gamma >= 1) {
    stop(paste0("gamma should be between (", alpha / 2, ", 1)"))
  } else {
    y <- 2 - 2 * pnorm(qnorm(1 - alpha / 2) / sqrt(t) -
                         qnorm(1 - gamma) * (1 - sqrt(t)) / sqrt(t))
    return(y)
  }
}

# Error spending function for the Pocock Lan-DeMets design
## Input
### alpha: one-sided significance level
### t: information time
## Output
### Vector of cumulative error spent
esf_Pocock_function <- function(alpha, t) {
  y <- alpha * log(1 + (exp(1) - 1) * t)
  return(y)
}

# Error spending function for the power family
## Input
### alpha: one-sided significance level
### t: information time
### rho: parameter rho
## Output
### Vector of cumulative error spent
esf_power_function <- function(alpha, t, rho) {
  y <- alpha * t^rho
  return(y)
}

#########################################################################
# Example call for Table 1
alpha <- 0.025
t <- c(0.25, 0.5, 0.75,1)
gamma <- 0.6
# Boundary approach
c_boundary <- uniroot(boundary_F1_function, interval = c(0.001, 10), alpha = alpha,
                      t = t, gamma = gamma)$root
a <- c_boundary / sqrt(t) - qnorm(1 - gamma) * sqrt(1 - t) / sqrt(t)
round(a, 3)
b <- CE_function(t, a)
round(b, 3)
# Error spending function approach
cumulative <- esf_F1_function(alpha, t, gamma)
a <- solver_boundary_esf_function(alpha, t, cumulative)
round(a, 3)
b <- CE_function(t, a)
round(b, 3)
