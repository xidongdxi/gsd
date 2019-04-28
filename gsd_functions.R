library(mvtnorm)

# Function to generate the correlation matrix
# t: information or information time
cr_function <- function(t){
  # Number of looks
  K <- length(t)
  # Correlation matrix based on planned information
  cr <- diag(K)
  for (i in 1:K){
    for (j in i:K){
      cr[i, j] <- sqrt(t[i] / t[j])
    }
  }
  cr <- cr + t(cr) - diag(K)
  return(cr)
}
