####
# Fisher scoring for CFA/SEM
# November 9, 2013
# Author: Paul Stey
####

## Read in Mardia, Kent, and Bibby (1979) "Open-Closed Book" data
dset <- read.table('http://www3.nd.edu/~kyuan/courses/CS/data/mardia.dat', header = FALSE, nrows = 88)


# One-factor model 
# Using last 3 variables

X <- dset[, 3:5]

# initialize arrays to be filled
f_dot <- rep(NA, length = 6)
sigma_dot <- array(NA, dim = c(3, 3, 6))
f_2dot <- matrix(NA, nrow = 6, ncol = 6)


# function to fill f_dot vector (first derivatives)
dx <- function(X, sigma, sigma_dot) {
    J <- 2*nrow(sigma)
    S <- cov(X)

    # Invert sigma outside loop
    sigma_inv <- solve(sigma)

    # Fill f_dot with this loop
    for (j in 1:J) {
        f_dot[j] <- -sum(diag(sigma_inv %*% (S - sigma) %*% sigma_inv %*% sigma_dot[, , j]))
    }
    return(f_dot)
}



# function to fill f_2dot matrix (second derivates)
d2x <- function(X, sigma, sigma_dot) {
    
    J <- 2*nrow(sigma)
    K <- J

    # Invert sigma outside loops
    sigma_inv <- solve(sigma)

    # Fill f_2dot with this loop
    for (j in 1:J) {
        for (k in 1:K) {
            f_2dot[j, k] <- sum(diag(sigma_inv %*% sigma_dot[, , j] %*% sigma_inv %*% sigma_dot[, , k]))
        }
    }
    return(f_2dot)
}





# Function to return delta theta, which is our step size
del_theta <- function(X, start) {
    
    theta <- start
    J <- length(theta)
    K <- J
    lam <- theta[1:3]
    
    # initialize e as identity matrix
    e <- array(c(1, 0, 0, 0, 1, 0, 0, 0, 1), dim = c(3, 1, 3))
    
    # define psi from theta values
    psi <- matrix(c(theta[4], 0, 0, 0, theta[5], 0, 0, 0, theta[6]), nrow = 3, ncol = 3)

    # Compute model-implied covariance matrix
    sigma <- lam %*% t(lam) + psi


    # loop to fill sigma_dot for j = 1, 2, 3
    for (i in 1:3) {
        sigma_dot[, , i] <- e[, , i] %*% t(lam) + lam %*% t(e[, , i])
    }

    # loop to fill sigma_dot for j = 4, 5, 6
    for (j in 4:6) {
        sigma_dot[, , j] <- e[, , j-3] %*% t(e[, , j-3])
    }

    # Apply function for first derivative of F
    f_dot <- dx(X, sigma, sigma_dot)

    # Apply function for second derivative of F
    f_2dot <- d2x(X, sigma, sigma_dot)

    ## compute step size 
    delta_theta <- solve(f_2dot) %*% f_dot

    return(delta_theta)
}


# Function that will iterate until sufficiently small step size
scoring <- function(X, start, cutoff, maxiter) {
    
    theta <- start
    count <- 1
    maxstep <- 1e5

    # while loop to minimize step
    while (maxstep > cutoff & count < maxiter) {
        
        maxstep <- max(abs(del_theta(X, theta)))
        theta <- theta - del_theta(X, theta)

        print(paste(count, maxstep))
        
        count <- count + 1
    }
    return(theta)
}



# set our start values and run the function
theta0 <- c(1, 1, 1, 1, 1, 1)

score_cfa <- scoring(X, start = theta0, cutoff = 1e-5, maxiter = 100)



# clean up and print result
score_cfa <- round(score_cfa, digits = 2)

print(paste0('The factor loadings are: ', score_cfa[1], ', ', score_cfa[2], ', ', 'and ', score_cfa[3], '.  ', 'The unique variances are: ', score_cfa[4], ', ', score_cfa[5], ', ', 'and ', score_cfa[6], '.'))





# Check results against lavaan
library(lavaan)

mod1 <- '
f1 =~ V3 + V4 + V5
'

fm1 <- cfa(mod1, data = X, std.lv = TRUE)
summary(fm1)





##########
# Two-factor model
# Using all five variables
##########


# read complete data set into X matrix
X <- dset



# initialize arrays to be filled
f_dot <- rep(NA, length = 11)
sigma_dot <- array(NA, dim = c(5, 5, 11))
f_2dot <- matrix(NA, nrow = 11, ncol = 11)

# define first derivative of lam
lam_dot <- array(0, dim = c(5, 2, 5))
lam_dot[1, 1, 1] <- 1
lam_dot[2, 1, 2] <- 1
lam_dot[3, 2, 3] <- 1
lam_dot[4, 2, 4] <- 1
lam_dot[5, 2, 5] <- 1

## define first derivative of uniquenesses
e <- array(0, dim = c(5, 5, 5))

for(i in 1:5){
    e[i, i, i] <- 1
}



## define first derivative of factor covariance matrix
phi_dot <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)




# The dx() function fills the f_dot vector
dx <- function(X, sigma, sigma_dot) {
    
    J <- (2 * nrow(sigma))+ 1
    S <- cov(X)

    ## Invert sigma outside loop
    sigma_inv <- solve(sigma)

    ## Fill f_dot with this loop
    for (j in 1:J) {
        f_dot[j] <- -sum(diag(sigma_inv %*% (S - sigma) %*% sigma_inv %*% sigma_dot[, , j]))
    }
    return(f_dot)
}





###
# function to fill f_2dot matrix
###
d2x <- function(X, sigma, sigma_dot){
    
    J <- (2 * nrow(sigma)) + 1
    K <- J

    ## Invert sigma outside loops
    sigma_inv <- solve(sigma)

    ## Fill f_2dot with this loop
    for (j in 1:J) {
        for (k in 1:K) {
            f_2dot[j, k] <- sum(diag(sigma_inv %*% sigma_dot[, , j] %*% sigma_inv %*% sigma_dot[, , k]))
        }
    }
    return(f_2dot)
}







###
# Function to return delta theta, which is our step size
###
del_theta <- function(X, start){
    
    theta <- start
    J <- length(theta)
    K <- J
    
    # define lam
    lam <- matrix(c(theta[1], theta[2], 0, 0, 0, 0, 0, theta[3], theta[4], theta[5]), nrow = 5, ncol = 2)
    
    # define psi from theta values
    psi <- matrix(c(theta[7], 0, 0, 0, 0, 0, theta[8], 0, 0, 0, 0, 0, theta[9], 0, 0, 0, 0, 0, theta[10], 0, 0, 0, 0, 0, theta[11]), nrow = 5, ncol = 5)

    # define Phi from theta values
    Phi <- matrix(c(1, theta[6], theta[6], 1), nrow = 2, ncol = 2)

    # Compute model-implied covariance matrix
    sigma <- lam %*% Phi %*% t(lam) + psi
    
    # loop to fill sigma_dot for j = 1, 2, 3, 4, 5
    for (j in 1:5) {
        sigma_dot[, , j] <- lam_dot[, , j] %*% Phi %*% t(lam) + lam %*% Phi %*% t(lam_dot[, , j])
    }


    # Compute sigma_dot for j = 6
    sigma_dot[, , 6] <- lam %*% phi_dot %*% t(lam)


    # Loop to fill sigma_dot for j = 7, 8, 9, 10, 11
    for (j in 7:11) {
        sigma_dot[, , j] <- e[, , j-6] %*% t(e[, , j-6])
    }

    f_dot <- dx(X, sigma, sigma_dot)                        # first deriv. of f
    f_2dot <- d2x(X, sigma, sigma_dot)                      # 2nd deriv. of f
    delta_theta <- solve(f_2dot) %*% f_dot                  # compute step size

    return(delta_theta)
}


# The scoring() function will iterate until 
# sufficiently small step size
scoring <- function(X, start, cutoff, maxiter){
    
    theta <- start
    count <- 1
    maxstep <- 1e5

    ## while loop to minimize step
    while (maxstep > cutoff & count < maxiter) {
        maxstep <- max(abs(del_theta(X, theta)))
        theta <- theta - del_theta(X, theta)
        
        ## progress monitor
        print(paste(count, maxstep))

        ## update counter
        count <- count + 1
    }
    return(theta)
}




###
# set our start values and run the function
###
theta0 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

score_cfa <- scoring(X, start = theta0, cutoff = 1e-5, maxiter = 100)




###
# clean up and print result
###
score_cfa <- round(score_cfa, digits = 2)

print(paste0('The factor loadings are: ', score_cfa[1], ', ', score_cfa[2], ', ', score_cfa[3], ', ', score_cfa[4], ', ', 'and ', score_cfa[5], '.  ', 'The unique variances are: ', score_cfa[7], ', ', score_cfa[8], ', ', score_cfa[9], ', ', score_cfa[10], ', ', 'and ', score_cfa[11],'.  ', 'The correlation between factors is ', score_cfa[6], '.'))



###
# Check results against lavaan
###

library(lavaan)

mod1 <- '
f1 =~ V1 + V2
f2 =~ V3 + V4 + V5
'

fm1 <- cfa(mod1, data = X, std.lv = TRUE)
summary(fm1)

