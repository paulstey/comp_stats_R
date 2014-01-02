####################################################
# Computational Stats.
# Homework 9
# Fisher Scoring for CFA/SEM
# November 9, 2013
# Author: Paul Stey
####################################################




## Read in Mardia, Kent, and Bibby (1979) "Open-Closed Book" data
dset <- read.table('http://www3.nd.edu/~kyuan/courses/CS/data/mardia.dat', header = FALSE, nrows = 88)






###########################################
# One-factor model 
# Using last 3 variables
###########################################

X <- dset[, 3:5]

## initialize arrays to be filled
F.dot <- rep(NA, length = 6)
Sigma.dot <- array(NA, dim = c(3, 3, 6))
F2.dot <- matrix(NA, nrow = 6, ncol = 6)


###
# function to fill F.dot vector (first derivatives)
###
dF <- function(X, Sigma, Sigma.dot){
	## bookkeeping
	J <- 2*nrow(Sigma)
	S <- cov(X)

	## Invert Sigma outside loop
	Sigma.inv <- solve(Sigma)

	## Fill F.dot with this loop
	for(j in 1:J){
		F.dot[j] <- -sum(diag(Sigma.inv %*% (S - Sigma) %*% Sigma.inv %*% Sigma.dot[, , j]))
	}
	return(F.dot)
}




###
# function to fill F2.dot matrix (second derivates)
###
d2F <- function(X, Sigma, Sigma.dot){
	## bookkeeping
	J <- 2*nrow(Sigma)
	K <- J

	## Invert Sigma outside loops
	Sigma.inv <- solve(Sigma)

	## Fill F2.dot with this loop
	for(j in 1:J){
		for(k in 1:K){
			F2.dot[j, k] <- sum(diag(Sigma.inv %*% Sigma.dot[, , j] %*% Sigma.inv %*% Sigma.dot[, , k]))
		}
	}
	return(F2.dot)
}








###
# Function to return delta theta, which is our step size
###
DelTheta <- function(X, start){
	## bookkeeping
	theta <- start
	J <- length(theta)
	K <- J
	lambda <- theta[1:3]
	
	## initialize e as identity matrix
	e <- array(c(1, 0, 0, 0, 1, 0, 0, 0, 1), dim = c(3, 1, 3))
	
	## define Psi from theta values
	Psi <- matrix(c(theta[4], 0, 0, 0, theta[5], 0, 0, 0, theta[6]), nrow = 3, ncol = 3)

	## Compute model-implied covariance matrix
	Sigma <- lambda %*% t(lambda) + Psi


	## loop to fill Sigma.dot for j = 1, 2, 3
	for(i in 1:3){
		Sigma.dot[, , i] <- e[, , i] %*% t(lambda) + lambda %*% t(e[, , i])
	}

	## Now, loop to fill Sigma.dot for j = 4, 5, 6
	for(j in 4:6){
		Sigma.dot[, , j] <- e[, , j-3] %*% t(e[, , j-3])
	}

	## Apply function for first derivative of F
	F.dot <- dF(X, Sigma, Sigma.dot)

	## Apply function for second derivative of F
	F2.dot <- d2F(X, Sigma, Sigma.dot)

	## compute step size 
	delta.theta <- solve(F2.dot) %*% F.dot

	return(delta.theta)
}




###
# Function that will iterate until sufficiently small step size
###

Scoring <- function(X, start, cutoff, max.iter){
	## bookkeeping
	theta <- start
	count <- 1
	max.step <- 1e5

	## while loop to minimize step
	while(max.step > cutoff & count < max.iter){
		max.step <- max(abs(DelTheta(X, theta)))
		theta <- theta - DelTheta(X, theta)

		## progress monitor
		print(paste(count, max.step))
		
		## update counter
		count <- count + 1
	}
	return(theta)
}




###
# set our start values and run the function
###
theta0 <- c(1, 1, 1, 1, 1, 1)

score.cfa <- Scoring(X, start = theta0, cutoff = 1e-5, max.iter = 100)


###
# clean up and print result
###
score.cfa <- round(score.cfa, digits = 2)

print(paste0('The factor loadings are: ', score.cfa[1], ', ', score.cfa[2], ', ', 'and ', score.cfa[3], '.  ', 'The unique variances are: ', score.cfa[4], ', ', score.cfa[5], ', ', 'and ', score.cfa[6], '.'))





###
# Check results against lavaan
###
library(lavaan)

mod1 <- '
f1 =~ V3 + V4 + V5
'

fm1 <- cfa(mod1, data = X, std.lv = TRUE)
summary(fm1)














################################################
# Two-factor model
# Using all five variables
################################################


## read complete data set into X matrix
X <- dset



## initialize arrays to be filled
F.dot <- rep(NA, length = 11)
Sigma.dot <- array(NA, dim = c(5, 5, 11))
F2.dot <- matrix(NA, nrow = 11, ncol = 11)

## define first derivative of Lambda
Lambda.dot <- array(0, dim = c(5, 2, 5))
Lambda.dot[1, 1, 1] <- 1
Lambda.dot[2, 1, 2] <- 1
Lambda.dot[3, 2, 3] <- 1
Lambda.dot[4, 2, 4] <- 1
Lambda.dot[5, 2, 5] <- 1

## define first derivative of uniquenesses
e <- array(0, dim = c(5, 5, 5))

for(i in 1:5){
	e[i, i, i] <- 1
}



## define first derivative of factor covariance matrix
Phi.dot <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)






###
# function to fill F.dot vector
###
dF <- function(X, Sigma, Sigma.dot){
	## bookkeeping
	J <- (2 * nrow(Sigma))+ 1
	S <- cov(X)

	## Invert Sigma outside loop
	Sigma.inv <- solve(Sigma)

	## Fill F.dot with this loop
	for(j in 1:J){
		F.dot[j] <- -sum(diag(Sigma.inv %*% (S - Sigma) %*% Sigma.inv %*% Sigma.dot[, , j]))
	}
	return(F.dot)
}





###
# function to fill F2.dot matrix
###
d2F <- function(X, Sigma, Sigma.dot){
	## bookkeeping
	J <- (2 * nrow(Sigma)) + 1
	K <- J

	## Invert Sigma outside loops
	Sigma.inv <- solve(Sigma)

	## Fill F2.dot with this loop
	for(j in 1:J){
		for(k in 1:K){
			F2.dot[j, k] <- sum(diag(Sigma.inv %*% Sigma.dot[, , j] %*% Sigma.inv %*% Sigma.dot[, , k]))
		}
	}
	return(F2.dot)
}







###
# Function to return delta theta, which is our step size
###
DelTheta <- function(X, start){
	## bookkeeping
	theta <- start
	J <- length(theta)
	K <- J
	
	
	## define Lambda
	Lambda <- matrix(c(theta[1], theta[2], 0, 0, 0, 0, 0, theta[3], theta[4], theta[5]), nrow = 5, ncol = 2)


	## define Psi from theta values
	Psi <- matrix(c(theta[7], 0, 0, 0, 0, 0, theta[8], 0, 0, 0, 0, 0, theta[9], 0, 0, 0, 0, 0, theta[10], 0, 0, 0, 0, 0, theta[11]), nrow = 5, ncol = 5)

	## define Phi from theta values
	Phi <- matrix(c(1, theta[6], theta[6], 1), nrow = 2, ncol = 2)


	## Compute model-implied covariance matrix
	Sigma <- Lambda %*% Phi %*% t(Lambda) + Psi
	


	## loop to fill Sigma.dot for j = 1, 2, 3, 4, 5
	for(j in 1:5){
		Sigma.dot[, , j] <- Lambda.dot[, , j] %*% Phi %*% t(Lambda) + Lambda %*% Phi %*% t(Lambda.dot[, , j])
	}


	## Compute Sigma.dot for j = 6
	Sigma.dot[, , 6] <- Lambda %*% Phi.dot %*% t(Lambda)


	## Loop to fill Sigma.dot for j = 7, 8, 9, 10, 11
	for(j in 7:11){
		Sigma.dot[, , j] <- e[, , j-6] %*% t(e[, , j-6])
	}


	## Apply function for first derivative of F
	F.dot <- dF(X, Sigma, Sigma.dot)

	## Apply function for second derivative of F
	F2.dot <- d2F(X, Sigma, Sigma.dot)

	## compute step size 
	delta.theta <- solve(F2.dot) %*% F.dot

	return(delta.theta)
}







###
# Function that will iterate until sufficiently small step size
###

Scoring <- function(X, start, cutoff, max.iter){
	## bookkeeping
	theta <- start
	count <- 1
	max.step <- 1e5

	## while loop to minimize step
	while(max.step > cutoff & count < max.iter){
		max.step <- max(abs(DelTheta(X, theta)))
		theta <- theta - DelTheta(X, theta)
		
		## progress monitor
		print(paste(count, max.step))

		## update counter
		count <- count + 1
	}
	return(theta)
}




###
# set our start values and run the function
###
theta0 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

score.cfa <- Scoring(X, start = theta0, cutoff = 1e-5, max.iter = 100)




###
# clean up and print result
###
score.cfa <- round(score.cfa, digits = 2)

print(paste0('The factor loadings are: ', score.cfa[1], ', ', score.cfa[2], ', ', score.cfa[3], ', ', score.cfa[4], ', ', 'and ', score.cfa[5], '.  ', 'The unique variances are: ', score.cfa[7], ', ', score.cfa[8], ', ', score.cfa[9], ', ', score.cfa[10], ', ', 'and ', score.cfa[11],'.  ', 'The correlation between factors is ', score.cfa[6], '.'))







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

