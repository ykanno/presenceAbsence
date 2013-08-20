# load libraries
library(Rcpp)
library(corpcor)
library(mcmc)
library(IMIS)
library(mvtnorm)
library(pso)

# set local wd
setwd("/users/eric.ward/documents/projects/presenceAbsence")
#setwd("/users/James.Thorson/Desktop/Project_Git/presenceAbsence")

# source local functions and C++ code
source("presenceAbsence functions.r")
sourceCpp("code.cpp", verbose=TRUE)

observedStages = 2 # juvenile and adult PA data observed, otherwise = 1

# For Ben: you'll need to create a vector 'y' of 0s and 1s that alternates between stage, then across sites. And a vector of covariates named 'temp'
nSite = 100
y = c(sample(c(0,1), size=200, replace=T))
temp = rnorm(nSite,0,1) # the covariate needs to be named 'temp'

# Fix the demographic rates 
fecundity = 0.9
juvSurvival = 0.8
adultSurvival = 0.99
migrationRate = 0.05

# Parameters: 1=log-Fecundity; 2=logit-junvenile survival; 3=logit-adult survival intercept; 4=logit-migration rate; 5=log-N; 6=survival slope
temp = scale(temp) # ~ N(0,1)

# UNCOMMENT TO FIT THE MAXIMUM LIKELIHOOD MODEL HERE
#Fixed = c(log(fecundity), qlogis(juvSurvival), qlogis(adultSurvival), qlogis(migrationRate), NA, NA) 
#mlEst = nlminb(start=runif(2), objective=FixedFn, Fixed=Fixed, control=list(trace=1))     
#Est = Fixed
#Est[5:6] = mlEst$par
#backCalculateParameters(Est)

# This is just another version of the function in "presenceAbsence functions.r" -- the difference is here, the first 4 parameters are set
calcLike = function(par) {
  pop = calcPop(c(log(fecundity), qlogis(juvSurvival), qlogis(adultSurvival), qlogis(migrationRate), par[5], par[6])) # calculate s.a.d.
  N = exp(par[5]) # this is local variable
  # calculate site-stage probabilities
  nSite2 = nSite*2
  if(observedStages == 1) popTrans = pop[seq(1,nSite2,2)]+pop[seq(2,nSite2,2)]
  if(observedStages == 2) popTrans = pop
  probs =1 - exp(-N * popTrans)
  #probs = ifelse(observedStages == 1, 1 - exp(-N * (pop[seq(1,nSite2,2)]+pop[seq(2,nSite2,2)])), 1 - exp(-N * pop)) 
  L = (sum(dbinom(y, size = 1, prob = probs, log=TRUE)))
  # combine stage probabilities, return neg log likelihood
  return(L)
}

# Estimate the Bayesian model here with SIR model. You should be able to let it run as is -- but you might want to adjust the maximum possible population size 'maxN' below -- for the ONP data, we used log(600). 
nDraws = 1000
iter = 0
cumLike = 0
thresh = 0
maxN = 600
savedDraws = matrix(0,nDraws,2)
while(iter < nDraws) {
	# draw from priors on log(N) and B.env
	logN = runif(1, 0, log(maxN))
	B.env = runif(1, -5, 5)
	
	# calculate likelihood
	like = exp(calcLike(c(log(fecundity), qlogis(juvSurvival), qlogis(adultSurvival), qlogis(migrationRate), logN, B.env)))
	if(is.na(like)==F) {
	  cumLike = cumLike + like
	  if(like> thresh) { # this if statement makes the algorithm reset if a better value is found
 	  	iter = 0
	  	cumLike = 0
	  	savedDraws = savedDraws*0
	  	thresh = like
	  }
	  if(cumLike > thresh) {
		while(cumLike > thresh) {
		  # save this draw
		  iter = iter+1
		  savedDraws[iter,] = c(logN,B.env)
		  cumLike = cumLike - thresh
		  print(paste("Iterations: ",iter,sep=""))
		}
	  }
	}
	#print(paste(cumLike,":",thresh))
}

par(mfrow=c(2,1), mai = c(1,1,0.2,0.2))
hist(savedDraws[,1], 40, xlab = "Log (N)", ylab="Posterior", col="grey",main="")
hist(savedDraws[,2], 40, xlab = "Environmental coefficient", ylab="Posterior", col="grey",main="")
