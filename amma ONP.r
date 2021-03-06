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

# load in data
load("/Users/Eric.Ward/Downloads/ONP_SoleDuc98_amphibs.RData")

# For Ben: you'll need to create a vector 'y' of 0s and 1s that alternates between stage, then across sites. For example 

# get data for amma
amma = dat[,c("AMMA","Northing","Easting","elev")]
# create distance matrix
amma.dist = dist(amma[,c("Northing","Easting")], diag=T, upper=T)
# convert AMMA data to 0s and 1s based on juvenile (2)/ adult (1)
# format is Site1-stage1, Site1-stage2, Site2-stage1, Site2-stage2
nSite=dim(amma)[1]
y = rep(0, nSite*2)
y[which(amma$AMMA==2)*2-1] = 1 # fill in juveniles
y[which(amma$AMMA==1)*2] = 1 # fill in adults

temp = log(amma$elev) # the covariate needs to be named 'temp'
# Parameters: 1=log-Fecundity; 2=logit-junvenile survival; 3=logit-adult survival intercept; 4=logit-migration rate; 5=log-N; 6=survival slope
temp = scale(temp) # ~ N(0,1)
# 1. FIT THE MAXIMUM LIKELIHOOD MODEL HERE
Fixed = c(log(0.9), qlogis(0.8), qlogis(0.99), qlogis(0.05), NA, NA) 
mlEst = nlminb(start=runif(2), objective=FixedFn, Fixed=Fixed, control=list(trace=1))     
Est = Fixed
Est[5:6] = mlEst$par
backCalculateParameters(Est)

# This is just another version of the function in "presenceAbsence functions.r" -- the difference is here, the first 4 parameters are set
calcLike = function(par) {
  pop = calcPop(c(log(0.9), qlogis(0.8), qlogis(0.8),qlogis(0.05), par[5], par[6])) # calculate s.a.d.
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

# Estimate the Bayesian model here with SIR model
nDraws = 1000
iter = 0
cumLike = 0
thresh = 4.0e-28
savedDraws = matrix(0,nDraws,2)
while(iter < nDraws) {
	# draw from priors on log(N) and B.env
	logN = runif(1, 0, log(600))
	B.env = runif(1, -5, 5)
	
	# calculate likelihood
	like = exp(calcLike(c(log(0.9), qlogis(0.8), qlogis(0.8),qlogis(0.05), logN, B.env)))
	if(is.na(like)==F) {
	  cumLike = cumLike + like
	  if(like> best) { # this if statement makes the algorithm reset if a better value is found
 	  	iter = 0
	  	cumLike = 0
	  	savedDraws = savedDraws*0
	  	best = like
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
hist(savedDraws[,1], 40, xlab = "Log (N)", ylab="Posterior", col="grey",main="AMMA ONP")
hist(savedDraws[,2], 40, xlab = "B.elevation", ylab="Posterior", col="grey",main="")
