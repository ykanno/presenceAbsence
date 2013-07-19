library(Rcpp)
library(corpcor)
library(mcmc)
library(IMIS)
library(mvtnorm)
 
sourceCpp("code.cpp", verbose=TRUE)

inv_logit <- plogis
#########################################################
# randPA generates a data frame of random presence-absence
# and a covariate
#########################################################
randPA <- function(nSites=10, lambda=1.5, sigma=1, B0 = 0.3, B1 = 7, distMat = NULL, temp=NULL) {
# 'nSites' is number of sites
# 'lambda' dampening coef for prob of occ by distance
# 'sigma' is the variance of the spatial process

# generate random (x,y) locations for sites
locs <- cbind(runif(nSites),runif(nSites))
 
# matrix of pairwise distances between sites
if(is.null(distMat)==T) {distMat <- as.matrix(dist(locs, upper=TRUE, diag=FALSE))}
 
# generate spatial var-cov matrix
vcSpace <- sigma*exp(-(lambda*distMat)^2)

if(is.null(temp)==T) {temp=mvrnorm(1, mu = rep(B0,nSites), Sigma = as.matrix(make.positive.definite(vcSpace)))}

PA = rbinom(nSite,1, prob=plogis(B1*temp))

return(list(PA=PA, distMat=distMat, vcSpace=vcSpace, temp = temp))
}


estimate = function(par) {
	fecund <- rep(exp(par[1]),nSite); # fecundity mean
	survJ <- rep(inv_logit(par[2] + par[6]*temp),nSite); # juv survival mean
	survA <- rep(inv_logit(par[3] + par[6]*temp),nSite); # adult survival mean
	migRate <- rep(inv_logit(par[4]),nSite); 
	N <- exp(par[5]);
	
	# fill in the diagonal elements
	bigMat <- bigM(nSite, migRate, fecund, survJ, survA)
    
    pop = Re(eigen(bigMat)$vectors[,1]);
    pop = pop/sum(pop)
	#stableAgeDist = eigen(bigMat)$vector[,1]
	
	# calculate site-stage probabilities
	probs = 1 - exp(-N * (pop[seq(1,nSite2,2)]+pop[seq(2,nSite2,2)]))

    NLL = sum(dbinom(y, size = 1, prob = probs, log=TRUE))	
	# combine stage probabilities, return NLL
	return(-NLL)
}

estimateIMIS = function(par) {
	fecund <- rep(exp(par[1]),nSite); # fecundity mean
	survJ <- rep(inv_logit(par[2] + par[6]*temp),nSite); # juv survival mean
	survA <- rep(inv_logit(par[3] + par[6]*temp),nSite); # adult survival mean
	migRate <- rep(inv_logit(par[4]),nSite); 
	N <- exp(par[5]);
	
	# fill in the diagonal elements
	bigMat <- bigM(nSite, migRate, fecund, survJ, survA)
    
    pop = Re(eigen(bigMat)$vectors[,1]);
    pop = pop/sum(pop)
	
	# calculate site-stage probabilities
	probs = 1 - exp(-N * (pop[seq(1,nSite2,2)]+pop[seq(2,nSite2,2)]))

    L = exp(sum(dbinom(y, size = 1, prob = probs, log=TRUE)))
	# combine stage probabilities, return NLL
	return(L)
}

likelihood <- function(theta) {
	if (is.matrix(theta)) {
		return(apply(X=theta, MARGIN=1, FUN=estimateIMIS))
	} else {
		return(estimateIMIS(theta))
	}
}

pri <- function(theta) {
	priorDF <- prod(c(
		dnorm(x=theta[1], mean=0, sd=0.4),
		dnorm(x=theta[2], mean=0, sd=0.5),
		dnorm(x=theta[3], mean=0, sd=0.5),
		dnorm(x=theta[4], mean=0, sd=1),
		dnorm(x=theta[5], mean=log(sum(y)), sd=1),
		dnorm(x=theta[6], mean=0, sd=1)
	))
	return(priorDF)
}

prior <- function(theta) {
	if (is.matrix(theta)) {
		return(apply(X=theta, MARGIN=1, FUN=pri))
	} else {
		return(pri(theta))
	}
}

sample.prior <- function(n) {
	theta <- matrix(data = c(	
		rnorm(n=n, mean=0, sd=0.4),
		rnorm(n=n, mean=0, sd=0.5),
		rnorm(n=n, mean=0, sd=0.5),
		rnorm(n=n, mean=0, sd=1),
		rnorm(n=n, mean=log(sum(y)), sd=1),
		rnorm(n=n, mean=0, sd=1)	
	), nrow=n, byrow=FALSE)
	return(theta)
}

####################################################################
# Simulation test # 1
# This simulation was done iterating over different values of spatial
# correlation and looking at how different degrees of spatial correlation
# change parameter estimates. 
####################################################################

nSite = 50
locs = cbind(runif(nSite),runif(nSite))

Ls = seq(0.1,2,0.05)# spatial correlation parameters
outputList = list()
nSite2 = nSite*2
for(i in 1:length(Ls)) {
  dat = randPA(nSites = nSite, lambda=Ls[i], sigma=0.1, distMat = as.matrix(dist(locs,upper=T,diag=T)))
  y = dat$PA
  temp = dat$temp
  outputList[[i]] = IMIS(B=500, B.re=500, number_k=100, D=10)
}

# make plot of the posteriors
q = matrix(0,length(Ls)[1],3)
for(i in 1:length(Ls)) {
	q[i,] = quantile(outputList[[i]]$resample[,6], c(0.025,0.5,0.975))
}
plot(Ls,q[,2],type="l",lwd=2,ylim=c(min(q),max(q)))
lines(Ls,q[,1])
lines(Ls,q[,3])


####################################################################
# Simulation test # 2
# This simulation was done using 20 replicate data sets of 0s and 1s, 
# with the exact same population structure & covariate temperature values. 
# The sensitivity question is whether if the 0s and 1s are resampled, 
# we can recover the coefficient values & whether migration / N estimates
# are similar
####################################################################
nSite = 50
outputList2 = list()
nSite2 = nSite*2
locs = cbind(runif(nSite),runif(nSite)) # randomize locations
# call randPA to generate the temp data
dat = randPA(nSites = nSite, lambda=1.5, sigma=0.1, distMat = as.matrix(dist(locs,upper=T,diag=T)))
temp = dat$temp

for(i in 1:20) {
  dat = randPA(nSites = nSite, lambda=1.5, sigma=0.1, distMat = as.matrix(dist(locs,upper=T,diag=T)), temp = temp)
  y = dat$PA
  o = IMIS(B=500, B.re=500, number_k=100, D=10)
  outputList2[[i]] = o
}
pdf("Simulation test 2.pdf")
# make boxplots of parameters
par(mfrow = c(3,2),mai=c(0.3,0.3,0.2,0.2))
names = c("fecund","sJuv","sAd","migRate","N","B")
for(param in 1:6) {
x = rep(1, 500)
y = outputList2[[1]]$resample[,param]
for(i in 2:20) {
	x = c(x, rep(i,500))
	y = c(y, outputList2[[i]]$resample[,param])
}
if(param %in%c(1,5)) y = exp(y)
if(param%in%c(2,3,4)) y = plogis(y)
boxplot(y~x,outline=F,main=names[param],col="grey70")
}
dev.off()
