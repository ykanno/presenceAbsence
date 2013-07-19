library(Rcpp)
library(IMIS)
library(corpcor)
sourceCpp("code.cpp", verbose=TRUE)
   
nSite = 100
locs = cbind(runif(nSite),runif(nSite))
distMat = as.matrix(dist(locs, upper=TRUE, diag=TRUE))
  
nSite2 = nSite*2
nProjections = 1

inv_logit <- plogis

randPA <- function(nSites=10, lambda=1.5, sigma=1, B0 = 0.3, B1 = 30) {
# 'nSites' is number of sites
# 'lambda' dampening coef for prob of occ by distance
# 'sigma' is the variance of the spatial process

# generate random (x,y) locations for sites
locs <- cbind(runif(nSites),runif(nSites))
 
# matrix of pairwise distances between sites
distMat <- as.matrix(dist(locs, upper=TRUE, diag=FALSE))
 
# pick random starting site
rStart <- sample(nSites)[1]

# convert distances to probs of occurences
probs <- exp(-lambda*distMat[rStart,])

# generate P/A for each site
PA <- rbinom(nSites,1,prob=probs)

# generate spatial var-cov matrix
vcSpace <- sigma*exp(-(lambda*distMat)^2)

    temp=rmvnorm(1, mean = rep(B0,nSites), make.positive.definite(vcSpace))
    PA = rbinom(nSite,1, prob=plogis(B1*rmvnorm(1, mean = rep(B0,nSite), make.positive.definite(vcSpace))))

return(list(PA=PA, distMat=distMat, vcSpace=vcSpace, temp = temp))
}

# generate real data that is spatially explicit 
  dat = randPA(nSites = nSite, lambda=1.5, sigma=0.1)
  y = dat$PA
  temp = dat$temp

estimate = function(par) {
	B <- par[6];
	fecund <- rep(exp(par[1]),nSite); # fecundity mean
	survJ <- rep(inv_logit(par[2] + B*temp),nSite); # juv survival mean
	survA <- rep(inv_logit(par[3] + B*temp),nSite); # adult survival mean
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
	return(NLL)
	
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

o = IMIS(B=500, B.re=5, number_k=100, D=10)




plot(exp(sample.prior(100)[,1]), main='fecund')
points(exp(o$resample[,1]), pch = '*', col='red')

plot(plogis(sample.prior(100)[,2]), main='b2')
points(plogis(o$resample[,2]), pch = '*', col='red')

plot(plogis(sample.prior(100)[,3]), main='s1')
points(plogis(o$resample[,3]), pch = '*', col='red')

plot(plogis(sample.prior(100)[,4]), main='s2')
points(plogis(o$resample[,4]), pch = '*', col='red')

plot(exp(sample.prior(100)[,5]), main='m')
points(exp(o$resample[,5]), pch = '*', col='red')

plot((sample.prior(100)[,6]), main='m')
points((o$resample[,6]), pch = '*', col='red')



