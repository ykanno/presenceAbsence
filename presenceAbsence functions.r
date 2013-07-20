
#########################################################
# randPA generates a data frame of random presence-absence
# and a covariate
#########################################################
randPA <- function(nSites=10, lambda=1.5, sigma=1, B1=2, locs = NULL, distMat = NULL, temp=NULL, fecund=0.9, survJ = 0.15, survA = 0.4, migRate = 0.1, N = 500) {
# 'nSites' is number of sites
# 'lambda' dampening coef for prob of occ by distance
# 'sigma' is the variance of the spatial process

# generate random (x,y) locations for sites
if(is.null(locs)==T) locs <- cbind(runif(nSites),runif(nSites))
 
# matrix of pairwise distances between sites
if(is.null(distMat)==T) {distMat <- as.matrix(dist(locs, upper=TRUE, diag=FALSE))}
 
# generate spatial var-cov matrix
vcSpace <- sigma*exp(-(lambda*distMat)^2)

# this line generates the gaussian random field, and is temporarily commented out
#if(is.null(temp)==T) {temp=mvrnorm(1, mu = rep(B0,nSites), Sigma = as.matrix(make.positive.definite(vcSpace)))}

if(is.null(temp)==T) temp = seq(0,1,length.out=nSites) # make temperature a uniform(0,1) RV
betaJ = log(survJ/(1-survJ)) # back transform from user-specified intercepts
betaA = log(survA/(1-survA))

logitJ = 0*temp*B1 + betaJ
survJ = plogis(logitJ); # turn survJ and survA into site-specific vectors of (survival)
logitA = temp*B1 + betaA
survA = plogis(logitA);
migRate = rep(migRate, nSites) # migrationRate and fecundity constant, but need to be vectors for the C call
fecund = rep(fecund,nSites)

# fill in the elements of the large transition matrix
bigMat <- bigM(nSites, migRate, fecund, survJ, survA)
    
pop = Re(eigen(bigMat)$vectors[,1]); # calculate dominant eigenvector, and rescale
pop = pop/sum(pop)
	
# calculate site-stage probabilities (Poisson). The order is [Stage1-Site1, Stage2-Site1, Stage1-Site2, Stage2-Site2, ...]
probs = 1 - exp(-N * (pop))

PA = rbinom(length(probs),1, prob=probs) # generate site-stage specific data

return(list(PA=PA, distMat=distMat, vcSpace=vcSpace, temp = temp, probs = probs, bigMat = bigMat, survA = survA, survJ = survJ, logitJ = logitJ, logitA = logitA))
}

estimateML = function(par) {
pop = calcPop(par) # calculate s.a.d.
N = exp(par[5]) # this is local variable
# calculate site-stage probabilities
probs = ifelse(observedStages == 1, 1 - exp(-N * (pop[seq(1,nSite2,2)]+pop[seq(2,nSite2,2)])), 1 - exp(-N * pop)) 
L = -(sum(dbinom(y, size = 1, prob = probs, log=TRUE)))
# combine stage probabilities, return neg log likelihood
return(L)
}

backCalculateParameters = function(par) {
	fecund <- rep(exp(par[1]),nSite); # fecundity mean
    survJ <- rep(plogis(par[2] + par[6]*temp),nSite); # juv survival mean
    survA <- rep(plogis(par[3] + par[6]*temp),nSite); # adult survival mean
    migRate <- rep(plogis(par[4]),nSite); 
    N <- exp(par[5]);
    return(list("fecund"=fecund,"survJ"=survJ,"survA"=survA,"migRate"=migRate,"N"=N))
}

# this function is used to calculate s.a.d
calcPop = function(par) {
transPars = backCalculateParameters(par)
fecund <- transPars$fecund; # fecundity mean
survJ <- transPars$survJ; # juv survival mean
survA <- transPars$survA; # adult survival mean
migRate <- transPars$migRate; 
N <- transPars$N;
	
# fill in the diagonal elements
bigMat <- bigM(nSite, migRate, fecund, survJ, survA)
    
pop = Re(eigen(bigMat)$vectors[,1]);
pop = pop/sum(pop)
return(pop)		
}

# this function is called by IMIS to generate posterior
estimateIMIS = function(par) {
pop = calcPop(par) # calculate s.a.d.
N = exp(par[5]) # this is local variable
# calculate site-stage probabilities
probs = ifelse(observedStages == 1, 1 - exp(-N * (pop[seq(1,nSite2,2)]+pop[seq(2,nSite2,2)])), 1 - exp(-N * pop)) 
L = exp(sum(dbinom(y, size = 1, prob = probs, log=TRUE)))
# combine stage probabilities, return likelihood in normal space
return(L)
}

# this function is called to calculate likelihood
likelihood <- function(theta) {
	if (is.matrix(theta)) {
		return(apply(X=theta, MARGIN=1, FUN=estimateIMIS))
	} else {
		return(estimateIMIS(theta))
	}
}

# this function is called to calculate prior density
pri <- function(theta) {
	priorDF <- prod(c(
		dnorm(x=theta[1], mean=0, sd=1),
		dnorm(x=theta[2], mean=0, sd=1),
		dnorm(x=theta[3], mean=0, sd=1),
		dnorm(x=theta[4], mean=0, sd=1),
		dunif(x=theta[5], log(sum(y)), log(1000)),
		dnorm(x=theta[6], mean=0, sd=1)
	))
	return(priorDF)
}

# this function returns a matrix of prior values
prior <- function(theta) {
	if (is.matrix(theta)) {
		return(apply(X=theta, MARGIN=1, FUN=pri))
	} else {
		return(pri(theta))
	}
}

# function to sample from priors
sample.prior <- function(n) {
	theta <- matrix(data = c(	
		rnorm(n=n, mean=0, sd=1),
		rnorm(n=n, mean=0, sd=1),
		rnorm(n=n, mean=0, sd=1),
		rnorm(n=n, mean=0, sd=1),
		runif(n=n, log(sum(y)), log(1000)),
		rnorm(n=n, mean=0, sd=1)	
	), nrow=n, byrow=FALSE)
	return(theta)
}

