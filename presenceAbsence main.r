# load libraries
library(Rcpp)
library(corpcor)
library(mcmc)
library(IMIS)
library(mvtnorm)

# set wd
setwd("/users/eric.ward/documents/projects/presenceAbsence")
#setwd("/users/James.Thorson/Desktop/Project_Git/presenceAbsence")

# source local functions and C++ code
source("presenceAbsence functions.r")
sourceCpp("code.cpp", verbose=TRUE)

observedStages = 2 # juvenile and adult PA data observed, otherwise = 1

####################################################################
# demonstrate how to simulate data
####################################################################
nSite=50
beta1 = 0.8
dat = randPA(nSites=nSite, lambda=1.5, sigma=1, B1=beta1, locs = NULL, distMat = NULL, temp=NULL, fecund=0.9, survJ = 0.1, survA = 0.2, migRate = 0.05, N = 100)
temp=dat$temp
pdf("Survival transformations.pdf")
# test the presence - absence function
par(mfrow=c(3,2), mai = c(0.6,0.6,0.1,0.1))
plot(seq(1,length(dat$temp)),dat$temp, xlab ="Temp", ylab="Uniform(0,1)",type="b",lwd=2)
legend('topleft',("1"))
plot(dat$temp, dat$logitA, xlab = "Temp", ylab = "Survival - logit space",type="b",lwd=2)
legend('topleft',("2"))
plot(dat$temp, dat$survA, xlab = "Temp", ylab = "Survival - normal space",type="b",lwd=2,ylim=c(0,1))
legend('topleft',("3"))
plot(dat$temp, dat$probs[seq(2,length(dat$PA),2)], xlab = "Temp", ylab = "Survival - postEigen")
legend('topleft',("4"))
plot(dat$temp,dat$PA[seq(2,length(dat$PA),2)], xlab="Temp",ylab="Data")
legend('topleft',("5"))
dev.off()

####################################################################
# demonstrate how to find MLEs with optim
####################################################################
nSite = 50

dat = randPA(nSites=50, lambda=1.5, sigma=1, B1=5, locs = NULL, distMat = NULL, temp=NULL, fecund=0.9, survJ = 0.1, survA = 0.2, migRate = 0.05, N = 100)
y = dat$PA
temp = dat$temp
mlEst = optim(runif(6), estimateML)
backCalculateParameters(mlEst$par)

####################################################################
# demonstrate how to find posteriors with ISIS
####################################################################
nSite = 50
dat = randPA(nSites=50, lambda=1.5, sigma=1, B1=5, locs = NULL, distMat = NULL, temp=NULL, fecund=0.9, survJ = 0.1, survA = 0.2, migRate = 0.05, N = 100)
y = dat$PA
temp = dat$temp
o = IMIS(B=500, B.re=500, number_k=100, D=10)

####################################################################
# Simulation test # 2
# This simulation was done using 20 replicate data sets of 0s and 1s, 
# with the exact same population structure & covariate temperature values. 
# The sensitivity question is whether if the 0s and 1s are resampled, 
# we can recover the coefficient values & whether migration / N estimates
# are similar
####################################################################
nSite = 100
outputList2 = list()
outputList = list()
nSite2 = nSite*2
locs = cbind(runif(nSite),runif(nSite)) # randomize locations
# call randPA to generate the temp data
dat = randPA(nSites = nSite, locs = locs, lambda=1.5, sigma=1, distMat = as.matrix(dist(locs,upper=T,diag=T)), N = 1000)
temp = dat$temp

for(i in 1:20) {
  dat = randPA(nSites = nSite, locs = locs, lambda=1.5, sigma=0.1, distMat = as.matrix(dist(locs,upper=T,diag=T)), temp = temp,N = 1000)
  y = dat$PA
  o = IMIS(B=500, B.re=500, number_k=100, D=10)
  outputList2[[i]] = o
  o = glm(y ~ c(temp), family = "binomial")
  outputList[[i]] = o  
}

pdf("Simulation test 2.pdf")
# make boxplots of parameters
par(mfrow = c(3,2),mai=c(0.3,0.3,0.2,0.2))
names = c("fecund","sJuv","sAd","migRate","N","B")
for(param in 1:6) {
x = rep(1, 500)
y = outputList2[[1]]$resample[,param]
for(i in 2:20) {
	x = c(x, rep(i,500))~/Desktop/Simulation test 2.pdf
	y = c(y, outputList2[[i]]$resample[,param])
}
if(param %in%c(1,5)) y = exp(y)
if(param%in%c(2,3,4)) y = plogis(y)
boxplot(y~x,outline=F,main=names[param],col="grey70")
}
dev.off()

m = matrix(0, 500, 100)
for(i in 1:500) {m[i,] = 1-exp(-2*exp(o$resample[i,6])*calcPop(o$resample[i,]))}

estimate(apply(o$resample,2,mean))
o2 = o$resample
o2[,5] = log(2) + (o2[,5])
estimate(apply(o2,2,mean))
