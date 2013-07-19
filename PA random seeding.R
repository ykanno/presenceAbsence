

randPA <- function(nSites=10,lambda=1.5) {
	# 'nSites' is number of sites
	# 'lambda' dampening coef for prob of occ by distance
	
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
	
	# plot present (solid) & absent (open) sites; solid blue is initial seed site
	plot(locs[PA==1,], pch=19, xlim=c(0,1), ylim=c(0,1),
		 ylab="", xlab="", xaxt="n", yaxt="n", asp=1)
	points(locs[PA==0,], pch=1)
	points(locs[rStart,1], locs[rStart,2], pch=19, col="blue")
		
	}
