
# load in data
load("/Users/Eric.Ward/Downloads/ONP_SoleDuc98_amphibs.RData")

# get data for amma
amma = dat[,c("AMMA","Northing","Easting","elev")]

# create distance matrix
amma.dist = dist(amma[,c("Northing","Easting")], diag=T, upper=T)

