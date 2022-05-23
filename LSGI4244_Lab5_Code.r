# LSGI4244 Lab 5 - Geostatistical Modelling
# R Programming Code
# Tang Justin Hayse Chi Wing G. - 20016345D 

# Assignment 1
install.packages("geoR",repos = "http://cran.r-project.org")
install.packages("maptools",repos = "http://cran.r-project.org")
install.packages("maps",repos = "http://cran.r-project.org")
install.packages("mvtnorm",repos = "http://cran.r-project.org")
library(geoR)
library(maps)
library(mvtnorm) # package for simulation
rm(list=ls()) # clear workspace

oh.range=map("state","ohio",plot=FALSE)$range
rx=oh.range[2]-oh.range[1] # range in x direction
ry=oh.range[4]-oh.range[3] # range in y direction

xg=seq(oh.range[1],oh.range[2], length.out=30) 
yg=seq(oh.range[3],oh.range[4], length.out=30)
ohgrid.locs=expand.grid(xg,yg)
map("state","ohio")
points(ohgrid.locs,pch=20)
ntot=dim(ohgrid.locs)[1] 

distmat=as.matrix(dist(ohgrid.locs))
max_dist=max(distmat)

theta1=2 # we are making these values in order to specify a true covariance matrix.
theta2=3
Sig=theta1*exp(-distmat/theta2) 

h = seq(0,max_dist,length.out=20) # lags where to estimate the variograms
plot( h, theta1*exp(-h/theta2), type="b", ylab="cov",
xlab="distance" )

y=rmvnorm(1,matrix(5,ntot,1),Sig) # takes some time
image(matrix(y,length(xg),length(yg)),x=xg,y=yg,
 col= rev(rainbow(100,start=0,end=.7)))
map("state","ohio",lwd=3,add=TRUE)
box() # add a boundary box to the plot

idxkeep=sort(sample(1:ntot,round(0.3*ntot)))
ymask=matrix(0,ntot,1)
ymask[idxkeep,]=1
image(matrix(y,length(xg),length(yg)),x=xg,y=yg,
 col=rev(rainbow(100,start=0,end=.7)))
image(matrix(ymask,length(xg),length(yg)),x=xg,y=yg,
 col=c("white","transparent"),add=TRUE)
map("state","ohio",add=TRUE,lwd=3)
box()

n=length(idxkeep)
ydat=matrix(y[idxkeep],n,1)
dat.locs=ohgrid.locs[idxkeep,]
ygeodat=as.geodata(cbind(dat.locs,ydat)) 

distsample=as.matrix(dist(dat.locs)) # distance between sample pairs
maxd = max(distsample)/2
y.v=variog(ygeodat,max.dist=maxd) # estimate empirical semi-variogram
y.v$beta.ols # estimated trend of the sample data

vt = theta1-theta1*exp(-y.v$u/theta2)

xmax = max(y.v$u,maxd) # xmax to set xlim for plot
ymax = max(y.v$v,vt) # ymax to set ylim for plot

plot(y.v, xlim = c(0,xmax), ylim = c(0,ymax) ) # empricial semivariogram
lines(y.v$u, vt,type="b",pch=20) # add true semivariogram


plot(variog4(ygeodat,max.dist=maxd)) 


y.v.wls=variofit(y.v, cov.model="exponential", wei="cressie")
y.v.wls
plot(y.v,pch=3, xlim=c(0,xmax), ylim=c(0,ymax)) # empricial semivariogram
lines(y.v.wls,col=4) # fitted semivariogram
lines(y.v$u,vt,lty=2, col=2) # add true semivariogram


pred.locs=ohgrid.locs[-idxkeep,] # predict un-sampled locations only will save time
y.krig.ok=krige.conv(ygeodat,loc=pred.locs,
krige=krige.control(type.krige="ok",obj.m=y.v.wls))
ypred.ok=matrix(0,ntot,1)
ypred.ok[idxkeep,]=ydat
ypred.ok[-idxkeep,]=y.krig.ok$predict
ypredvar.ok=matrix(0,ntot,1)
ypredvar.ok[-idxkeep,]=y.krig.ok$krige.var 

windows(6,6)
par(mar=c(4,3,2,0.5)) # set plot margin
par(mfrow=c(2,2), pty="s")

image(matrix(y,length(xg),length(yg)),
 x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),
 main="Full Dataset")
map("state","ohio",lwd=2,add=TRUE)
box()


image(matrix(y,length(xg),length(yg)),x=xg,y=yg,
 col=rev(rainbow(100,start=0,end=.7)),main="Sample Data") 

image(matrix(ymask,length(xg),length(yg)),
 x=xg,y=yg,col=c("white","transparent"),add=TRUE) 
map("state","ohio",lwd=2,add=TRUE)
box() 

image(matrix(ypred.ok,length(xg),length(yg)),
 x=xg,y=yg, col=rev(rainbow(100,start=0,end=.7)),
 main="Prediction(OK)")
map("state","ohio",lwd=2,add=TRUE)
box()
image(matrix(sqrt(ypredvar.ok),length(xg),length(yg)),
 x=xg,y=yg,col=rev(rainbow(100,start=0, end=.7)),
 main="Prediction Standard Errors (OK)")
image(matrix(1-ymask,length(xg),length(yg)),
 x=xg,y=yg,col=c("white","transparent"), add=TRUE)
map("state","ohio",lwd=2,add=TRUE)
box()

par(mfrow = c(1,1)) 
plot(ypred.ok,y, xlab="Predicted Data", ylab="True Data",main="Scatterplot of Prediction vs True Data")

# Assignment 2 
library(geoR)
library(maptools)
rm(list=ls())
data(meuse) 

coords = cbind(meuse$x,meuse$y) # read coordinates
coordinates(meuse) = coords
bubble(meuse, "copper") 
data(meuse.grid)

Cd.df = SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")],
 data = meuse.grid)
coords.grid = cbind(meuse.grid$x,meuse.grid$y) # read coordinates
coordinates(meuse.grid) = coords.grid
plot(meuse.grid) # grid locations to predict


Cdgeodat=as.geodata(cbind(meuse$x,meuse$y,meuse$copper))
distmat=as.matrix( dist(cbind(meuse$x,meuse$y)) )
maxd = max(distmat)/2
maxd
y.v=variog(Cdgeodat,max.dist= maxd)
y.v$beta.ols

xmax = max(y.v$u,maxd)
ymax = max(y.v$v)
plot(y.v,xlim = c(0,xmax), ylim = c(0,ymax) )

y.v.wls=variofit(y.v, cov.model="exponential", wei="cressie")
y.v.wls
plot(y.v,xlim = c(0,xmax), ylim = c(0,ymax))
lines(y.v.wls,col=2)

y.krig.ok=krige.conv(Cdgeodat,loc=coords.grid,
krige=krige.control(type.krige="ok",obj.m=y.v.wls))
ypred.ok = y.krig.ok$pred
yVar.ok = y.krig.ok$krige.var
Cd.df$OK_pred = ypred.ok
Cd.df$OK_se = sqrt(yVar.ok)

(maxCd = ceiling(max(ypred.ok))) # maxCd = 13
(minCd = floor(min(ypred.ok))) # minCd = 0
maxCd - minCd # range = 13

bluepal = colorRampPalette(c("azure3", "blue"))
brks = seq(minCd,maxCd,2)
cols = bluepal(length(brks) - 1)
maxSe = ceiling(max(Cd.df$OK_se)) # maxSe = 4
minSe = floor(min(Cd.df$OK_se)) # minSe = 2
maxSe-minSe # range = 2
brks.se = seq(minSe,maxSe,0.3)
cols.se = bluepal(length(brks.se) - 1)

dev.new(width=8, height=8) # create new plot space with specific size
par(mfrow = c(1, 1), mar=c(1,1,2,1), pty = "s")
image(Cd.df, "OK_pred", col = cols)
symbols(coords,circles=meuse$copper*5,fg="red",inches=F,add=T)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)
title(main = "Copper (Cu) - Ordinary Kriging")
box()

image(Cd.df, "OK_se", col = cols.se)
symbols(coords,circles=meuse$copper*0+20,fg="red",inches=F,add=T)
legend("topleft", fill=cols.se,
legend=leglabs(brks.se),bty="n",cex=0.8)
title(main = "Copper (Cu) - Kriging Error(OK)")
box()
