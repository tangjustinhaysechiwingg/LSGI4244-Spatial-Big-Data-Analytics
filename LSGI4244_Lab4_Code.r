LSGI4244 Lab 4 R Programming code

install.packages("spdep",repos = "http://cran.r-project.org")
install.packages('RSEIS') # matrix manipulation
install.packages('RColorBrewer')

library(spdep)
library(rgdal) 
library('RSEIS')
library(RColorBrewer)


setwd("C:/Users/20016345d/Documents/R/Lab4/Lab4_Area_pattern_analysis/data")
getwd()

ColData = readOGR(".", "columbus") # columbus.shp is in folder: data
coords = coordinates(ColData) # coordinates of the census tract centroids
IDs = row.names(as(ColData, "data.frame"))

col_nb = poly2nb(ColData, queen = TRUE)
summary(col_nb)
class(col_nb) # check the class of col_nb
help(poly2nb) # check help to see detailed information about poly2nb function

nb = unlist(col_nb)
class(nb) # check the class of nb
nb
nbc1 = table(nb) #table shows the number of neighbors for each polygon
nbc2 = table(nbc1) #table shows the distribution of neighbor connections
 #i.e., the number of cases that have certain number of neighbors
nbc1
nbc2
barplot(nbc2,col=rainbow(20), xlab="number of neighbors",ylab="number
of cases")

plot(ColData, border = "grey60")
plot(col_nb, coords, pch = 19, cex = 0.6, add = TRUE)
title("LSGI4244 Lab 4 - Adjacency Connection")
box()

nbk1 = knn2nb(knearneigh(coords, k=1), row.names=IDs)
nbk2 = knn2nb(knearneigh(coords, k=2), row.names=IDs)

plot(ColData, border="grey60")
plot(nbk1, coords, add=TRUE, pch=19, cex=0.6)
text(bbox(ColData)[1,1] + 0.5, bbox(ColData)[2,2], labels="k=1")
box()

plot(ColData, border="grey60")
plot(nbk2, coords, add=TRUE, pch=19, cex=0.6,)
text(bbox(ColData)[1,1] + 0.5, bbox(ColData)[2,2], labels="k=2")
box()


dsts = unlist(nbdists(col_nb, coords))
summary(dsts)
hist(dsts)

dist1=2000
dist2=3000
nb1 = dnearneigh(coords, d1=0, d2=dist1, row.names=IDs)
nb2 = dnearneigh(coords, d1=0, d2=dist2, row.names=IDs)


dev.off()
dev.new(width=8, height=4)
par(mfrow = c(1, 2), mar=c(1,1,1,1), pty = "s")

# Band width=2000 m
plot(ColData, border="grey60")
plot(nb1, coords, add=TRUE, pch=19, cex=0.5)
text(bbox(ColData)[1,1] + 3000, bbox(ColData)[2,2], labels="0.5")
box()



plot(ColData, border="grey60")
plot(nb2, coords, add=TRUE, pch=19, cex=0.8)
text(bbox(ColData)[1,1] + 3000, bbox(ColData)[2,2], labels="0.8")
box()







# Assignment 1 of Lab 4


col.W = nb2listw(col_nb, style="W") # weight list
crime = ColData$CRIME

#create a thematic map of crime to check the spatial pattern visually
 source('C:/Users/20016345d/Documents/R/Lab4/Lab4_Area_pattern_analysis/code/themeMaps3.r')
 n = 7 # the number of classes in the thematic map
 wh = c(5,4) # plot window size
 mrgn=c(0.5,0.5,2,0.5) # set plot margin
 lgdloc = 'topleft' # location of the legend
 windows(wh[1],wh[2])
 par(mfrow=c(1,1), mar=mrgn)
 themeMaps3(ColData, crime, n, lgdloc)
 title(main = 'LSGI4244 Lab 4 - Thematic Map of Crime in Columbus, 1980')
 box()


windows(10,10)
 par(mar=c(4,4,1,0.1)) # set plot margin
 moran.plot(crime, col.W, pch=19) # Moran scatterplot


moran.test(crime, col.W,zero.policy=T)
 windows(6,5)
 par(mar=c(4,4,1,0.1)) # set plot margin


# plot the permutation test results
 windows(4,4)
 par(mar=c(4,4,1,0.1)) # set plot margin
 graph999 = hist(mI_perm999$res, breaks=100,col="light blue",
xlab="Moran's I",main="Permutation Test for Moran's I - 999
permutations")


# add the Moran's I of the data
 lines(mI_perm999$statistic,max(graph999$counts),
type="h",col="red",lwd=2)



# Assignment 2 of Lab 4


mI_local = localmoran(crime, col.W)

quadrant = vector(mode="numeric",length=nrow(mI_local))
 cCrime = crime - mean(crime) # centers the variable around its mean
 C_mI = mI_local[,1] - mean(mI_local[,1]) # centers the local Moran's around the
mean
 signif = 0.05 # set a statistical significance level
quadrant[cCrime >0 & C_mI>0] = 1 # high-high cluster
quadrant[cCrime <0 & C_mI>0] = 2 # low-low cluster
quadrant[cCrime >0 & C_mI<0] = 3 # high-low cluster
quadrant[cCrime <0 & C_mI<0] = 4 # low-high cluster
quadrant[mI_local[,5]>signif] = 0 # non-significant Moran's in the category '0'
cols.c = c('ghostwhite','red','blue','pink','lightblue')
clr.c = cols.c[quadrant+1]

windows(5,4)
par(mar=c(0.5,0.5,2,0.5))
plot(ColData, col=clr.c)
leglabs = c("Not Significant","High-High","Low-Low","High-Low","LowHigh")
legend('topleft', fill = cols.c,legend = leglabs,bty = 'n', cex = 1.1)
title(main = 'LISA Cluster Map (CRIME)')
box()

dev.copy(png,"C:/Users/20016345d/Pictures/Saved Pictures/lisaCluster.png",width=5,height=4,units="in",res=500
)
dev.off()







# Assignment 3 of Lab 4

Col.lm=lm(CRIME ~ INC + HOVAL, data=ColData)
summary(Col.lm)
ols.resd=resid(Col.lm) # get residuals
ols.res=(ols.resd-min(ols.resd))/diff(range(ols.resd)) #normalized the residuals to 0-1 to better show the difference among area units

source('C:/Users/20016345d/Documents/R/Lab4/Lab4_Area_pattern_analysis/code/themeMaps3.r') #call an external function
n = 7 # the number of classes in the thematic map
wh = c(5,4) # plot window size
mrgn=c(0.5,0.5,2,0.5) # set plot margin
lgdloc = 'topleft' # location of the legend
windows(wh[1],wh[2])
par(mfrow=c(1,1), mar=mrgn)
themeMaps3(ColData,ols.res,n,lgdloc)
title('Residuals of OLS Model')
box()
dev.copy(png,"C:/Users/20016345d/Pictures/Saved Pictures/ols.residuals.map.png",width=8,height=8,
units="in",res=500)
dev.off()

 ols.resd=resid(Col.lm) # get residuals
 windows(8,8)
 par(mfrow=c(2,2), mar=c(4,4,2,0.5))
 plot(ColData$HOVAL,ols.resd,ylab="Residuals",xlab="HOVAL")
 abline(h=0)
 title('OLS Model: Residuals vs. HOVAL')
 plot(ColData$INC,ols.resd,ylab="Residuals",xlab="INC")
 abline(h=0)
 title('OLS Model: Residuals vs. INC')
 hist(ols.resd, xlab="Residuals")
 qqnorm(ols.resd, ylab="Residuals")
 qqline(ols.resd)
 dev.copy(png,"C:/Users/20016345d/Pictures/Saved Pictures/ols.png",width=8,height=8,units="in",res=500)
 dev.off()

 # Moran plot
 windows(5,5)
 par(mar=c(4,4,0.1,0.1)) # set plot margin
 # Moran's I scatterplot: points with large influence is labelled.
 moran.plot(ols.resd, col.W,pch=19)

 dev.copy(png,"C:/Users/20016345d/Pictures/Saved Pictures/moran.ols.resd.png",width=5,height=5,
 units="in",res=500)
 dev.off ()

lm.morantest(Col.lm,col.W)
 ols.mI_perm999 = moran.mc(ols.resd, col.W, 999)
 windows(5,5)
 par(mar=c(4,4,1,0.1)) # set plot margin
 ols.graph999=hist(ols.mI_perm999$res, breaks=100,col="lightblue", xlab="Moran's I", main="Morans'I permutation test - 999 permutations")
 lines(ols.mI_perm999$statistic,max(ols.graph999$counts),
type="h",col="red",lwd=2)
 dev.copy(png,"C:/Users/20016345d/Pictures/Saved Pictures/ols_mI_Permutation_2.png",width=8,height=8,units="in",res=500)

