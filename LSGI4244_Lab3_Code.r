# LSGI4244 Lab 3 Code

install.packages("spatstat",repos = "http://cran.r-project.org")

data(japanesepines) # load Japanese pine data (Diggle, p.1)
jp = japanesepines # save it with another name
class(jp)

summary(jp)
plot(jp,axes=T, main=" 65 Japanese black pine saplings") # plot


Assignment 1 


jp.Z.05 = density.ppp(jp, 0.05)
par(mar=c(0,0,1,1)) # set plot margin
tl.05 = expression(paste("Kernel Estimation of JP: ", sigma, " = 0.05"))
plot(jp.Z.05, main = tl.05)
plot(jp.Z.05, main="Kernel Estimation of JP: sigma = 0.05")
points(jp,pch="+",col="6")

jp.Z.1 = density.ppp(jp, 0.1)
par(mar=c(0,0,1,1)) # set plot margin
tl.1 = expression(paste("Kernel Estimation of JP: ", sigma, " = 0.1"))
plot(jp.Z.1, main=tl.1)
points(jp,pch="+",col="6")


Assignment 2 

G Function:

jp.ghat = Gest(jp)
g.max = max(jp.ghat$r)
plot(jp.ghat,cbind(rs,theo)~r, main="Ghat for Japanese Pine Sapling", 
xlim=c(0,g.max), xlab="r", ylab="G(r)")



F Function:

jp.fhat = Fest(jp)
f.max = max(jp.fhat$r)
plot(jp.fhat, cbind(rs,theo)~r, main="Fhat for Japanese Pine Sapling", 
xlim=c(0,f.max), xlab="r", ylab="F(r)")



K Function:

jp.khat = Kest(jp)
plot(jp.khat, cbind(border, theo)~r, main="Khat for Japanese Pine Sapling")



L Function:

plot(jp.khat, sqrt(border/pi)-r ~ r, ylab="L(r)", 
 main="Lhat for Japanese Pine Sapling",ylim=c(-0.025,0.025))
abline(h=0,lty=2,col='red')

or

plot(jp.khat, sqrt(cbind(border,theo)/pi)-r ~ r, ylab="L(r)", 
 main="Lhat for Japanese Pine Sapling",ylim=c(-0.025,0.025))



Without speciying xlim argument:


G Function:

jp.ghat = Gest(jp)
g.max = max(jp.ghat$r)
plot(jp.ghat,cbind(rs,theo)~r, main="Ghat for Japanese Pine Sapling", 
xlab="r", ylab="G(r)")



F Function:

jp.fhat = Fest(jp)
f.max = max(jp.fhat$r)
plot(jp.fhat, cbind(rs,theo)~r, main="Fhat for Japanese Pine Sapling", 
xlab="r", ylab="F(r)")



K Function:

jp.khat = Kest(jp)
plot(jp.khat, cbind(border, theo)~r, main="Khat for Japanese Pine Sapling")



L Function:

plot(jp.khat, sqrt(border/pi)-r ~ r, ylab="L(r)", 
 main="Lhat for Japanese Pine Sapling",ylim=c(-0.025,0.025))
abline(h=0,lty=2,col='red')

or

plot(jp.khat, sqrt(cbind(border,theo)/pi)-r ~ r, ylab="L(r)", 
 main="Lhat for Japanese Pine Sapling",ylim=c(-0.025,0.025))











Ghat and Fhat with simulating bounds for Japanese pine sapling



-----------------------------Ghat-----------------------------------
?runifpoint
N = 65 # number of points to generate
r1 = runifpoint(N) #Generate N uniform random points
 r2 = runifpoint(N)
par(mar=c(1,1,1,1)) # set plot margin
plot(r1,pch="+",main="65 points under CSR")


r1.ghat=Gest(r1)
r1.ghat$rs 

par(mfrow=c(1,2))
plot(r1.ghat)
plot(r1.ghat$r, r1.ghat$rs,type="l",xlim=c(0,0.1))
r2.ghat=Gest(r2)
plot(r2.ghat)
plot(r2.ghat$r, r2.ghat$rs,type="l",xlim=c(0,0.1))
par(mfrow=c(1,1))

hold = matrix(0, nrow=100, ncol=length(r1.ghat$r)) #To save results
dim(hold)
for (i in 1:100) {
 rp =runifpoint(65)
 rp.ghat = Gest(rp, r1.ghat$r)
 rp.ghat.rs = rp.ghat$rs
 hold[i,] = rp.ghat.rs
}
r1.ghat$r[100] #the 100th distance point
summary(hold[,100]) #summary of Ghat at the 100th distance point
max(hold[,100])

min(hold[,100]) #lower bound of 100 simulations at the 100th point
apply(hold,2,max)[100]
# get the upper bound and lower bound at every point 
ubnd = apply(hold,2,max)
lbnd = apply(hold,2,min)

# plot the results with JP data (Figure 5)
plot(jp.ghat,rs~r,xlim=c(0,max(r1.ghat$r)))
par(mar=c(4,4.5,1.5,.5))
plot(jp.ghat,rs~r,xlim=c(0,max(r1.ghat$r)))
lines(r1.ghat$r,ubnd,lty=2,col=2)
lines(r1.ghat$r,lbnd,lty=2,col=2)
dev.copy(jpeg,"Gplot.jpg",width=5,height=5,units="in",res=300)
dev.off()







ghat.env = function(n, s, r, win=owin(c(0,1),c(0,1)))
{
 hold = matrix(0, s, length(r))
 for(i in 1:s)
 {
 hold[i,] = Gest(runifpoint(n, win=win), r=r)$rs
 }
 mn = apply(hold, 2, mean)
 Up = apply(hold, 2, max)
 Down = apply(hold, 2, min)
 return(data.frame(mn, Up, Down))
}
jp.ghat = Gest(jp)
jp.win = window(jp)
plot(jp.ghat,rs~r, main="G estimates")
jp.genv = ghat.env(n=jp$n, s=100, r=jp.ghat$r, win=jp.win)

# upper and lower envelopes
lines(jp.ghat$r, jp.genv$Up, lty=5, col=2)
lines(jp.ghat$r, jp.genv$Down, lty=5, col=3)


---------------------------Fhat------------------------------------

N = 65 # number of points to generate
r1 = runifpoint(N) #Generate N uniform random points
 r2 = runifpoint(N)
par(mar=c(1,1,1,1)) # set plot margin
plot(r1,pch="+",main="65 points under CSR")


r1.fhat=Fest(r1)
r1.fhat$rs 

par(mfrow=c(1,2))
plot(r1.fhat)
plot(r1.fhat$r, r1.fhat$rs,type="l",xlim=c(0,0.1))
r2.fhat=Fest(r2)
plot(r2.fhat)
plot(r2.fhat$r, r2.fhat$rs,type="l",xlim=c(0,0.1))
par(mfrow=c(1,1))

hold = matrix(0, nrow=100, ncol=length(r1.ghat$r)) #To save results
dim(hold)
for (i in 1:100) {
 rp =runifpoint(65)
 rp.fhat = Gest(rp, r1.fhat$r)
 rp.fhat.rs = rp.fhat$rs
 hold[i,] = rp.fhat.rs
}
r1.fhat$r[100] #the 100th distance point
summary(hold[,100]) #summary of Ghat at the 100th distance point
max(hold[,100])

min(hold[,100]) #lower bound of 100 simulations at the 100th point
apply(hold,2,max)[100]
# get the upper bound and lower bound at every point 
ubnd = apply(hold,2,max)
lbnd = apply(hold,2,min)

# plot the results with JP data (Figure 5)
plot(jp.fhat,rs~r,xlim=c(0,max(r1.fhat$r)))
par(mar=c(4,4.5,1.5,.5))
plot(jp.fhat,rs~r,xlim=c(0,max(r1.fhat$r)))
lines(r1.fhat$r,ubnd,lty=2,col=2)
lines(r1.fhat$r,lbnd,lty=2,col=2)
dev.copy(jpeg,"Fplot.jpg",width=5,height=5,units="in",res=300)
dev.off()



fhat.env = function(n, s, r, win=owin(c(0,1),c(0,1)))

{
 hold = matrix(0, s, length(r))
 for(i in 1:s)
 {
 hold[i,] = Fest(runifpoint(n, win=win), r=r)$rs
 }
 mn = apply(hold, 2, mean)
 Up = apply(hold, 2, max)
 Down = apply(hold, 2, min)
 return(data.frame(mn, Up, Down))
}

jp.fhat = Fest(jp)
jp.win = window(jp)
jp.fenv = fhat.env(n=jp$n, s=100, r=jp.fhat$r, jp.win)
plot(jp.fhat,rs~r, main="F estimates")

# upper and lower envelopes
lines(jp.fhat$r, jp.fenv$Up, lty=5, col=2)
lines(jp.fhat$r, jp.fenv$Down, lty=5, col=3)



Khat:

 Kest(runifpoint(n, win=win), r=r)$border


khat.env = function(n, s, r, win=owin(c(0,1),c(0,1)))

{
 hold = matrix(0, s, length(r))
 for(i in 1:s)
 {
 hold[i,] = Kest(runifpoint(n, win=win), r=r)$border
 }
 mn = apply(hold, 2, mean)
 Up = apply(hold, 2, max)
 Down = apply(hold, 2, min)
 return(data.frame(mn, Up, Down))
}

jp.khat = Kest(jp)
jp.win = window(jp)
jp.fenv = khat.env(n=jp$n, s=100, r=jp.khat$r, jp.win)
plot(jp.khat,border~r, main="K estimates")

# upper and lower envelopes
lines(jp.khat$r, jp.fenv$Up, lty=5, col=2)
lines(jp.khat$r, jp.fenv$Down, lty=5, col=3)


