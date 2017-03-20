# Some (but not all) R code that may be useful. 
# Problem 3
require(spatstat)

blackoak<-split(lansing)$blackoak
maple<-split(lansing)$maple

Lblackoak<-envelope(blackoak,fun=Lest,correction="iso")
plot(Lblackoak,.-r~r,legend=F)
Lmaple<-envelope(maple,fun=Lest,correctin="iso")
plot(Lmaple,.-r~r,legend=F)

require(splancs)
# specify radii
h<-seq(0,.5,l=100)
# get coordinates
tree.poly<-list(x=c(blackoak$x,maple$x),y=c(blackoak$y,maple$y))
# recompute the K functions
kblackoak<-khat(as.points(blackoak),bboxx(bbox(as.points(tree.poly))),h)
kmaple<-khat(as.points(maple),bboxx(bbox(as.points(tree.poly))),h)
# get the differences
k.diff<-kblackoak - kmaple
# generate the envelope
env<-Kenv.label(as.points(blackoak),as.points(maple),
bboxx(bbox(as.points(tree.poly))),nsim=99,s=h)
# plot the results
plot(h,seq(-0.15,0.05,l=length(h)),type="n",ylab="Kdiff",
main="Envelopes for Kdiff")
lines(h,k.diff)
lines(h,env$low,lty=2)
lines(h,env$up,lty=2)
abline(h=0)

Kenv<-envelope(lansing,Kcross, i="maple",j="blackoak")
plot(Kenv,sqrt(./pi)-r~r,ylab="Lij - h",main="Cross L Function - Maple
and Black Oak",legend=F)

plot(density(blackoak))
plot(density(maples))

# Problem 5

# get the CN data
names(CN.dat)<-c("x","y","tn","tc","cn")
attach(CN.dat)
CN.lm<-lm(tc~tn)
resids<-residuals(CN.lm)
# convert to a geodata object 
require(geoR)
resids.dat<-cbind(CN.dat$x,CN.dat$y,resids)
resids.dat<-data.frame(resids.dat)
names(resids.dat)<-c("x","y","resids")
resids.geodat<-as.geodata(resids.dat,coords.col=1:2,data.col=3)
# now fit the models. 


