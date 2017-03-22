require(sp)
require(gstat)
require(lattice)

# use Co
# dont use Ni, Ni, Ni, Ni
# Ni is most promising

# use jura dataset
data("jura")
head(jura.val)
jura.dat <- NULL
jura.dat <- jura.val[,c("long", "lat", "Co", "Ni")]
dim(jura.dat)

jura.max <- apply(jura.dat,2,max)
jura.min <- apply(jura.dat,2,min)
lat.range <- c(seq(jura.min[2], jura.max[2], length.out = 1000))
long.range <- c(seq(jura.min[1], jura.max[2], length.out = 1000))
jura.new <- NULL
jura.new$lat <- sample(lat.range, 10, replace = FALSE)
jura.new$long <- sample(long.range, 10, replace = FALSE)

dist.seq <- c(seq(0,0.005,10))
jura.new$pairs.lat <- jura.new$lat
jura.new$pairs.long <- jura.new$long + dist.seq
jura.new$co <- c(seq(min(jura.Co$Co), max(jura.Co$Co), length=10))
k <- c(seq(1:10))
jura.new$pairs.co <- sqrt(k) + jura.new$co

lt <- c(jura.dat$lat, jura.new$lat, jura.new$pairs.lat)
lg <- c(jura.dat$long, jura.new$long, jura.new$pairs.long)
out <- c(jura.dat$Co, jura.new$co, jura.new$pairs.co)
out.ni <- c(jura.dat$Ni, rep(NA, 10))

jura.dat <- NULL
jura.dat <- data.frame(cbind(lat=lt, long=lg, Co=out, Ni=out.ni))

jura.dat <- jura.dat[-which(duplicated(jura.dat[,c("Co")]) == TRUE),]
jura.dat <- jura.dat[-which(duplicated(jura.dat[,c("Ni")]) == TRUE),]
dim(jura.dat)

jura.Ni <- jura.val[,c("long", "lat", "Ni")]
jura.Co <- jura.val[,c("long", "lat", "Co")]

jura.all <- jura.dat



par(mfrow=c(1,1))
histogram(jura.dat$Ni)
histogram(log10(jura.dat$Ni))


histogram(jura.dat$Co, Nint=12)
histogram(log10(jura.dat$Co), Nint=12)


# try to predict Co given Ni
par(mfrow=c(1,1))
xyplot(Ni ~ Co, as.data.frame(jura.dat))
cor(jura.dat$Ni,jura.dat$Co)^2
#51.9% of variation in Nickel explained by Copper

# Ni == pull low little lower and medium little higher
##jura.Ni$Ni[which(jura.Ni$Ni < quantile(jura.Ni$Ni, 0.25))] <- 
  #jura.Ni$Ni[which(jura.Ni$Ni < quantile(jura.Ni$Ni, 0.25))] - 3
#jura.Ni$Ni[which(jura.Ni$Ni > quantile(jura.Ni$Ni,0.5))] <- 
  #jura.Ni$Ni[which(jura.Ni$Ni > quantile(jura.Ni$Ni,0.5))] + 5

summary(jura.dat)
# no missing values
# make every fourth Ni value missing

jura.Ni$Ni[sample(1:dim(jura.Ni)[1],25, replace = FALSE)] <- NA

#1. The variogram method Computes the experimental variogram as a
#variogram object;
#2. The vgm method Nieates a variogram model object;
#3. The fit.variogram method adjusts a variogram model object to
#a variogram object.

jura.Ni <- jura.Ni[-which(is.na(jura.Ni$Ni)),]
jura.extra <- jura.all[setdiff(rownames(jura.all), rownames(jura.Ni)),
                       c("long", "lat", "Ni")]



coordinates(jura.dat) <- ~ long + lat
coordinates(jura.Ni) <- ~ long + lat
coordinates(jura.grid) <- ~ long + lat
coordinates(jura.Co) <- ~ long + lat

coordinates(jura.extra) <- ~ long+lat


# empirical variogram
par(mfrow=c(1,2))
v.Co <- variogram(Co ~ 1, data = jura.Co)
plot(v.Co, main = "Copper Empirical Variogram")
v.Ni <- variogram(Ni ~ 1, data=jura.Ni)
plot(v.Ni, main = "Ni Empirical Variogram")

# estimate variogram model form and parameters by eye
# options = psill, "Exp", range, nugget
m.Co <- vgm(0.8-0.2,"Sph",0.04,0.2)
plot(v.Co, pl=T, model=m.Co)
# fit model parameters by weighted least-squares
(m.Co.f <- fit.variogram(v.Co, m.Co))
plot(v.Co, pl=T, model=m.Co.f)

m.Ni <- vgm(77-10,"Sph",0.015, 10)
plot(v.Ni, pl=T, model=m.Ni)
# fit model parameters by weighted least-squares
(m.Ni.f <- fit.variogram(v.Ni, m.Ni))
plot(v.Ni, pl=T, model=m.Ni.f)


# idea: went from intial estimates and then automatically fit by least squares

k.o <- krige(Ni ~1, locations=jura.Ni, newdata=jura.grid, model=m.Ni.f)
# summary statistics
summary(k.o)

#  Display a map of the predictions and their errors.
#plot.kresults(k.o, "var1", jura.dat, jura.Ni, f=1,"OK")


  # predict at the extra points
  k <- krige(Ni ~ 1, jura.Ni, jura.extra, m.Ni.f)

# Nioss.valid
# Compute and summarize validation errors
summary(k)
diff <- k$var1.pred - jura.extra$Ni
summary(diff)
sqrt(sum(diff^2)/length(diff)) # RMSE (precision)
sum(diff)/length(diff) # mean error (bias)
median(jura.extra$Ni) # median error

diff <- as.data.frame(diff)
coordinates(diff) <- coordinates(jura.extra)
bubble(diff, zCol="diff", main="OK validation errors at undersampled points, Ni")

cv.o <- krige.cv(Ni ~ 1, jura.Ni, model=m.Ni.f, nfold=nrow(jura.Ni))
summary(cv.o)
res <- as.data.frame(cv.o)$residual
sqrt(mean(res^2))
mean(res)
mean(res^2/as.data.frame(cv.o)$var1.var)
rm(res)

plot.valids(k, "var1", jura.extra, "Ni", cv.o, "OK")

#  Interpolate at the extra points, using OK from the sample
#   points, and determine the actual prediction bias and precision.


# section 6

# 0. Compare variogram structure to target variable
# psill = nugget  partial sill
m.Co.f$range[2]; m.Ni.f$range[2]
round(m.Co.f$psill[1]/sum(m.Co.f$psill),2)
round(m.Ni.f$psill[1]/sum(m.Ni.f$psill),2)

# build gstat model piece by piece

# 1. Nieate gstat object with two frames
# both direct and Nioss variograms can be Computed and plotted from these
g <- NULL

g <- gstat(NULL, id = "Ni", form = Ni~1, data=jura.Ni)
g <- gstat(g, id = "Co", form = Co~1, data=jura.Co)

v.Nioss <- variogram(g)
str(v.Nioss)
plot(v.Nioss, pl=T)


# 2.  Add variogram models to the gstat object and fit a them
# using the linear model of Co-regionalisation.
g <- gstat(g, id = "Ni", model = m.Ni.f, fill.all=T)
#We use the fitted model for the target value as a starting point for all
#three variogram models. This is because we will use the linear model of
#Co-regionalisation, which requires a single range and structure. By filling
#all the frames with one model (using the fill.all = T argument), these
#Conditions are automatically met.

(g <- fit.lmc(v=v.Nioss, g=g))

plot(variogram(g), model=g$model)
# The gstat object now Contains both the data and the models
# Note: Could specify individually for each direct and Nioss-variogram
# however, diffiColt to ensure the resulting CK system is positive defiNite

# 3. Fit all three variograms together
# This takes the iNitial estimate,
# fits all the variograms, and then each of the partial sills is adjusted (by
# least squares) to the closest value that will result in a positive defiNite
# matrices.

(g <- fit.lmc(v.Nioss, g, fit.method = 1))
plot(variogram(g), model=g$model)
#  Note that the ranges are not adjusted by
##fit.lmc, however all the partial sills (both for the spherical model and# for the nugget) of the Co-variable and Nioss-variogram were adjusted. So,
#the fitted variograms have the same range but different sills and nuggets.

# Comparing fits for direct variograms from Nioss (Coregulization)
# and individually fit (regulization)

g$model$Co$psill - m.Co.f$psill
sum(g$model$Co$psill) - sum(m.Co.f$psill)
sum(g$model$Co$psill)
g$model$Ni$psill - m.Ni.f$psill
sum(g$model$Ni$psill) - sum(m.Ni.f$psill)
sum(g$model$Ni$psill)

## at k.c endpoint, make plot using 
## spplot()
## works instead

## in their plot function, change quotes to make their functions work

## try demo(cokriging)
## use fit.method=6

## steve wants to see more theory
## so what is a co-variogram
## see the text he gave me

# email author and see why the default isn't working
