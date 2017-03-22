require(pscl) #inverse gamma function
require(mvtnorm)
require(geoR)
require(MCMCpack) # this was in the mh code, but i can't remember what it was for
require(MTS) # to square root a matrix


# set data
#setwd("~/Documents/alm0823_hub/writing_project_s17/writing_project/data")
#http://leg.ufpr.br/geoR/geoRdoc/geoRintro.html


se <- read.csv("SEDataForAnalysisSpSurveyFinal.csv")
head(se)
colnames(se)

se.list <- c("xcoord", "ycoord", "SeWater", "SeBugs" ,"SeSed","SeRoots")
require(psych)
cor.plot(se[,se.list[-c(1,2)]])
#The fields you are likely interested in for now are the xcoord, ycoord 
#(locations for the sample points); SiteID the unique id for a survey location; 

#mdcaty.x = refers to the water zones within a management unit: 
#flooded, saturated, intermittent; stratum.x = the unit name on the refuge so the different management units;

#Acres and Hectares are the size of the areal frame within each of those water zones within a unit; 
#Date_ = survey date when the data was collected;
#SeSed is selenium concentration within a sediment sample;
#SeRoot is selenium concentration within a root sample; 
#SeWater is selenium concentration within a water sample; and 
#SeBugs is selenium concentration within a bug sample.



se.cols <- c("xcoord", "ycoord", "SeBugs", "SeWater")

se_sub <- se[,se.cols]

require(sp)
require(gstat)
require(lattice)

dim(se_sub)
row.na <- apply(se_sub[,3:4], 1, function(t){all(is.na(t))})

se_sub <- se_sub[-which(row.na == TRUE),]
dim(se_sub)

length(which(is.na(se_sub$SeWater)))

par(mfrow=c(1,2))
hist(se_sub$SeBugs)
hist(se_sub$bugs.log <- log10(se_sub$SeBugs))


histogram(se_sub$SeWater)
histogram(se_sub$water.log <- log10(se_sub$SeWater), Nint=12)


# try to predict water given sed
par(mfrow=c(1,1))
xyplot(SeWater ~ SeBugs, as.data.frame(se_sub))
cor(se_sub$SeBugs[-which(is.na(se_sub$SeWater))], se_sub$SeWater[-which(is.na(se_sub$SeWater))])^2

se_sub.Bugs <- se_sub[,c("xcoord", "ycoord", "SeBugs")]
se_sub.Bugs <- se_sub.Bugs[complete.cases(se_sub.Bugs$SeBugs),]
se_sub.Water <- se_sub[,c("xcoord", "ycoord", "SeWater")]
se_sub.Water <- se_sub[complete.cases(se_sub.Water$SeWater),]
se_sub <- data.frame(se_sub[,c("xcoord", "ycoord", "SeBugs", "SeWater")])

se.extra <- (se_sub[setdiff(rownames(data.frame(se_sub)), rownames(se_sub.Water)),])

se_pred.Water <- data.frame(se_sub.Water[-c(1,2,8,9,10,11),])
se_pred.extra <- se_sub[setdiff(rownames(data.frame(se_sub.Water)), rownames(data.frame(se_pred.Water))),]

coordinates(se.extra) <- ~ xcoord + ycoord
coordinates(se_sub.Bugs) <- ~ xcoord + ycoord
coordinates(se_sub.Water) <- ~ xcoord + ycoord
coordinates(se_sub) <- ~ xcoord + ycoord

coordinates(se_pred.extra) <- ~ xcoord + ycoord
coordinates(se_pred.Water) <- ~ xcoord + ycoord



# empirical variogram
par(mfrow=c(1,2))
v.water <- variogram(SeWater ~ 1, data = se_sub.Water)
plot(v.water, main = "Water [SE] Empirical Variogram")
v.bugs <- variogram(SeBugs ~ 1, data=se_sub.Bugs)
plot(v.bugs, main = "Bugs [SE] Empirical Variogram")

par(mfrow=c(1,1))
v.pred.water <- variogram(SeWater ~ 1, data = se_pred.Water)
plot(v.pred.water, main = "Subsetted Water [SE] Empirical Variogram")

par(mfrow=c(1,2))
# estimate variogram model form and parameters by eye
# options = psill, "Exp", range, nugget
m.water <- vgm(1-0.5,"Sph",1500,0.5)
plot(v.water, pl=T, model=m.water)
# fit model parameters by weighted least-squares
(m.water.f <- fit.variogram(v.water, m.water))
plot(v.water, pl=T, model=m.water.f)

m.pred.water <- vgm(1-0.5,"Sph",1500,0.5)
plot(v.pred.water, pl=T, model=m.pred.water)
# fit model parameters by weighted least-squares
(m.pred.water.f <- fit.variogram(v.pred.water, m.pred.water))
plot(v.pred.water, pl=T, model=m.pred.water.f)

m.bugs <- vgm(1.7-0.5,"Sph",1500, 0.5)
plot(v.bugs, pl=T, model=m.bugs)
# fit model parameters by weighted least-squares
(m.bugs.f <- fit.variogram(v.bugs, m.bugs))
plot(v.bugs, pl=T, model=m.bugs.f)


# idea: went from intial estimates and then automatically fit by least squares

se.grid <- se_sub[,1:2]
coordinates(se.grid) <- ~ xcoord + ycoord

k.o <- krige(SeWater ~1, locations=se_sub.Water, newdata=se.grid, model=m.water.f)
# summary statistics
summary(k.o)

k.pred.o <- krige(SeWater~1, locations = se_pred.Water, newdata=se.grid, model=m.pred.water.f)

#  Display a map of the predictions and their errors.
#plot.kresults(k.o, "var1", jura.dat, jura.Ni, f=1,"OK")

# predict at the extra points
k <- krige(SeWater ~ 1, se_sub.Water, se.extra, m.water.f)

# predict at extra points using observed water
k.pred <- krige(SeWater~1, se_pred.Water, se_pred.extra, m.pred.water.f)
# Nioss.valid
# Compute and summarize validation errors
summary(k)


diff <- k.pred$var1.pred - se_pred.extra$SeWater
summary(diff)
sqrt(sum(diff^2)/length(diff)) # RMSE (precision)
sum(diff)/length(diff) # mean error (bias)
median(se_pred.extra$SeWater) # median error

diff <- as.data.frame(diff)
coordinates(diff) <- coordinates(se_pred.extra)
bubble(diff, zCol="diff", main="OK validation errors at undersampled points, Ni")

cv.o <- krige.cv(SeWater ~ 1, se_sub.Water, model=m.water.f, nfold=nrow(se_sub.Water))
#nfold = cross validation once for each point
summary(cv.o)
res <- as.data.frame(cv.o)$residual
sqrt(mean(res^2))
mean(res)
mean(res^2/as.data.frame(cv.o)$var1.var)
rm(res)

#plot.valids(k, "var1", se_pred.extra, "Water", cv.o, "OK")


diff <- (k.pred$var1.pred) - (se_pred.extra$SeWater)
summary(diff)
histogram(diff, col="darkseagreen2", nint=12, type="count",
          main="Validation errors", xlab="Water, [SE]")
diff <- (cv.o$var1.pred) - (se_sub.Water$SeWater)
summary(diff)
histogram(diff, col="lightblue2", nint=12, type="count",
          main="Cross-Validation errors", xlab="Water, [SE]")
rm(diff)

#  Interpolate at the extra points, using OK from the sample
#   points, and determine the actual prediction bias and precision.


# section 6

# 0. Compare variogram structure to target variable
# psill = nugget  partial sill
m.water.f$range[2]; m.bugs.f$range[2]
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

(g <- fit.lmc(v.Nioss, g))
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


