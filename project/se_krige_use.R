require(pscl) #inverse gamma function
require(mvtnorm)
require(geoR)
require(MCMCpack) # this was in the mh code, but i can't remember what it was for
require(MTS) # to square root a matrix


# set data
#setwd("~/Documents/alm0823_hub/writing_project_s17/writing_project/data")
#http://leg.ufpr.br/geoR/geoRdoc/geoRintro.html


se <- read.csv("BugsataForAnalysisSpSurveyFinal.csv")
head(se)
colnames(se)

se.list <- c("xcoord", "ycoord", "SeWater", "SeBugs" ,"SeSed","SeRoots")
require(psych)
cor.plot(se[,se.list[-c(1,2)]])
#The fields you are likely interested in for now are the xcoord, ycoord 
#(locations for the sample points); SiteID the unique id for a survey location; 

#mdcaty.x = refers to the Bugs zones within a management unit: 
#flooded, saturated, intermittent; stratum.x = the unit name on the refuge so the different management units;

#Acres and Hectares are the size of the areal frame within each of those Bugs zones within a unit; 
#Date_ = survey date when the data was collected;
#SeSed is selenium concentration within a Sediment sample;
#SeRoot is selenium concentration within a root sample; 
#SeBugs is selenium concentration within a Bugs sample; and 
#SeSed is selenium concentration within a bug sample.



se.cols <- c("xcoord", "ycoord", "SeSed", "SeBugs")

se_sub <- se[,se.cols]

require(sp)
require(gstat)
require(lattice)

dim(se_sub)
row.na <- apply(se_sub[,3:4], 1, function(t){all(is.na(t))})

se_sub <- se_sub[-which(row.na == TRUE),]
dim(se_sub)

length(which(is.na(se_sub$SeBugs)))

par(mfrow=c(1,2))
hist(se_sub$SeSed)
hist(se_sub$Sed.log <- log10(se_sub$SeSed))


histogram(se_sub$SeBugs)
histogram(se_sub$Bugs.log <- log10(se_sub$SeBugs), Nint=12)


# try to predict Bugs given Sed
par(mfrow=c(1,1))
xyplot(SeBugs ~ SeSed, as.data.frame(se_sub))
cor(se_sub$SeSed[-which(is.na(se_sub$SeBugs))], se_sub$SeBugs[-which(is.na(se_sub$SeBugs))])^2

se_sub <- data.frame(se_sub)[,-c(6,7)]

se_sub.Sed <- se_sub[,c("xcoord", "ycoord", "SeSed")]
se_sub.Sed <- se_sub.Sed[complete.cases(se_sub.Sed$SeSed),]
se_sub.Bugs <- se_sub[,c("xcoord", "ycoord", "SeBugs")]
se_sub.Bugs <- se_sub[complete.cases(se_sub.Bugs$SeBugs),]
se_sub <- data.frame(se_sub[,c("xcoord", "ycoord", "SeSed", "SeBugs")])

se.extra <- (se_sub[setdiff(rownames(data.frame(se_sub)), rownames(se_sub.Bugs)),])

pred.vals <- c(1,2,8,9,10,11)

se_pred.Bugs <- se_sub.Bugs[pred.vals,]
se_sub.Bugs <- data.frame(se_sub.Bugs[-pred.vals,])

se_all <- data.frame(rbind(se_sub.Bugs[,1:2], se_sub.Sed[,1:2]))
#coordinates(se_all) <- ~ xcoord + ycoord

coordinates(se.extra) <- ~ xcoord + ycoord
coordinates(se_sub.Sed) <- ~ xcoord + ycoord
coordinates(se_sub.Bugs) <- ~ xcoord + ycoord
coordinates(se_sub) <- ~ xcoord + ycoord
coordinates(se_sub1.Bugs) <- ~ xcoord + ycoord

coordinates(se_pred.extra) <- ~ xcoord + ycoord
coordinates(se_pred.Bugs) <- ~ xcoord + ycoord


#se_sub1.Bugs = locations of bugs to use
#se_pred.Bugs = locations of bugs to predict at


# empirical variogram
par(mfrow=c(1,2))
v.Bugs <- variogram(SeBugs ~ 1, data = se_sub.Bugs)
plot(v.Bugs, main = "Bugs [SE] Empirical Variogram")
v.Sed <- variogram(SeSed ~ 1, data=se_sub.Sed)
plot(v.Sed, main = "Sed [SE] Empirical Variogram")

par(mfrow=c(1,1))
v.pred.Bugs <- variogram(SeBugs ~ 1, data = se_pred.Bugs)
plot(v.pred.Bugs, main = "Subsetted Bugs [SE] Empirical Variogram")

m.pred.Bugs <- vgm(1-0.2,"Exp",600,0.2)
plot(v.pred.Bugs, pl=T, model=m.pred.Bugs)
# fit model parameters by weighted least-squares
(m.pred.Bugs.f <- fit.variogram(v.pred.Bugs, m.pred.Bugs))
plot(v.pred.Bugs, pl=T, model=m.pred.Bugs.f)


par(mfrow=c(1,2))
# estimate variogram model form and parameters by eye
# options = psill, "Exp", range, nugget
m.Bugs <- vgm(1-0.5,"Exp",1500,0.5)
plot(v.Bugs, pl=T, model=m.Bugs)
# fit model parameters by weighted least-squares
(m.Bugs.f <- fit.variogram(v.Bugs, m.Bugs))
plot(v.Bugs, pl=T, model=m.Bugs.f)

m.Sed <- vgm(1.7-0.5,"Exp",1500, 0.5)
plot(v.Sed, pl=T, model=m.Sed)
# fit model parameters by weighted least-squares
(m.Sed.f <- fit.variogram(v.Sed, m.Sed))
plot(v.Sed, pl=T, model=m.Sed.f)


# idea: went from intial estimates and then automatically fit by least squares

se.grid <- se_sub[,1:2]
coordinates(se.grid) <- ~ xcoord + ycoord

k.o <- krige(SeBugs ~1, locations=se_sub.Bugs, newdata=se.grid, model=m.Bugs.f)
# summary statistics
summary(k.o)

k.pred.o <- krige(SeBugs~1, locations = se_pred.Bugs, newdata=se.grid, model=m.pred.Bugs.f)

#  Display a map of the predictions and their errors.
#plot.kresults(k.o, "var1", jura.dat, jura.Ni, f=1,"OK")

# predict at the extra points
k <- krige(SeBugs ~ 1, se_sub.Bugs, se.extra, m.Bugs.f)

# predict at extra points using observed Bugs
k.pred <- krige(SeBugs~1, se_pred.Bugs, se_pred.extra, m.pred.Bugs.f)
# Nioss.valid
# Compute and summarize validation errors
summary(k)


diff <- k.pred$var1.pred - se_pred.extra$SeBugs
summary(diff)
sqrt(sum(diff^2)/length(diff)) # RMSE (precision)
sum(diff)/length(diff) # mean error (bias)
median(se_pred.extra$SeBugs) # median error

diff <- as.data.frame(diff)
coordinates(diff) <- coordinates(se_pred.extra)
bubble(diff, zCol="diff", main="OK validation errors at undersampled points, Ni")

cv.o <- krige.cv(SeBugs ~ 1, se_sub.Bugs, model=m.Bugs.f, nfold=nrow(se_sub.Bugs))
#nfold = cross validation once for each point
summary(cv.o)
res <- as.data.frame(cv.o)$residual
sqrt(mean(res^2))
mean(res)
mean(res^2/as.data.frame(cv.o)$var1.var)
rm(res)

#plot.valids(k, "var1", se_pred.extra, "Bugs", cv.o, "OK")


diff <- (k.pred$var1.pred) - (se_pred.extra$SeBugs)
summary(diff)
histogram(diff, col="darkseagreen2", nint=12, type="count",
          main="Validation errors", xlab="Bugs, [SE]")
diff <- (cv.o$var1.pred) - (se_sub.Bugs$SeBugs)
summary(diff)
histogram(diff, col="lightblue2", nint=12, type="count",
          main="Cross-Validation errors", xlab="Bugs, [SE]")
rm(diff)

#  Interpolate at the extra points, using OK from the sample
#   points, and determine the actual prediction bias and precision.


# section 6
m.Bugs.f$range[2]; m.Sed.f$range[2]
# can't use other variables because different ranges
# compare variogram structure to target variable
m.Bugs.f$range[2]; m.Sed.f$range[2]
#ratio of nugget to total sill
round(m.Bugs.f$psill[1]/sum(m.Bugs.f$psill),2)
round(m.Sed.f$psill[1]/sum(m.Sed.f$psill),2)


#Now comes the hardest part of co-kriging: modelling the co-regionalisation.
##We have to fit models to both the direct and cross-variograms
#simultaneously, and these models must lead to a positive definite cokriging
#system. The easiest way to ensure this is to fit a linear model of
#co-regionalisation: all models (direct and cross) have the same shape and
#range, but may have different partial sills and nuggets.

# build gstat object with subsample of and complete sample of 
(g <- gstat(NULL, id = "SeBugs", form = SeBugs ~ 1, data=se_sub.Bugs))
(g <- gstat(g, id = "SeSed", form = SeSed ~ 1, data=se_sub.Sed))

duplicated(g$SeBugs$data@coords)
#Compute and display the two direct variograms and one
#cross-variogram.

v.cross <- variogram(g)
str(v.cross)
plot(v.cross, pl=T)

# compare cross variogram to direct variograms
# Add variogram models to the gstat object and fit a them
# using the linear model of co-regionalisation.

# requires similar range and model structure
# fit using that of target variable (bugs)

(g <- gstat(g, id = "SeBugs", model = m.Bugs.f, fill.all=T))


#Now we fit all three variograms together, ensuring they lead to a positive
#definite co-kriging system. For this we use the fit.lmc method
#("fit linear model of co-regionalization"). This takes the initial estimate,
#fits all the variograms, and then each of the partial sills is adjusted (by
#least squares) to the closest value that will result in a positive definite
#matrices
zerodist2(se_sub.Bugs, se_sub.Sed)

# note default fit.lmc(fit.method) is fit.method=7, which doesn't allow zero
# distances because it weights based on count within distance h divided by distance h
# individually, I did not have duplicated coordinates for sebugs and sesed, but
# did for the combined either v.cross object or g object
# so i got an error
# it seems that I would want to have duplicated at some points, but not all,
# so i don't understand this
# so use unweighted = OLS

(g <- fit.lmc(v.cross, g, fit.method = 6))
plot(variogram(g), model=g$model)

#any(duplicated(coordinates(v.cross)))

# ERRORS SKIP TO SECTION 7
k.c <- predict(g, se_pred.extra)


###
a <- read.table("http://www.stat.ucla.edu/~nchristo/statistics_c173_c273/soil.txt", header=TRUE) 

# note this data set has all co-located variables
# use log sed but not log bugs
# 1 = variable want to predict

#CO-KRIGING:
#Create a gstat object:
g1 <- NULL
g1 <- gstat(id="bugs", formula = SeBugs~1, data = se_sub.Bugs) 
#Append:
g1 <- gstat(g1,id="log_sed", formula = log(SeSed)~1, data = se_sub) 

#Plot the variograms and cross-variograms:
plot(variogram(g1))

#Fit a model variogram to the target variogram:
g <- gstat(id="bugs", formula = SeBugs~1, data = se_sub.Bugs) 

v.fit <- fit.variogram(variogram(g), vgm(0.5,"Sph",1000,0.1)) 

#Fit a model variogram to all the variograms:
g1 <- gstat(id="bugs", formula = SeBugs~1, data = se_sub.Bugs) 
g1 <- gstat(g1,id="log_sed", formula = log(SeSed)~1, data = se_sub) 

vm <- variogram(g1) 

vm.fit <- fit.lmc(vm, g1, model=v.fit, fit.method = 6) 

#Plot the fitted variograms to all the sample variograms:
plot(variogram(g1),vm.fit)

#Create the grid for predictions:
x.range <- as.integer(range(coordinates(se_sub)[1])) 
y.range <- as.integer(range(coordinates(se_sub)[2])) 
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=50), 
                   y=seq(from=y.range[1], to=y.range[2], by=50)) 

#Perform co-kriging predictions:
ck <- predict.gstat(vm.fit, grd) 

#Perform cross-validation:
cv_ck <- gstat.cv(vm.fit)


#COMPARE ORDINARY KRIGING WITH CO-KRIGING:
sum(cv_ok$residual^2)

sum(cv_ck$residual^2)


