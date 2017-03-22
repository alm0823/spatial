y1 <- se_sub.Bugs$SeBugs
y2 <- se_sub$SeSed

x1 <- data.frame(cbind(c(rep(1,length(y1))), 
                       coordinates(se_sub.Bugs)[1], coordinates(se_sub.Bugs)[2]))
x2 <- data.frame(cbind(c(rep(1,length(y2))),
                       coordinates(se_sub.Sed)[1], coordinates(se_sub.Sed)[2]))
x0 <- data.frame(cbind(c(rep(1,length(se_pred.Bugs))), coordinates(se_pred.Bugs)[1],
                       coordinates(se_pred.Bugs)[2]))
dim(x1)
dim(x2)

x1 <- c(rep(1,length(y1)))
x2 <- c(rep(1,length(y2)))
X <- data.frame(c(x1,x2))

#C11 == covariance matrix of Bugs
#we assume it is known
#i estimated it using OLS
#not accounting for the fact that i estimated it
# nugget = 0.2107182
# sill = 3.86876
# range = 4239.226

# for h put in gstat object with coordinates of interest

exp.fn <- function(tau2,sigma2,phi,h){
t <- as.matrix(dist(rbind(c(0,0),coordinates(h)[1:nrow(h),])))[-1,-1]
diag(t) <- 1
out <- ifelse(t>0,sigma2*exp(-t/phi),tau2+sigma2)
return(as.matrix(out, nrow=nrow(t), ncol=ncol(t)))
}

C11 <- exp.fn(0.2107182, 3.86876-0.2107182, 4239.226, se_sub.Bugs)
dim(C11)
C11[1:5,1:5]

# nugget = 0
# sill = 50=.55
# range = 222.551
C22 <- exp.fn(0,5.55,222.551,se_sub.Sed)
dim(C22)
C22[1:5,1:5]

# nugget = 0.15
# sill = 1.45877
# range = 180.0092
cbind(se_sub.Bugs, se_sub.Sed)
CALL <- exp.fn(0,1.45877,180.0092, se_all)
C12 <- CALL[36:266,1:35]
C21 <- CALL[1:35,36:266]

C.row <- data.frame(rbind(C11,C12))
C.col <- data.frame(rbind(C21,C22))
C <- data.frame(cbind(C.row,C.col))
dim(C)

# note I had to do this the long way to use the three different covariance functions
X <- as.matrix(X)
C <- as.matrix(C)
CALL <- as.matrix(CALL)

V <- solve(t(X)%*%solve(CALL)%*%X)
dim(C)
dim(X)

############ ASK STEVE, BECAUSE X CONTAINS DUPLICATE LOCATIONS, THIS IS SOLVE(t(x)%*%x)
# is singular