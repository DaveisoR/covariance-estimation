rm(list = ls())
Returns <- read.csv("D:/Research/R data/ret no dates.csv", header = TRUE)
Relative.volatility <- read.csv("D:/Research/R data/relative volatility no dates.csv", header = TRUE)
Volume <- read.csv("D:/Research/R data/volume no dates.csv", header = TRUE)

stderr <- Relative.volatility/(Volume^0.5)

Relative.range <- read.csv("D:/Research/R data/relative range.csv", header = TRUE)
Relative.volatility2 <- Relative.volatility[-2709, ]
rm(Relative.volatility)
holding <- 21
lookback <- 252
periods <- trunc((nrow(Returns)-lookback)/holding)
sundry <- as.list(c(holding, lookback, periods))
names(sundry) <- c("holding", "lookback", "periods")

###################SORTING FUNCTIONS############################

lookback.bucket <- function(datasource, period, sundry){
  datasource[(1+(sundry$holding*(period-1))):(sundry$lookback+(sundry$holding*(period-1))), ]
}

holding.bucket <- function(datasource, period, sundry){
  datasource[(1+sundry$lookback+(sundry$holding*(period-1))):(sundry$lookback+(sundry$holding*(period))), ]
}

####################SORTING##############################

lb.returns <- vector(mode = "list", length = periods)
hld.returns <- vector(mode = "list", length = periods)

for(x in 1:periods){
  lb.returns[[x]] <- lookback.bucket(Returns, x, sundry)
  hld.returns[[x]] <- holding.bucket(Returns, x, sundry)
}

lb.stderr <- vector(mode = "list", length = periods)
hld.stderr <- vector(mode = "list", length = periods)

for(x in 1:periods){
  lb.stderr[[x]] <- lookback.bucket(stderr, x, sundry)
  hld.stderr[[x]] <- holding.bucket(stderr, x, sundry)
}

###################COVARIANCE CALCS#######################

co.var <- cov.wt(lb.returns[[1]], wt = 1/(lb.stderr[[1]][ ,1]), method = "unbiased")
#Note it only allows one weight vector, not an entire weight matrix like we would like.

#manual attempt

#PREPARING THE DATA
#standardised weights function **??????Do we divide by the average so that the returns arent shrunk prematurely? Otherwise divide by sum. standardising formula gives negative weights.
std.weights <- function(object){
  (1/object)/mean(1/object)
}

#apply function "std.weights" to each column of each matrix in the list.
lb.weights <- lapply(lb.stderr, function(x){apply(x, 2, std.weights)})
hld.weights <- lapply(hld.stderr, function(x){apply(x, 2, std.weights)})

#weighted returns (this is unecessary if you use the weighted.mean function)
#lb.w.ret <- Map("*", lb.weights, lb.returns)

#deviations function

# lb.weighted.devs <- Map("-", lb.returns, lapply(lb.w.ret, function(a){apply(a, 2, mean)}))
#the above subtracts the mean of the weighted returns from the unweighted returns but it doesnt seem to pass the test below

#set the structure of the list and the matrices contained in each list element beforehand
element.matrix <- matrix(0, 252, 24)
lb.weighted.devs <- rep(list(element.matrix), 116)
hld.weighted.devs <- rep(list(element.matrix), 116)
for (a in 1:116){
  for(b in 1:24){
    lb.weighted.devs[[a]][ , b] <- (lb.returns[[a]][ , b]-weighted.mean(lb.returns[[a]][ , b], lb.weights[[a]][ , b]))
    hld.weighted.devs[[a]][ , b] <- (hld.returns[[a]][ , b]-weighted.mean(hld.returns[[a]][ , b], hld.weights[[a]][ , b]))
  }
}
#this method works for below test!!!

#test (dont match after first 1.). why...? Probably because the list of means is unequal to the list of unweighted returns
b <- (lb.weighted.devs[[1]][ , 1])
c <- (lb.returns[[1]][ , 1]-weighted.mean(lb.returns[[1]][ , 1], lb.weights[[1]][ , 1]))
#comment below sorted by code above

#Need to calculate joint weights for all combinations *probably have to use loops
##***????*** Do i need to weight the deviations seperately if i have already weighted the returns and not just the means to calc devs...? #but then the deviation can actually be upweighted when the weight is lower if the mean is positive and the return is less than the mean value. ie downweight makes deviation greater...? which means the "mean" should be weighted but subtracted from actual ret and not weighted ret. (ie adjust deviation function)

weighted.mean(lb.weighted.devs[[1]][ , 1], lb.weights[[1]][ , 1])
#shows that mean of deviations is basically zero, as it should be
#Need to calculate joint weights. require triple nested loop

cov.matrix <- matrix(0,24, 24)
lb.example <- rep(list(cov.matrix), 116)
hld.example <- rep(list(cov.matrix), 116)
for(l in 1:116){
  for(f in 1:24){
    for(s in 1:24){
      lb.example[[l]][f, s] <- sum(((lb.weights[[l]][ , f]*lb.weights[[l]][ , s])^0.5)*(lb.weighted.devs[[l]][ , f]*lb.weighted.devs[[l]][ , s]))/251
      hld.example[[l]][f, s] <- sum(((hld.weights[[l]][ , f]*hld.weights[[l]][ , s])^0.5)*(hld.weighted.devs[[l]][ , f]*hld.weighted.devs[[l]][ , s]))/251
    }
  }
}
# example.1 <- sum(((lb.weights[[1]][ , 1]*lb.weights[[1]][ , 2])^0.5)*(lb.weighted.devs[[1]][ , 1]*lb.weighted.devs[[1]][ , 2]))
# 
# 
# tester <- matrix(cbind(lb.weighted.devs[[1]][ , 1], lb.weighted.devs[[1]][ , 2]), nrow = 252, ncol = 2)
# test.weight <- (lb.weights[[1]][ , 1]*lb.weights[[1]][ , 2])^0.5
# example[[1]]
# cov.wt(tester, wt = test.weight/252)
# sum((lb.weights[[l]][ , f]*lb.weights[[l]][ , s])^0.5)/251

#compare to unweighted
lb.example.uw <- lapply(lb.returns, cov)
hld.example.uw <- lapply(hld.returns, cov)

lb.example[[2]][1, 2]
lb.example.uw[[1]][1, 1]
#quite different

#note the number of covariance terms is (n^2-n)/2 + n. which is ((24^2)-24)/2 + 24 = 300
lb.cov.ts <- matrix(0, 116, 300)
lb.cov.ts.uw <- matrix(0, 116, 300)
hld.cov.ts <- matrix(0, 116, 300)
hld.cov.ts.uw <- matrix(0, 116, 300)
for(x in 1:300){
  lb.cov.ts[ , x] <- sapply(lb.example, "[[", x)
  lb.cov.ts.uw[ , x] <- sapply(lb.example.uw, "[[", x)
  hld.cov.ts[ , x] <- sapply(hld.example, "[[", x)
  hld.cov.ts.uw[ , x] <- sapply(hld.example.uw, "[[", x)
}

# the four matrices represent the lookback and holding covariances for the unweighted and weighted methods each row represents all 300 of the covariances for each period. each column contains all of the same covariances over time. 
# lb.single.cov.ts <- sapply(lb.example, "[[", 1)
# lb.single.cov.ts.uw <- sapply(lb.example.uw, "[[", 2)
# 
# hld.single.cov.ts <- sapply(hld.example, "[[", 2)
# hld.single.cov.ts.uw <- sapply(hld.example.uw, "[[", 2)

#find residuals between weighted and unweighted covaraince matrices
lb.cov.residual <- lb.cov.ts - lb.cov.ts.uw
hld.cov.residual <- hld.cov.ts - hld.cov.ts.uw

#to plot all columns in a matrix as individual variables
matplot(lb.cov.residual, type = "l")
matplot(hld.cov.residual, type = "l")
#the fact that the residuals are generally negative shows that the weighted covariances tend to be smaller.
#then test hld to hld prediction errors
hld.predicted.cov.residual <- hld.cov.ts[-116, ]-hld.cov.ts[-1, ]
hld.predicted.cov.residual.uw <- hld.cov.ts.uw[-116, ]-hld.cov.ts.uw[-1, ]
hld.prediction.improvement <- hld.predicted.cov.residual - hld.predicted.cov.residual.uw

matplot(hld.prediction.improvement, type = "l")

#and lb to hld prediction errors
predicted.cov.residual <- lb.cov.ts[-116, ]-hld.cov.ts[-1, ]
predicted.cov.residual.uw <- lb.cov.ts.uw[-116, ]-hld.cov.ts.uw[-1, ]
prediction.improvement <- predicted.cov.residual - predicted.cov.residual.uw
matplot(prediction.improvement, type = "l")

x.axis <- seq(as.Date("2008/2/1"), by = "month", length.out = 115)

matplot(x.axis, prediction.improvement, type = "l")
matplot(x.axis, prediction.improvement, type = "l")
matplot(x = as.Date(x.axis), y = prediction.improvement, type = "l")

matplot(prediction.improvement, type = "l")
mean(prediction.improvement)
#RMSD

RMSE.weighted <- (mean((lb.cov.ts - hld.cov.ts)^2))^0.5
RMSE.unweighted <- (mean((lb.cov.ts.uw - hld.cov.ts.uw)^2))^0.5
improvement <- (RMSE.weighted-RMSE.unweighted)/RMSE.unweighted

x.axis[115]

plot(hld.cov.ts.uw[ , 1])
lines(hld.cov.ts[ , 1])
sd(hld.single.cov.ts.uw)
sd(hld.single.cov.ts)
#mean(lb.returns[[1]][ , 1]*(1/lb.stderr[[1]][ , 1]))

