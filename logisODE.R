## ----include=FALSE----------------------------------------------------------------------------------------
source("/Users/FJRubio/Dropbox/ODESurv/ODESurv/Codes/routines/routines.R")


## ----eval=FALSE-------------------------------------------------------------------------------------------
## source("routines.R")


## ---------------------------------------------------------------------------------------------------------
# Parameter values
l <- 1
k <- 10
h0s <- c(1/10,2,1)*k

haz1 <- Vectorize(function(t) hlogisode(t, l, k, h0s[1]))
curve(haz1, 0, 10, ylim = c(0,20), xlab = "time", ylab = "Hazard", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 1, n = 250)

haz2 <- Vectorize(function(t) hlogisode(t, l, k, h0s[2]))
curve(haz2, 0, 10, ylim = c(0,20), xlab = "time", ylab = "Hazard", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE, n = 250)

haz3 <- Vectorize(function(t) hlogisode(t, l, k, h0s[3]))
curve(haz3, 0, 10, ylim = c(0,20), xlab = "time", ylab = "Hazard", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 3, add = TRUE, n = 250)

legend("topright", legend = c(expression(paste(h[0], " = ", kappa/10)),
                              expression(paste(h[0], " = ", 2*kappa)),
                              expression(paste(h[0], " = ", kappa))),
       lty = c(1,2,3), lwd = c(2,2,2))



## ---------------------------------------------------------------------------------------------------------
# Parameter values
l <- 1
k <- 10
h0s <- c(1/10,2,1)*k

surv1 <- Vectorize(function(t) exp(-chlogisode(t, l, k, h0s[1])))
curve(surv1, 0, 2, ylim = c(0,1), xlab = "time", ylab = "Survival", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 1, n = 250)

surv2 <- Vectorize(function(t) exp(-chlogisode(t, l, k, h0s[2])))
curve(surv2, 0, 10, ylim = c(0,1), xlab = "time", ylab = "Survival", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 2, add = TRUE, n = 250)

surv3 <- Vectorize(function(t) exp(-chlogisode(t, l, k, h0s[3])))
curve(surv3, 0, 10, ylim = c(0,1), xlab = "time", ylab = "Survival", 
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, lty = 3, add = TRUE, n = 250)

legend("topright", legend = c(expression(paste(h[0], " = ", kappa/10)),
                              expression(paste(h[0], " = ", 2*kappa)),
                              expression(paste(h[0], " = ", kappa))),
       lty = c(1,2,3), lwd = c(2,2,2))



## ---------------------------------------------------------------------------------------------------------
# Parameter values
l <- 1
k <- 10
h0s <- c(1/10,2,1)*k
n <- 10000

set.seed(1234)
dat1 <- rlogisode(n, l, k, h0s[1])
den1 <- Vectorize(function(t) dlogisode(t, l, k, h0s[1]))
hist(dat1, breaks = 50, probability = TRUE, xlab = "time", ylab = "Density", 
      cex.axis = 1.5, cex.lab = 1.5, main = "", ylim = c(0,h0s[1]) )
curve(den1, 0, 3, lwd = 2, lty = 1, n = 250, add = TRUE)
box()

set.seed(1234)
dat2 <- rlogisode(n, l, k, h0s[2])
den2 <- Vectorize(function(t) dlogisode(t, l, k, h0s[2]))
hist(dat2, breaks = 50, probability = TRUE, xlab = "time", ylab = "Density", 
      cex.axis = 1.5, cex.lab = 1.5, main = "", ylim = c(0,h0s[2]) )
curve(den2, 0, 0.5, lwd = 2, lty = 1, n = 250, add = TRUE)
box()

set.seed(1234)
dat3 <- rlogisode(n, l, k, h0s[3])
den3 <- Vectorize(function(t) dlogisode(t, l, k, h0s[3]))
hist(dat3, breaks = 50, probability = TRUE, xlab = "time", ylab = "Density", 
      cex.axis = 1.5, cex.lab = 1.5, main = "", ylim = c(0,h0s[3]) )
curve(den3, 0, 1, lwd = 2, lty = 1, n = 250, add = TRUE)
box()




