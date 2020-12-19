## ----eval=FALSE---------------------------------------------------------------
#  set.seed(4)
#  N <- 1000
#  theta <- c(5,8,1,0.5,1)
#  mu <- theta[1:2]
#  sigma <- matrix(c(theta[3],theta[4],theta[4],theta[5]),2,2)
#  prob <- 0.6
#  #  prob observe y1
#  psi1 <- c(-1,0.3)
#  psi2 <- c(-1,0.1)
#  x <- NULL
#  r <- NULL
#  
#  for (i in 1:N){
#    obs <- rbinom(1,1,prob)
#    #  observation indicator
#    y <- mu + t(chol(sigma))%*%rnorm(2,0,1)
#    #  Generate full data
#    ep1 <- exp(psi1[1]+psi1[2]*y[1])
#    p1 <- ep1/(1+ep1)
#    ep2 <- exp(psi2[1]+psi2[2]*y[2])
#    p2 <- ep2/(1+ep2)
#    r1 <- obs + (1-obs)*rbinom(1,1,p2)
#    r2 <- (1-obs) + obs*rbinom(1,1,p1)
#    x <- rbind(x,t(y))
#    r <- rbind(r,c(r1,r2))
#  }
#  x <- round(x,3)
#  #  Round to 3 decimal to make it like real data
#  x[r==0] <- NA
#  #  Replacing missing with NA
#  bvn <- x
#  save(bvn,file="bvn.rda")

## ----eval=FALSE---------------------------------------------------------------
#  my.em <- function(x,max=100,tol=1e-5,pr){
#    r <- 1-is.na(x)
#    #  missing index
#    N <- nrow(x)
#    iter <- exist <- 1
#    xc <- x[complete.cases(x), ]
#    #  complete case
#    Nc <- nrow(xc)
#    x[is.na(x)] <- -999
#    #  Change NAs to nonsensical value
#    mu0 <- apply(xc,2,mean)
#    sigma0 <- t(sweep(data.matrix(xc),2,mu0))%*%sweep(data.matrix(xc),2,mu0)/(Nc-1)
#    #  complete case sample mean and sample covariance as initial
#    theta <- c(mu0,sigma0[1,1],sigma0[1,2],sigma0[2,2])
#    record <- c(0,theta)
#    #  record
#    while (! criterion(iter,max,tol,exist)){
#      mu <- theta[1:2]
#      sigma <- matrix(c(theta[3],theta[4],theta[4],theta[5]),2,2)
#      ey1 <- (r[,1]==1)*x[,1] + (r[,1]==0)*(mu[1]+ sigma[1,2]*(x[,2]-mu[2])/sigma[2,2])
#      ey2 <- (r[,2]==1)*x[,2] + (r[,2]==0)*(mu[2]+ sigma[1,2]*(x[,1]-mu[1])/sigma[1,1])
#      ey1.2 <- (r[,1]==1)*x[,1]^2 + (r[,1]==0)*(sigma[1,1]-sigma[1,2]^2/sigma[2,2]+ ey1^2)
#      ey2.2 <- (r[,2]==1)*x[,2]^2 + (r[,2]==0)*(sigma[2,2]-sigma[1,2]^2/sigma[1,1]+ ey2^2)
#      ey1y2 <- (r[,1]*r[,2]*x[,1]*x[,2] + r[,1]*(1-r[,2])*x[,1]*ey2
#                + (1-r[,1])*r[,2]*ey1*x[,2])
#      #  E-Step
#      munew <- c(mean(ey1),mean(ey2))
#      s11 <- mean(ey1.2)-munew[1]^2
#      s22 <- mean(ey2.2)-munew[2]^2
#      s12 <- mean(ey1y2)-munew[1]*munew[2]
#      #  M-Step
#      thetanew <- c(munew,s11,s12,s22)
#      if (pr==1) print(c(iter,thetanew))
#      #  print the results of each iteration
#      record <- rbind(record,c(iter,thetanew))
#      #  record
#      exist <- (thetanew-theta)/theta
#      #  relative change
#      theta <- thetanew
#      iter <- iter+1
#    }
#    results <- list(theta.em=theta,record=record)
#    return(results)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  my.boot <- function(x,max=100,tol=1e-5,B){
#    N <- nrow(x)
#    theta.boot <- NULL
#    for (b in 1:B){
#      i <- sample(N,N,replace=TRUE)
#      #  Sample with replacement
#      xc.b <- x[i,]
#      #  Call EM
#      em.b <- my.em(xc.b,max,tol,0)
#      theta.b <- em.b$theta.em
#      theta.boot <- rbind(theta.boot,theta.b)
#  
#    }
#    boot.mean <- apply(theta.boot,2,mean)
#    boot.se <- apply(theta.boot,2,sd)
#    return(list(boot.mean=boot.mean,boot.se=boot.se))
#  }

## -----------------------------------------------------------------------------
set.seed(4) 
library("StatComp20011")
data(bvn)
pr <- 1
max <- 100
tol <- 1e-5
a <- my.em(bvn,max,tol,pr)$theta.em
b <- my.boot(bvn,max,tol,250)$boot.se
df <- data.frame(a,b)
names(df) <- c("mean", "standard error")
rownames(df) <- c("mu1","mu2","sigma[1,1]","sigma[1,2]","sigma[2,2]")
knitr::kable(df)

## -----------------------------------------------------------------------------
library(norm)
prelim.em <- prelim.norm(bvn)
#  set up the data
theta <- em.norm(prelim.em,showits=FALSE)
theta <- getparam.norm(prelim.em,theta,corr=TRUE)
theta$mu <- c(5.042471,7.952805)
theta$sdv <- sqrt(c(1.026614,1.033791))
theta$r[1,2] <- 0.494137/sqrt(1.026614*1.033791)
#  Use the sample complete case standard deviations
#  Use the sample complete case correlation
theta <- makeparam.norm(prelim.em,theta)
theta.final <-
  em.norm(prelim.em,start=theta,showits=TRUE,maxits=100,criterion=1e-5)        
theta.final <- getparam.norm(prelim.em,theta.final,corr=TRUE)
mu.final <- theta.final$mu
sd.final <- theta.final$sdv
var.final <- sd.final^2
covar.final <- theta.final$r[1,2]*sd.final[1]*sd.final[2]
round(c(mu.final,var.final[1],covar.final,var.final[2]),6)

