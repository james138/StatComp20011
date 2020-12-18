#' @title Bivariate normal dataset with missing.
#' @name bvn
#' @description A dataset simulated from "randomized monotone missingness process" and is used to carry out \code{my.em} and \code{my.boot}.
#' @importFrom stats complete.cases sd
NULL

#' @title Stopping criterion.
#' @description It is used to decide whether to stop EM algorithm or not.
#' @param iter iteration of EM algorithm (numeric)
#' @param max maximum iteration of EM algorithm (numeric)
#' @param tol convergence tolerance (numeric)
#' @param exist existing relative change
#' @return Return \code{TRUE} if the iteration gets the maximum or existing relative change reach the convergence tolerance else return \code{FALSE}.
#' @examples
#' \dontrun{
#' iter <- 10
#' max <- 100
#' tol <- 1e-5
#' eist <- 9.76e-6
#' criterion()
#' }
#' @export
criterion <- function(iter, max ,tol , exist) max(abs(exist))<tol|iter == max

#' @title My EM algorithm.
#' @description EM algorithm for bivariate normal data
#' @param x dataset 
#' @param max maximum iteration of EM algorithm 
#' @param tol convergence tolerance
#' @param pr print the results of each iteration or not
#' @return If \code{pr==TRUE}, print the results of each iteration, return the estimation of theta as \code{theta.em} and the results of each iteration in \code{record}. If \code{pr==FALSE}, just return the estimation of theta as \code{theta.em} and the results of each iteration in \code{record}.
#' @examples
#' \dontrun{
#' data(bvn)
#' attach(bvn)
#' pr <- 1
#' my.em(bvn,max,tol,pr)
#' }
#' @export
my.em <- function(x,max=100,tol=1e-5,pr){
  r <- 1-is.na(x)
  #  missing index
  N <- nrow(x)
  iter <- exist <- 1
  xc <- x[complete.cases(x), ]
  #  complete case
  Nc <- nrow(xc)
  x[is.na(x)] <- -999
  #  Change NAs to nonsensical value
  mu0 <- apply(xc,2,mean)
  sigma0 <- t(sweep(data.matrix(xc),2,mu0))%*%sweep(data.matrix(xc),2,mu0)/(Nc-1)
  #  complete case sample mean and sample covariance as initial
  theta <- c(mu0,sigma0[1,1],sigma0[1,2],sigma0[2,2])
  record <- c(0,theta)
  #  record
  while (! criterion(iter,max,tol,exist)){
    mu <- theta[1:2]
    sigma <- matrix(c(theta[3],theta[4],theta[4],theta[5]),2,2)
    ey1 <- (r[,1]==1)*x[,1] + (r[,1]==0)*(mu[1]+ sigma[1,2]*(x[,2]-mu[2])/sigma[2,2])
    ey2 <- (r[,2]==1)*x[,2] + (r[,2]==0)*(mu[2]+ sigma[1,2]*(x[,1]-mu[1])/sigma[1,1])
    ey1.2 <- (r[,1]==1)*x[,1]^2 + (r[,1]==0)*(sigma[1,1]-sigma[1,2]^2/sigma[2,2]+ ey1^2)
    ey2.2 <- (r[,2]==1)*x[,2]^2 + (r[,2]==0)*(sigma[2,2]-sigma[1,2]^2/sigma[1,1]+ ey2^2)
    ey1y2 <- (r[,1]*r[,2]*x[,1]*x[,2] + r[,1]*(1-r[,2])*x[,1]*ey2
              + (1-r[,1])*r[,2]*ey1*x[,2])
    #  E-Step
    munew <- c(mean(ey1),mean(ey2))
    s11 <- mean(ey1.2)-munew[1]^2
    s22 <- mean(ey2.2)-munew[2]^2
    s12 <- mean(ey1y2)-munew[1]*munew[2]
    #  M-Step
    thetanew <- c(munew,s11,s12,s22)
    if (pr==1) print(c(iter,thetanew))
    #  print the results of each iteration
    record <- rbind(record,c(iter,thetanew))
    #  record
    exist <- (thetanew-theta)/theta
    #  relative change
    theta <- thetanew
    iter <- iter+1
  }
  results <- list(theta.em=theta,record=record)
  return(results)
}

#' @title My bootstrap algorithm.
#' @description Bootstrap algorithm to compute the variance in EM algorithm
#' @param x dataset 
#' @param max maximum iteration of EM algorithm 
#' @param tol convergence tolerance
#' @param B number of bootstrap
#' @return Return the Bootstrap means and the Bootstrap standard errors.
#' @examples
#' \dontrun{
#' data(bvn)
#' attach(bvn)
#' B <- 250
#' my.boot(x=bvn,B=B)
#' }
#' @export
my.boot <- function(x,max=100,tol=1e-5,B){
  N <- nrow(x)
  theta.boot <- NULL
  for (b in 1:B){
    i <- sample(N,N,replace=TRUE)
    #  Sample with replacement
    xc.b <- x[i,]
    #  Call EM
    em.b <- my.em(xc.b,max,tol,0)
    theta.b <- em.b$theta.em
    theta.boot <- rbind(theta.boot,theta.b)
    
  }
  boot.mean <- apply(theta.boot,2,mean)
  boot.se <- apply(theta.boot,2,sd)
  return(list(boot.mean=boot.mean,boot.se=boot.se))
}