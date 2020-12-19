## -----------------------------------------------------------------------------
pairs(USArrests, pch=23, bg='orange',                                          
      cex.labels=1.5) 

## ---- results = 'asis'--------------------------------------------------------
knitr::kable(head(USArrests))

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 10000
u <- runif(n)
x <- 2/sqrt(u)
hist(x, prob=TRUE, main = 'Pareto(2,2)', xlim = c(0,50), breaks=500)
y <- seq(2, 50, .1)
lines(y, 8/(y^3),lwd = 1.5, col = "red")

## -----------------------------------------------------------------------------
repan <- function(n){
   u <- runif(3*n, min = -1, max = 1)
for(i in 1:n){
   if(abs(u[3*i]) >= abs(u[3*i-1]) && abs(u[3*i]) >= abs(u[3*i-2])) x[i] <- u[3*i-1] else x[i] <- u[3*i]
}
   return(x)
}
set.seed(12345)
x <- repan(10000)
hist(x, prob=TRUE, main = 'Epanechnikov kernel', breaks = 50)

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 1000
u <- runif(n)
x <- 2*(1-u^(1/4))/u^(1/4)
hist(x, prob=TRUE, main = 'Pareto', xlim = c(0,10), breaks=100)
y <- seq(0, 10, .1)
lines(y, 64/(2+y)^5,lwd = 1.5, col = "red")

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 1e4
x <- runif(n, min = 0, max = pi/3)
theta.hat <- mean(sin(x))*pi/3
print(theta.hat)
print(1/2)

## -----------------------------------------------------------------------------
set.seed(12345)
m <- 100000
U <- runif(m)
T1 <- exp(U)
T2 <- 1/2*exp(U)+1/2*exp(1-U)
mean(T1)
mean(T2)
(var(T1)-var(T2))/var(T1)

## -----------------------------------------------------------------------------
n <- 10000
set.seed(12345)
theta_hat <- numeric(2)
se <- numeric(2)
g <- function(x){
  x^2*exp(-x^2/2)/(sqrt(2*pi))* (x > 1)
}

u <- runif(n) #f1,inverse transform method
x <- 1-log(u)
fg <- g(x)/exp(1-x)
theta_hat[1] <- mean(fg)
se[1] <- sd(fg)

u <- runif(n) #f2,inverse transform method
x <- tan(pi/4*(u+1))
fg <- g(x)/(4/((1+x^2)*pi))
theta_hat[2] <- mean(fg)
se[2] <- sd(fg)

rbind(theta_hat, se)

## -----------------------------------------------------------------------------
x <- seq(1, 10, 0.01)
w <- 2
f1 <- exp(1-x)
f2 <- 4/((1+x^2)*pi)
g <- x^2/(sqrt(2*pi)) * exp(-x^2/2)

plot(x, g, type = "l", ylim=c(0,1.5), lwd = w)
lines(x, f1, lty = 2, lwd = w)
lines(x, f2, lty = 3, lwd = w) 
legend("topright", legend = c("g", 1:2), lty = 1:3, lwd = w, inset = 0.02)

plot(x, g/f1, type = "l", ylim = c(0,3), lwd = w, lty = 2)
lines(x, g/f2, lty = 3, lwd = w)
legend("topright", legend = c(1:2), lty = 2:3, lwd = w, inset = 0.02)

## ----5.15---------------------------------------------------------------------
n <- 10000
k <- 5
r <- n/k
N <- 50
set.seed(12345)

T2 <- numeric(k)
est <- matrix(0, N, 2)
g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)
for (i in 1:N){
  u <- runif(n)
  x <- - log(1 - u * (1 - exp(-1)))
  fg <- g(x) / (exp(-x) / (1 - exp(-1)))
  est[i, 1] <- mean(fg)
  for(j in 1:k){ 
    u <- runif(n/k)
    t <- u * exp(-j/5)+(1 - u) * exp(-(j-1)/5)
    x <- - log(t)
    fg <- 5*g(x) / (exp(-x) / (exp(-(j-1)/5) - exp(-j/5)))
    T2[j] <- mean(fg)
  }
  est[i, 2] <- mean(T2)
} 
theta_hat <- apply(est,2,mean)
se <- apply(est,2,sd)
rbind(theta_hat, se)

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
set.seed(1)
UCL <- replicate(1000, expr = {
  x <- rlnorm(n, 0, 2)
  sqrt(n)*mean(log(x))/sqrt(var(log(x)))/ qt(1-alpha/2,n-1)
  })
mean(abs(UCL) < 1)

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
UCL <- replicate(1000, expr = {
  x <- rchisq(n, df = 2)
  sqrt(n)*(mean(x)-2)/sqrt(var(x))/ qt(1-alpha/2,n-1)
  })
mean(abs(UCL) < 1)

UCL <- replicate(1000, expr = {
  x <- rchisq(n, df = 2)
  (n-1)*var(x)/qchisq(alpha,n-1)
  })
mean(UCL > 4)

## -----------------------------------------------------------------------------
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

set.seed(12345)
a <- .1              #significance level
n <- 30 
m <- 3000
alpha <- seq(.1, 20.1, 1)
N <- length(alpha)
pwr_Beta <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1 - a / 2, 0, sqrt(6 * (n - 2) / ((n + 1) * (n + 3))))

for (i in 1:N) {     #for each alpha 
  sktests_Beta <- numeric(m)
  for (j in 1:m) {   #for each replicate
    x_Beta <- rbeta(n, alpha[i], alpha[i])
    sktests_Beta[j] <- as.integer(abs(sk(x_Beta)) >= cv)
  }
  pwr_Beta[i] <- mean(sktests_Beta)
}
#plot power vs alpha 
plot(alpha, pwr_Beta, xlab = bquote(alpha), ylab = "power",
      type = "b", ylim = c(0, 0.1))
abline(h = .1, lty = 3)                         #significance level
se_Beta <- sqrt(pwr_Beta * (1 - pwr_Beta) / m)  #standard error
#plot confidence interval
lines(alpha, pwr_Beta + se_Beta, lty = 3)
lines(alpha, pwr_Beta - se_Beta, lty = 3)

## -----------------------------------------------------------------------------
nu <- seq(.1, 20.1, 1)
N2 <- length(nu)
pwr_t <- numeric(N2)

for (i in 1:N2) {        #for each nu 
  skt_t <- numeric(m)
  for (j in 1:m) {       #for each replicate
    x_t <- rt(n, nu[i])
    skt_t[j] <- as.integer(abs(sk(x_t)) >= cv)
  }
  pwr_t[i] <- mean(skt_t)   #power
}

#plot power vs nu
plot(nu, pwr_t, xlab = bquote(nu), ylab = "power", 
     type = "b", ylim = c(0, 1))
abline(h = .1, lty = 3)                #significance level
se_t <- sqrt(pwr_t * (1 - pwr_t) / m)  #standard error
#plot confidence interval
lines(nu, pwr_t + se_t, lty = 3)
lines(nu, pwr_t - se_t, lty = 3)

## -----------------------------------------------------------------------------
sigma1 <- 1
sigma2 <- 1.5
m <- 10000

count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}

power_count <- function(m,n,sigma1,sigma2){
  tests <- replicate(m, expr = {
  x <- rnorm(n, 0, sigma1)
  y <- rnorm(n, 0, sigma2)
  count5test(x, y)
  } )
  return(mean(tests))
}

power_f <- function(m,n,sigma1,sigma2){
  tests <- replicate(m, expr = {
  x <- rnorm(n, 0, sigma1)
  y <- rnorm(n, 0, sigma2)
  return(as.integer(var.test(x,y,conf.level=0.945)$p.value<0.055))
  } )
  return(mean(tests))
}

n <- c(20, 100, 500)
k = length(n)
power_1 = numeric(k)
power_2 = numeric(k)

for(i in 1:k){
  power_1[i] <- power_count(m,n[i],sigma1,sigma2)
  power_2[i] <- power_f(m,n[i],sigma1,sigma2)
}
df<-data.frame(power_1, power_2)
names(df) <- c("Count Five test","F-test")
row.names(df) <- c("n=20", "n=100", "n=500")
knitr::kable(df)

## -----------------------------------------------------------------------------
set.seed(1234)
m <- 1000
library(MASS)
#compute the skewness function
sk2 <- function(X) {
        k <- nrow(X)
        Y <- matrix(rep(colMeans(X),k), nrow=k, byrow=TRUE)     #colmean 
        Sigma <- t(X-Y)%*%(X-Y)/k                               #covariance         
        b <- sum(((X-Y)%*%solve(Sigma)%*%t(X-Y))^3)
        return(b/k^2)
}
sktests=function(n,cv,d){
  tests <- replicate(m, expr = {                #tests
          x=mvrnorm(n,rep(0,d),diag(d))         #generation
          return(as.integer(sk2(x)>=cv)) 
          } )
mean(tests)
}
f <- function(d){
         n <- c(10,20,30,50,100,500)
         #critical value for the skewness test
         cv <- qchisq(0.95, d*(d+1)*(d+2)/6)*6/n
         p.reject <- numeric(length(n))
        for(i in 1:length(n)){
             p.reject[i] <- sktests(n[i],cv[i],d)
        }
 p.reject          #percent of reject
}
df2 <- data.frame(f(1), f(2))
names(df2) <- c("d=1","d=2")
row.names(df2) <- c("n=10", "n=20", "n=30", "n=50", "n=100", "n=500")
knitr::kable(df2)

## -----------------------------------------------------------------------------
alpha=0.1            #significance level
g <- function(d,n){
    epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
    N <- length(epsilon)
    pwr <- numeric(N)
    #critical value for the skewness test
    cv <- qchisq(.95, d*(d+1)*(d+2)/6)*6/n  
    for (j in 1:N) {              #for each epsilon
        e <- epsilon[j]
        sktests <- numeric(m)
        for (i in 1:m) {           #for each replicate
                x <- matrix(0, nrow = n, ncol = d)
                for(k in 1:n){
                   sigma <- sample(c(1, 10), replace = TRUE,
                                size = 1, prob = c(1-e, e))        
                   x[k,] <- mvrnorm(1, rep(0, d), (sigma^2)*diag(d))     
                }
              sktests[i] <- as.integer(sk2(x)>=cv)  
        }
        pwr[j] <- mean(sktests)      #power
    }  
     #plot power vs epsilon
     plot(epsilon, pwr,xlab = bquote(epsilon), ylab = "power",
          type = "b", ylim = c(0,1))
     abline(h = .1, lty = 3)       #significance level
     se <- sqrt(pwr * (1-pwr) / m) #standard errors
     #plot confidence interval
     lines(epsilon, pwr+se, lty = 3)
     lines(epsilon, pwr-se, lty = 3) 
}
g(1,30) #d=1,n=30
g(2,30) #d=2,n=30

## -----------------------------------------------------------------------------
data(law, package = "bootstrap") #import data
n <- nrow(law)
x <- law$LSAT
y <- law$GPA
theta_hat <- cor(x,y)
#compute the leave-one-out estimates
theta_jack <- rep(0, n)
for (i in 1:n)
  theta_jack[i] <- cor(x[-i],y[-i])
bias <- (n - 1) * (mean(theta_jack) - theta_hat)
se <- (n - 1) * sqrt(var(theta_jack) / n)
df <- data.frame(theta_hat, bias, se) 
knitr::kable(df)


## -----------------------------------------------------------------------------
set.seed(12345)
options(warn = -1)
library(boot)
B <- 2000                                                         #bootstrap number
stat <- function(dat, index){mean(dat[index])}                    #function
bootstrap_result <- boot(data = aircondit$hours, stat, R = B)     #boot result
boot.ci(bootstrap_result, conf = 0.95, type = c("norm", "basic", "perc", "bca")) 

## -----------------------------------------------------------------------------
data(scor, package = "bootstrap") #import data
n <- nrow(scor)
lambda_hat <- eigen((n-1)*cov(scor)/n)$values
theta_hat <- lambda_hat[1]/sum(lambda_hat)
theta_j <- rep(0, n) 
for (i in 1:n) { 
  x <- scor [-i,]
  m <- n-1
  lambda <- eigen((m-1)*cov(x)/m)$values 
  theta_jack[i] <- lambda[1] / sum(lambda)   
  }
bias <- (n - 1) * (mean(theta_jack) - theta_hat)
se <- (n - 1) * sqrt(var(theta_jack) / n)
df2 <- data.frame(theta_hat, bias, se)
knitr::kable(df2)


## -----------------------------------------------------------------------------
data(ironslag, package = "DAAG")
magnetic <- ironslag$magnetic
chemical <- ironslag$chemical
n <- nrow(ironslag) 
e1 <- e2 <- e3 <- e4 <- matrix(0,n,n)
# for n-fold cross validation
# fit models on leave-two-out samples
for (k in 1:n) {
  for(l in 1:n){
    y <- magnetic[c(-k,-l)]
    x <- chemical[c(-k,-l)]
    
    J1 <- lm(y ~ x)
    yhat11 <- J1$coef[1] + J1$coef[2] * chemical[k]
    yhat12 <- J1$coef[1] + J1$coef[2] * chemical[l]
    e1[k,l] <- ((magnetic[k] - yhat11)^2 + (magnetic[l] - yhat12)^2)/2       #mean error
    
    J2 <- lm(y ~ x + I(x^2))
    yhat21 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
    yhat22 <- J2$coef[1] + J2$coef[2] * chemical[l] + J2$coef[3] * chemical[l]^2
    e2[k,l] <- ((magnetic[k] - yhat21)^2 + (magnetic[l] - yhat22)^2)/2      #mean error
     
    J3 <- lm(log(y) ~ x)
    logyhat31 <- J3$coef[1] + J3$coef[2] * chemical[k]
    logyhat32 <- J3$coef[1] + J3$coef[2] * chemical[l]
    yhat31 <- exp(logyhat31)
    yhat32 <- exp(logyhat32)
    e3[k,l] <- ((magnetic[k] - yhat31)^2 + (magnetic[l] - yhat32)^2)/2       #mean error
    
    J4 <- lm(log(y) ~ log(x))
    logyhat41 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    logyhat42 <- J4$coef[1] + J4$coef[2] * log(chemical[l])
    yhat41 <- exp(logyhat41)
    yhat42 <- exp(logyhat42)
    e4[k,l] <- ((magnetic[k] - yhat41)^2 + (magnetic[l] - yhat42)^2)/2       #mean error
  }
}
df3 <- data.frame(mean(e1), mean(e2), mean(e3), mean(e4))
names(df3)=c("mean1", "mean2", "mean3", "mean4")
knitr::kable(df3)

## -----------------------------------------------------------------------------
L2 <- lm(magnetic ~ chemical + I(chemical^2))
L2
plot(L2$fit, L2$res)            #residuals vs fitted values
abline(0, 0)                    #reference line 
qqnorm(L2$res)                  #qq plot
qqline(L2$res)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))       
  return(as.integer(max(c(outx, outy)) > 5))        # return 1 (reject) or 0 (do not reject H0)
}

set.seed(12345)                                     # seed
n1 <- 20 
n2 <- 30
mu1 <- mu2 <-0
sigma1 <- sigma2 <-1
x <- rnorm(n1, mu1, sigma1)                         # sample generation
y <- rnorm(n2, mu2, sigma2)
z <- c(x, y)
R <- 1e4                                            # cycle times
n <- (n1 + n2) / 2
m <- c(1 : (n1 + n2))

res <- numeric(R)                                   # storage 
for (i in 1:R) {                                    # loop
  k <- sample(m, size = n, replace = FALSE)         # split
  x <- z[k]
  y <- z[-k]
  res[i] <- count5test(x, y)
}
mean(res)

## -----------------------------------------------------------------------------
options(warn = -1)
library(RANN)
library(boot)
library(energy)
library(Ball)


m <- 50; k<-3; p<-2; 
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)

Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,mean = 1,sd = 0.6),ncol=p);
  y <- matrix(rnorm(n2*p,mean = 1,sd = 0.85),ncol=p);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
print(pow)
barplot(pow,col = "lightgreen",main = "Unequal variances and unequal expectations",
        names.arg=c("NN","energy","ball"))

## -----------------------------------------------------------------------------
m <- 50; k<-3; p<-2; 
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,mean = 0.4,sd = 0.6),ncol=p);
  y <- matrix(rnorm(n2*p,mean = 0.5,sd = 0.85),ncol=p);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
print(pow)
barplot(pow,col = "lightblue",main = "Unequal variances and equal expectations",
        names.arg=c("NN","energy","ball"))

## -----------------------------------------------------------------------------
m <- 50; k<-3; p<-2; 
n1 <- n2 <- 20; R<-999; n <- n1+n2; N = c(n1,n2)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  # t distribution with 1 df (heavy-tailed distribution)
  x <- matrix(rt(n1*p,df = 1),ncol=p); 
  #bimodel distribution (mixture of two normal distributions)
  y <- cbind(rnorm(n2,mean = 0.4),rnorm(n2,mean = 0.5));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
print(pow)
barplot(pow,col = "pink",main = "Non-normal distributions",
        names.arg=c("NN","energy","ball"))

## -----------------------------------------------------------------------------
m <- 50; k<-3; p<-2; 
n1 <- 10;n2 <- 100;R<-999; n <- n1+n2; N = c(n1,n2)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- c(rnorm(n1,mean = 1,sd = 1)); # n1 = 10
  y <- c(rnorm(n2,mean = 2,sd = 2)); # n2 = 100
  z <- c(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
print(pow)
barplot(pow,col = "lightgreen",main = "Unbalanced sample",names.arg=c("NN","energy","ball"))

## -----------------------------------------------------------------------------
options(warn = -1)
library(GeneralizedHyperbolic)
Laplace.Metropolis <- function(sigma, N){
  x <- numeric(N)
  x[1] <- rnorm(1,0,sigma)
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if(u[i] <= dskewlap(y)/dskewlap(x[i-1]))
      x[i] <- y
    else{
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x = x, k = k))
}
N <- 1e4
sigma <- c(1, 4, 9, 16)
set.seed(1234)
laplace1<-Laplace.Metropolis(sigma[1], N)
laplace2<-Laplace.Metropolis(sigma[2], N)
laplace3<-Laplace.Metropolis(sigma[3], N)
laplace4<-Laplace.Metropolis(sigma[4], N)
df <- data.frame(c(1-laplace1$k/N,1-laplace2$k/N,1-laplace3$k/N,1-laplace4$k/N))
names(df) <- c("acceptance rate")
rownames(df) <- c("sigma=1", "sigma=4", "sigma=9", "sigma=16")
knitr::kable(df)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}
Laplace.Metropolis <- function(sigma, N, X1){
  x <- numeric(N)
  x[1] <- X1
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1,x[i-1],sigma)
    if(u[i] <= dskewlap(y)/dskewlap(x[i-1]))
      x[i] <- y
    else{
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(x)
}
k <- 4 #number of chains to generate
n <- 15000
b <- 1000
sigma <- 1
x0 <- c(-10, -5, 5, 10)
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- Laplace.Metropolis(sigma, n, x0[i])

psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
#plot psi for the four chains
for (i in 1:k)
plot(psi[i, (b+1):n], type="l",
xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0,n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
u <- c(seq(4,25),100,500,1000)
b <- numeric(length(u))
c <- numeric(length(u))
s1 <- function(a){   #S_k-1
  1-pt(sqrt(a^2*(k-1)/(k-a^2)),k-1)
}
s2 <- function(a){   #S_k
  1-pt(sqrt(a^2*k/(k+1-a^2)),k)
}
s <- function(a){   #s1-s2
  1-pt(sqrt(a^2*(k-1)/(k-a^2)),k-1)-(1-pt(sqrt(a^2*k/(k+1-a^2)),k))
}
for(i in 1:length(u)){   
  k <- u[i]
  b[i] <- uniroot(s,c(0.01,2))$root    #get root
  c[i] <- s1(b[i])
}

df2 <- data.frame(u,b,c)
names(df2) <- c("k", "abscissa", "Y-axis")
knitr::kable(df2)

## ----echo=FALSE---------------------------------------------------------------
options(warn=-1)
dat <- rbind(Genotype=c('AA','BB','OO','AO','BO','AB','Sum'),
             Frequency=c('p^2','q^2','r^2','2pr','2qr','2pq',1),
             Count=c('nAA','nBB','nOO','nAO','nBO','nAB','n'))
knitr::kable(dat)

## -----------------------------------------------------------------------------
p0 <- 0.5; q0 <- 0.3                                   # initial value
nA <- 444; nB <- 132; nOO <- 361; nAB <- 63            # data
tol <- .Machine$double.eps^0.5                         # tolerance
N <- 10000                                                            
p <- p0; q <- q0
logl <- numeric(N)                                     # log-likelihood                       
k <- 1 
L.old <- c(p0, q0)
for (i in 1:N) {
  r <- 1-p-q
  a <- p/(2-p-2*q)
  b <- q/(2-2*p-q)
  Hp <- nA*(1+a)+nAB
  Hq <- nB*(1+b)+nAB
  Hr <- nA*(1-a)+nB*(1-b)+2*nOO
  H2 <- nA*(1-a)+nB*(1-b)+nAB
  logl[i] <- Hp*log(p)+Hq*log(q)+Hr*log(r)+H2*log(2)    # log-likelihood
  p.old <- p; q.old <- q                
  p <- Hp/(Hp+Hq+Hr); q <- Hq/(Hp+Hq+Hr)                # update
  L <- c(p, q)
  if (sum(abs(L-L.old)/L.old)<tol) break
  L.old <- L
  k <- k+1
}  
p.hat<-p; q.hat<-q
logl<-logl[1:k]
df <- data.frame(c(p.hat,q.hat))
names(df) <- "value"
rownames(df) <- c("p.hat", "q.hat")
knitr::kable(df)
plot(c(1:length(logl)), logl, type = "b", xlab = "iteration", ylab = " log-likelihood")

## -----------------------------------------------------------------------------

attach(mtcars)

formulas = list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
#1 for loops
f3 = vector("list", length(formulas))
for (i in seq_along(formulas)){
  f3[[i]] = lm(formulas[[i]], data = mtcars)
}
f3
#2 lapply
la3 = lapply(formulas, function(x) lm(formula = x, data = mtcars))
la3


## -----------------------------------------------------------------------------
trials <- replicate(
  100, 
  t.test(rpois(10, 10), rpois(7, 10)), 
  simplify = FALSE 
) 

## -----------------------------------------------------------------------------
# anonymous
sapply(trials, function(x) x$p.value)
# without anonymous function
sapply(trials, '[[', "p.value")

## -----------------------------------------------------------------------------
Mvap <- function(x, fun, value){
  out <- Map(function(y) vapply(y, fun, value), x)
  unlist(out)
}
options(warn = -1)
testlist <- list(cars, faithful, iris)
Mvap(testlist, mean, numeric(1))

## -----------------------------------------------------------------------------
library(Rcpp)
library(GeneralizedHyperbolic)
library(microbenchmark)
f <- function(x) 0.5*exp(-abs(x))
rwMetropolosR <- function(sigma, x0, N){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= f(y)/f(x[i-1]))
      x[i] <- y 
    else {
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x = x,k = k))
}
cppFunction('NumericVector rwMetropolosC(double sigma, double x0, int N){
  NumericVector x(N+1);
  double k=0;
  x[0]=x0;
  double y=0;
  NumericVector u(N);
  for (int i=0;i<N;i++){
    u[i]=runif(1)[0];
  }
  for (int j=1;j<N-1;j++){
    y=rnorm(1,x[j-1],sigma)[0];
    if (u[j]<=0.5*exp(-abs(y))/(0.5*exp(-abs(x[j-1]))))
      x[j]=y;
    else {
      x[j]=x[j-1];
      k++;
    }
  }
  x[N]=k;
  return x;
}')
N <- 2e3
set.seed(1234)
sigma <- c(.05, .5, 2, 16)
x0 <- 25
#rwMetropolosR() returns a list with x and k, with the sample vector and the rejection number. 
rwR1 <- rwMetropolosR(sigma[1],x0,N)
rwR2 <- rwMetropolosR(sigma[2],x0,N)
rwR3 <- rwMetropolosR(sigma[3],x0,N)
rwR4 <- rwMetropolosR(sigma[4],x0,N)
#rwMetropolosC() returns a vector x, with the first N the sample and the last the rejection number.
rwC1 <- rwMetropolosC(sigma[1],x0,N)
rwC2 <- rwMetropolosC(sigma[2],x0,N)
rwC3 <- rwMetropolosC(sigma[3],x0,N)
rwC4 <- rwMetropolosC(sigma[4],x0,N)
refline <- qskewlap(c(.025, .975))
rwR <- list(rwR1,rwR2,rwR3,rwR4)
rwC <- list(rwC1,rwC2,rwC3,rwC4)
for (i in 1:4){
  plot(1:2000, rwR[[i]]$x, type = "l", xlab = paste("sigma=",sigma[i]), ylab = "x", main = "Using R codes")
  abline(h=refline)
}
for (i in 1:4) {
  plot(1:2000, rwC[[i]][1:2000], type = "l", xlab = paste("sigma=",sigma[i]), ylab = "x", main = "Using Rcpp")
  abline(h=refline)
}
#acceptance rate with rwMetropolosR
print(c(1-rwR1$k/N, 1-rwR2$k/N, 1-rwR3$k/N, 1-rwR4$k/N))
#acceptance rate with rwMetropolosC
print(c(1-rwC1[N+1]/N, 1-rwC2[N+1]/N, 1-rwC3[N+1]/N, 1-rwC4[N+1]/N))

## -----------------------------------------------------------------------------
#qqplot
a <- ppoints(100)
for (i in 1:4) {
  QR <- qskewlap(a)
  Q <- quantile(rwR[[i]]$x,a)
  qqplot(QR, Q, main=paste("Using R codes with sigma=",sigma[i]),xlab = "standard Laplace Quantiles", ylab = "Sample Quantiles")
  lines(QR,QR)
}
for (i in 1:4) {
  QR <- qskewlap(a)
  Q <- quantile(rwC[[i]][1:2000],a)
  qqplot(QR, Q, main=paste("Using Rcpp with sigma=",sigma[i]), xlab = "standard Laplace Quantiles", ylab = "Sample Quantiles")
  lines(QR,QR)
}

## -----------------------------------------------------------------------------
summary <- list(4)
for (i in 1:4) {
  ts <- microbenchmark(rwR[[i]] <- rwMetropolosR(sigma[i],x0,N), rwC[[i]] <- rwMetropolosC(sigma[i],x0,N))
  summary[[i]]<-list(paste0('sigma=',sigma[i]),summary(ts)[,c(1,3,5,6)])
}
print(summary)

