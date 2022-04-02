library(MASS)
library(GIGrvg)
library(truncnorm)
library(tibble)
library(rmutil)
library(msm)
library(invgamma)

# Sample beta
sampleBeta <- function(X,z,w,tau2,theta,kappa,beta0,sigma=1) {
  term1 <- 0
  for (i in 1:length(z)) {
    term1 <- term1+X[i,]%*%t(X[i,])/(kappa^2*sigma^2*w[i])
  }
  cov <- solve(diag(1/tau2)+term1)
  term2 <- 0
  for (i in 1:length(z)) {
    term2 <- term2+X[i,]*(z[i]-theta*sigma*w[i]-beta0)/(kappa^2*sigma^2*w[i])
  }
  mean <- cov%*%term2
  return(mvrnorm(1,mean,cov))
}

# Sample w
sampleW <- function(X,z,beta,beta0,theta,kappa,sigma=1) {
  lambda <- 0.5
  chi <- NULL
  for (i in 1:length(z)) {
    chi <- c(chi,(z[i]-beta0-array(X[i,],dim=c(1,length(beta)))%*%beta)^2/(sigma^2*kappa^2))
  }
  psi <- (theta^2/(sigma^2*kappa^2) + 2)
  return(sapply(1:length(z),function(x) rgig(1,lambda,chi[x],psi)))
}

# Sample z
sampleZ <- function(y,X,beta,beta0,delta,w,kappa,theta,sigma=1) {
  mean <- NULL
  for (i in 1:length(w)) {
    mean <- c(mean,beta0+array(X[i,],dim=c(1,length(beta)))%*%beta+theta*sigma*w[i])
  }
  sd <- kappa*sigma*sqrt(w)
  delta.new <- c(-10000000,delta,10000000)
  return(rtruncnorm(length(y),delta.new[y],delta.new[y+1],mean,sd))
}

# Sample beta0
sampleBeta0 <- function(X,z,w,beta,theta,kappa,alpha,sigma=1) {
  alpha.tilde <- 1/((1/alpha)+sum(1/(kappa^2*sigma^2*w)))
  phi <- 0
  for (i in 1:length(z)) {
    phi <- phi+(1/(kappa^2*sigma^2*w[i]))*(z[i]-array(X[i,],dim=c(1,length(beta)))%*%beta-theta*sigma*w[i])
  }
  return(rnorm(1,phi*alpha.tilde,sqrt(alpha.tilde)))
}

# Sample tau
sampleTau2 <- function(lambda2,beta,sigma=1) {
  lambda <- 0.5
  chi <- beta^2/sigma^2
  psi <- lambda2
  return(sapply(1:length(beta),function(x) rgig(1,lambda,chi[x],psi[x])))
}

# Sample lambda
sampleLambda2 <- function(r,gamma,tau2,sigma=1) {
  return(rgamma(length(tau2),r+1,(gamma+(tau2/(2)))))
}

# Sample delta
sampleDelta <- function(y,z,delta) {
  delta.new <- c(-10000000,delta,10000000)
  for (i in 3:(length(delta.new)-1)) {
    lower <- min(c(max(z[y==(i-1)]),delta.new[i],10000000))
    upper <- max(c(min(z[y==i]),delta.new[i-1],-10000000))
    delta.new[i] <- runif(1,lower,upper)
  }
  return(delta.new[2:(length(delta.new)-1)])
}

# Sample sigma
sampleSigma <- function(z,kappa,w,beta,theta,tau2,lambda2) {
  term1 <- (3*length(z)/2) + (3*length(beta))/2
  term2 <- sum((z-rowSums(X%*%beta)-theta*w)/(2*kappa^2*w))+sum(w)+
    sum(lambda2*tau2/2)+sum(beta^2/(2*tau2))
  return(rinvgamma(1,term1,term2))
}

# Sampler
BOQR_BAL <- function(data,q,iters=10000) {
  # Initialize parameters
  y <- data[,1]
  X <- as.matrix(data[,2:ncol(data)])
  beta0 <- 0
  beta <- rep(1,ncol(X))
  theta <- (1-2*q)/(q*(1-q))
  kappa <- sqrt(2/(q*(1-q)))
  w <- rexp(length(y),1)
  delta <- 0:(max(y)-2)
  delta.new <- c(-10000000,delta,10000000)
  z <- rtruncnorm(length(y),delta.new[y],delta.new[y+1],0,1)
  r <- 1
  gamma <- 0.1
  alpha <- 10
  lambda2 <- rep(sqrt(0.2),length(beta))
  tau2 <- rep(0.001,length(beta))
  sigma <- 1
  
  # Store estimates
  beta.chain <- array(NA,dim=c(iters,length(beta)))
  delta.chain <- array(NA,dim=c(iters,length(delta)))
  
  # Run sampler
  for (i in 1:iters) {
    beta <- sampleBeta(X,z,w,tau2,theta,kappa,beta0,sigma)
    w <- sampleW(X,z,beta,beta0,theta,kappa,sigma)
    z <- sampleZ(y,X,beta,beta0,delta,w,kappa,theta,sigma)
    beta0 <- sampleBeta0(X,z,w,beta,theta,kappa,alpha,sigma)
    tau2 <- sampleTau2(lambda2,beta,sigma)
    lambda2 <- sampleLambda2(r,gamma,tau2,sigma)
    delta <- sampleDelta(y,z,delta)
    
    beta.chain[i,] <- beta
    delta.chain[i,] <- delta 
  }
  
  return(list(beta.chain,delta.chain))
}
