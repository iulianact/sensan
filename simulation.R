#----------------------------------------------------------------------------------
# Simulate nrep populations of size N and estimate target causal parameters.
#----------------------------------------------------------------------------------

rm(list=ls())
set.seed(1)

library(MASS)
library(sandwich)

source("functions.R")

nrep <- 10000

#----------------------------------------------------------------------------------
# Simulation parameters
#----------------------------------------------------------------------------------

# sample size
N <- 500
#-----
# U,Z ~ MultivariateNormal(Mu,Sig)
muU <- 0
muZ <- 0
sigU <- 1
sigZ <- 1
rho <- 0.2
Mu <- c(muU,muZ)
Sig <- matrix(c(sigU^2, rho*sigU*sigZ, rho*sigU*sigZ, sigZ^2),nrow=2,ncol=2)
#-----
# U|Z ~ Normal(muUZ,sigUZ)
muUZ <- function(z){muU + rho * sigU/sigZ * (z-muZ)}
sigUZ <- sqrt(1-rho^2) * sigU

#----------------------------------------------------------------------------------
# Simulation function
#----------------------------------------------------------------------------------

#default scenario is Scenario I
sim_binary <- function(a0=0.5, a1=-0.5, a2=-2, a3=0, a4=0, #X ~ Bernoulli(expit(a0+a1*Z+a2*U+a3*Z*U+a4*Z^2))
                       b0=-0.5, b1=2, b2=2, b3=2, b4=0, b5=0 # Y ~ Bernoulli(expit(b0+b1*Z+b2*U+b3*X+b4*X*Z+b5*U*Z))
                       ){
  
  a0<<-a0; a1<<-a1; a2<<-a2; a3<<-a3; a4<<-a4
  b0<<-b0; b1<<-b1; b2<<-b2; b3<<-b3; b4<<-b4; b5<<-b5
  
  UZ <- mvrnorm(N,Mu,Sig)
  U <- UZ[,1]
  Z <- UZ[,2]
  X <- rbinom(N,1,plogis(a0 + a1*Z + a2*U + a3*Z*U + a4*Z^2))
  Y <- rbinom(N,1,plogis(b0 + b1*Z + b2*U + b3*X + b4*X*Z + b5*U*Z))
  
  dat <- data.frame(Y=Y,X=X,Z=Z)

  alpha01 <- alphafun(0,1)
  alpha10 <- alphafun(1,0)

  realYvals <- realY(alpha01,alpha10)
  compYvals <- compY(alpha01,alpha10)

  list(
       REALYest = c(realYvals[1],realYvals[2]),
       COMPYest = c(compYvals[1],compYvals[2]),

       ORYest=ORY(dat,alpha01,alpha10,level=0.95),
       IPWYest=IPWY(dat,alpha01,alpha10,level=0.95),
       DRYest=DRY(dat,alpha01,alpha10,level=0.95),

       stdYest=ORY(dat,alpha01=0,alpha10=0,level=0.95),
       ipwYest=IPWY(dat,alpha01=0,alpha10=0,level=0.95),
       drYest=DRY(dat,alpha01=0,alpha10=0,level=0.95),
       
       alpha=c(alpha01,alpha10)
       )
}

set.seed(4) 
resSIM_binary_I <- lapply(1:nrep,FUN=function(x){sim_binary()})

set.seed(51) 
resSIM_binary_II <- lapply(1:nrep,FUN=function(x){sim_binary(a2=1,b2=0,a4=-2)})

set.seed(61) 
resSIM_binary_III <- lapply(1:nrep,FUN=function(x){sim_binary(a2=-0.5,b4=-2)})
