#----------------------------------------------------------------------------------
# Auxiliary functions
#----------------------------------------------------------------------------------

gfun <- qlogis; ginvfun <- plogis

#--------------------------

exapproxfun <- function(mu,sigma,fn,qtol=10){
  fnint <- function(z){fn(z) * 1/sqrt(2*pi*sigma^2) * exp(- 1/(2*sigma^2) * (z-mu)^2)}
  integrate(fnint,lower=qnorm(10^(-qtol),mu,sigma),upper=qnorm(10^(-qtol),mu,sigma,lower.tail = FALSE))$value
}

#--------------------------

RDfun <- function(est){
  p0 <- est[1]
  p1 <- est[2]
  return(p1-p0)
}
RRfun <- function(est){
  p0 <- est[1]
  p1 <- est[2]
  return(p1/p0)
}
ORfun <- function(est){
  p0 <- est[1]
  p1 <- est[2]
  return((p1/(1-p1))/(p0/(1-p0)))
}

#--------------------------

derRDfun <- function(est){
  return(c(1,-1))
}

derRRfun <- function(est){
  p0 <- est[1]
  p1 <- est[2]
  return(c(1/p0,-p1/(p0^2)))
}

derORfun <- function(est){
  p0 <- est[1]
  p1 <- est[2]
  return(c((1-p0)/(p0*(1-p1)^2),-p1/(p0^2*(1-p1))))
}

#--------------------------

rhatfun_or <- function(x1,x2,dat){
  N <- nrow(dat)
  fitr <- glm(Y ~ .,family="binomial", data=dat)
  newdat <- dat
  newdat$X <- rep(x1,N)
  predr <- predict(fitr,newdata=newdat,type="response")
  sum(predr * as.numeric(dat$X==x2)) / sum(as.numeric(dat$X==x2))
}

rhatfun_ipw <- function(x1,x2,dat){
  N <- nrow(dat)
  datx <- dat[,-1]
  fitrx <- glm(X ~ .,family="binomial", data=datx)
  predrx <- predict(fitrx,type="response")
  if(x1==x2){
    ww <- rep(1,N)
  } else {
    if(x2==1){
      ww <- predrx/(1-predrx)
    } else {
      ww <- (1-predrx)/predrx
    }
  }
  sum(ww * as.numeric(dat$X==x1) * dat$Y) / sum(as.numeric(dat$X==x2))
}

rhatfun_dr <- function(x1,x2,dat){
  N <- nrow(dat)
  newdat <- dat
  newdat$X <- rep(x1,N)
  fitr <- glm(Y ~ .,family="binomial", data=dat)
  predr <- predict(fitr,newdata=newdat,type="response")
  
  datx <- dat[,-1]
  fitrx <- glm(X ~ .,family="binomial", data=datx)
  predrx <- predict(fitrx,type="response")
  
  if(x1==x2){
    ww <- as.numeric(dat$X==x1)/mean(dat$X==x2)
  } else {
    if(x2==1){
      ww <- (predrx/(1-predrx)) * (as.numeric(dat$X==x1)/mean(dat$X==x2))
    } else { 
      ww <- ((1-predrx)/predrx) * (as.numeric(dat$X==x1)/mean(dat$X==x2))
    }
  }
  sum(ww * dat$Y - (ww - as.numeric(dat$X==x2)/mean(dat$X==x2))*predr) / N
}

#--------------------------

alphafun <- function(x1,x2){
  
  #Ax = int_{U|Z}{P(X=x|Z,U)dF(U|z)} (a function of Z)
  A1 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2)}) }) }
  A0 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){1-plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2)}) }) }
  
  #Bx = int_Z{Ax dF(z)}
  B1 <- exapproxfun(muZ,sigZ,fn=A1)
  B0 <- exapproxfun(muZ,sigZ,fn=A0)
  
  if(x1==x2){
    0
  } else {
    if(x1==1) {
      #s(1,0,0)
      inners100 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){plogis(b0+b1*w+b2*u+b3+b4*w+b5*w*u)*(1-plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2))}) }) }
      s100 <- (1/B0) * exapproxfun(muZ,sigZ,fn=inners100)
      #s(1,1,0)
      inners110 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){(A0(w)/A1(w))*plogis(b0+b1*w+b2*u+b3+b4*w+b5*w*u)*plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2)}) }) }
      s110 <- (1/B0) * exapproxfun(muZ,sigZ,fn=inners110)
      #alpha(1,0)
      alpha10 <- gfun(s100) - gfun(s110)
      
      alpha10
      
    } else {
      #s(0,1,1)
      inners011 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){plogis(b0+b1*w+b2*u+b5*w*u)*plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2)}) }) }
      s011 <- (1/B1) * exapproxfun(muZ,sigZ,fn=inners011)
      #s(0,0,1)
      inners001 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){(A1(w)/A0(w))*plogis(b0+b1*w+b2*u+b5*w*u)*(1-plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2))}) }) }
      s001 <- (1/B1) * exapproxfun(muZ,sigZ,fn=inners001)
      #alpha(0,1)
      alpha01 <- gfun(s011)-gfun(s001)
      
      alpha01
      
    }
  }
}

#--------------------------
realY <- function(alpha01,alpha10){
  #P(Y(1)=1)
  inner1 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){plogis(b0+b1*w+b2*u+b3+b4*w+b5*w*u)}) })}
  rY11 <- exapproxfun(muZ,sigZ,fn=inner1)
  #P(Y(0)=1)
  inner0 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){plogis(b0+b1*w+b2*u+b5*w*u)}) })}
  rY01 <- exapproxfun(muZ,sigZ,fn=inner0)
  
  c(rY01,rY11)
}

compY <- function(alpha01,alpha10){
  #Ax = int_{U|Z}{P(X=x|Z,U)dF(U|z)} (a function of Z)
  A1 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2)}) }) }
  A0 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){1-plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2)}) }) }
  #Bx = int_Z{Ax dF(z)}
  B1 <- exapproxfun(muZ,sigZ,fn=A1)
  B0 <- exapproxfun(muZ,sigZ,fn=A0)

  pX <- B1
  
  #r(1,1)
  innerr11 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){plogis(b0+b1*w+b2*u+b3+b4*w+b5*w*u)*plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2)}) }) }
  r11 <- (1/B1) * exapproxfun(muZ,sigZ,fn=innerr11)
  #r(0,0)
  innerr00 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){plogis(b0+b1*w+b2*u+b5*w*u)*(1-plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2))}) }) }
  r00 <- (1/B0) * exapproxfun(muZ,sigZ,fn=innerr00)
  
  #s(1,1,0)
  inners110 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){(A0(w)/A1(w))*plogis(b0+b1*w+b2*u+b3+b4*w+b5*w*u)*plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2)}) }) }
  s110 <- (1/B0) * exapproxfun(muZ,sigZ,fn=inners110)
  #s(0,0,1)
  inners001 <- function(z) { sapply(z,function(w) { exapproxfun(muUZ(w),sigUZ,fn=function(u){(A1(w)/A0(w))*plogis(b0+b1*w+b2*u+b5*w*u)*(1-plogis(a0+a1*w+a2*u+a3*w*u+a4*w^2))}) }) }
  s001 <- (1/B1) * exapproxfun(muZ,sigZ,fn=inners001)
  
  #P(Y(1)=1) = r(1,1)P(X=1) + g^{-1}( g(r(1,0)) + alpha(1,0) )P(X=0)
  cY11 <- r11 * pX + ginvfun(gfun(s110) + alpha10) * (1-pX)
  #P(Y(0)=1) = r(0,0)P(X=0) + g^{-1}( g(r(0,1)) + alpha(0,1) )P(X=1)
  cY01 <- r00 * (1-pX) + ginvfun(gfun(s001) + alpha01) * pX
  
  c(cY01,cY11)
}


#----------------------------------------------------------------------------------
# Estimate p(Y(1)=1), p(Y(0)=1) and contrasts (RD, RR, OR) using 
# outcome regression, inverse probability weighting and doubly robust estimators
#----------------------------------------------------------------------------------
#ESTIMATE using OR 
ORY <- function(dat,alpha01,alpha10,level=0.95){
  N <- nrow(dat)
  X <- dat$X
  nz <- ncol(dat)-2
  
  rhat00_or <- rhatfun_or(0,0,dat)
  rhat01_or <- rhatfun_or(0,1,dat)
  rhat10_or <- rhatfun_or(1,0,dat)
  rhat11_or <- rhatfun_or(1,1,dat)
  
  eY11_or <- sum( ginvfun( gfun(rhat11_or * X + rhat10_or *(1-X)) + alpha10 * as.numeric(X==0) ) )/N
  eY01_or <- sum( ginvfun( gfun(rhat01_or * X + rhat00_or *(1-X)) + alpha01 * as.numeric(X==1) ) )/N
  
  est <- c(eY01_or,eY11_or)
  
  AOR0 <- AOR(0,dat,alpha01,alpha10)
  AOR1 <- AOR(1,dat,alpha01,alpha10)
    
  AOR_ext <- cbind(AOR1,matrix(0,ncol=3,nrow=nrow(AOR1)))
  AOR_ext <- rbind(AOR_ext, matrix(0,nrow=3,ncol=ncol(AOR_ext)))
  rownames(AOR_ext) <- colnames(AOR_ext) <- c(colnames(AOR1),"tx2","tx1","qx2")
  diag(AOR_ext[(nrow(AOR_ext)-2):nrow(AOR_ext),(ncol(AOR_ext)-2):ncol(AOR_ext)]) <- 1
  AOR_ext["tx2",1:(2+nz+2)] <- AOR0["rx1",c(2,1,3:(2+nz+2))] 
  AOR_ext["tx1",1:(2+nz+2)] <- AOR0["rx2",c(2,1,3:(2+nz+2))]
  AOR_ext["qx2",(ncol(AOR_ext)-2):ncol(AOR_ext)] <- AOR0["qx1",(ncol(AOR0)-2):ncol(AOR0)]
    
  newdat1 <- dat
  newdat1$X <- rep(1,N)
  newdat0 <- dat
  newdat0$X <- rep(0,N)
  fitr <- glm(Y ~ .,family="binomial", data=dat)
  predr1 <- predict(fitr,newdata=newdat1,type="response")
  predr0 <- predict(fitr,newdata=newdat0,type="response")
      
  estfunfitr <- estfun(fitr)
      
  alphaeasy <- function(x1,x2){
    if(x1==x2){
      0
    } else {
      if(x1==1) {
        alpha10
      } else {
        alpha01
      }
    }
  }
  rhatfun <- rhatfun_or
      
  mean1 <- sum(X==1)/N
  mean0 <- sum(X==0)/N
      
  rhat1_or_all <- sapply(X,FUN=function(x){ifelse(x==1,rhat11_or,rhat10_or)})
  rhat0_or_all <- sapply(X,FUN=function(x){ifelse(x==1,rhat01_or,rhat00_or)})
  alphaeasy1 <- sapply(X,FUN=function(x){ifelse(x==1,0,alpha10)})
  alphaeasy0 <- sapply(X,FUN=function(x){ifelse(x==1,alpha01,0)})
      
  U <- cbind((X==1) - mean1,
             (X==0) - mean0,
             estfunfitr,
             (X==1) / mean1 * predr1 - rhat11_or,
             (X==0) / mean0 * predr1 - rhat10_or,
             ginvfun(gfun(rhat1_or_all)+alphaeasy1) - eY11_or,
             (X==0) / mean0 * predr0 - rhat00_or,
             (X==1) / mean1 * predr0 - rhat01_or,
             ginvfun(gfun(rhat0_or_all)+alphaeasy0) - eY01_or
            )
  colnames(U) <- c(colnames(AOR1),"tx2","tx1","qx2")
    
  BOR_ext <- var(U)
    
  ss <- solve(AOR_ext)
  SWOR_ext <- ss %*% BOR_ext %*% t(ss) 
    
  idq <- c(which(colnames(SWOR_ext)=="qx1"),which(colnames(SWOR_ext)=="qx2"))
  covq <- data.frame(qx1 = SWOR_ext["qx1",idq] ,qx2 = SWOR_ext["qx2",idq])

  asvar_eY01_or <- covq[2,2] 
  asvar_eY11_or <- covq[1,1] 
  se <- c(sqrt(asvar_eY01_or)/sqrt(N),sqrt(asvar_eY11_or)/sqrt(N))
    
  z <- qnorm(1-(1-level)/2)
  ciest <- matrix(NA,nrow=2,ncol=2)
  ciest[1,] <- c(eY01_or - z*se[1], eY01_or + z*se[1])
  ciest[2,] <- c(eY11_or - z*se[2], eY11_or + z*se[2])
  rownames(ciest) <- c("P(Y(0)=1)","P(Y(1)=1)")
  colnames(ciest) <- c("lower","upper")
    
  fnest <- list()
  fnci <- list()
  fnse <- list()
    
  fnest[[1]] <- RDfun(est)
  derfun <- derRDfun
  var <- t(derfun(est)) %*% as.matrix(covq) %*% derfun(est)
  fnse[[1]] <- sqrt(var)/sqrt(N)
  fnci[[1]] <- c(fnest[[1]] - z * fnse[[1]], fnest[[1]] + z * fnse[[1]])
  
  fnest[[2]] <- RRfun(est)
  derfun <- derRRfun
  var <- t(derfun(est)) %*% as.matrix(covq) %*% derfun(est)
  fnse[[2]] <- sqrt(var)/sqrt(N)
  fnci[[2]] <- c(fnest[[2]] - z * fnse[[2]], fnest[[2]] + z * fnse[[2]])
    
  fnest[[3]] <- ORfun(est)
  derfun <- derORfun
  var <- t(derfun(est)) %*% as.matrix(covq) %*% derfun(est)
  fnse[[3]] <- sqrt(var)/sqrt(N)
  fnci[[3]] <- c(fnest[[3]] - z * fnse[[3]], fnest[[3]] + z * fnse[[3]])
    
  names(fnest) <- names(fnci) <- names(fnse) <- c("RD","RR","OR")

  list(est=est,
       se=se,
       ciest=ciest,
       fnest=fnest,
       fnse=fnse,
       fnci=fnci
       )
}

AOR <- function(x,dat,alpha01,alpha10){
  x1 <- x ; x2 <- 1-x 
  nz <- ncol(dat)-2 
  N <- nrow(dat)
  X <- dat$X
  
  rhatfun <- rhatfun_or
  
  newdat1 <- dat
  newdat1$X <- rep(x1,N)
  fitr <- glm(Y ~ .,family="binomial", data=dat)
  predr1 <- predict(fitr,newdata=newdat1,type="response")/c(1+exp(predict(fitr,newdata=newdat1)))
  
  alphaeasy <- function(x1,x2){
    if(x1==x2){
      0
    } else {
      if(x1==1) {
        alpha10
      } else {
        alpha01
      }
    }
  }
  
  Amat <- diag(1,7+nz) 
  Amat[3:(4+nz),3:(4+nz)] <- solve(bread(fitr)) 
  Amat[(4+nz+1),1] <- rhatfun(x1,x1,dat)/(sum(X==x1)/N)
  
  Amat[(4+nz+1),3] <- -(1/sum(X==x1)) * sum(as.numeric(X==x1)  * predr1)
  Amat[(4+nz+1),4] <- -(1/sum(X==x1)) * sum(as.numeric(X==x1) * x1 * predr1)
  for(k in 1:nz){
    Amat[(4+nz+1),(4+k)] <- -(1/sum(X==x1)) * sum(as.numeric(X==x1) * dat[,2+k] * predr1)
  }
  
  Amat[(4+nz+2),2] <- rhatfun(x1,x2,dat)/(sum(X==x2)/N)
  
  Amat[(4+nz+2),3] <- -(1/sum(X==x2)) * sum(as.numeric(X==x2) * predr1)
  Amat[(4+nz+2),4] <- -(1/sum(X==x2)) * sum(as.numeric(X==x2) * x1 * predr1)
  for(k in 1:nz){
    Amat[(4+nz+2),(4+k)] <- -(1/sum(X==x2)) * sum(as.numeric(X==x2) * dat[,2+k] * predr1)
  }
  
  Amat[(4+nz+3),(4+nz+1)] <- -sum(X==x1)/N
  
  dqdr2 <- exp(alphaeasy(x1,x2))/((1+rhatfun(x1,x2,dat)*(exp(alphaeasy(x1,x2))-1))^2)
  Amat[(4+nz+3),(4+nz+2)] <- -sum(X==x2)/N * dqdr2
  
  rownames(Amat) <- colnames(Amat) <- c("px1","px2",paste0("beta",seq(0,nz+1)),"rx1","rx2","qx1")
  Amat
}

#-----------------------------------------
#ESTIMATE using IPW 
IPWY <- function(dat,alpha01,alpha10,level=0.95){
  N <- nrow(dat)
  X <- dat$X
  Y <- dat$Y
  nz <- ncol(dat)-2
  
  rhat00_ipw <- rhatfun_ipw(0,0,dat)
  rhat01_ipw <- rhatfun_ipw(0,1,dat)
  rhat10_ipw <- rhatfun_ipw(1,0,dat)
  rhat11_ipw <- rhatfun_ipw(1,1,dat)
  
  eY11_ipw <- sum( ginvfun( gfun(rhat11_ipw * X + rhat10_ipw *(1-X)) + alpha10 * as.numeric(X==0) ) )/N
  eY01_ipw <- sum( ginvfun( gfun(rhat01_ipw * X + rhat00_ipw *(1-X)) + alpha01 * as.numeric(X==1) ) )/N
  
  est <- c(eY01_ipw,eY11_ipw)

  AIPW0 <- AIPW(0,dat,alpha01,alpha10)
  AIPW1 <- AIPW(1,dat,alpha01,alpha10)
    
  AIPW_ext <- cbind(AIPW1,matrix(0,ncol=3,nrow=nrow(AIPW1)))
  AIPW_ext <- rbind(AIPW_ext, matrix(0,nrow=3,ncol=ncol(AIPW_ext)))
  rownames(AIPW_ext) <- colnames(AIPW_ext) <- c(colnames(AIPW1),"tx2","tx1","qx2")
  diag(AIPW_ext[(nrow(AIPW_ext)-2):nrow(AIPW_ext),(ncol(AIPW_ext)-2):ncol(AIPW_ext)]) <- 1
  AIPW_ext["tx2",1:(2+nz+1)] <- AIPW0["rx1",c(2,1,3:(2+nz+1))] 
  AIPW_ext["tx1",1:(2+nz+1)] <- AIPW0["rx2",c(2,1,3:(2+nz+1))] 
  AIPW_ext["qx2",(ncol(AIPW_ext)-2):ncol(AIPW_ext)] <- AIPW0["qx1",(ncol(AIPW0)-2):ncol(AIPW0)]

  datx <- dat[,-1]
  fitrx <- glm(X ~ .,family="binomial", data=datx)
    
  estfunfitrx <- estfun(fitrx)
    
  predpi0 <- predict(fitrx,type="response")/(1-predict(fitrx,type="response")) 
  predpi1 <- (predict(fitrx,type="response")/(1-predict(fitrx,type="response")))^(-1)  

  rhatfun <- rhatfun_ipw
    
  alphaeasy <- function(x1,x2){
    if(x1==x2){
      0
    } else {
      if(x1==1) {
        alpha10
      } else {
        alpha01
      }
    }
  }
    
  rhat1_ipw_all <- sapply(X,FUN=function(x){ifelse(x==1,rhat11_ipw,rhat10_ipw)})
  rhat0_ipw_all <- sapply(X,FUN=function(x){ifelse(x==1,rhat01_ipw,rhat00_ipw)})
  alphaeasy1 <- sapply(X,FUN=function(x){ifelse(x==1,0,alpha10)})
  alphaeasy0 <- sapply(X,FUN=function(x){ifelse(x==1,alpha01,0)})
    
  mean1 <- mean(X==1)
  mean0 <- mean(X==0)
    
  U <- cbind((X==1) - mean1,
             (X==0) - mean0,
             estfunfitrx,
             Y * (X==1) / (mean1) - rhat11_ipw,
             Y * (X==1) / (mean0) * predpi1 - rhat10_ipw,
             ginvfun(gfun(rhat1_ipw_all)+alphaeasy1) - eY11_ipw,
             Y * (X==0) / (mean0) - rhat00_ipw,
             Y * (X==0) / (mean1) * predpi0 - rhat01_ipw,
             ginvfun(gfun(rhat0_ipw_all)+alphaeasy0) - eY01_ipw
            )
  colnames(U) <- c(colnames(AIPW1),"tx2","tx1","qx2")
    
  BIPW_ext <- var(U)
    
  ss <- solve(AIPW_ext)
  SIPW_ext <- ss %*% BIPW_ext %*% t(ss) 
    
  idq <- c(which(colnames(SIPW_ext)=="qx1"),which(colnames(SIPW_ext)=="qx2"))
  covq <- data.frame(qx1 = SIPW_ext["qx1",idq] ,qx2 = SIPW_ext["qx2",idq])

  asvar_eY01_ipw <- covq[2,2] 
  asvar_eY11_ipw <- covq[1,1] 
    
  se <- c(sqrt(asvar_eY01_ipw)/sqrt(N),sqrt(asvar_eY11_ipw)/sqrt(N))
    
  z <- qnorm(1-(1-level)/2)
  ciest <- matrix(NA,nrow=2,ncol=2)
  ciest[1,] <- c(eY01_ipw - z*se[1], eY01_ipw + z*se[1])
  ciest[2,] <- c(eY11_ipw - z*se[2], eY11_ipw + z*se[2])
  rownames(ciest) <- c("P(Y(0)=1)","P(Y(1)=1)")
  colnames(ciest) <- c("lower","upper")
    
  fnest <- list()
  fnci <- list()
  fnse <- list()
    
  fnest[[1]] <- RDfun(est)
  derfun <- derRDfun
  var <- t(derfun(est)) %*% as.matrix(covq) %*% derfun(est)
  fnse[[1]] <- sqrt(var)/sqrt(N)
  fnci[[1]] <- c(fnest[[1]] - z * fnse[[1]], fnest[[1]] + z * fnse[[1]])
    
  fnest[[2]] <- RRfun(est)
  derfun <- derRRfun
  var <- t(derfun(est)) %*% as.matrix(covq) %*% derfun(est)
  fnse[[2]] <- sqrt(var)/sqrt(N)
  fnci[[2]] <- c(fnest[[2]] - z * fnse[[2]], fnest[[2]] + z * fnse[[2]])
    
  fnest[[3]] <- ORfun(est)
  derfun <- derORfun
  var <- t(derfun(est)) %*% as.matrix(covq) %*% derfun(est)
  fnse[[3]] <- sqrt(var)/sqrt(N)
  fnci[[3]] <- c(fnest[[3]] - z * fnse[[3]], fnest[[3]] + z * fnse[[3]])
    
  names(fnest) <- names(fnci) <-  names(fnse) <- c("RD","RR","OR")

  list(est=est,
       se=se,
       ciest=ciest,
       fnest=fnest,
       fnse=fnse, 
       fnci=fnci)
}

AIPW <- function(x,dat,alpha01,alpha10){
  x1 <- x ; x2 <- 1-x 
  nz <- ncol(dat)-2 
  N <- nrow(dat)
  Y <- dat$Y
  X <- dat$X
  
  rhatfun <- rhatfun_ipw
  
  datx <- dat[,-1]
  fitrx <- glm(X ~ .,family="binomial", data=datx)

  if(x2==1){
    preddpi <- exp(predict(fitrx)) 
  } else {
    preddpi <- -1/exp(predict(fitrx))  
  }
  
  alphaeasy <- function(x1,x2){
    if(x1==x2){
      0
    } else {
      if(x1==1) {
        alpha10
      } else {
        alpha01
      }
    }
  }
  
  Amat <- diag(1,6+nz)
  Amat[3:(3+nz),3:(3+nz)] <- solve(bread(fitrx))
  
  Amat[4+nz,1] <- rhatfun(x1,x1,dat)/(sum(X==x1)/N)
  
  Amat[4+nz+1,2] <- rhatfun(x1,x2,dat)/(sum(X==x2)/N)
  
  Amat[4+nz+1,3] <- -(1/sum(X==x2)) * sum(as.numeric(X==x1) * preddpi * Y)
  for(k in 1:nz){
    Amat[4+nz+1,3+k] <- -(1/sum(X==x2)) * sum(as.numeric(X==x1) * dat[,2+k] * preddpi * Y)
  }
  
  Amat[4+nz+2,4+nz] <- -sum(X==x1)/N
  
  dqdr2 <- exp(alphaeasy(x1,x2))/((1+rhatfun(x1,x2,dat)*(exp(alphaeasy(x1,x2))-1))^2)
  Amat[4+nz+2,4+nz+1] <- -sum(X==x2)/N * dqdr2
  
  rownames(Amat) <- colnames(Amat) <- c("px1","px2",paste0("beta",seq(0,nz)),"rx1","rx2","qx1")
  Amat
}

#-----------------------------------------
#ESTIMATE using DR 
DRY <- function(dat,alpha01,alpha10,level=0.95){
  X <- dat$X 
  Y <- dat$Y
  N <- nrow(dat)
  nz <- ncol(dat)-2
  
  rhat00_dr <- rhatfun_dr(0,0,dat)
  rhat01_dr <- rhatfun_dr(0,1,dat)
  rhat10_dr <- rhatfun_dr(1,0,dat)
  rhat11_dr <- rhatfun_dr(1,1,dat)
  
  eY11_dr <- sum( ginvfun( gfun(rhat11_dr * X + rhat10_dr *(1-X)) + alpha10 * as.numeric(X==0) ) )/N
  eY01_dr <- sum( ginvfun( gfun(rhat01_dr * X + rhat00_dr *(1-X)) + alpha01 * as.numeric(X==1) ) )/N
  
  est <- c(eY01_dr,eY11_dr)

  ADR0 <- ADR(0,dat,alpha01,alpha10)
  ADR1 <- ADR(1,dat,alpha01,alpha10)
    
  ADR_ext <- cbind(ADR1,matrix(0,ncol=3,nrow=nrow(ADR1)))
  ADR_ext <- rbind(ADR_ext, matrix(0,nrow=3,ncol=ncol(ADR_ext)))
  rownames(ADR_ext) <- colnames(ADR_ext) <- c(colnames(ADR1),"tx2","tx1","qx2")
  diag(ADR_ext[(nrow(ADR_ext)-2):nrow(ADR_ext),(ncol(ADR_ext)-2):ncol(ADR_ext)]) <- 1
  ADR_ext["tx2",1:(2+2*nz+2+1)] <- ADR0["rx1",c(2,1,3:(2+2*nz+2+1))]
  ADR_ext["tx1",1:(2+2*nz+2+1)] <- ADR0["rx2",c(2,1,3:(2+2*nz+2+1))]
  ADR_ext["qx2",(ncol(ADR_ext)-2):ncol(ADR_ext)] <- ADR0["qx1",(ncol(ADR0)-2):ncol(ADR0)]
    
  newdat1 <- dat
  newdat1$X <- rep(1,N)
  newdat0 <- dat
  newdat0$X <- rep(0,N)
  fitr <- glm(Y ~ .,family="binomial", data=dat)
  predr1 <- predict(fitr,newdata=newdat1,type="response")
  predr0 <- predict(fitr,newdata=newdat0,type="response")
    
  datx <- dat[,-1]
  fitrx <- glm(X ~ .,family="binomial", data=datx)
    
  predpi0 <- predict(fitrx,type="response")/(1-predict(fitrx,type="response"))
  predpi1 <- (predict(fitrx,type="response")/(1-predict(fitrx,type="response")))^(-1)

  estfunfitrx <- estfun(fitrx)
  estfunfitr <- estfun(fitr)
      
  alphaeasy <- function(x1,x2){
    if(x1==x2){
      0
    } else {
      if(x1==1) {
        alpha10
      } else {
        alpha01
      }
    }
  }
    
  rhatfun <- rhatfun_dr
    
  mean1 <- mean(X==1)
  mean0 <- mean(X==0)
    
  rhat1_dr_all <- sapply(X,FUN=function(x){ifelse(x==1,rhat11_dr,rhat10_dr)})
  rhat0_dr_all <- sapply(X,FUN=function(x){ifelse(x==1,rhat01_dr,rhat00_dr)})
  alphaeasy1 <- sapply(X,FUN=function(x){ifelse(x==1,0,alpha10)})
  alphaeasy0 <- sapply(X,FUN=function(x){ifelse(x==1,alpha01,0)})
    
  U <- cbind((X==1) - mean1,
             (X==0) - mean0,
             estfunfitrx,
             estfunfitr,
             ((X==1) / (mean1)) * Y - rhat11_dr,
             ((X==1) / (mean0)) * predpi1 * (Y - predr1) + ((X==0) / (mean0)) * predr1 - rhat10_dr,
             ginvfun(gfun(rhat1_dr_all)+alphaeasy1) - eY11_dr,
             ((X==0) / (mean0)) * Y - rhat00_dr,
             ((X==0) / (mean1)) * predpi0 * (Y - predr0) + ((X==1) / (mean1)) * predr0 - rhat01_dr,
             ginvfun(gfun(rhat0_dr_all)+alphaeasy0) - eY01_dr
            )
  colnames(U) <- c("px1","px2",paste0("gamma",seq(0,nz)),paste0("beta",seq(0,nz+1)),"rx1","rx2","qx1","tx2","tx1","qx2")

  BDR_ext <- var(U)
    
  ss <- solve(ADR_ext)
  SDR_ext <- ss %*% BDR_ext %*% t(ss) 
    
  idq <- c(which(colnames(SDR_ext)=="qx1"),which(colnames(SDR_ext)=="qx2"))
  covq <- data.frame(qx1 = SDR_ext["qx1",idq] ,qx2 = SDR_ext["qx2",idq])

  asvar_eY01_dr <- covq[2,2] 
  asvar_eY11_dr <- covq[1,1] 
    
  se <- c(sqrt(asvar_eY01_dr)/sqrt(N),sqrt(asvar_eY11_dr)/sqrt(N))
    
  z <- qnorm(1-(1-level)/2)
  ciest <- matrix(NA,nrow=2,ncol=2)
  ciest[1,] <- c(eY01_dr - z*se[1], eY01_dr + z*se[1])
  ciest[2,] <- c(eY11_dr - z*se[2], eY11_dr + z*se[2])
  rownames(ciest) <- c("P(Y(0)=1)","P(Y(1)=1)")
  colnames(ciest) <- c("lower","upper")
    
  fnest <- list()
  fnci <- list()
  fnse <- list()
    
  fnest[[1]] <- RDfun(est)
  derfun <- derRDfun
  var <- t(derfun(est)) %*% as.matrix(covq) %*% derfun(est)
  fnse[[1]] <- sqrt(var)/sqrt(N)
  fnci[[1]] <- c(fnest[[1]] - z * fnse[[1]], fnest[[1]] + z * fnse[[1]])
    
  fnest[[2]] <- RRfun(est)
  derfun <- derRRfun
  var <- t(derfun(est)) %*% as.matrix(covq) %*% derfun(est)
  fnse[[2]] <- sqrt(var)/sqrt(N)
  fnci[[2]] <- c(fnest[[2]] - z * fnse[[2]], fnest[[2]] + z * fnse[[2]])
    
  fnest[[3]] <- ORfun(est)
  derfun <- derORfun
  var <- t(derfun(est)) %*% as.matrix(covq) %*% derfun(est)
  fnse[[3]] <- sqrt(var)/sqrt(N)
  fnci[[3]] <- c(fnest[[3]] - z * fnse[[3]], fnest[[3]] + z * fnse[[3]])
    
  names(fnest) <- names(fnci) <- names(fnse) <-c("RD","RR","OR")

  list(est=est,
       se=se,
       ciest=ciest,
       fnest=fnest,
       fnse=fnse, 
       fnci=fnci)
}

ADR <- function(x,dat,alpha01,alpha10){
  x1 <- x ; x2 <- 1-x 
  N <- nrow(dat)
  nz <- ncol(dat)-2
  X <- dat$X
  Y <- dat$Y
  
  rhatfun <- rhatfun_dr
  
  newdat1 <- dat
  newdat1$X <- rep(x1,N)
  fitr <- glm(Y ~ .,family="binomial", data=dat)
  predr <- predict(fitr,newdata=newdat1,type="response")
  predr1 <- predict(fitr,newdata=newdat1,type="response")/c(1+exp(predict(fitr,newdata=newdat1)))
  
  datx <- dat[,-1]
  fitrx <- glm(X ~ .,family="binomial", data=datx)
  
  if(x2==1){
    preddpi <- exp(predict(fitrx)) 
    predpi <- predict(fitrx,type="response")/(1-predict(fitrx,type="response"))
  } else {
    preddpi <- -1/exp(predict(fitrx)) 
    predpi <- (predict(fitrx,type="response")/(1-predict(fitrx,type="response")))^(-1)
  }
  
  alphaeasy <- function(x1,x2){
    if(x1==x2){
      0
    } else {
      if(x1==1) {
        alpha10
      } else {
        alpha01
      }
    }
  }

  Amat <- diag(1,5+1+2*(nz+1))
  
  Amat[3:(3+nz),3:(3+nz)] <- solve(bread(fitrx))
  Amat[(4+nz):(5+2*nz),(4+nz):(5+2*nz)] <- solve(bread(fitr))
  
  Amat[6+2*nz,1] <- rhatfun(x1,x1,dat)/(sum(X==x1)/N) 
  
  Amat[7+2*nz,2] <- rhatfun(x1,x2,dat)/(sum(X==x2)/N) 
  
  Amat[7+2*nz,3] <- -(1/sum(X==x2)) * sum(as.numeric(X==x1) * preddpi * (Y-predr))
  for(k in 1:nz){
    Amat[7+2*nz,3+k] <- -(1/sum(X==x2)) * sum(as.numeric(X==x1) * dat[,2+k] * preddpi * (Y-predr))
  }
  
  Amat[7+2*nz,4+nz] <- -(1/sum(X==x2)) * sum(predr1 * (as.numeric(X==x2) - predpi * as.numeric(X==x1)))
  Amat[7+2*nz,5+nz] <- -(1/sum(X==x2)) * sum(predr1 * x1 *(as.numeric(X==x2) - predpi * as.numeric(X==x1)))
  for(k in 1:nz){
    Amat[7+2*nz,5+nz+k] <- -(1/sum(X==x2)) * sum(predr1 * dat[,2+k] * (as.numeric(X==x2) - predpi * as.numeric(X==x1)))
  }
  
  Amat[8+2*nz,6+2*nz] <- -sum(X==x1)/N
  dqdr2 <- exp(alphaeasy(x1,x2))/((1+rhatfun(x1,x2,dat)*(exp(alphaeasy(x1,x2))-1))^2)
  Amat[8+2*nz,7+2*nz] <- -sum(X==x2)/N * dqdr2
  
  rownames(Amat) <- colnames(Amat) <- c("px1","px2",paste0("gamma",seq(0,nz)),paste0("beta",seq(0,nz+1)),"rx1","rx2","qx1")
  Amat
}
