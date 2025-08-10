# First, install the 'mixtools' package if it is not already installed
if (!requireNamespace("mixtools", quietly = TRUE)) {
  install.packages("mixtools")
}
# Load the 'mixtools' package
library(mixtools)
# %% [code]
Zrand2sidedalternative = function(K,propvec, cvvec, sigmavec){ #sample from a mixture model
  uvec = runif(K)
  J = length(propvec)
  ret = numeric(K)
  
  Hnull = which(uvec<propvec[1])
  ret[Hnull] = rnorm(length(Hnull), mean = cvvec[1], sd = sigmavec[1]) 
  
  for (j in 2:J){
    if (j<J){
      Hj = which (uvec>= sum(propvec[1:(j-1)]) & uvec<sum(propvec[1:j])) 
    }
    if (j==J){
      Hj = which (uvec>=sum(propvec[1:(j-1)]))
    }
    ret[Hj] = rnorm(length(Hj), mean = cvvec[j], sd = sigmavec[j]) 
  }
  ret
}



marg.dense2sidedalternative = function (z,propvec, cvvec, sigmavec){ #compute the marginal density of a  mixture model 
  J = length(propvec)
  out=rep(0, length(z))
  for (j in 1:J){
    out = out + propvec[j]*dnorm(z, mean = cvvec[j], sd = sigmavec[j])  
  }  
  return(out)
  
}

scaled.null.dens = function(z, propvec, cvvec, sigmavec, alternative = c("min", "two.sided", "less", "greater")){ 
  J = length(propvec)
  emM<-matrix(0,3,J)
  emM[1,]<-propvec
  emM[2,]<-cvvec
  emM[3,]<-sigmavec
  emD<-emM
  rankmu<-rank(emM[2,])
  emM[,rankmu]<-emD
  out=rep(0, length(z))
  outup=rep(0, length(z))
  outdown=rep(0, length(z))
  if (alternative =="min"){
    outup = emM[1,3]*dnorm(z, mean = emM[2,3], sd = emM[3,3])
    outdown = emM[1,1]*dnorm(z, mean = emM[2,1], sd = emM[3,1])  
    tmin=pmin(outup,outdown)+emM[1,2]*dnorm(z, mean = emM[2,2], sd = emM[3,2])
    minindex=(outup<=outdown)*1
    out<-list(out=tmin,minindex=minindex)
    return(out)
  }
  if (alternative =="less"){
    out = emM[1,3]*dnorm(z, mean = emM[2,3], sd = emM[3,3]) + emM[1,2]*dnorm(z, mean = emM[2,2], sd = emM[3,2])
    minindex=rep(0,length(z))
    out<-list(out=out,minindex=minindex)
    return(out)
  }  
  
  if (alternative =="greater"){
    out = emM[1,1]*dnorm(z, mean = emM[2,1], sd = emM[3,1]) + emM[1,2]*dnorm(z, mean = emM[2,2], sd = emM[3,2])
    minindex=rep(0,length(z))
    out<-list(out=out,minindex=minindex)
    return(out)
  }  
  if (alternative =="two.sided"){
    out = emM[1,2]*dnorm(z, mean = emM[2,2], sd = emM[3,2])
    minindex=rep(0,length(z))
    out<-list(out=out,minindex=minindex)
    return(out)
  }
}


fZV = function(K,  mus, Rreg, Cpp, propvec, cvvec,sigmavec, alternative = c("min", "two.sided", "less", "greater"), maxiter=10000){#in propvec and cvvec, respectively the mean and  proportion of zero, pos, neg 
  lev = pow=ev=maxr=numeric(length(mus))
  maxK = ifelse(K>200, max(round(K*(1-max(propvec))),200),K)
  for (iter in 1:maxiter){
    z = Zrand2sidedalternative(K, propvec, cvvec, sigmavec) #if sample from the mixture density with parameters propvec and cvvec and sigmavec
    Pz = marg.dense2sidedalternative(z,propvec, cvvec, sigmavec)
    Tz = scaled.null.dens(z,propvec, cvvec, sigmavec, alternative =alternative )$out/Pz 
    Tz = sort(Tz)
    az = 1-Tz[1:maxK]
    bz = BZCpp(Tz[1:maxK])
    for (mui in 1:length(mus)){
      mu=mus[mui]
      Rz = az-mu*bz
      Dz = DCpp(Rz)
      lev[mui] = lev[mui] + sum(bz[Dz==1])
      pow[mui] = pow[mui] + sum(az[Dz==1])
      maxr[mui] = max(maxr[mui],sum(Dz==1))
    } 
  }
  return(list(lev = lev/maxiter,pow = pow/maxiter))
  
}
library(Rcpp)
cppFunction('NumericVector BZCpp(NumericVector Tz){
            int K=Tz.size();
            NumericMatrix S(K,K);
            std::fill(S.begin(), S.end(), 0.);
            NumericVector b(K,0.);
            S(0,0)=1;
            for (int i=1;i<K;i++){
            S(i,0) = S(i-1,0)*(1-Tz[i-1]);
            for (int v=1;v<=i;v++){
            S(i,v) = S(i-1,v)*(1-Tz[i-1])+S(i-1,v-1)*Tz[i-1];
            }
            }
            
            b[0] = S(0,0)*Tz[0];
            for (int i=1;i<K;i++){
            for (int j = 0; j <= i;j++){
            b[i] += S(i,j)* ((i-j)*Tz[i] - j*(1-Tz[i]))/(i*(i+1));
            }
            }
            
            return b;
            }')



cppFunction('IntegerVector DCpp(NumericVector R) {
            int K=R.size();
            IntegerVector D(K,0);
            double cR=0.;
            
            for (int i =0;i<K;i++){
            cR=0.;
            for (int j = i;j<K;j++){
            cR += R[j];
            if (cR>0) {D[i]=1;break;}
            }
            if (D[i]==0) break;
            } 
            return D;
            }')  



findmuOMT = function(K, propvec, cvvec, sigmavec, alternative = c("min", "two.sided", "less", "greater"), alpha = 0.05, mus = seq(100,6100,500), murange=1000, maxiter= 10000){
  
  out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative, maxiter=maxiter)
  
  mbest = which.min(abs(out$lev-alpha))
  #if mbest equals 1 start again with all mu's below min(mus), if mbest equals length(mus) start again with all mu's above max(mus)
  while(mbest==1 | mbest== length(mus)){
    if (mbest==1){mus = seq(max(0, min(mus)-murange), min(mus)+500,100)}
    if (mbest==length(mus)){mus = seq(max(0,max(mus)-500), max(mus)+murange,100)}
    out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
    mbest = which.min(abs(out$lev-alpha))
  }
  print(cbind(mbest, mus[mbest]))
  
  mus = seq(max(mus[mbest]-500,0),mus[mbest]+500,100)
  out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
  mbest = which.min(abs(out$lev-alpha))
  print(cbind(mbest, mus[mbest], out$lev[mbest]))
  
  mus = seq(max(mus[mbest]-100,0),mus[mbest]+100,10)
  out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
  mbest = which.min(abs(out$lev-alpha))
  print(cbind(mbest, mus[mbest], out$lev[mbest]))
  
  mus = seq(max(mus[mbest]-10,0),mus[mbest]+10,1)
  out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
  mbest = which.min(abs(out$lev-alpha))
  print(cbind(mbest, mus[mbest], out$lev[mbest]))
  
  # mus = seq(max(mus[mbest]-1,0),mus[mbest]+1,0.1)
  # out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
  # mbest = which.min(abs(out$lev-alpha))
  # print(cbind(mbest, mus[mbest], out$lev[mbest]))
  # 
  # mus = seq(max(mus[mbest]-0.1,0),mus[mbest]+0.1,0.01)
  # out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
  # mbest = which.min(abs(out$lev-alpha))
  # print(cbind(mbest, mus[mbest], out$lev[mbest]))
  
  return(mus[mbest])
  
}

TFDR.func<-function(zv,out,alpha){
  K<-length(zv)
  ###################################################
  propvec = out$pi
  print(propvec)
  cvvec = out$mu
  print(cvvec)
  sigmavec = out$sigma
  print(sigmavec)
  #For est-OMT-FDR: search for the optimal mu given the mixture parameter
  mus = seq(100,6100,500)
  murange = 1000
  muOMT = findmuOMT(K, propvec, cvvec, sigmavec, alternative = "min", alpha = alpha, mus = mus, murange=murange, maxiter= 10000) 
  #find est-OMT-FDR number of discoveries
  Pz = marg.dense2sidedalternative(zv,propvec, cvvec, sigmavec)
  lfdr_1<-scaled.null.dens(zv,propvec, cvvec, sigmavec, alternative = "min")$out/Pz 
  D<-matrix(0,3,K)
  t_1<-lfdr_1
  order_1<-order(t_1)
  t_1<-sort(t_1)
  a<-1-t_1
  b<-BZCpp(t_1)
  r_1<-a-muOMT*b
  m<-c(rep(0,K))
  m[K]<-max(0,r_1[K])
  for (i in (K-1):1) {
    m[i]<-max(0,m[i+1]+r_1[i])
  }
  delta<-(cumsum(sign(m))==seq(1,K,1))
  p<-c(rep(0,K))
  p[order_1]<-delta
  delta<-p
  minindex<-scaled.null.dens(zv,propvec, cvvec, sigmavec, alternative = "min")$minindex
  D[1,]<-delta*(1-minindex)###+++
  D[3,]<-delta*minindex###---
  D[2,]<-1-D[1,]-D[3,]
  y<-list(Decision=D)
  return (y)
}
TFDR.fun<-function(zv,out,muOMT){
  K<-length(zv)
  ###################################################
  propvec = out$pi
  cvvec = out$mu
  sigmavec = out$sigma
  #find est-OMT-FDR number of discoveries
  Pz = marg.dense2sidedalternative(zv,propvec, cvvec, sigmavec)
  lfdr_1<-scaled.null.dens(zv,propvec, cvvec, sigmavec, alternative = "min")$out/Pz 
  D<-matrix(0,3,K)
  t_1<-lfdr_1
  order_1<-order(t_1)
  t_1<-sort(t_1)
  a<-1-t_1
  b<-BZCpp(t_1)
  r_1<-a-muOMT*b
  m<-c(rep(0,K))
  m[K]<-max(0,r_1[K])
  for (i in (K-1):1) {
    m[i]<-max(0,m[i+1]+r_1[i])
  }
  delta<-(cumsum(sign(m))==seq(1,K,1))
  p<-c(rep(0,K))
  p[order_1]<-delta
  delta<-p
  minindex<-scaled.null.dens(zv,propvec, cvvec, sigmavec, alternative = "min")$minindex
  D[1,]<-delta*(1-minindex)###+++
  D[3,]<-delta*minindex###---
  D[2,]<-1-D[1,]-D[3,]
  y<-list(Decision=D)
  return (y)
}
mu.vec<-seq(1.8,3,0.2)
# muOMT.vec<-c(1486, 2340,3218,3698,3837,3685,3321)
ndots<-length(mu.vec)
n.repeat<-100
##
DDTP<-matrix(0,ndots,n.repeat)
DDFDP<-matrix(0,ndots,n.repeat)
##
t4TP<-matrix(0,ndots,n.repeat)
t4FDP<-matrix(0,ndots,n.repeat)
##
t10TP<-matrix(0,ndots,n.repeat)
t10FDP<-matrix(0,ndots,n.repeat)
##
chiTP<-matrix(0,ndots,n.repeat)
chiFDP<-matrix(0,ndots,n.repeat)
for (i.parameter in 1:ndots) {
  K<-5000
  p1<-0.1
  p_1<-0.1
  p0<-1-p1-p_1
  mu0<-0
  mu1<-mu.vec[i.parameter]
  mu_1<--mu1
  sigma0<-1
  sigma1<-1
  sigma_1<-1
  alpha<-0.05
  ###################################################
  propvec = c(0.1,0.8,0.1)
  cvvec = c(mu_1,0,mu1)
  sigmavec = c(1,1,1)
  em<-list(pi=propvec,mu=cvvec,sigma=sigmavec)
  #For est-OMT-FDR: search for the optimal mu given the mixture parameter
  mus = seq(100,6100,500)
  murange = 1000
  muOMT = findmuOMT(K, propvec, cvvec, sigmavec, alternative = "min", alpha = alpha, mus = mus, murange=murange, maxiter= 10000)
  # muOMT<-muOMT.vec[i.parameter]
  for (j.repeat in 1:n.repeat) {
    z<-c(rep(0,K))
    components.random <- sample(1:3,prob=c(p1,p0,p_1),size=K,replace=TRUE)
    mus <- c(mu1,mu0,mu_1)
    sds <- c(sigma1,sigma0,sigma_1)
    z <- rnorm(n=K,mean=mus[components.random],sd=sds[components.random])
    #####################################
    index<-matrix(0,3,K)
    for (i in 1:K) {
      if(components.random[i]==1){index[1,i]<-1}
      if(components.random[i]==2){index[2,i]<-1}
      if(components.random[i]==3){index[3,i]<-1}
    }
    theta1<-index[1,]
    theta_1<-index[3,]
    ###################################
    TFDR<-TFDR.fun(z,em,muOMT)
    delta1<-TFDR$Decision[1,]
    delta_1<-TFDR$Decision[3,]
    DDTP[i.parameter,j.repeat]<-sum(theta1*delta1)+sum(theta_1*delta_1)
    DDFDP[i.parameter,j.repeat]<-(sum((1-theta1)*delta1)+sum((1-theta_1)*delta_1))/(max(1,(sum(delta1)+sum(delta_1))))
    ###################################
    z<-c(rep(0,K))
    for (i in 1:K) {
      if(components.random[i]==1){
        z[i]<-mu1+rt(1,4)/sqrt(2)
      }
      if(components.random[i]==2){
        z[i]<-rnorm(1)
      }
      if(components.random[i]==3){
        z[i]<-mu_1+rt(1,4)/sqrt(2)
      }
    }
    #####################################
    TFDR<-TFDR.fun(z,em,muOMT)
    delta1<-TFDR$Decision[1,]
    delta_1<-TFDR$Decision[3,]
    t4TP[i.parameter,j.repeat]<-sum(theta1*delta1)+sum(theta_1*delta_1)
    t4FDP[i.parameter,j.repeat]<-(sum((1-theta1)*delta1)+sum((1-theta_1)*delta_1))/(max(1,(sum(delta1)+sum(delta_1))))
    ###################################
    z<-c(rep(0,K))
    for (i in 1:K) {
      if(components.random[i]==1){
        z[i]<-mu1+rt(1,10)/sqrt(5/4)
      }
      if(components.random[i]==2){
        z[i]<-rnorm(1)
      }
      if(components.random[i]==3){
        z[i]<-mu_1+rt(1,10)/sqrt(5/4)
      }
    }
    ###################################
    mTFDR<-TFDR.fun(z,em,muOMT)
    delta1<-mTFDR$Decision[1,]
    delta_1<-mTFDR$Decision[3,]
    t10TP[i.parameter,j.repeat]<-sum(theta1*delta1)+sum(theta_1*delta_1)
    t10FDP[i.parameter,j.repeat]<-(sum((1-theta1)*delta1)+sum((1-theta_1)*delta_1))/(max(1,(sum(delta1)+sum(delta_1))))
    ###################################
    z<-c(rep(0,K))
    for (i in 1:K) {
      if(components.random[i]==1){
        z[i]<-mu1+(rchisq(1,4)-4)/sqrt(8)
      }
      if(components.random[i]==2){
        z[i]<-rnorm(1)
      }
      if(components.random[i]==3){
        z[i]<-mu_1-(rchisq(1,4)-4)/sqrt(8)
      }
    }
    ###################################
    dEBFDR<-TFDR.fun(z,em,muOMT)
    delta1<-dEBFDR$Decision[1,]
    delta_1<-dEBFDR$Decision[3,]
    chiTP[i.parameter,j.repeat]<-sum(theta1*delta1)+sum(theta_1*delta_1)
    chiFDP[i.parameter,j.repeat]<-(sum((1-theta1)*delta1)+sum((1-theta_1)*delta_1))/(max(1,(sum(delta1)+sum(delta_1))))
  }
}
TFDR.ETP<-rowMeans(DDTP)
TFDR.alpha<-rowMeans(DDFDP)
###
t4.ETP<-rowMeans(t4TP,na.rm = T)
t4.alpha<-rowMeans(t4FDP,na.rm = T)
###
t10.ETP<-rowMeans(t10TP)
t10.alpha<-rowMeans(t10FDP)
###
chi.ETP<-rowMeans(chiTP)
chi.alpha<-rowMeans(chiFDP)
result<-matrix(0,8,ndots)
result[1,]<-TFDR.ETP
result[2,]<-TFDR.alpha
result[3,]<-t4.ETP
result[4,]<-t4.alpha
result[5,]<-t10.ETP
result[6,]<-t10.alpha
result[7,]<-chi.ETP
result[8,]<-chi.alpha
result
write.csv(result,"onnt.txt",row.names=FALSE)