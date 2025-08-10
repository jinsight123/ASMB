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
EBFDR.func<-function(z,em,alpha1,alpha_1){
  emM<-matrix(0,3,3)
  emM[1,]<-em$pi
  emM[2,]<-em$mu
  emM[3,]<-em$sigma
  emD<-emM
  rankmu<-rank(emM[2,])
  emM[,rankmu]<-emD
  p_1<-emM[1,1]
  p0<-emM[1,2]
  p1<-emM[1,3]
  mu_1<-emM[2,1]
  mu0<-emM[2,2]
  mu1<-emM[2,3]
  sigma_1<-emM[3,1]
  sigma0<-emM[3,2]
  sigma1<-emM[3,3]
  f0<-function(x){dnorm(x,mu0,sigma0)}
  f1<-function(x){dnorm(x,mu1,sigma1)}
  f_1<-function(x){dnorm(x,mu_1,sigma_1)}
  funv_1<-function(x){(p_1*f_1(x))/(p_1*f_1(x)+p0*f0(x)+p1*f1(x))}
  funv1<-function(x){(p1*f1(x))/(p_1*f_1(x)+p0*f0(x)+p1*f1(x))}
  funv0<-function(x){(p0*f0(x))/(p_1*f_1(x)+p0*f0(x)+p1*f1(x))}
  D<-matrix(0,3,length(z))
  v0<-funv0(z)
  v1<-funv1(z)
  v_1<-funv_1(z)
  v01<-v0+v1
  v0_1<-v0+v_1
  Dn<-c(v_1>v0)*c(v_1>v1)+0
  Dp<-c(v1>v0)*c(v1>v_1)+0
  #rejn
  rejn<-v01[which(Dn==1)]
  ordern<-order(rejn)
  temp<-sort(rejn)
  for (ni in 1:length(rejn)) {
    rejn[ni]<-mean(temp[1:ni])
  }
  rejn<-c(rejn<=alpha_1)+0
  q<-c(rep(0,length(rejn)))
  q[ordern]<-rejn
  rejn<-q
  #rejp
  rejp<-v0_1[which(Dp==1)]
  orderp<-order(rejp)
  temp<-sort(rejp)
  for (pi in 1:length(rejp)) {
    rejp[pi]<-mean(temp[1:pi])
  }
  rejp<-c(rejp<=alpha1)+0
  q<-c(rep(0,length(rejp)))
  q[orderp]<-rejp
  rejp<-q
  for (di in 1:length(z)) {
    if(Dn[di]==1){D[3,di]<-rejn[1]
    rejn<-rejn[-1]}
    if(Dp[di]==1){D[1,di]<-rejp[1]
    rejp<-rejp[-1]}
  }
  D[2,]<-1-D[1,]-D[3,]
  y<-list(Decision=D)
  return (y)
}
G3mFDR.func<-function(z,em,alpha1,alpha_1){
  K<-length(z)
  iter.max<-1000
  a<-c(rep(0,K))
  b<-c(rep(0,K))
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  emM<-matrix(0,3,3)
  emM[1,]<-em$pi
  emM[2,]<-em$mu
  emM[3,]<-em$sigma
  emD<-emM
  rankmu<-rank(emM[2,])
  emM[,rankmu]<-emD
  p_1<-emM[1,1]
  p0<-emM[1,2]
  p1<-emM[1,3]
  mu_1<-emM[2,1]
  mu0<-emM[2,2]
  mu1<-emM[2,3]
  sigma_1<-emM[3,1]
  sigma0<-emM[3,2]
  sigma1<-emM[3,3]
  f0<-function(x){dnorm(x,mu0,sigma0)}
  f1<-function(x){dnorm(x,mu1,sigma1)}
  f_1<-function(x){dnorm(x,mu_1,sigma_1)}
  f<-function(x){p0*f0(x)+p1*f1(x)+p_1*f_1(x)}
  t1<-function(x){1-(p1*f1(x))/(p0*f0(x)+p1*f1(x)+p_1*f_1(x))}
  t_1<-function(x){1-(p_1*f_1(x))/(p0*f0(x)+p1*f1(x)+p_1*f_1(x))}
  N<-function(lambda1,lambda_1,z){
    T1<-t1(z)
    T_1<-t_1(z)
    ######################### from lambda to D
    cost.M<-matrix(0,3,K)
    cost.M[1,]<-((1-T_1)+lambda1*(T1-alpha1))
    cost.M[2,]<-1-T1+1-T_1
    cost.M[3,]<-((1-T1)+lambda_1*(T_1-alpha_1))
    for (i in 1:K) {
      cost.M[,i]<-(cost.M[,i]==min(cost.M[,i]))
    }
    D<-(cost.M==1)
    for (i in 1:K) {
      x<-c(D[,i])
      D[,i]<-c(0,0,0)
      r<-resample(which(x==1),1)
      D[r,i]<-1
    }
    ######################### from D to N
    N1<-(D[1,]%*%matrix(T1-alpha1))/K
    N_1<-(D[3,]%*%matrix(T_1-alpha_1))/K
    e<-matrix(c(N1,N_1),2,1)
    return(e)
  }
  T1<-t1(z)
  T_1<-t_1(z)
  T1<-sort(T1)
  T_1<-sort(T_1)
  for (i in 1:K) {
    a[i]<-mean(T1[1:i])-alpha1
    b[i]<-mean(T_1[1:i])-alpha_1
  }
  r1<-max(length(which(a<=0)),1) 
  r_1<-max(length(which(b<=0)),1)
  #the orbit of lambda
  lambda1.vec<-(1-alpha1)/(T1-alpha1)-1
  lambda_1.vec<-(1-alpha_1)/(T_1-alpha_1)-1
  lambda1<-lambda1.vec[r1]
  lambda_1<-lambda_1.vec[r_1]
  #########################
  counter<-0
  lambda1up<-lambda1.vec[r1+1]
  lambda_1up<-lambda_1.vec[r_1+1]
  while ((counter<iter.max)&(1-(N(lambda1,lambda_1,z)[1]<=0&N(lambda1,lambda_1,z)[2]<=0&N(lambda1up,lambda_1,z)[1]>0&N(lambda1,lambda_1up,z)[2]>0))) {
    counter<-counter+1
    # update r1 r_1 lambda1 lambda_1 lambda1up lambda_1up
    while((N(lambda1,lambda_1,z)[1]<=0)&(r1+1<=K)){
      r1<-r1+1
      lambda1<-lambda1.vec[r1]
    }
    r1<-max((r1-1),1)
    while((N(lambda1,lambda_1,z)[2]<=0)&(r_1+1<=K)){
      r_1<-r_1+1
      lambda_1<-lambda_1.vec[r_1]
    }
    r_1<-max((r_1-1),1)
    lambda1<-lambda1.vec[r1]
    lambda1up<-lambda1.vec[r1+1]
    lambda_1<-lambda_1.vec[r_1]
    lambda_1up<-lambda_1.vec[r_1+1]
  }
  T1<-t1(z)
  T_1<-t_1(z)
  ######################### from lambda to D
  cost.M<-matrix(0,3,K)
  cost.M[1,]<-((1-T_1)+lambda1*(T1-alpha1))
  cost.M[2,]<-1-T1+1-T_1
  cost.M[3,]<-((1-T1)+lambda_1*(T_1-alpha_1))
  for (i in 1:K) {
    cost.M[,i]<-(cost.M[,i]==min(cost.M[,i]))
  }
  D<-(cost.M==1)
  for (i in 1:K) {
    x<-c(D[,i])
    D[,i]<-c(0,0,0)
    r<-resample(which(x==1),1)
    D[r,i]<-1
  }
  y<-list(Decision=D)
  return (y)
}
findmunuOMT = function(K, propvec, cvvec, sigmavec, alpha1 = 0.05, alpha_1 = 0.05, muOMT=1000, nuOMT=1000, maxiter= 10000){
  out =   fZVmunu(K,  mus=seq(0,muOMT+500,500), nus=seq(0,nuOMT+500,500), Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec, maxiter=maxiter)
  net<-(out$nlev<=alpha_1)*(out$plev<=alpha1)
  mubest.index<-min(which(rowSums(net)>=1))
  nubest.index<-min(which(colSums(net)>=1))
  print(cbind(mubest.index, mus[mubest.index]))
  print(cbind(nubest.index, nus[nubest.index]))
  
  muOMT<-mus[mubest.index]
  nuOMT<-nus[nubest.index]
  
  mus = seq(max(muOMT-500,0),muOMT+500,100)
  nus = seq(max(nuOMT-500,0),nuOMT+500,100)
  out =   fZVmunu(K,  mus, nus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,maxiter=maxiter)
  net<-(out$nlev<=alpha_1)*(out$plev<=alpha1)
  mubest.index<-min(which(rowSums(net)>=1))
  nubest.index<-min(which(colSums(net)>=1))
  print(cbind(mubest.index, mus[mubest.index]))
  print(cbind(nubest.index, nus[nubest.index]))
  
  muOMT<-mus[mubest.index]
  nuOMT<-nus[nubest.index]
  
  mus = seq(max(muOMT-100,0),muOMT+100,10)
  nus = seq(max(nuOMT-100,0),nuOMT+100,10)
  out =   fZVmunu(K,  mus, nus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,maxiter=maxiter)
  net<-(out$nlev<=alpha_1)*(out$plev<=alpha1)
  mubest.index<-min(which(rowSums(net)>=1))
  nubest.index<-min(which(colSums(net)>=1))
  print(cbind(mubest.index, mus[mubest.index]))
  print(cbind(nubest.index, nus[nubest.index]))
  
  muOMT<-mus[mubest.index]
  nuOMT<-nus[nubest.index]
  
  mus = seq(max(muOMT-10,0),muOMT+10,1)
  nus = seq(max(nuOMT-10,0),nuOMT+10,1)
  out =   fZVmunu(K,  mus, nus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
  net<-(out$nlev<=alpha_1)*(out$plev<=alpha1)
  mubest.index<-min(which(rowSums(net)>=1))
  nubest.index<-min(which(colSums(net)>=1))
  print(cbind(mubest.index, mus[mubest.index]))
  print(cbind(nubest.index, nus[nubest.index]))
  
  muOMT<-mus[mubest.index]
  nuOMT<-nus[nubest.index]
  
  # mus = seq(max(muOMT-1,0),muOMT+1,0.1)
  # nus = seq(max(nuOMT-1,0),nuOMT+1,0.1)
  # out =   fZVmunu(K,  mus, nus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
  # net<-(out$nlev<=alpha_1)*(out$plev<=alpha1)
  # mubest.index<-min(which(rowSums(net)>=1))
  # nubest.index<-min(which(colSums(net)>=1))
  # print(cbind(mubest.index, mus[mubest.index]))
  # print(cbind(nubest.index, nus[nubest.index]))
  # 
  # muOMT<-mus[mubest.index]
  # nuOMT<-nus[nubest.index]
  # 
  # mus = seq(max(muOMT-0.1,0),muOMT+0.1,0.01)
  # nus = seq(max(nuOMT-0.1,0),nuOMT+0.1,0.01)
  # out =   fZVmunu(K,  mus, nus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
  # net<-(out$nlev<=alpha_1)*(out$plev<=alpha1)
  # mubest.index<-min(which(rowSums(net)>=1))
  # nubest.index<-min(which(colSums(net)>=1))
  # print(cbind(mubest.index, mus[mubest.index]))
  # print(cbind(nubest.index, nus[nubest.index]))
  # 
  # muOMT<-mus[mubest.index]
  # nuOMT<-nus[nubest.index]
  return(list(mu=muOMT,nu=nuOMT))
}
EBFDR.func<-function(z,em,alpha1,alpha_1){
  emM<-matrix(0,3,3)
  emM[1,]<-em$pi
  emM[2,]<-em$mu
  emM[3,]<-em$sigma
  emD<-emM
  rankmu<-rank(emM[2,])
  emM[,rankmu]<-emD
  p_1<-emM[1,1]
  p0<-emM[1,2]
  p1<-emM[1,3]
  mu_1<-emM[2,1]
  mu0<-emM[2,2]
  mu1<-emM[2,3]
  sigma_1<-emM[3,1]
  sigma0<-emM[3,2]
  sigma1<-emM[3,3]
  f0<-function(x){dnorm(x,mu0,sigma0)}
  f1<-function(x){dnorm(x,mu1,sigma1)}
  f_1<-function(x){dnorm(x,mu_1,sigma_1)}
  funv_1<-function(x){(p_1*f_1(x))/(p_1*f_1(x)+p0*f0(x)+p1*f1(x))}
  funv1<-function(x){(p1*f1(x))/(p_1*f_1(x)+p0*f0(x)+p1*f1(x))}
  funv0<-function(x){(p0*f0(x))/(p_1*f_1(x)+p0*f0(x)+p1*f1(x))}
  D<-matrix(0,3,length(z))
  v0<-funv0(z)
  v1<-funv1(z)
  v_1<-funv_1(z)
  v01<-v0+v1
  v0_1<-v0+v_1
  Dn<-c(v_1>v0)*c(v_1>v1)+0
  Dp<-c(v1>v0)*c(v1>v_1)+0
  #rejn
  rejn<-v01[which(Dn==1)]
  ordern<-order(rejn)
  temp<-sort(rejn)
  for (ni in 1:length(rejn)) {
    rejn[ni]<-mean(temp[1:ni])
  }
  rejn<-c(rejn<=alpha_1)+0
  q<-c(rep(0,length(rejn)))
  q[ordern]<-rejn
  rejn<-q
  #rejp
  rejp<-v0_1[which(Dp==1)]
  orderp<-order(rejp)
  temp<-sort(rejp)
  for (pi in 1:length(rejp)) {
    rejp[pi]<-mean(temp[1:pi])
  }
  rejp<-c(rejp<=alpha1)+0
  q<-c(rep(0,length(rejp)))
  q[orderp]<-rejp
  rejp<-q
  for (di in 1:length(z)) {
    if(Dn[di]==1){D[3,di]<-rejn[1]
    rejn<-rejn[-1]}
    if(Dp[di]==1){D[1,di]<-rejp[1]
    rejp<-rejp[-1]}
  }
  D[2,]<-1-D[1,]-D[3,]
  y<-list(Decision=D)
  return (y)
}



G3mFDR.func<-function(z,em,alpha1,alpha_1){
  K<-length(z)
  iter.max<-1000
  a<-c(rep(0,K))
  b<-c(rep(0,K))
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  emM<-matrix(0,3,3)
  emM[1,]<-em$pi
  emM[2,]<-em$mu
  emM[3,]<-em$sigma
  emD<-emM
  rankmu<-rank(emM[2,])
  emM[,rankmu]<-emD
  p_1<-emM[1,1]
  p0<-emM[1,2]
  p1<-emM[1,3]
  mu_1<-emM[2,1]
  mu0<-emM[2,2]
  mu1<-emM[2,3]
  sigma_1<-emM[3,1]
  sigma0<-emM[3,2]
  sigma1<-emM[3,3]
  f0<-function(x){dnorm(x,mu0,sigma0)}
  f1<-function(x){dnorm(x,mu1,sigma1)}
  f_1<-function(x){dnorm(x,mu_1,sigma_1)}
  f<-function(x){p0*f0(x)+p1*f1(x)+p_1*f_1(x)}
  t1<-function(x){1-(p1*f1(x))/(p0*f0(x)+p1*f1(x)+p_1*f_1(x))}
  t_1<-function(x){1-(p_1*f_1(x))/(p0*f0(x)+p1*f1(x)+p_1*f_1(x))}
  N<-function(lambda1,lambda_1,z){
    T1<-t1(z)
    T_1<-t_1(z)
    ######################### from lambda to D
    cost.M<-matrix(0,3,K)
    cost.M[1,]<-((1-T_1)+lambda1*(T1-alpha1))
    cost.M[2,]<-1-T1+1-T_1
    cost.M[3,]<-((1-T1)+lambda_1*(T_1-alpha_1))
    for (i in 1:K) {
      cost.M[,i]<-(cost.M[,i]==min(cost.M[,i]))
    }
    D<-(cost.M==1)
    for (i in 1:K) {
      x<-c(D[,i])
      D[,i]<-c(0,0,0)
      r<-resample(which(x==1),1)
      D[r,i]<-1
    }
    ######################### from D to N
    N1<-(D[1,]%*%matrix(T1-alpha1))/K
    N_1<-(D[3,]%*%matrix(T_1-alpha_1))/K
    e<-matrix(c(N1,N_1),2,1)
    return(e)
  }
  T1<-t1(z)
  T_1<-t_1(z)
  T1<-sort(T1)
  T_1<-sort(T_1)
  for (i in 1:K) {
    a[i]<-mean(T1[1:i])-alpha1
    b[i]<-mean(T_1[1:i])-alpha_1
  }
  r1<-max(length(which(a<=0)),1) 
  r_1<-max(length(which(b<=0)),1)
  #the orbit of lambda
  lambda1.vec<-(1-alpha1)/(T1-alpha1)-1
  lambda_1.vec<-(1-alpha_1)/(T_1-alpha_1)-1
  lambda1<-lambda1.vec[r1]
  lambda_1<-lambda_1.vec[r_1]
  #########################
  counter<-0
  lambda1up<-lambda1.vec[r1+1]
  lambda_1up<-lambda_1.vec[r_1+1]
  while ((counter<iter.max)&(1-(N(lambda1,lambda_1,z)[1]<=0&N(lambda1,lambda_1,z)[2]<=0&N(lambda1up,lambda_1,z)[1]>0&N(lambda1,lambda_1up,z)[2]>0))) {
    counter<-counter+1
    # update r1 r_1 lambda1 lambda_1 lambda1up lambda_1up
    while((N(lambda1,lambda_1,z)[1]<=0)&(r1+1<=K)){
      r1<-r1+1
      lambda1<-lambda1.vec[r1]
    }
    r1<-max((r1-1),1)
    while((N(lambda1,lambda_1,z)[2]<=0)&(r_1+1<=K)){
      r_1<-r_1+1
      lambda_1<-lambda_1.vec[r_1]
    }
    r_1<-max((r_1-1),1)
    lambda1<-lambda1.vec[r1]
    lambda1up<-lambda1.vec[r1+1]
    lambda_1<-lambda_1.vec[r_1]
    lambda_1up<-lambda_1.vec[r_1+1]
  }
  T1<-t1(z)
  T_1<-t_1(z)
  ######################### from lambda to D
  cost.M<-matrix(0,3,K)
  cost.M[1,]<-((1-T_1)+lambda1*(T1-alpha1))
  cost.M[2,]<-1-T1+1-T_1
  cost.M[3,]<-((1-T1)+lambda_1*(T_1-alpha_1))
  for (i in 1:K) {
    cost.M[,i]<-(cost.M[,i]==min(cost.M[,i]))
  }
  D<-(cost.M==1)
  for (i in 1:K) {
    x<-c(D[,i])
    D[,i]<-c(0,0,0)
    r<-resample(which(x==1),1)
    D[r,i]<-1
  }
  y<-list(Decision=D)
  return (y)
}
G3FDR.func<-function(zv,out,alpha1,alpha_1){
  K<-length(zv)
  ###################################################
  propvec = out$pi
  print(propvec)
  cvvec = out$mu
  print(cvvec)
  sigmavec = out$sigma
  print(sigmavec)
  #For est-OMT-FDR: search for the optimal mu given the mixture parameter
  mus = seq(4500,15000,500)
  murange = 4000
  muOMT = findmuOMT(K, propvec, cvvec, sigmavec, alternative = "less", alpha = alpha_1, mus = mus, murange=murange, maxiter= 10000) 
  nus = seq(4500,15000,500)
  nurange = 4000
  nuOMT = findmuOMT(K, propvec, cvvec, sigmavec, alternative = "greater", alpha = alpha1, mus = nus, murange=nurange, maxiter= 10000)
  maxK = ifelse(K>200, max(round(K*(1-max(propvec))),200),K)
  if(maxK*2>=K){
    munu<-findmunuOMT(K, propvec, cvvec, sigmavec, alpha1, alpha_1, muOMT, nuOMT, maxiter= 10000)
    muOMT<-munu$mu
    nuOMT<-munu$nu
  }
  #find est-OMT-FDR number of discoveries
  Pz = marg.dense2sidedalternative(zv,propvec, cvvec, sigmavec)
  Tnz<-scaled.null.dens(zv,propvec, cvvec, sigmavec, alternative = "less")$out/Pz 
  Tpz<-scaled.null.dens(zv,propvec, cvvec, sigmavec, alternative = "greater")$out/Pz
  Tnz<-sort(Tnz)
  Tpz<-sort(Tpz)
  az<-1-Tnz
  bz<-BZCpp(Tnz)
  cz<-1-Tpz
  dz<-BZCpp(Tpz)
  Rz<-az-muOMT*bz
  Lz<-cz-nuOMT*dz
  D<-matrix(0,3,K)
  ###
  Dnz = DCpp(Rz)
  Dpz = DCpp(Lz)
  if((sum(Dnz)+sum(Dpz))>K){
    nlap<-(sum(Dnz)+sum(Dpz))-K
    temp<-rep(0,nlap)
    for (i in 1:nlap) {
      temp[i]<-sum(Rz[1:(sum(Dnz)-c)])+sum(Lz[1:(sum(Dpz)+c-nlap)])
    }
    c<-which.max(temp)
    Dnz<-c(rep(1,(sum(Dnz)-c)),rep(0,(K-sum(Dnz)+c)))
    Dpz<-c(rep(1,(sum(Dpz)+c-nlap)),rep(0,(K-sum(Dpz)-c+nlap)))
  }
  ###
  nr<-sum(Dnz)
  if(nr!=0){
    nc<-sort(zv)[nr]
    D[3,]<-as.numeric(zv<=nc)
  }
  ###
  pr<-sum(Dpz)
  if(pr!=0){
    pc<-sort(zv,decreasing = T)[pr]
    D[1,]<-as.numeric(zv>=pc) 
  }
  D[2,]<-1-D[1,]-D[3,]
  y<-list(Decision=D)
  return (y)
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

zdFDR.func<-function(z,out,alpha){
  K<-length(z)
  emM<-matrix(0,3,3)
  emM[1,]<-out$pi
  emM[2,]<-out$mu
  emM[3,]<-out$sigma
  emD<-emM
  rankmu<-rank(emM[2,])
  emM[,rankmu]<-emD
  p_1<-emM[1,1]
  p0<-emM[1,2]
  p1<-emM[1,3]
  mu_1<-emM[2,1]
  mu0<-emM[2,2]
  mu1<-emM[2,3]
  sigma_1<-emM[3,1]
  sigma0<-emM[3,2]
  sigma1<-emM[3,3]
  ###################################################
  propvec = out$pi
  print(propvec)
  cvvec = out$mu
  print(cvvec)
  sigmavec = out$sigma
  print(sigmavec)
  u <- pnorm((z-mu0)/sigma0)
  A.set<-(u>0.25&u<0.75)*1
  R.set<-(u<0.25|u>0.75)*1
  mask.set<-A.set+R.set
  hat.dfdr<-1
  t<-1
  reflection<-function(u){
    (1.5-u)*(u>0.5)+(0.5-u)*(u<=0.5)
  }
  leaveone<-rep(1,K)
  while (hat.dfdr>alpha&t<=K) {
    t<-t+1
    u.mask<-mask.set*reflection(u)+(1-mask.set)*u
    z.mask<-qnorm(u.mask,mean = mu0, sd= sigma0)
    z.mix<-(z<0)*(pmin(z,z.mask))+(z>=0)*(pmax(z,z.mask))
    tmin<-scaled.null.dens(z.mix,propvec, cvvec, sigmavec, alternative = "min")$out
    tminle<-tmin*leaveone
    mask.set[which.max(tminle)]<-0
    A.set[which.max(tminle)]<-0
    R.set[which.max(tminle)]<-0
    mask.set<-A.set+R.set
    A.numb<-sum(A.set)
    R.numb<-sum(R.set)
    hat.dfdr<-(1+A.numb)/max(R.numb,1)
    # hat.dfdr
    leaveone[which.max(tminle)]<-0
    # sum(leaveone)
    # t
  }
  D<-matrix(0,3,K)
  D[1,]<-R.set*(z>0)
  D[3,]<-R.set*(z<0)
  D[2,]<-1-D[1,]-D[3,]
  y<-list(Decision=D,iter=t, A=A.numb, R=R.numb, L=sum(leaveone))
  return (y)
}
###
mTFDR.func<-function(z,em,alpha){
  K<-length(z)
  ###################################################
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  emM<-matrix(0,3,3)
  emM[1,]<-em$pi
  emM[2,]<-em$mu
  emM[3,]<-em$sigma
  emD<-emM
  rankmu<-rank(emM[2,])
  emM[,rankmu]<-emD
  p_1<-emM[1,1]
  p0<-emM[1,2]
  p1<-emM[1,3]
  mu_1<-emM[2,1]
  mu0<-emM[2,2]
  mu1<-emM[2,3]
  sigma_1<-emM[3,1]
  sigma0<-emM[3,2]
  sigma1<-emM[3,3]
  ###
  f0<-function(x){dnorm(x,mu0,sigma0)}
  f1<-function(x){dnorm(x,mu1,sigma1)}
  f_1<-function(x){dnorm(x,mu_1,sigma_1)}
  f<-function(x){p0*f0(x)+p1*f1(x)+p_1*f_1(x)}
  t1<-function(x){1-(p1*f1(x))/(p0*f0(x)+p1*f1(x)+p_1*f_1(x))}
  t_1<-function(x){1-(p_1*f_1(x))/(p0*f0(x)+p1*f1(x)+p_1*f_1(x))}
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  N<-function(lambda,z){
    T1<-t1(z)
    T_1<-t_1(z)
    ######################### from lambda to D
    cost.M<-matrix(0,3,K)
    cost.M[1,]<-((1-T_1)+lambda*(T1-alpha))
    cost.M[2,]<-1-T1+1-T_1
    cost.M[3,]<-((1-T1)+lambda*(T_1-alpha))
    for (i in 1:K) {
      cost.M[,i]<-(cost.M[,i]==min(cost.M[,i]))
    }
    D<-(cost.M==1)
    for (i in 1:K) {
      x<-c(D[,i])
      D[,i]<-c(0,0,0)
      r<-resample(which(x==1),1)
      D[r,i]<-1
    }
    ######################### from D to N
    N1<-(D[1,]%*%matrix(T1-alpha))/K
    N_1<-(D[3,]%*%matrix(T_1-alpha))/K
    e<-N1+N_1
    return(e)
  }
  iter.max<-1000
  D<-matrix(0,3,K)
  ###
  T1<-t1(z)
  T_1<-t_1(z)
  Tmin<-pmin(T1, T_1)
  Tmin<-sort(Tmin)
  a<-rep(1,K)
  for (i in 1:K) {
    a[i]<-mean(Tmin[1:i])-alpha
  }
  r1<-max(length(which(a<=0)),1) 
  #the orbit of lambda
  lambda.vec<-(1-alpha)/(Tmin-alpha)-1
  lambda<-lambda.vec[r1]
  #########################
  counter<-0
  lambdaup<-lambda.vec[r1+1]
  while ((counter<iter.max)&(1-(N(lambda,z)<=0&N(lambdaup,z)>0))) {
    counter<-counter+1
    # update r1 r_1 lambda1 lambda_1 lambda1up lambda_1up
    while((N(lambda,z)<=0)&(r1+1<=K)){
      r1<-r1+1
      lambda<-lambda.vec[r1]
    }
    r1<-max((r1-1),1)
    lambda<-lambda.vec[r1]
    lambdaup<-lambda.vec[r1+1]
  }
  T1<-t1(z)
  T_1<-t_1(z)
  ######################### from lambda to D
  cost.M<-matrix(0,3,K)
  cost.M[1,]<-((1-T_1)+lambda*(T1-alpha))
  cost.M[2,]<-1-T1+1-T_1
  cost.M[3,]<-((1-T1)+lambda*(T_1-alpha))
  for (i in 1:K) {
    cost.M[,i]<-(cost.M[,i]==min(cost.M[,i]))
  }
  D<-(cost.M==1)
  for (i in 1:K) {
    x<-c(D[,i])
    D[,i]<-c(0,0,0)
    r<-resample(which(x==1),1)
    D[r,i]<-1
  }
  y<-list(Decision=D)
  return (y)
}
###
dEBFDR.func<-function(z,em,alpha){
  K<-length(z)
  emM<-matrix(0,3,3)
  emM[1,]<-em$pi
  emM[2,]<-em$mu
  emM[3,]<-em$sigma
  emD<-emM
  rankmu<-rank(emM[2,])
  emM[,rankmu]<-emD
  p_1<-emM[1,1]
  p0<-emM[1,2]
  p1<-emM[1,3]
  mu_1<-emM[2,1]
  mu0<-emM[2,2]
  mu1<-emM[2,3]
  sigma_1<-emM[3,1]
  sigma0<-emM[3,2]
  sigma1<-emM[3,3]
  ###
  f0<-function(x){dnorm(x,mu0,sigma0)}
  f1<-function(x){dnorm(x,mu1,sigma1)}
  f_1<-function(x){dnorm(x,mu_1,sigma_1)}
  funv_1<-function(x){(p_1*f_1(x))/(p_1*f_1(x)+p0*f0(x)+p1*f1(x))}
  funv1<-function(x){(p1*f1(x))/(p_1*f_1(x)+p0*f0(x)+p1*f1(x))}
  funv0<-function(x){(p0*f0(x))/(p_1*f_1(x)+p0*f0(x)+p1*f1(x))}
  ###################################################
  D<-matrix(0,3,K)
  v0<-funv0(z)
  v1<-funv1(z)
  v_1<-funv_1(z)
  v01<-v0+v1
  v0_1<-v0+v_1
  Dn<-c(v_1>v0)*c(v_1>v1)+0
  Dp<-c(v1>v0)*c(v1>v_1)+0
  maxrej<-sum(Dn)+sum(Dp)
  ###
  psi<-pmin(v01,v0_1)
  orderpsi<-order(psi)
  temp<-sort(psi)
  rejpsi<-rep(1,length(psi))
  for (ni in 1:maxrej) {
    rejpsi[ni]<-mean(temp[1:ni])
  }
  rejpsi<-c(rejpsi<=alpha)+0
  q<-c(rep(0,length(rejpsi)))
  q[orderpsi]<-rejpsi
  rejpsi<-q
  D[1,]<-rejpsi*Dp
  D[3,]<-rejpsi*Dn
  D[2,]<-1-D[1,]-D[3,]
  y<-list(Decision=D)
  return (y)
}