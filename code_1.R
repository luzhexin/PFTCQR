rm(list=ls())
#install.packages("fda")
#install.packages('truncnorm')
#install.packages('actuar')
#install.packages('statmod')
#install.packages("ghyp")
library(fda)
library(truncnorm)
library(actuar)
library(statmod)
library(ghyp)


####################################
## Start Simulation 1 ##################
####################################

#set.seed(0)
N <-500
n<-200;  tau=0.25;
nobs <- 101;  d<-101; 
tobs <- seq(0,1,length=nobs);
tgrid <- seq(0,1,length=d)
M <- 150   

#--------------------Generate data--------------------------# 
##X(t)
for(iter in 1:N){
  xiM <- matrix(0,n,M)
  for(j in 1:M){
     xiM[,j] <- rnorm(n,0,1/j)
  }
  w <- 1:M;
  X <- matrix(0,n,nobs)
  for(i in 1:n){
     for(j in 1:nobs){
        X[i,j] <- sqrt(1/2)*sum(xiM[i,]*cos(w*pi*tobs[j])) 
        }
  }

##beta
   betafun=function(t) {
    betafun=0.2*(1+tau)*sin(pi*t/2)+sqrt(2)*cos(4*pi*t/2)
  }

  truebeta<-betafun(tgrid)
  muB <- (tobs[2]-tobs[1])*c(X%*%truebeta)  

##Z
  Z<-rnorm(n,0,0.5)
  truealpha=0.5+tau/5
  muA <- truealpha*Z

##Error
  E <- rnorm(n,0,1)-qnorm(tau,0,1)

##Response
  Y <- c(muA+muB)+E;
  T=Y;
  c=runif(n,0,15); 
  T[T>c]=0;

#--------------------Performing FPCA --------------------------#  
  W <- t(X);  
  nbasis<-150
  xbasis<-create.bspline.basis(rangeval=c(0,1),nbasis=nbasis)    
  xfd<-Data2fd(tgrid,W,xbasis)
  x.fd<-eval.fd(tgrid,xfd) 
  pcafd<-pca.fd(xfd,nharm=nbasis)
  eve1<-eval.fd(tgrid,pcafd$harmonics)
  U1<-inprod(xfd,pcafd$harmonics)
  m<-2;
  Um1<-U1[,1:m]; 
  xi<-Um1;  

#------------------Gibbs-----------------------#
chains=10000;burn_in=2000;
b_k=matrix(0,chains+1,m)  
sigma=matrix(0,chains+1,1)  
vn=matrix(0,chains+1,n)
alphait=matrix(0,chains+1,1) 

theta_tau0=rep(0,1); A_tau0=diag(1,1,1)
B_tau0=rep(0,m); B1_tau0=diag(1,m,m)
n0=0.5;s0=0.1;
sigma[1,]=1
vn[1,]=1
b_k[1,]=0.5
alphait[1,]=truealpha 
gamma=(1-2*tau)/(tau*(1-tau))
lambd2=2/(tau*(1-tau))


for(it in 2:(chains+1)){
#sample Y 
    Y=T
    for(i in 1:n){
    Y1=rtruncnorm(n,a=c,b=Inf,mean=Z[i]*alphait[it-1,]+xi[i,]%*%b_k[it-1,]+gamma*vn[it-1,i],sd=sqrt(sigma[it-1,]*lambd2*vn[it-1,i]));
    }
    id=(T==0);
    Y[id]=Y1[id];
  
#sample alpha
    A_tau1=0
    for(i in 1:n){
    A_tau1=A_tau1+Z[i]%*%t(Z[i])/(lambd2*sigma[it-1,]*vn[it-1,i])
    }  
    A_tau4=A_tau1+ginv(A_tau0)

    theta_tau1=0
    for(i in 1:n){
    theta_tau1=theta_tau1+(Y[i]-xi[i,]%*%b_k[it-1,]-gamma*vn[it-1,i])*Z[i]/(lambd2*sigma[it-1,]*vn[it-1,i])
    }   
    theta_tau4=ginv(A_tau4)%*%(theta_tau1+ginv(A_tau0)%*%theta_tau0)
    alphait[it,]=rnorm(1,theta_tau4,ginv(A_tau4))

#sample beta
    B_tau1=0
    for(i in 1:n){
    B_tau1=B_tau1+xi[i,]%*%t(xi[i,])/(lambd2*sigma[it-1,]*vn[it-1,i])
    }  
    B_tau4=B_tau1+ginv(B1_tau0)

    B1_tau1=0
    for(i in 1:n){
    B1_tau1=B1_tau1+(Y[i]-Z[i]*alphait[it,]-gamma*vn[it-1,i])*xi[i,]/(lambd2*sigma[it-1,]*vn[it-1,i])
    }   
    B1_tau4=ginv(B_tau4)%*%(B1_tau1+ginv(B1_tau0)%*%B_tau0)
    b_k[it,]=mvrnorm(1, B1_tau4,ginv(B_tau4))
  
#sample sigma
    lat=2*sum(vn[it-1,])+s0+sum((Y-Z*alphait[it,]-xi%*%b_k[it,]-gamma*vn[it-1,])^2/(lambd2*vn[it-1,]))
    sigma[it,]=sqrt(rinvgamma(1,shape=(3*n+n0)/2,rate=lat/2))

#sample vn   
   for(i in 1:n){
   chi=(Y[i]-Z[i]*alphait[it,]-xi[i,]%*%b_k[it,])^2/(lambd2*sigma[it,])
   psi=2/sigma[it,]+(gamma^2)/(lambd2*sigma[it,])
   } 
   vn[it,]=rgig(n,lambda=0.5,chi=chi,psi=psi)

#print(it)      
}
print(iter)
}



                                                                                                                                                                                                   