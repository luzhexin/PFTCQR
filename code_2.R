rm(list=ls())
#install.packages("fda")
#install.packages('truncnorm')
#install.packages('actuar')
#install.packages('statmod')
#install.packages("ghyp")
#install.packages("corpcor")
library(fda)
library(truncnorm)
library(actuar)
library(statmod)
library(ghyp)
library(corpcor) 

####################################
## Start Simulation 2 ##################
####################################

#set.seed(0)
N <-500
N_subj=n<-200;  
tau=0.25;
dim_image <- 300
K_x = 3

for(iter in 1:N){
eigen_vector = array(0,c(K_x,dim_image*dim_image))
eigen_vector[3,] = c(rep(c(rep(1,dim_image/3),rep(0,dim_image*2/3)),dim_image/3),rep(0,(dim_image*2/3)*dim_image))
eigen_vector[2,] = c(rep(0,(dim_image*1/3)*dim_image),rep(c(rep(0,dim_image/3),rep(1,dim_image/3),rep(0,dim_image/3)),dim_image/3),rep(0,(dim_image*1/3)*dim_image))
eigen_vector[1,] = c(rep(0,(dim_image*2/3)*dim_image),rep(c(rep(0,dim_image*2/3),rep(1,dim_image*1/3)),dim_image*1/3))

####standarize the eigenimage
for(jj in 1:K_x){
  eigen_vector[jj,] = eigen_vector[jj,]/sqrt(sum(eigen_vector[jj,]^2))
}

####generate the eigen_score 
eigen_score = array(0,c(N_subj,K_x))
eigen_var = rep(0,K_x)
eigen_var[1:K_x] = 0.5^((1:K_x) - 1)
eigen_sd = sqrt(eigen_var)
for(jj in 1:K_x){
  eigen_score[,jj] = rnorm(N_subj, mean = 0, sd = eigen_sd[jj])
}

####generate the images, X
true.funcs = eigen_score%*%eigen_vector  
mean.funcs = apply(true.funcs,2,mean)
for(ii in 1:N_subj){
  true.funcs[ii,] = true.funcs[ii,] - mean.funcs
}

####fast svd to get eigenvalues and eigenvectors 
eigen_svd = fast.svd(t(true.funcs),tol = 0.0001)
eigenimage_est = t(eigen_svd$u)
ind_rev = rep(0,K_x)
for(jj in 1:K_x){
  if(sum(eigenimage_est[jj,]*eigen_vector[jj,])<0){
    eigenimage_est[jj,] = - eigenimage_est[jj,]
    ind_rev[jj] = 1
  }
  
}

#####Response
beta_true <- c(0.5+tau/8,0.5+tau/8,-0.5+tau/8)
true_Beta = beta_true[1]*eigen_vector[1,] + beta_true[2]*eigen_vector[2,] + beta_true[3]*eigen_vector[3,] 
E <- rnorm(N_subj,0,0.5)-qnorm(tau,0,0.5)
truealpha=1+tau/8
Z<-rnorm(n,0,1)
xi=eigen_score
outcomes <- sapply(1:N_subj, function(u) sum(true.funcs[u,]*true_Beta))+truealpha*Z
Y <- outcomes+E;
T=Y;
c=runif(n,0,15); 
T[T>c]=0;
eve1=eigen_svd$u 
m=3

#-------------------  Gibbs  -----------------------#
chains=10000;burn_in=2000;
b_k1=matrix(0,chains+1,1) 
b_k2=matrix(0,chains+1,1) 
b_k3=matrix(0,chains+1,1) 
b_k=cbind(b_k1, b_k2, b_k3) 
alphait=matrix(0,chains+1,1) 
sigma=matrix(0,chains+1,1)  
vn=matrix(0,chains+1,n)
  
theta_tau1=rep(0,1); A1_tau=diag(10,1,1)
theta_tau2=rep(0,1); A2_tau=diag(10,1,1)
theta_tau3=rep(0,1); A3_tau=diag(10,1,1)
B_tau0=rep(0,1); B1_tau0=diag(10,1,1)
n0=0.5;s0=1 
sigma[1,]=0.001
vn[1,]=0.5
b_k1[1,]=beta_true[1]
b_k2[1,]=beta_true[2]
b_k3[1,]=beta_true[3]
b_k[1,]=cbind(b_k1[1,], b_k2[1,], b_k3[1,]) 
alphait[1,]=truealpha 
gamma=(1-2*tau)/(tau*(1-tau))
lambd2=2/(tau*(1-tau))
  

for(it in 2:(chains+1)){
#sample Y
   Y=T;
    for(i in 1:n){
      Y1=rtruncnorm(n,a=c,b=Inf,mean=Z[i]*alphait[it-1,]+xi[i,]%*%b_k[it-1,]+gamma*vn[it-1,i],sd=sqrt(sigma[it-1,]*lambd2*vn[it-1,i]));
    }
    id=(T==0);
    Y[id]=Y1[id];

#sample alpha
    B_tau1=0
    for(i in 1:n){
    B_tau1=B_tau1+Z[i]%*%t(Z[i])/(lambd2*sigma[it-1,]*vn[it-1,i])
    }  
    B_tau4=B_tau1+ginv(B1_tau0)

    B1_tau1=0
    for(i in 1:n){
    B1_tau1=B1_tau1+((Y[i]-xi[i,]%*%b_k[it-1,]-gamma*vn[it-1,i])*Z[i]/(lambd2*sigma[it-1,]*vn[it-1,i]))
    }   
    B1_tau4=ginv(B_tau4)%*%(B1_tau1+ginv(B1_tau0)%*%B_tau0)
    alphait[it,]=rnorm(1, B1_tau4,ginv(B_tau4))
   
#sample b1
    A1_tau1=0
    for(i in 1:n){
    A1_tau1=A1_tau1+xi[i,1]%*%t(xi[i,1])/(lambd2*sigma[it-1,]*vn[it-1,i])
    }  
    A1_tau4=A1_tau1+ginv(A1_tau)
    
    theta_tau1=0
    for(i in 1:n){
    theta_tau1=theta_tau1+((Y[i]-Z[i]*alphait[it,]-gamma*vn[it-1,i])*xi[i,1]/(lambd2*sigma[it-1,]*vn[it-1,i]))
    }   
    theta_tau4=ginv(A1_tau4)%*%(theta_tau1+ginv(A1_tau)%*%theta_tau1)
    b_k1[it,]=rnorm(1,theta_tau4,ginv(A1_tau4))

#sample b2
    A2_tau1=0
    for(i in 1:n){
    A2_tau1=A2_tau1+xi[i,2]%*%t(xi[i,2])/(lambd2*sigma[it-1,]*vn[it-1,i])
    }  
    A2_tau4=A2_tau1+ginv(A2_tau)
    
    theta_tau2=0
    for(i in 1:n){
    theta_tau2=theta_tau2+((Y[i]-Z[i]*alphait[it,]-gamma*vn[it-1,i])*xi[i,2]/(lambd2*sigma[it-1,]*vn[it-1,i]))
    }   
    theta_tau24=ginv(A2_tau4)%*%(theta_tau2+ginv(A2_tau)%*%theta_tau2)
    b_k2[it,]=rnorm(1,theta_tau24,ginv(A2_tau4))

#sample b3
    A3_tau1=0
    for(i in 1:n){
    A3_tau1=A3_tau1+xi[i,3]%*%t(xi[i,3])/(lambd2*sigma[it-1,]*vn[it-1,i])
    }  
    A3_tau4=A3_tau1+ginv(A3_tau)
    
    theta_tau3=0
    for(i in 1:n){
    theta_tau3=theta_tau3+((Y[i]-Z[i]*alphait[it,]-gamma*vn[it-1,i])*xi[i,3]/(lambd2*sigma[it-1,]*vn[it-1,i]))
    }   
    theta_tau34=ginv(A3_tau4)%*%(theta_tau3+ginv(A3_tau)%*%theta_tau3)
    b_k3[it,]=rnorm(1,theta_tau34,ginv(A3_tau4))
    b_k[it,]=cbind(b_k1[it,], b_k2[it,], b_k3[it,]) 

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

