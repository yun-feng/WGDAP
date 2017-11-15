
g<-rep(c(1,2),each=10)


a_m<-c(rep(-1,40),rep(0,20),rep(-1,40))
b_m<-c(rep(0,40),rep(1,10),rep(-1,10),rep(0,40))

a_p<-c(rep(-1,20),rep(1,20),rep(0,60))
b_p<-c(rep(1,20),rep(0,20),rep(1,10),rep(-1,10),rep(0,40))
C_m=(1+b_m)%*%t(g)+a_m%*%(t(g)-t(rep(1,length(g))))
C_p=(1+b_p)%*%t(g)+a_p%*%t(rep(1,length(g)))
C=C_m+C_p
N=nrow(C)
K=ncol(C)

E_m=matrix(nrow=N,ncol=K)
E_p=matrix(nrow=N,ncol=K)
for(i in 1:N){
  for (j in 1:K){
    E_m[i,j]=(-1)^rbinom(1,1,0.5)*rpois(1,0.5)#0.5 simu
    E_p[i,j]=(-1)^rbinom(1,1,0.5)*rpois(1,0.5)
  }
}

real_C_m=C_m+E_m;
real_C_p=C_p+E_p
for(i in 1:N){
  for(j in 1:K){
    real_C_m[i,j]=max(0,real_C_m[i,j])
    real_C_p[i,j]=max(0,real_C_p[i,j])
  }
}

coordinate=(1:100)*1000
coordinate=data.frame(rep(1,100),coordinate)

wkdata=list(nsample=K,nloci=N,minor=real_C_m,major=real_C_p,loci=coordinate)