source('anova_ada_high.R')
source('subfun.R')
library(matrixcalc)
library(Matrix)
library(quadprog)
library(glmnet)
library(assist)
library(gss)
library(matrixcalc)

phi1x<-function(x)
{
 temp=sqrt(12)*(x-0.5)
 temp
}

R0<-function(x,z)
{ len=length(x)
temp=NULL
for(i in 1:len)
  temp=rbind(temp,1+12*(x[i]-0.5)*(z-0.5))

temp
}

R01<-function(x,z)
{ len=length(x)
temp=NULL
for(i in 1:len)
  temp=rbind(temp,12*(x[i]-0.5)*(z-0.5))

temp
}

funanova.general.parallel.ho.l1<-function(x,datax,datay,lambda1_list=seq(-3,3,by=0.5),tau1_list=seq(-8,2,by=0.5),M_list=seq(0.5,0.9,by=0.05),ny=10,nbasis=10,norder=4,lamb=0.0001,bf.lamb=0.001,maxit=100,xmin=0,xmax=1,lam_range=c(-3,3),inclu=TRUE)
{ 
  n=dim(datay)[2]
  m=dim(datay)[1]
  # Fitting the PCA on response data
  timebasis = create.bspline.basis(c(xmin,xmax),nbasis=nbasis,norder=norder)
  harmLfd = int2Lfd(m=2)
  tempfdPar = fdPar(timebasis,harmLfd,lamb)
  tempfd = smooth.basis(x,datay,tempfdPar)
  
  # Evlauate  fPCA's
  temppca = pca.fd(tempfd$fd,nharm=ny)
  #### The First ny PC's
  harmfd = temppca$harmonics
  harmvals = eval.fd(x,harmfd)
  ################  compute y_ij based on the EFPC ##############
  
  tmp1=matrix(rep(1,ny))
  tmp2=matrix(rep(1,n))
  tmp3=matrix(rep(1,n*ny))
  temp=cbind(kronecker(t(datay),tmp1),kronecker(tmp2,t(harmvals)),kronecker(matrix(x,1,m),tmp3))
  Yij=apply(temp,1,numinterxy)
  ################  compute X_ij based on the EFPC ##############
  #g(t,s,x(t),x(s)): 1 t; 2 s; 3 x(t); 4 x(s)
  ################# create basis functions #####################
  Xij=NULL
  for(i in 1:n)
  { xi=datax[,i]
  for(j in 1:ny)
  { y1=harmvals[,j]
  ##################### 0  1 #####################
  tmp=y1
  temp=numIntegrate(tmp,x)
  
  ##################### 1  t #####################
  x1=phi1x(x) 
  tmp=x1*y1
  temp=c(temp,numIntegrate(tmp,x))
  
  #################### 3  xt ####################
  x1=phi1x(xi) 
  tmp=x1*y1
  temp=c(temp,numIntegrate(tmp,x))
  #################### 4  xs ####################
  x1=phi1x(xi)  
  tmp=numIntegrate(x1,x)*numIntegrate(y1,x)
  temp=c(temp,tmp)
  
  #################### (1  3)  t xt #############
  x1=phi1x(xi)*phi1x(x) 
  tmp=x1*y1
  temp=c(temp,numIntegrate(tmp,x))
  
  #################### (1  4)  t xs #############
  x1=phi1x(x)
  x2=phi1x(xi)  
  tmp=x1*y1
  tmp=numIntegrate(x2,x)*numIntegrate(tmp,x)
  temp=c(temp,tmp)
  #################### (3  4)  xt xs #############
  x1=phi1x(xi)  
  tmp=x1*y1
  tmp=numIntegrate(tmp,x)*numIntegrate(x1,x)
  temp=c(temp,tmp)
  #################### (2  4)  s xs #############
  x1=phi1x(xi)*phi1x(x) 
  tmp=numIntegrate(y1,x)*numIntegrate(x1,x)
  temp=c(temp,tmp)
  
  if(inclu){
    #################### (1 2  4)  t s xs #############
    x1=phi1x(xi)*phi1x(x) #(datax[,i]-0.5)*(x-0.5)
    x2=y1*phi1x(x) #(x-0.5)
    tmp=numIntegrate(x2,x)*numIntegrate(x1,x)
    temp=c(temp,tmp)
    
    #################### (2 3 4)  s xt  xs #############
    x1=phi1x(xi)*phi1x(x)  #(datax[,i]-0.5)*(x-0.5)
    tmp=phi1x(xi)*y1
    tmp=numIntegrate(tmp,x)*numIntegrate(x1,x)
    temp=c(temp,tmp)
    
    #################### (1 2 3 4)  t  s xt  xs #############
    x1=phi1x(xi)*phi1x(x)  #(datax[,i]-0.5)*(x-0.5)
    tmp=x1*y1
    tmp=numIntegrate(tmp,x)*numIntegrate(x1,x)
    temp=c(temp,tmp)
  }
  
  Xij=rbind(Xij,temp)
  }
  }
  ################  compute Sigma_ij based on the EFPC ##############
  Rmat.x=R1(x,x)
  parind=cbind(rep(c(1:n),each=n),rep(c(1:n),n))
  #dopar
  tempres=foreach(i=1:dim(parind)[1]) %dopar%  {
    
    lab=parind[i,]
    x1=datax[,lab[1]]
    x2=datax[,lab[2]]
    tmp=R1(x1,x2)
    tmp
  }
  
  Rmat=array(0,c(n,n,m,m))
  for(i in 1:dim(parind)[1])
  {
    lab=parind[i,]
    Rmat[lab[1],lab[2],,]=tempres[[i]]
  }
  
  rm(tempres)
  
  Rmat.xt=Rmat
  
  Rmat.xxt12=array(0,c(n,n,m,m))
  Rmat.xxt21=Rmat.xxt12
  Rmat.xxt22=Rmat.xxt12
  
  parind=cbind(rep(c(1:n),each=n),rep(c(1:n),n))
  tempres=foreach(i=1:dim(parind)[1]) %dopar%  {
    lab=parind[i,]
    x1=datax[,lab[1]]
    x2=datax[,lab[2]]
    tmp=matrix(rep(1,m))
    tmp1=cbind(x,x1)
    tmp2=cbind(x,x2)
    tmp=cbind(kronecker(tmp1,tmp),kronecker(tmp,tmp2))
    
    tmp12=matrix(apply(tmp,1,R12),m,m, byrow=TRUE)
    tmp21=matrix(apply(tmp,1,R21),m,m, byrow=TRUE)
    tmp22=matrix(apply(tmp,1,R22),m,m, byrow=TRUE)
    
    tmp=list(r12=tmp12,r21=tmp21,r22=tmp22)
    tmp
  }
  
  for(i in 1:dim(parind)[1])
  {
    lab=parind[i,]
    Rmat.xxt12[lab[1],lab[2],,]=tempres[[i]][[1]]
    Rmat.xxt21[lab[1],lab[2],,]=tempres[[i]][[2]]
    Rmat.xxt22[lab[1],lab[2],,]=tempres[[i]][[3]]
    
  }
  
  rm(tempres)
  
  ######### create RK matrices ################
  
  Sigma1s=matrix(0,n*ny,n*ny)
  Sigma3s=matrix(0,n*ny,n*ny)
  Sigma4s=matrix(0,n*ny,n*ny)
  Sigma13=matrix(0,n*ny,n*ny)
  Sigma13ls=matrix(0,n*ny,n*ny)
  Sigma13sl=matrix(0,n*ny,n*ny)
  Sigma13ss=matrix(0,n*ny,n*ny)
  Sigma14=matrix(0,n*ny,n*ny)
  Sigma14ls=matrix(0,n*ny,n*ny)
  Sigma14sl=matrix(0,n*ny,n*ny)
  Sigma14ss=matrix(0,n*ny,n*ny)
  Sigma24=matrix(0,n*ny,n*ny)
  Sigma24ls=matrix(0,n*ny,n*ny)
  Sigma24sl=matrix(0,n*ny,n*ny)
  Sigma24ss=matrix(0,n*ny,n*ny)
  Sigma34=matrix(0,n*ny,n*ny)
  Sigma34ls=matrix(0,n*ny,n*ny)
  Sigma34sl=matrix(0,n*ny,n*ny)
  Sigma34ss=matrix(0,n*ny,n*ny)
  Sigma124=matrix(0,n*ny,n*ny)
  Sigma234=matrix(0,n*ny,n*ny)
  Sigma1234=matrix(0,n*ny,n*ny)
  
  
  parind=cbind(rep(c(1:n),each=ny),rep(c(1:ny),n))
  tempres=foreach(ij=1:dim(parind)[1]) %dopar%  {
    lab=parind[ij,]
    i=lab[1]
    j=lab[2]
    
    x1=datax[,i]
    y1=harmvals[,j]
    tempres2=NULL
    for(i1 in 1:n)
      for(j1 in 1:ny)
      {
        x2=datax[,i1]
        y2=harmvals[,j1]
        tempres1=NULL
        
        tmp=matrix(rep(1,m))
        y2xkro=kronecker(matrix(c(y2,x),1,),tmp)
        xkro=kronecker(matrix(c(x),1,),tmp)
        
        temp=cbind(Rmat.x,y2xkro)
        tmp=apply(temp,1,numinterxy)
        tmp=tmp*y1
        tmp1s=numIntegrate(tmp,x)
        tempres1=c(tempres1,tmp1s)  #Sigma1s
        
        
        temp=cbind(Rmat.xt[i,i1,,],y2xkro)
        tmp=apply(temp,1,numinterxy)
        tmp=tmp*y1
        tmp3s=numIntegrate(tmp,x)
        tempres1=c(tempres1,tmp3s)  #Sigma3s
        
        temp=cbind(Rmat.xt[i,i1,,],xkro)
        tmp=apply(temp,1,numinterx)
        tmp4s=numIntegrate(tmp,x)
        tmp1=numIntegrate(y1,x)*numIntegrate(y2,x)
        tmp=tmp4s*tmp1
        tempres1=c(tempres1,tmp)  #Sigma4s
        
        ################ (1  3) ########################
        temp=cbind(Rmat.xxt12[i,i1,,],y2xkro)
        tmp=apply(temp,1,numinterxy)
        tmp=tmp*y1
        tmp1312=numIntegrate(tmp,x)
        tempres1=c(tempres1,tmp1312)  #Sigma13ls
        
        temp=cbind(Rmat.xxt21[i,i1,,],y2xkro)
        tmp=apply(temp,1,numinterxy)
        tmp=tmp*y1
        tmp1321=numIntegrate(tmp,x)
        tempres1=c(tempres1,tmp1321)  #Sigma13sl
        
        temp=cbind(Rmat.xxt22[i,i1,,],y2xkro)
        tmp=apply(temp,1,numinterxy)
        tmp=tmp*y1
        tmp1322=numIntegrate(tmp,x)
        tempres1=c(tempres1,tmp1322)  #Sigma13ss
        
        ##########################  (1  4)  #################
        
        tmp1=tmp4s 
        tmp=numIntegrate(phi1x(x)*y1,x)*numIntegrate(phi1x(x)*y2,x)
        tmp=tmp*tmp1
        tempres1=c(tempres1,tmp)  #Sigma14ls
        
        tmp1=numIntegrate(phi1x(x1),x)*numIntegrate(phi1x(x2),x)
        tmp=tmp1s*tmp1
        tempres1=c(tempres1,tmp)  #Sigma14sl
		
        tmp1=tmp1s  
        tmp=tmp4s 
        tmp=tmp*tmp1
        tempres1=c(tempres1,tmp)  #Sigma14ss
        
        ##########################(2  4) #######################
        
        temp=cbind(Rmat.xxt12[i,i1,,],xkro)
        tmp=apply(temp,1,numinterx)
        tmp2412=numIntegrate(tmp,x)
        tmpy1y2=numIntegrate(y1,x)*numIntegrate(y2,x)
        tmp=tmp2412*tmpy1y2
        tempres1=c(tempres1,tmp)  #Sigma24ls
        
        temp=cbind(Rmat.xxt21[i,i1,,],xkro)
        tmp=apply(temp,1,numinterx)
        tmp2421=numIntegrate(tmp,x)
        tmp=tmp2421*tmpy1y2
        tempres1=c(tempres1,tmp)  #Sigma24sl
        
        temp=cbind(Rmat.xxt22[i,i1,,],xkro)
        tmp=apply(temp,1,numinterx)
        tmp2422=numIntegrate(tmp,x)
        tmp=tmp2422*tmpy1y2
        tempres1=c(tempres1,tmp)  # Sigma24ss
        #######################(3   4) ############################
        
        tmp1=numIntegrate(phi1x(x1)*y1,x)*numIntegrate(phi1x(x2)*y2,x)
        tmp=tmp4s 
        tmp=tmp*tmp1
        tempres1=c(tempres1,tmp)  #Sigma34ls
        
        
        tmp1=tmp3s 
        tmp=numIntegrate(phi1x(x1),x)*numIntegrate(phi1x(x2),x)
        tmp=tmp*tmp1
        tempres1=c(tempres1,tmp)  #Sigma34sl
        
        
        tmp1=tmp3s 
        tmp=tmp4s 
        tmp=tmp*tmp1
        tempres1=c(tempres1,tmp)  # Sigma34ss
        
        
        ##########################(1  2  4)  t s xs  #######################
        tmp124=NULL
        
        xy1=phi1x(x)*y1
        xy2=phi1x(x)*y2
        tmpy1y2=numIntegrate(xy1,x)*numIntegrate(xy2,x)
        tmp=tmp2412*tmpy1y2
        tmp124=c(tmp124,tmp)  #Sigma124lls
        
        tmp=tmp2421*tmpy1y2
        tmp124=c(tmp124,tmp)  #Sigma124lsl
        
        tmp=tmp2422*tmpy1y2
        tmp124=c(tmp124,tmp)  #Sigma124lss
        
        xx1=phi1x(x1)*phi1x(x) 
        xx2=phi1x(x2)*phi1x(x) 
        tmp2411=numIntegrate(xx1,x)*numIntegrate(xx2,x)
        tmp=tmp1s*tmp2411
        tmp124=c(tmp124,tmp)  #Sigma124sll
        
        tmp=tmp2412*tmp1s
        tmp124=c(tmp124,tmp)  #Sigma124sls
        
        tmp=tmp2421*tmp1s
        tmp124=c(tmp124,tmp)  #Sigma124ssl
        
        tmp=tmp2422*tmp1s
        tmp124=c(tmp124,tmp)  #Sigma124sss
        
        tempres1=c(tempres1,sum(tmp124))
        #######################(2  3   4)  s xt  xs ############################
        tmp234=NULL
        
        xy1=phi1x(x1)*y1
        xy2=phi1x(x2)*y2
        tmpy1y2=numIntegrate(xy1,x)*numIntegrate(xy2,x)
        tmp=tmp2412*tmpy1y2
        tmp234=c(tmp234,tmp)  #Sigma234lls
        
        tmp=tmp2421*tmpy1y2
        tmp234=c(tmp234,tmp)  #Sigma234sll
        
        tmp=tmp2422*tmpy1y2
        tmp234=c(tmp234,tmp)  #Sigma234sls
        
        tmp=tmp3s*tmp2411
        tmp234=c(tmp234,tmp)  #Sigma234lsl
        
        tmp=tmp2412*tmp3s
        tmp234=c(tmp234,tmp)  #Sigma234lss
        
        tmp=tmp2421*tmp3s
        tmp234=c(tmp234,tmp)  #Sigma234ssl
        
        tmp=tmp2422*tmp3s
        tmp234=c(tmp234,tmp)  #Sigma234sss
        
        tempres1=c(tempres1,sum(tmp234))
        
        #######################(1  2  3   4) t s xt  xs ############################
        tmp1234=NULL
        
        xy1=phi1x(x1)*phi1x(x)*y1
        xy2=phi1x(x2)*phi1x(x)*y2
        tmp1311=numIntegrate(xy1,x)*numIntegrate(xy2,x)
        tmp=tmp2412*tmp1311
        tmp1234=c(tmp1234,tmp)  #Sigma1234llls
        
        tmp=tmp2421*tmp1311
        tmp1234=c(tmp1234,tmp)  #Sigma1234lsll
        
        tmp=tmp2422*tmp1311
        tmp1234=c(tmp1234,tmp)  #Sigma1234lsls
        
        tmp=tmp2411*tmp1312
        tmp1234=c(tmp1234,tmp)  #Sigma1234llsl
        
        tmp=tmp2412*tmp1312
        tmp1234=c(tmp1234,tmp)  #Sigma1234llss
        
        tmp=tmp2421*tmp1312
        tmp1234=c(tmp1234,tmp)  #Sigma1234lssl
        
        tmp=tmp2422*tmp1312
        tmp1234=c(tmp1234,tmp)  #Sigma1234lsss
        
        
        tmp=tmp2411*tmp1321
        tmp1234=c(tmp1234,tmp)  #Sigma1234slll
        
        
        tmp=tmp2412*tmp1321
        tmp1234=c(tmp1234,tmp)  #Sigma1234slls
        
        tmp=tmp2421*tmp1321
        tmp1234=c(tmp1234,tmp)  #Sigma1234ssll
        
        tmp=tmp2422*tmp1321
        tmp1234=c(tmp1234,tmp)  #Sigma1234ssls
        
        tmp=tmp2411*tmp1322
        tmp1234=c(tmp1234,tmp)  #Sigma1234slsl
        
        tmp=tmp2412*tmp1322
        tmp1234=c(tmp1234,tmp)  #Sigma1234slss
        
        tmp=tmp2421*tmp1322
        tmp1234=c(tmp1234,tmp)  #Sigma1234sssl
        
        tmp=tmp2422*tmp1322
        tmp1234=c(tmp1234,tmp)  #Sigma1234ssss
        
        tempres1=c(tempres1,sum(tmp1234))
        
        
        tempres2=cbind(tempres2,tempres1)
      }
    tempres2
  }
  
  for(ij in 1:dim(parind)[1])
  {
    lab=parind[ij,]
    i=lab[1]
    j=lab[2]
    
    Sigma1s[(i-1)*ny+j,] = tempres[[ij]][1,]
    Sigma3s[(i-1)*ny+j,] = tempres[[ij]][2,]
    Sigma4s[(i-1)*ny+j,] = tempres[[ij]][3,]
    
    Sigma13ls[(i-1)*ny+j,] = tempres[[ij]][4,]
    Sigma13sl[(i-1)*ny+j,] = tempres[[ij]][5,]
    Sigma13ss[(i-1)*ny+j,] = tempres[[ij]][6,]
    
    Sigma14ls[(i-1)*ny+j,] = tempres[[ij]][7,]
    Sigma14sl[(i-1)*ny+j,] = tempres[[ij]][8,]
    Sigma14ss[(i-1)*ny+j,] = tempres[[ij]][9,]
    
    Sigma24ls[(i-1)*ny+j,] = tempres[[ij]][10,]
    Sigma24sl[(i-1)*ny+j,] = tempres[[ij]][11,]
    Sigma24ss[(i-1)*ny+j,] = tempres[[ij]][12,]
    
    Sigma34ls[(i-1)*ny+j,] = tempres[[ij]][13,]
    Sigma34sl[(i-1)*ny+j,] = tempres[[ij]][14,]
    Sigma34ss[(i-1)*ny+j,] = tempres[[ij]][15,]
    
    Sigma124[(i-1)*ny+j,] = tempres[[ij]][16,]
    Sigma234[(i-1)*ny+j,] = tempres[[ij]][17,]
    Sigma1234[(i-1)*ny+j,] = tempres[[ij]][18,]
    
  }
  
  rm(tempres)
  
  #####################################################
  Yijvec=c(Yij)
  
  Sigma13=Sigma13ls+Sigma13sl+Sigma13ss
  Sigma14=Sigma14ls+Sigma14sl+Sigma14ss
  Sigma24=Sigma24ls+Sigma24sl+Sigma24ss
  Sigma34=Sigma34ls+Sigma34sl+Sigma34ss
  
  if(inclu){
    Q1=list(Sigma1s,Sigma3s,Sigma4s,Sigma13,Sigma14,Sigma24,Sigma34,Sigma124,Sigma234,Sigma1234)
  }
  else{
    Q1=list(Sigma1s,Sigma3s,Sigma4s,Sigma13,Sigma14,Sigma24,Sigma34)
  }
  
  
  rm(Sigma1s,Sigma3s,Sigma4s,Sigma13ls,Sigma13sl,Sigma13ss,Sigma14ls,Sigma14sl,Sigma14ss,
     Sigma24ls,Sigma24sl,Sigma24ss,Sigma34ls,Sigma34sl,Sigma34ss,Sigma124,Sigma234,Sigma1234)
  rm(Rmat.xxt21,Rmat.xxt12,Rmat.xxt22)
  RK=NULL
  nu=dim(Q1[[1]])[1]
  for (i in 1:length(Q1)) {
    RK<-array(c(RK,Q1[[i]]),c(nu,nu,i))
  }
  w1=rep(1,dim(Xij)[2])
  w2=rep(1,dim(RK)[3])
  obj=anova_l1(Yijvec,Xij,RK,lambda1_list,M_list,w1,w2,tol1=0.001,tol2=0.001,maxit=100,k=5,lam_range=lam_range)
  d0=obj$d
  theta0=obj$theta
  #cross-validation to set the threshold
  thre_list1=c(0.5,0.3,0.2) 
  thre_list2=c(0.5,0.3,0.2)
  count=0
  ise_thre=rep(0,length(thre_list1)*length(thre_list2))
  for (thre1 in thre_list1) {
    for (thre2 in thre_list2) {
      count=count+1
      c=obj$c
      d=obj$d
      d[abs(d)<thre1]=0
      theta=obj$theta
      theta[theta<thre2]=0
      ghat=Xij%*%d+obj$R%*%theta
      ise_thre[count]=mean((Yijvec-ghat)^2)
    }
  }
  min_ise_pos=which.min(ise_thre)
  if (min_ise_pos%%length(thre_list2)==0){
    thre1_opt=thre_list1[min_ise_pos%/%length(thre_list2)]
    thre2_opt=thre_list2[length(thre_list2)]
  }
  else{
    thre1_opt=thre_list1[min_ise_pos%/%length(thre_list2)+1]
    thre2_opt=thre_list2[min_ise_pos%%length(thre_list2)]
  }
  d=obj$d
  d[abs(d)<thre1_opt]=0
  theta=obj$theta
  theta[theta<thre2_opt]=0
  result=list(cmat=obj$c,dmat=d,theta=theta,lambda=obj$tau0,lambda1=obj$lambda1,M=obj$M,Psi=harmvals,harmfd=harmfd,w2=obj$w2)
  rm(Q1)
  gc()
  result
}


