source("subfun.R")
sourceCpp("Funanova.cpp")
library(fda)
library(assist)
library(Rcpp)
library(Matrix)
library(parallel)

phi1x<-function(x)
{
 temp= x-0.5
 temp
}

R0<-function(x,z)
{ len=length(x)
temp=NULL
for(i in 1:len)
  temp=rbind(temp,1+(x[i]-0.5)*(z-0.5))

temp
}

R01<-function(x,z)
{ len=length(x)
temp=NULL
for(i in 1:len)
  temp=rbind(temp, (x[i]-0.5)*(z-0.5))

temp
}

#################### The proposed model L2 fit 
# ds: a vector of length 8 with component value 0 or 1, where 1 means the corresponding parametric component used in model fit. 
# thetas: a vector of length 7 with component value 0 or 1, where 1 means the corresponding nonparametric component used in model fit.

funanova.general.parallel.l2<-function(x,datax,datay,ny=10,nbasis=10,norder=4,lamb=0.0001,bf.lamb=0.001,maxit=100,xmin=0,xmax=1,limnla=c(-0.1,0.1),ds=NULL,thetas=NULL)
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
  ################  compute y_ij based on the EFPC 
  
  tmp1=matrix(rep(1,ny))
  tmp2=matrix(rep(1,n))
  tmp3=matrix(rep(1,n*ny))
  temp=cbind(kronecker(t(datay),tmp1),kronecker(tmp2,t(harmvals)),kronecker(matrix(x,1,m),tmp3))
  Yij=apply(temp,1,numinterxy)
  
  ################  compute X_ij based on the EFPC 
  #g(t,s,x(t),x(s)): 1 t; 2 s; 3 x(t); 4 x(s)
  
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
  
  Xij=rbind(Xij,temp)
  }
  }
  
  if(is.null(ds)){ds=rep(1,dim(Xij)[2])}
  Xij=Xij[,abs(ds)>0.0001,drop=FALSE]
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
    
    
  }
  
  rm(tempres)
  #####################################################
  Yijvec=c(Yij)
  Sigma13=Sigma13ls+Sigma13sl+Sigma13ss
  Sigma14=Sigma14ls+Sigma14sl+Sigma14ss
  Sigma24=Sigma24ls+Sigma24sl+Sigma24ss
  Sigma34=Sigma34ls+Sigma34sl+Sigma34ss
  
  Q1tmp=list(Sigma1s,Sigma3s,Sigma4s,Sigma13,Sigma14,Sigma24,Sigma34)
  qlen=length(Q1tmp)
  if(is.null(thetas)){thetas=rep(1,qlen)}
  Q1=list()
  ij=1
  for(i in 1:qlen){
    if(thetas[i]>0.00001){
      Q1[[ij]]<-Q1tmp[[i]]
      ij=ij+1
    }
  }
  
  rm(Sigma1s,Sigma3s,Sigma4s,Sigma13ls,Sigma13sl,Sigma13ss,Sigma14ls,Sigma14sl,Sigma14ss,
     Sigma24ls,Sigma24sl,Sigma24ss,Sigma34ls,Sigma34sl,Sigma34ss,Sigma13,Sigma14,Sigma24,Sigma34)
  rm(Rmat.xxt21,Rmat.xxt12,Rmat.xxt22)
  
  
  if(length(Q1)>0){
    lambda1=NULL
    for(i in 1:length(Q1))
    {q=Q1[[i]]
    lambdatmp=sum(diag(q))
    lambda1=c(lambda1,log10(length(Yijvec)*lambdatmp))
    }
    lambda1[lambda1<0]=0
    
    prec=1
    n.bf=length(Q1)+1
    f=matrix(0,n*ny,n.bf)
    theta=rep(0,length(Q1))
    j=0
    while((prec>bf.lamb) &  (j<maxit)){
      fold=f
      yw=Yijvec-apply(f[,-1,drop=FALSE],1,sum)
      tmp=lm(yw~Xij-1)
      f[,1]=tmp$fitted.values
      for(i in 2:n.bf)
      {
        yw=Yijvec-apply(f[,-i,drop=FALSE],1,sum)
        q=Q1[[i-1]]
        lambda0=lambda1[i-1]+limnla
        res=dsidr(y = yw, q = q, limnla=lambda0)
        f[,i]=res$fit
        lambda=(10^res$nlaht)/length(yw)
        theta[i-1]=1/lambda
      }
      prec=mean(abs(f-fold))
      j=j+1
    }
    
    q=matrix(0,n*ny,n*ny)
    for(i in 1:(n.bf-1))
    {
      q=q+theta[i]*Q1[[i]]
    }
    
    tmp=eigen(q)
    vec=tmp$vectors
    val=tmp$values
    if(any(val<0.00001)){q=NPD(q)}
    
    result=1000
    tmp=eigen(q)
    val=tmp$values
    if(all(val>=0) &  all(val!=1/0) ) 
    {
      res=dsidr(y = Yijvec, q = q, s = Xij)
      cmat=res$c
      dmat=res$d
      lambda=(10^res$nlaht)/length(Yijvec)
     
      ds[abs(ds)>0.00001]=dmat
      thetas[abs(thetas)>0.00001]=theta
      result=list(cmat=cmat,dmat=ds,Psi=harmvals,lambda=lambda,theta=thetas,harmfd=harmfd)
    }
    rm(Q1)
  } else {
    
    res=lm(Yijvec~Xij-1)
    
    cmat=rep(0,n*ny)
    dmat=res$coefficients
    lambda=1 
    ds[abs(ds)>0.00001]=dmat
    result=list(cmat=cmat,dmat=ds,Psi=harmvals,lambda=lambda,theta=thetas,harmfd=harmfd)
    
    
  }
  
  gc()
  result
}

######################## Extended  FLM model L2 fit 
funanova.EFLM.parallel.l2<-function(x,datax,datay,ny=10,nbasis=10,norder=4,lamb=0.0001,bf.lamb=0.001,maxit=100,xmin=0,xmax=1,limnla=c(-0.1,0.1))
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
  ##############################################################
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
  #################### 4  xs ####################
  x1=phi1x(xi)  
  tmp=numIntegrate(x1,x)*numIntegrate(y1,x)
  temp=c(temp,tmp)
  
  #################### (1  4)  t xs #############
  x1=phi1x(x)
  x2=phi1x(xi)  
  tmp=x1*y1
  tmp=numIntegrate(x2,x)*numIntegrate(tmp,x)
  temp=c(temp,tmp)
  #################### (2  4)  s xs #############
  x1=phi1x(xi)*phi1x(x) 
  tmp=numIntegrate(y1,x)*numIntegrate(x1,x)
  temp=c(temp,tmp)
  
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
  #Rmat.xxt11=Rmat.xxt12
  
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
    # tmp11=matrix(apply(tmp,1,R11),m,m, byrow=TRUE)
    
    tmp=list(r12=tmp12,r21=tmp21,r22=tmp22)
    tmp
  }
  
  for(i in 1:dim(parind)[1])
  {
    lab=parind[i,]
    Rmat.xxt12[lab[1],lab[2],,]=tempres[[i]][[1]]
    Rmat.xxt21[lab[1],lab[2],,]=tempres[[i]][[2]]
    Rmat.xxt22[lab[1],lab[2],,]=tempres[[i]][[3]]
    #  Rmat.xxt11[lab[1],lab[2],,]=tempres[[i]][[4]]
    
  }
  
  rm(tempres)
  
  ###############################################
  
  Sigma1s=matrix(0,n*ny,n*ny)
  Sigma4s=matrix(0,n*ny,n*ny)
  Sigma14=matrix(0,n*ny,n*ny)
  Sigma14ls=matrix(0,n*ny,n*ny)
  Sigma14sl=matrix(0,n*ny,n*ny)
  Sigma14ss=matrix(0,n*ny,n*ny)
  Sigma24=matrix(0,n*ny,n*ny)
  Sigma24ls=matrix(0,n*ny,n*ny)
  Sigma24sl=matrix(0,n*ny,n*ny)
  Sigma24ss=matrix(0,n*ny,n*ny)
  
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
        
        
        temp=cbind(Rmat.xt[i,i1,,],xkro)
        tmp=apply(temp,1,numinterx)
        tmp4s=numIntegrate(tmp,x)
        tmp1=numIntegrate(y1,x)*numIntegrate(y2,x)
        tmp=tmp4s*tmp1
        tempres1=c(tempres1,tmp)  #Sigma4s
        
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
    Sigma4s[(i-1)*ny+j,] = tempres[[ij]][2,]
      
    Sigma14ls[(i-1)*ny+j,] = tempres[[ij]][3,]
    Sigma14sl[(i-1)*ny+j,] = tempres[[ij]][4,]
    Sigma14ss[(i-1)*ny+j,] = tempres[[ij]][5,]
    
    Sigma24ls[(i-1)*ny+j,] = tempres[[ij]][6,]
    Sigma24sl[(i-1)*ny+j,] = tempres[[ij]][7,]
    Sigma24ss[(i-1)*ny+j,] = tempres[[ij]][8,]
    
  }
  
  rm(tempres)
  
  #####################################################
  Yijvec=c(Yij)
  Sigma14=Sigma14ls+Sigma14sl+Sigma14ss
  Sigma24=Sigma24ls+Sigma24sl+Sigma24ss
  
  Q1=list(Sigma1s,Sigma4s,Sigma14,Sigma24)
  
  
  rm(Sigma1s,Sigma4s,Sigma14ls,Sigma14sl,Sigma14ss,Sigma24ls,Sigma24sl,Sigma24ss)
  rm(Rmat.xxt21,Rmat.xxt12,Rmat.xxt22)
  
  lambda1=NULL
  for(i in 1:length(Q1))
  {q=Q1[[i]]
  lambdatmp=sum(diag(q))
  lambda1=c(lambda1,log10(length(Yijvec)*lambdatmp))
  }
  
  
  prec=1
  n.bf=length(Q1)+1
  f=matrix(0,n*ny,n.bf)
  theta=rep(0,length(Q1))
  j=0
  while((prec>bf.lamb) &  (j<maxit)){
    fold=f
    yw=Yijvec-apply(f[,-1,drop=FALSE],1,sum)
    tmp=lm(yw~Xij-1)
    f[,1]=tmp$fitted.values
    for(i in 2:n.bf)
    {
      yw=Yijvec-apply(f[,-i,drop=FALSE],1,sum)
      q=Q1[[i-1]]
      lambda0=lambda1[i-1]+limnla
      res=dsidr(y = yw, q = q, limnla=lambda0)
      f[,i]=res$fit
      lambda=(10^res$nlaht)/length(yw)
      theta[i-1]=1/lambda
    }
    prec=mean(abs(f-fold))
    j=j+1
  }
  
  q=matrix(0,n*ny,n*ny)
  for(i in 1:(n.bf-1))
  {
    q=q+theta[i]*Q1[[i]]
  }
  
  tmp=eigen(q)
  vec=tmp$vectors
  val=tmp$values
  if(any(val<0.00001)){q=NPD(q)}
 
  result=1000
  tmp=eigen(q)
  val=tmp$values
  if(all(val>=0) &  all(val!=1/0) ) 
  {
    res=dsidr(y = Yijvec, q = q, s = Xij)
    #res=dmudr(y = Yijvec, q = Q1, s = Xij)
    
    cmat=res$c
    dmat=res$d
    lambda=(10^res$nlaht)/length(Yijvec)
	
    result=list(cmat=cmat,dmat=dmat,Psi=harmvals,lambda=lambda,theta=theta,harmfd=harmfd)
  }
  rm(Q1)
  gc()
  result
}


#################### NCM model L2 fit 
funanova.NCM.parallel.l2<-function(x,datax,datay,ny=10,nbasis=10,norder=4,lamb=0.0001,bf.lamb=0.001,maxit=100,xmin=0,xmax=1,limnla=c(-0.1,0.1))
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
  ##############################################################
  Xij=NULL
  for(i in 1:n)
  { 
    for(j in 1:ny)
    { y1=harmvals[,j]
    xi=datax[,i]
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
    #################### (1  3)  t xt #############
    x1=phi1x(xi)*phi1x(x) 
    tmp=x1*y1
    temp=c(temp,numIntegrate(tmp,x))
    
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
#Rmat.xxt11=Rmat.xxt12

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

  ###############################################
  
  Sigma1s=matrix(0,n*ny,n*ny)
  Sigma3s=matrix(0,n*ny,n*ny)
  Sigma13ls=matrix(0,n*ny,n*ny)
  Sigma13sl=matrix(0,n*ny,n*ny)
  Sigma13ss=matrix(0,n*ny,n*ny)
  
  
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
    ################ (1  3) ########################
    temp=cbind(Rmat.xxt12[i,i1,,],y2xkro)
    tmp=apply(temp,1,numinterxy)
    tmp=tmp*y1
    tmp=numIntegrate(tmp,x)
    tempres1=c(tempres1,tmp)  #Sigma13ls
    
    temp=cbind(Rmat.xxt21[i,i1,,],y2xkro)
    tmp=apply(temp,1,numinterxy)
    tmp=tmp*y1
    tmp=numIntegrate(tmp,x)
    tempres1=c(tempres1,tmp)  #Sigma13sl
    
    temp=cbind(Rmat.xxt22[i,i1,,],y2xkro)
    tmp=apply(temp,1,numinterxy)
    tmp=tmp*y1
    tmp=numIntegrate(tmp,x)
    tempres1=c(tempres1,tmp)  #Sigma13ss
    tempres2=cbind(tempres2,tempres1)
      }
    return(tempres2)
  }
 
  for(ij in 1:dim(parind)[1])
  {
    lab=parind[ij,]
    i=lab[1]
    j=lab[2]
    
    Sigma1s[(i-1)*ny+j,] = tempres[[ij]][1,]
    Sigma3s[(i-1)*ny+j,] = tempres[[ij]][2,]
    
    Sigma13ls[(i-1)*ny+j,] = tempres[[ij]][3,]
    Sigma13sl[(i-1)*ny+j,] = tempres[[ij]][4,]
    Sigma13ss[(i-1)*ny+j,] = tempres[[ij]][5,]
    
  }
  
  rm(tempres)
  
  
  ######################################################
  Yijvec=c(Yij)
  #Sigma13ss=Sigma13sl+Sigma13ss
  Sigma13=Sigma13ls+Sigma13sl+Sigma13ss
  Q1=list(Sigma1s,Sigma3s,Sigma13)
  rm(Sigma1s,Sigma3s,Sigma13ls,Sigma13sl,Sigma13ss)
  rm(Rmat.xxt21,Rmat.xxt12,Rmat.xxt22)
  
  lambda1=NULL
  for(i in 1:length(Q1))
  {q=Q1[[i]]
  lambdatmp=sum(diag(q))
  lambda1=c(lambda1,log10(length(Yijvec)*lambdatmp))
  }
  
  
  prec=1
  n.bf=length(Q1)+1
  f=matrix(0,n*ny,n.bf)
  theta=rep(0,length(Q1))
  j=0
  while((prec>bf.lamb) &  (j<maxit)){
    fold=f
    yw=Yijvec-apply(f[,-1,drop=FALSE],1,sum)
    tmp=lm(yw~Xij-1)
    f[,1]=tmp$fitted.values
    for(i in 2:n.bf)
    {
      yw=Yijvec-apply(f[,-i,drop=FALSE],1,sum)
      q=Q1[[i-1]]
	    lambda0=lambda1[i-1]+limnla
      res=dsidr(y = yw, q = q, limnla=lambda0)
      f[,i]=res$fit
      lambda=(10^res$nlaht)/length(yw)
      theta[i-1]=1/lambda
    }
    prec=mean(abs(f-fold))
    j=j+1
  }
  
  q=matrix(0,n*ny,n*ny)
  for(i in 1:(n.bf-1))
  {
    q=q+theta[i]*Q1[[i]]
  }
  
  tmp=eigen(q)
  vec=tmp$vectors
  val=tmp$values
  if(any(val<0.00001)){q=NPD(q)}
  
  result=1000
  tmp=eigen(q)
  val=tmp$values
  if(all(val>=0) &  all(val!=1/0)) 
  {
    res=dsidr(y = Yijvec, q = q, s = Xij)
   
    cmat=res$c
    dmat=res$d
    lambda=(10^res$nlaht)/length(Yijvec)
    
    result=list(cmat=cmat,dmat=dmat,Psi=harmvals,lambda=lambda,theta=theta,harmfd=harmfd)
    
  }
  rm(Q1)
  gc()
  result
}


####################  NCM model L2 fit 

funanova.NCM.parallel.l2<-function(x,datax,datay,ny=10,nbasis=10,norder=4,lamb=0.0001,bf.lamb=0.001,maxit=100,xmin=0,xmax=1,limnla=c(-0.1,0.1))
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
  #g(t, x(t)): 1 t;  3 x(t)
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
    #################### (1  3)  t xt #############
    x1=phi1x(xi)*phi1x(x) 
    tmp=x1*y1
    temp=c(temp,numIntegrate(tmp,x))
    
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

  ###############################################
  
  Sigma1s=matrix(0,n*ny,n*ny)
  Sigma3s=matrix(0,n*ny,n*ny)
  Sigma13ls=matrix(0,n*ny,n*ny)
  Sigma13sl=matrix(0,n*ny,n*ny)
  Sigma13ss=matrix(0,n*ny,n*ny)
  
  
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
    ################ (1  3) ########################
    temp=cbind(Rmat.xxt12[i,i1,,],y2xkro)
    tmp=apply(temp,1,numinterxy)
    tmp=tmp*y1
    tmp=numIntegrate(tmp,x)
    tempres1=c(tempres1,tmp)  #Sigma13ls
    
    temp=cbind(Rmat.xxt21[i,i1,,],y2xkro)
    tmp=apply(temp,1,numinterxy)
    tmp=tmp*y1
    tmp=numIntegrate(tmp,x)
    tempres1=c(tempres1,tmp)  #Sigma13sl
    
    temp=cbind(Rmat.xxt22[i,i1,,],y2xkro)
    tmp=apply(temp,1,numinterxy)
    tmp=tmp*y1
    tmp=numIntegrate(tmp,x)
    tempres1=c(tempres1,tmp)  #Sigma13ss
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
    
    Sigma13ls[(i-1)*ny+j,] = tempres[[ij]][3,]
    Sigma13sl[(i-1)*ny+j,] = tempres[[ij]][4,]
    Sigma13ss[(i-1)*ny+j,] = tempres[[ij]][5,]
    
  }
  
  rm(tempres)
  
  
  ######################################################
  Yijvec=c(Yij)
  #Sigma13ss=Sigma13sl+Sigma13ss
  Sigma13=Sigma13ls+Sigma13sl+Sigma13ss
  Q1=list(Sigma1s,Sigma3s,Sigma13)
  rm(Sigma1s,Sigma3s,Sigma13ls,Sigma13sl,Sigma13ss)
  rm(Rmat.xxt21,Rmat.xxt12,Rmat.xxt22)
  
  lambda1=NULL
  for(i in 1:length(Q1))
  {q=Q1[[i]]
  lambdatmp=sum(diag(q))
  lambda1=c(lambda1,log10(length(Yijvec)*lambdatmp))
  }
  
  
  prec=1
  n.bf=length(Q1)+1
  f=matrix(0,n*ny,n.bf)
  theta=rep(0,length(Q1))
  j=0
  while((prec>bf.lamb) &  (j<maxit)){
    fold=f
    yw=Yijvec-apply(f[,-1,drop=FALSE],1,sum)
    tmp=lm(yw~Xij-1)
    f[,1]=tmp$fitted.values
    for(i in 2:n.bf)
    {
      yw=Yijvec-apply(f[,-i,drop=FALSE],1,sum)
      q=Q1[[i-1]]
	  lambda0=lambda1[i-1]+limnla
      res=dsidr(y = yw, q = q, limnla=lambda0)
      f[,i]=res$fit
      lambda=(10^res$nlaht)/length(yw)
      theta[i-1]=1/lambda
    }
    prec=mean(abs(f-fold))
    j=j+1
  }
  
  q=matrix(0,n*ny,n*ny)
  for(i in 1:(n.bf-1))
  {
    q=q+theta[i]*Q1[[i]]
  }
  
  tmp=eigen(q)
  vec=tmp$vectors
  val=tmp$values
  if(any(val<0.00001)){q=NPD(q)}
  
  result=1000
  tmp=eigen(q)
  val=tmp$values
  if(all(val>=0) &  all(val!=1/0)) 
  {
    res=dsidr(y = Yijvec, q = q, s = Xij)
   
    cmat=res$c
    dmat=res$d
    lambda=(10^res$nlaht)/length(Yijvec)
    
    result=list(cmat=cmat,dmat=dmat,Psi=harmvals,lambda=lambda,theta=theta,harmfd=harmfd)
    
  }
  rm(Q1)
  gc()
  result
}

