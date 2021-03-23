library(MASS)

###########################  subfunctions #########################

##### Generation of Gaussian process 
G.pro<-function(x,mat)  
{ if(is.matrix(x)){ tmp=mvrnorm(n=1,mu=rep(0,dim(x)[1]),Sigma=mat)} else
{tmp=mvrnorm(n=1,mu=x-x,Sigma=mat)}
  list(x=x,y=tmp)  
}

########## subfunctions to  compute reproducing kernels
k0x<-function(x){rep(1,length(x))}

k1x<-function(x){x-0.5}

k2x<-function(x){(k1x(x)^2-1/12)/2}

k4x<-function(x){(k1x(x)^4-k1x(x)^2/2+7/240)/24}

R1<-function(x,z)
{ len=length(x)
temp=NULL
for(i in 1:len)
{xz=abs(x[i]-z)
temp=rbind(temp,k2x(x[i])*k2x(z)-k4x(xz))
}
temp
}

R11<-function(x)
{ 
  temp=R01(x[1],x[3])*R01(x[2],x[4])
  c(temp)
}

R12<-function(x)
{ 
  temp=R01(x[1],x[3])*R1(x[2],x[4])
  c(temp)
}

R21<-function(x)
{ 
  temp=R1(x[1],x[3])*R01(x[2],x[4])
  c(temp)
}

R22<-function(x)
{ 
  temp=R1(x[1],x[3])*R1(x[2],x[4])
  c(temp)
}

R112<-function(x)
{ 
  temp=R01(x[1],x[4])*R01(x[2],x[5])*R1(x[3],x[6])
  temp
}

R121<-function(x)
{ 
  temp=R01(x[1],x[4])*R1(x[2],x[5])*R01(x[3],x[6])
  temp
}

R211<-function(x)
{ 
  temp=R1(x[1],x[4])*R01(x[2],x[5])*R01(x[3],x[6])
  temp
}

R122<-function(x)
{ 
  temp=R01(x[1],x[4])*R1(x[2],x[5])*R1(x[3],x[6])
  temp
}

R221<-function(x)
{ 
  temp=R1(x[1],x[4])*R1(x[2],x[5])*R01(x[3],x[6])
  temp
}

R212<-function(x)
{ 
  temp=R1(x[1],x[4])*R01(x[2],x[5])*R1(x[3],x[6])
  temp
}

R222<-function(x)
{ 
  temp=R1(x[1],x[4])*R1(x[2],x[5])*R1(x[3],x[6])
  temp
}

############## Functions to compute numerical integration 

num.integrate<-function(x,z)
{ nlen=length(x)
tmp=sum((x[-1]+x[-nlen])/2*(z[-1]-z[-nlen]))
tmp 
}


numinterxy<-function(z)
{m1=length(z)/3; tmp=z[1:m1]*z[(m1+1):(2*m1)];tmp1=z[-(1:(2*m1))];
tmp2=numIntegrate(tmp,tmp1); 
tmp2}

numinterx<-function(z)
{m1=length(z)/2; tmp=z[1:m1];tmp1=z[-(1:m1)];
tmp2=numIntegrate(tmp,tmp1); 
tmp2}


numinterxyR<-function(z)
{m1=length(z)/3; tmp=z[1:m1]*z[(m1+1):(2*m1)];tmp1=z[-(1:(2*m1))];
tmp2=num.integrate(tmp,tmp1); 
tmp2}

numinterxR<-function(z)
{m1=length(z)/2; tmp=z[1:m1];tmp1=z[-(1:m1)];
tmp2=num.integrate(tmp,tmp1); 
tmp2}


############## Functions to revise non-positive definite matrix 
 
matVec.prod<-function(mat, vec, left = TRUE) 
{ 
  ## to multiply mat with vec 
  ## v%*%M or M%*%v  
  if(left && (length(vec) != nrow(mat))) stop( 
    "length of Vector does not match row of Matrix") 
  if(!left && (length(vec) != ncol(mat))) 
    stop("length of Vector does not match column of Matrix") 
  as.vector(apply(mat, 1 + left, function(x, z) 
    sum(x * z), z = vec)) 
} 

NPD<-function(q)
{
  tmp=eigen(q)
  vec=tmp$vectors
  val=tmp$values
  val[val<0]=0
  q=vec%*%diag(val)%*%t(vec)
  
  tmp=eigen(q)
  vec=tmp$vectors
  val=tmp$values
  if(any(val<0.00001)) {
    val[val<0.00001]=0.00001
    q=vec%*%diag(val)%*%t(vec)}
  
  tmp=eigen(q)
  vec=tmp$vectors
  val=tmp$values
  if(any(val<0.0000001)) {
    val[val<0.00001]=0
    q=vec%*%diag(val)%*%t(vec)+0.00001*diag(rep(1,length(val)))}
  
  #if(any(eigen(q)$values<0)) {q=matrix(c(matrix(nearPD(q)$mat)),n*ny,n*ny)}
  
  q
  
}

############# functions to compute kernel functions in Gaussian process

#Data.f, Data.new.f: column 1st, t; row, s; column, covariate
#cov.lin.f=xixj.f(Data.f,Data.new.f) 
#v.power.f=xixj_sta.f(Data.f,Data.new.f,power=1)


cov.linear=function(hyper,Data=NULL,Data.new=NULL,cov.lin.f=NULL){
  #k(x^1,x^2)=sum_{q=1}^{Q} a_q * x^1_q * x^2_q + sum_{q=1}^{Q1} b_q *\int z^1_q(s) * z^2_q(s) ds
  #hyper is a list of hyper-parameters
  ##cov.lin.f: list of covariance for functional part
  if(is.null(Data.new)) Data.new=Data
  hyper=lapply(hyper,exp)
  cov.lin=0
  if(!is.null(Data)) {n.hyper=length(hyper$linear.a)
                      cov.lin=cov.lin+xixj(Data.new,Data,hyper$linear.a)}
  
  if(!is.null(cov.lin.f)){
  hyp.f=hyper$linear.f
  tmp=mapply(function(x,y) x*y, cov.lin.f,hyp.f,SIMPLIFY=FALSE)
  tmp=Reduce("+",tmp)
  cov.lin=cov.lin+tmp
  }
  return(cov.lin)
}

cov.pow.ex=function(hyper,Data=NULL,Data.new=NULL,v.power.f=NULL,gamma=1){
  #hyper is a list of hyper-parameters
  ##cv.power.f: list of covariance for functional part
  if(is.null(gamma)) gamma=1
  hyper=lapply(hyper,exp)
  v.power=0
if(!is.null(Data)) {v.power=v.power+xixj_sta(Data,Data.new,hyper$pow.ex.w,power=gamma)}
if(!is.null(v.power.f)){
  hyp.f=hyper$pow.ex.f
  tmp=mapply(function(x,y) x*y, v.power.f,hyp.f,SIMPLIFY=FALSE)
  tmp=(Reduce("+",tmp))
  v.power=v.power+tmp
  }
  ## tmp is too big!!
  exp.v.power=hyper$pow.ex.v*exp(-v.power/2)  #exp(-v.power/2-tmp/2)
  return(exp.v.power)
}

cov.rat.qu=function(hyper,Data,Data.new=NULL){
  #hyper is a list of hyper-parameters
  hyper=lapply(hyper,exp);
  datadim=dim(as.matrix(Data))
  
  v.power=xixj_sta(Data,Data.new,hyper$rat.qu.w,power=1)
  sa=20^(1/hyper$rat.qu.a)-1
  rat.qu.v=(1+sa*v.power)^(-hyper$rat.qu.a)
  return(rat.qu.v)
}

cov.von.mi=function(hyper,Data,Data.new=NULL){
  #hyper is a list of hyper-parameters
  hyper=lapply(hyper,exp);
  datadim=dim(as.matrix(Data))
  
  v.power=xixj_sta2(Data,Data.new,power=1)
  von.mi.v=hyper$von.mi.v*exp(v.power*hyper$von.mi.w)
  return(von.mi.v)
}

cov.matern=function(hyper,Data,Data.new=NULL){
  #hyper is a list of hyper-parameters
  hyper=lapply(hyper,exp);
  datadim=dim(as.matrix(Data))
  w=rep(1,datadim[2])
  v.power=xixj_sta(Data,Data.new,w,power=1)
  v.power=hyper$matern.w*((v.power)^0.5)
  matern.v=(1+v.power)*exp(-v.power)
  return(matern.v)
}


xixj=function(mat,mat.new=NULL,a=NULL){
  mat=as.matrix(mat)
  mdim=dim(mat)
  #   err=1
  
  if(is.null(mat.new)){
    #     err=0
    mat.new=mat
  }
  mat.new=as.matrix(mat.new)
  if(is.null(a))  a=rep(1,mdim[2])
  if(length(a)<mdim[2]) {
    a1=rep(1,mdim[2])
    a1[1:length(a)]=a
    a=a1;rm(a1)
    warning('number of "a" is less than the number of columns, use 1 as the missing "a"')
  }
  if(length(a)>mdim[2]) {
    a=a[1:mdim[2]]
    warning('number of "a" is more than the number of columns, omit the extra "a"')
  }
  
  aa=matrix(rep(a,mdim[1]),ncol=mdim[2],byrow=T)
  out=(aa*mat)%*%t(mat.new)
  return(out)
}

xixj.f=function(mat,mat.new=NULL,w=NULL){
## cov.linear part
# mat, mat.new: column 1st, t; row, s; column, covariate; t: observed times, s: curve times
  if(is.null(mat.new)) mat.new=mat
  n.s=length(table(mat[,1]))
  n.t=table(mat[,1])[[1]]
  tims=mat[1:n.s,1]
  if(is.null(w)) w=c(tims[1],diff(tims))
  aa=matrix(rep(w,n.t),n.t,n.s,byrow=T)

  n.s.new=length(table(mat.new[,1]))
  n.t.new=table(mat.new[,1])[[1]]
  p=dim(mat)[2]-1

  out=list()
  for(i in 1:p)
  {
    tmp=mat[,i+1]
    tmp1=mat.new[,i+1]

    mat1=matrix(tmp,n.t,n.s,byrow=TRUE)  ## t times s
    mat2=matrix(tmp1,n.t.new,n.s,byrow=TRUE)

    out[[i]]=mat2%*%t(aa*mat1)
  }

  return(out)
}

xixj_sta=function(mat,mat.new=NULL,w=NULL,power=NULL){
  mat=as.matrix(mat)
  if(is.null(mat.new)) mat.new=mat
  mat.new=as.matrix(mat.new)
  mdim=dim(mat);mdim.new=dim(mat.new)
  cov.=matrix(sapply(1:mdim[1],function(i) matrix(rep(mat[i,],mdim.new[1]),nrow=mdim.new[1],byrow=T)-mat.new),ncol=mdim[1])
  #mn<<-mat.new;mmm<<-mat

  if(is.null(power)) power=1
  cov.=((cov.)^2)^power;
  if(is.null(w)) {
    w=rep(1,mdim[2])
    warning('missing "weight", use 1 instead')
  }
  if(length(w)==1&mdim[2]>1){
    w=rep(w,mdim[2])
    warning('only one "weight" found, applied to all columns')
  }
  
  if(length(w)>1&length(w)<mdim[2]){
    w1=rep(1,mdim[2])
    w1[1:length(w)]=w
    w=w1;rm(w1)
    warning('number of "weight" is less than the number of columns, use 1 as the missing "weight"')
  }
  if(length(w)>mdim[2]){
    w=w[1:mdim[2]]
    warning('number of "weight" is more than the number of columns, omit the extra "weight"')
  }

  wmat=matrix(rep(w,each=dim(cov.)[1]*dim(cov.)[2]/mdim[2]),ncol=dim(cov.)[2],byrow=T)

  cov.=wmat*cov.

  cov..=matrix(0,ncol=mdim[1],nrow=mdim.new[1])
  if(mdim[2]>1){
    for(i in 1:(mdim[2]-1)){
      cov..=cov..+cov.[1:mdim.new[1],,drop=F];cov.=cov.[-(1:mdim.new[1]),,drop=F]}
    cov.=cov..+cov.
  }
  return(cov.)  
}

xixj_sta.f=function(mat,mat.new=NULL,w=NULL,power=NULL){
  ## where is the power
  if(is.null(mat.new)) mat.new=mat
  n.s=length(table(mat[,1]))## timelen of functional variable
  n.t=table(mat[,1])[[1]]  ## timelen of longitudial observation
  tims=mat[1:n.s,1]
  if(is.null(w)) w=c(tims[1],diff(tims))
  aa=matrix(rep(w,n.t),n.t,n.s,byrow=T)
  
  n.s.new=length(table(mat.new[,1]))
  n.t.new=table(mat.new[,1])[[1]]
  p=dim(mat)[2]-1

  out=list()
  for(i in 1:p) {
    tmp=mat[,i+1]
    tmp1=mat.new[,i+1]

    mat1=matrix(tmp,n.t,n.s,byrow=TRUE)
    mat2=matrix(tmp1,n.t.new,n.s,byrow=TRUE)
    cmat=matrix(0,n.t.new,n.t)
    for(k in 1:n.t.new)
        for(j in 1:n.t)
         cmat[k,j]=sum(w*(mat2[k,]-mat1[j,])^2) 

    out[[i]]=cmat
  }

  return(out)  
}


xixj_sta2=function(mat,mat.new=NULL,w=NULL,power=NULL){
  mat=as.matrix(mat)
  if(is.null(mat.new)) mat.new=mat
  mdim=dim(mat);mdim.new=dim(mat.new)
  cov.=matrix(sapply(1:mdim[1],function(i) matrix(rep(mat[i,],mdim.new[1]),nrow=mdim.new[1],byrow=T)-mat.new),ncol=mdim[1])
  #mn<<-mat.new;mmm<<-mat
  
  if(is.null(power)) power=1
  cov.=(cos(cov.))^power;
  if(is.null(w)) {
    w=rep(1,mdim[2])
   # warning('missing "weight", use 1 instead')
  }
  if(length(w)==1&mdim[2]>1){
    w=rep(w,mdim[2])
   # warning('only one "weight" found, applied to all columns')
  }
  
  if(length(w)>1&length(w)<mdim[2]){
    w1=rep(1,mdim[2])
    w1[1:length(w)]=w
    w=w1;rm(w1)
   # warning('number of "weight" is less than the number of columns, use 1 as the missing "weight"')
  }
  if(length(w)>mdim[2]){
    w=w[1:mdim[2]]
  #  warning('number of "weight" is more than the number of columns, omit the extra "weight"')
  }
  
  wmat=matrix(rep(w,each=dim(cov.)[1]*dim(cov.)[2]/mdim[2]),ncol=dim(cov.)[2],byrow=T)
  
  cov.=wmat*cov.
  
  cov..=matrix(0,ncol=mdim[1],nrow=mdim.new[1])
  if(mdim[2]>1){
    for(i in 1:(mdim[2]-1)){
      cov..=cov..+cov.[1:mdim.new[1],,drop=F];cov.=cov.[-(1:mdim.new[1]),,drop=F]}
    cov.=cov..+cov.
  }
  cov.=cov.-mdim[2]
  return(cov.)  
}

