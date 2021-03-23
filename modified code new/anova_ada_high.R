library(quadprog)
library(gss)
source('subfun.R')
# My_solve.QP is a modified version of solve.QP for constrained quadratic programming
My_solve.QP=function(Dmat,dvec,Amat,bvec)
{
  solution=tryCatch(solve.QP(Dmat,dvec,Amat,bvec)$solution,error=function(x) NA)
  if(is.na(solution[1]))
  {  M <- solve(Dmat) 
  Dmat <- t(M)%*%M
  sc <- norm(Dmat,"2")
  solution=tryCatch(solve.QP(Dmat = Dmat/sc, dvec=dvec/sc, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution,error=function(x) NA)
  if(is.na(solution[1]))
  {
    Dmat=diag(diag(Dmat))
    solution=tryCatch(solve.QP(Dmat,dvec,Amat,bvec)$solution,error=function(x) NA)
  }
  }
  if(is.na(solution[1]))
  {  M <- solve(Dmat) 
  Dmat <- t(M)%*%M
  nn2 = sqrt(norm(dvec,"2"))
  solution= tryCatch(solve.QP(Dmat = Dmat*nn2, dvec = dvec/(nn2^2), Amat = Amat, bvec = bvec, factorized=TRUE)$solution,error=function(x) NA)
  }
  return(solution)
}

#anova_l2 is the function to do L2 estimation 
anova_l2<-function(y,X,RK,w1=NULL,w2=NULL,maxit=100){
  qt=dim(RK)[3]
  n=length(y)
  mt=dim(X)[2]
  nm=dim(RK)[1]
  if (is.null(w1)){
    w1=rep(1,mt)
  }
  if (is.null(w2)){
    w2=rep(1,qt)
  }
  lambda1=NULL
  for(i in 1:qt)
  {q=RK[,,i]
  lambdatmp=sum(diag(q))
  lambda1=c(lambda1,log10(length(y)*lambdatmp))
  }
  
  limnla=c(-3,3)
  prec=1
  n.bf=qt+1
  ny=10
  f=matrix(0,n,n.bf)
  theta=rep(0,qt)
  j=0
  bf.lamb=0.001
  #backfitting procedure to estimate theta
  while((prec>bf.lamb) &  (j<maxit)){
    fold=f
    yw=y-apply(f[,-1],1,sum)
    tmp=lm(yw~X-1)
    f[,1]=tmp$fitted.values
    for(i in 2:n.bf)
    {
      yw=y-apply(f[,-i],1,sum)
      q=RK[,,(i-1)]
      lambda0=lambda1[i-1]+limnla
      res=dsidr(y = yw, q = q,limnla = lambda0)
      f[,i]=res$fit
      lambda=(10^res$nlaht)/length(yw)
      theta[i-1]=1/lambda
    }
    prec=mean(abs(f-fold))
    j=j+1
  }
  q=matrix(0,n,n)
  for(i in 1:(n.bf-1))
  {
    q=q+theta[i]*RK[,,i]
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
    res=dsidr(y = y, q = q, s = X)
    cmat=res$c
    dmat=res$d
    lambda=(10^res$nlaht)/length(y)
  }
  result=list(cmat=cmat,dmat=dmat,lambda=lambda,theta=theta)
  result
}

#anova_l1 is the function to perform L1 model selection and ada is the option to apply adaptive lasso or not
anova_l1<-function(y,X,RK,lambda1_list,M_list,w1=NULL,w2=NULL,tol1=0.001,tol2=0.001,maxit=100,k=10,r1=1,r2=1,ada=TRUE,lam_range=c(-3,3))
{ 
qt=dim(RK)[3]
n=length(y)
mt=dim(X)[2]
nm=dim(RK)[1]
# fit l2 estimation first to get good initial values
l2est<-anova_l2(y,X,RK,w1=NULL,w2=NULL,maxit=100)
d0=l2est$dmat
cmat=l2est$cmat
theta=l2est$theta
lambda=l2est$lambda
theta0=rep(1,qt)
for (i in 1:qt) {
  theta0[i]=sqrt(sum((theta[i]*RK[,,i]%*%cmat)^2))
}
M_cand=max(2*length(theta0),sum(theta0))
if (ada==TRUE){
  w0=rep(1,qt)
  for (i in 1:qt) {
    w0[i]=sqrt(sum((theta[i]*RK[,,i]%*%cmat)^2))
  }
  #set weights w1 and w2 based on L2 estimates
  for (i in 1:qt) {
    w2[i]=(1/w0[i])^r2
  }
  for (i in 1:length(d0)) {
    w1[i]=(1/abs(d0[i]))^r1
  }
  d0=rep(1,mt)
  theta0=rep(1,qt)
  M_cand=length(theta0)
}
# do not penalize cosntant
w1[1]=0
# Pre-select d, c and theta before iteration update
#Step 1 use c as in the L2 estimate
c_updated=cmat
tau0=lambda
#Step 2 update d using glmnet lasso

R=matrix(0,nm,qt)
for (i in 1:qt){
  R[,i]=RK[,,i]%*%c_updated
}
y_new=y-R%*%theta0
lasso<-cv.glmnet(X,y_new,type.measure = "mse", nfolds = 20,intercept = FALSE,parallel = TRUE,penalty.factor=w1)
d_updated=as.vector(coef(lasso,s = "lambda.1se"))[-1]
lambda1=lasso$lambda.min
# Step 3 update theta

M=tune.M(y,X,d_updated,R,w2,tau0,c_updated,M_cand,M_list,mt,qt,k)[[1]]
Hmatrix=2*t(R)%*%R
w=tau0*t(R)%*%c_updated
Dmat=Hmatrix
if (is.positive.definite(Dmat)==FALSE){
  Dmat=Dmat+10^(-6)*diag(dim(Dmat)[1])
}
dvec=2*t(R)%*%(y-X%*%d_updated)-w
Amat=rbind(diag(qt),-w2)
bvec=c(rep(0,qt),-M)
theta_updated=My_solve.QP(Dmat,dvec,t(Amat),bvec)
theta_updated[theta_updated<1e-8]=0

theta_mat=matrix(0,qt,maxit+1)
theta_mat[,1]=theta_updated
dmat=matrix(0,length(d_updated),maxit+1)
dmat[,1]=as.vector(d_updated)
diff_c=1
diff_d=1
diff_theta=1
iter=1
zero_d=1
zero_theta=1
while ((diff_d>tol1 |diff_theta>tol1)&(iter<maxit)&(zero_d>=1|zero_theta>=1)) {
  #Step 1 update c
  
  z=rep(0,n)
  z=y-X%*%d_updated
  Sigma=matrix(0,nm,nm)
  for (v in 1:qt) {
    Sigma=Sigma+(1/w2[v])*theta_updated[v]*RK[,,v]
  }
  res=dsidr(y=z,q=Sigma,limnla = lam_range)
  c_updated=res$c
  tau0=(10^res$nlaht)/length(z)
  # Step 2 update d using glmnet lasso
  
  R=matrix(0,nm,qt)
  for (i in 1:qt){
    R[,i]=(1/w2[i])*RK[,,i]%*%c_updated
  }
  y_new=y-R%*%theta_updated
  lasso<-cv.glmnet(X,y_new,type.measure = "mse", nfolds = 20,intercept = FALSE,parallel = TRUE,penalty.factor=w1)
  d_updated=as.vector(coef(lasso,s = "lambda.1se"))[-1]
  lambda1=lasso$lambda.min
  # Step 3 update theta
  
  # cross-validation to select M
  M=tune.M(y,X,d_updated,R,w2,tau0,c_updated,M_cand,M_list,mt,qt,k)[[1]]
  Hmatrix=2*t(R)%*%R
  w=tau0*t(R)%*%c_updated
  Dmat=Hmatrix
  if (is.positive.definite(Dmat)==FALSE){
    Dmat=Dmat+10^(-6)*diag(dim(Dmat)[1])
  }
  dvec=2*t(R)%*%(y-X%*%d_updated)-w
  Amat=rbind(diag(qt),-w2)
  bvec=c(rep(0,qt),-M)
  theta_updated=My_solve.QP(Dmat,dvec,t(Amat),bvec)
  iter=iter+1
  theta_mat[,iter]=theta_updated
  dmat[,iter]=as.vector(d_updated)
  # compute the criteria for convergence checking
  if (iter<=5){
    zero_d=1
    zero_theta=1
  }
  else {
    zero_d=sum(d_updated>1e-6)-sum(dmat[,iter-5]>1e-6)
    zero_theta=sum(theta_updated>1e-6)-sum(theta_mat[,iter-5]>1e-6)
  }
  diff_d=sqrt(sum((dmat[iter]-dmat[iter-1])^2))/(sqrt(sum((dmat[iter-1])^2))+1e-6)
  diff_theta=sqrt(sum((theta_mat[iter]-theta_mat[iter-1])^2))/(sqrt(sum((theta_mat[iter-1])^2))+1e-6)
}
theta_updated[theta_updated<1e-8]=0
fit_values=X%*%d_updated+R%*%theta_updated
l2error=mean((y-fit_values)^2)
list(fit=fit_values,l2error=l2error,R=R,c=c_updated,d=d_updated,theta=theta_updated,tau0=tau0,lambda1=lambda1,M=M,iter=iter,w1=w1,w2=w2)
}

# tune.M is the function to select M using cross-validation
tune.M<-function(y,X,d,R,w2,tau0,c_updated,M_cand,M_list,mt,qt,k){
  n=length(y)
  M_list1=M_cand*M_list
  l2_error_mat=rep(0,length(M_list1))
  count=0
  for (M in M_list1) {
    count=count+1
    w=tau0*t(c_updated)%*%R
    w=t(w)
    score=rep(0,k)
    grps=cut(1 : n, k, labels = FALSE)[sample(n)]
    for (kth in 1:k) {
      omit=which(grps == kth)
      R_omit=R[omit,]
      y_omit=y[omit]
      R_keep=R[-omit,]
      y_keep=y[-omit]
      Dmat=2*t(R_keep)%*%R_keep
      if (is.positive.definite(Dmat)==FALSE){
        Dmat=Dmat+10^(-6)*diag(dim(Dmat)[1])
      }
      dvec=-(2*t(-R_keep)%*%(y_keep-X[-omit,]%*%d)+w)
      Amat=rbind(diag(qt),-w2)
      bvec=c(rep(0,qt),-M)
      theta_updated=My_solve.QP(Dmat,dvec,t(Amat),bvec)
      theta_updated[theta_updated<1e-8]=0
      fit_omit=X[omit,]%*%d-R_omit%*%theta_updated
      score[kth]=mean((y_omit-fit_omit)^2)
    }
    l2_error_mat[count]=mean(score)
    
  }
  min_error_pos=which.min(l2_error_mat)
  M_opt=M_list1[min_error_pos]
  list(M_opt=M_opt)
}
# compute the TPR and FPR for ROC
ROCval0<-function(truePara, estimatedPara,m) {
  p <- length(truePara)
  if (m==1){
    estimatedPara[abs(estimatedPara)<1e-6]=0
    estimatedPara[abs(estimatedPara)>1e-6]=1
  }else{
    estimatedPara[estimatedPara<1e-6]=0
    estimatedPara[estimatedPara>1e-6]=1
  }
  # for Theta
  P <- sum(truePara)
  TP <- sum(truePara * estimatedPara)
  N <- p- P
  allOne <- rep(1,p)
  FP <- sum(estimatedPara * (allOne - truePara))
  # calculate TPR & FPR
  if (P > 0) {
    TPR <- TP / P
  } else {
    TPR <- 0
  }
  if (N > 0) {
    FPR <- FP / N
  } else {
    FPR <- 0
  }
  
  list(P=P, N=N, FPR=FPR, TPR=TPR)
}

# compute the SPE, SEN and F1 score
Calmetric<-function(perf)
{
  FPR <- perf$FPR; TPR <- perf$TPR; P <- perf$P; N <- perf$N;
  FP <- N * FPR
  TP <- P * TPR
  if(P == 0) P <- P + 1e-5
  if(TP == 0) TP <- TP + 1e-5
  FN <- P - TP
  TN <- N - FP
  mcc <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
  Sp <- TN / N
  Sn <- TP / P
  F1 <- 2* TP / (2*TP+FN+FP)
  list(mcc=mcc, Sp=Sp, Sn=Sn, F1=F1)
}
# compute the mean and standard deviation
getMeanAndSd <- function(data) {
  meanSp <- round(mean(data$Spe,na.rm=TRUE),3)
  meanSn <- round(mean(data$Sen,na.rm=TRUE),3)
  meanF1 <- round(mean(data$F1,na.rm=TRUE),3)
  sdSp <- round(sd(data$Spe,na.rm=TRUE),3)
  sdSn <- round(sd(data$Sen,na.rm=TRUE),3)
  sdF1 <- round(sd(data$F1,na.rm=TRUE),3)
  list(meanSp = meanSp, meanSn = meanSn, meanF1 = meanF1,
       sdSp = sdSp, sdSn = sdSn, sdF1 = sdF1)
}

