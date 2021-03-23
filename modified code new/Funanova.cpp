#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double numIntegrate(NumericVector x, NumericVector z) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < (n-1); ++i) {
    total += (x[i]+x[i+1])/2*(z[i+1]-z[i]);
  }
  return total;
}

// [[Rcpp::export]]
NumericVector k0xc(NumericVector x) {
  int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i]= 1;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector k1xc(NumericVector x) {
  int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i]= x[i]-0.5;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector k2xc(NumericVector x) {
  int n = x.size();
  NumericVector out(n);
  NumericVector tmp(1);
  for(int i = 0; i < n; ++i) {
    tmp[0] = x[i];
    tmp = k1xc(tmp);
    out[i]= (tmp[0] * tmp[0] - 1.0/12.0)/2;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector k4xc(NumericVector x) {
  int n = x.size();
  NumericVector out(n);
  NumericVector tmp(1);
  double tmp1;
  double tmp2;
  for(int i = 0; i < n; ++i) {
    tmp[0] = x[i];
    tmp = k1xc(tmp);
	tmp1 = pow(tmp[0],4);
	tmp2 = pow(tmp[0],2);
    out[i]= (tmp1 - tmp2/2.0 + 7.0/240.0)/24.0;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector phi1xc(NumericVector x) {
  int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i]= x[i]-0.5;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector phi2xc(NumericVector x) {
  NumericVector out(1);
  out[0]=(x[0]-0.5)*(x[1]-0.5);
  return out;
}

// [[Rcpp::export]]
NumericMatrix R01c(NumericVector x, NumericVector z) {
  int nx = x.size();
  int nz = z.size();
  NumericMatrix out(nx,nz);
  for(int i = 0; i < nx; ++i) {
    for(int j = 0; j < nz; ++j) {
    out(i, j)= (x[i]-0.5) * (z[j]-0.5);
   }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix R1c(NumericVector x, NumericVector z) {
  int nx = x.size();
  int nz = z.size();
  double tmp;
  NumericMatrix out(nx,nz);
  NumericVector xi(1), zj(1);

  for(int i = 0; i < nx; ++i) {
    for(int j = 0; j < nz; ++j) {
    xi[0] = x[i];
	  zj[0] = z[j];
	  tmp=xi[0]-zj[0];
	  xi[0] = k2xc(xi)[0] * k2xc(zj)[0];
	  zj[0] = fabs(tmp);
    zj[0] = k4xc(zj)[0];
    out(i, j)= xi[0] - zj[0];
   }
  }
  return out;
}

// [[Rcpp::export]]
double R12c(NumericVector x) {
  double total = 0;
  NumericVector tmp0(1), tmp1(1), tmp2(1), tmp3(1);
  tmp0[0]=x[0];
  tmp1[0]=x[1];
  tmp2[0]=x[2];
  tmp3[0]=x[3];
  
  tmp0[0]=R01c(tmp0, tmp2)(0,0);
  tmp1[0]=R1c(tmp1, tmp3)(0,0);
  total=tmp0[0]*tmp1[0];
  return total;
}



// [[Rcpp::export]]
double R21c(NumericVector x) {
  double total = 0;
  NumericVector tmp0(1), tmp1(1), tmp2(1), tmp3(1);
  tmp0[0]=x[0];
  tmp1[0]=x[1];
  tmp2[0]=x[2];
  tmp3[0]=x[3];
  
  tmp0[0]=R1c(tmp0, tmp2)(0,0);
  tmp1[0]=R01c(tmp1, tmp3)(0,0);
  total=tmp0[0]*tmp1[0];
  return total;
}

// [[Rcpp::export]]
double R22c(NumericVector x) {
  double total = 0;
  NumericVector tmp0(1), tmp1(1), tmp2(1), tmp3(1);
  tmp0[0]=x[0];
  tmp1[0]=x[1];
  tmp2[0]=x[2];
  tmp3[0]=x[3];
  
  tmp0[0]=R1c(tmp0, tmp2)(0,0);
  tmp1[0]=R1c(tmp1, tmp3)(0,0);
  total=tmp0[0]*tmp1[0];
  return total;
}

// L1 estimation for simulation
// [[Rcpp::export]]
double gestimation124cl1(NumericVector xtest,  NumericVector x, NumericMatrix Xt, NumericMatrix Psi, NumericVector theta, NumericVector dmat, NumericVector cmat, NumericVector w2) {
  int m = Xt.nrow(), n = Xt.ncol(), ny = Psi.ncol();
  
  double total, temp, xij, y1x, xy1x, xix1;
  NumericVector u(1), v(1), xu(1), xv(1), y1(m), xi(m);
  NumericVector tmp0(2), tmp(m), tmp1(4), tmp2(m), tmp3(m);
  NumericMatrix x1(1,m);
  
  u[0]=xtest[0];
  v[0]=xtest[1];
  xu[0]=xtest[2];
  xv[0]=xtest[3];
  
  total=0;
  total=total+dmat[0];
  total=total+phi1xc(u)[0]*dmat[1];
  total=total+phi1xc(xu)[0]*dmat[2];
  total=total+phi1xc(xv)[0]*dmat[3];
  tmp0[0]=u[0];
  tmp0[1]=xu[0];
  total=total+phi2xc(tmp0)[0]*dmat[4];
  tmp0[0]=u[0];
  tmp0[1]=xv[0];
  total=total+phi2xc(tmp0)[0]*dmat[5];
  tmp0[0]=xu[0];
  tmp0[1]=xv[0];
  total=total+phi2xc(tmp0)[0]*dmat[6];
  tmp0[0]=v[0];
  tmp0[1]=xv[0];
  total=total+phi2xc(tmp0)[0]*dmat[7];
  
  xij=0;
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < ny; ++j) {
      for(int k=0; k < m; ++k){
        y1[k]=Psi(k,j);
        xi[k]=Xt(k,i);
      }
      x1=R1c(u,x);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      temp=numIntegrate(tmp,x)*theta[0]/w2[0]; // 1s
      
      x1=R1c(xu,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      temp=temp+numIntegrate(tmp,x)*theta[1]/w2[1];  //3s
      
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      y1x=numIntegrate(y1,x);
      temp=temp+numIntegrate(tmp,x)*y1x*theta[2]/w2[2]; //4s
      
      tmp1[0]=u[0];
      tmp1[1]=xu[0];
      for(int k=0; k < m; ++k){  
        tmp1[2]=x[k];
        tmp1[3]=xi[k];
        tmp[k]=R12c(tmp1)*y1[k];
        tmp2[k]=R21c(tmp1)*y1[k];
        tmp3[k]=R22c(tmp1)*y1[k];
      }
      temp=temp+(numIntegrate(tmp,x)+numIntegrate(tmp2,x)+numIntegrate(tmp3,x))*theta[3]/w2[3]; //1 3
      
      tmp=phi1xc(x);
      for(int k=0; k < m; ++k){  
        tmp[k]=tmp[k]*y1[k];
      }
      xy1x=numIntegrate(tmp,x);
      tmp1[0]=phi1xc(u)[0]*xy1x;
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[0]=numIntegrate(tmp,x)*tmp1[0]; // 1 4 ls
      
      x1=R1c(u,x);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      tmp1[3]=numIntegrate(tmp,x);
      tmp=phi1xc(xi);
      xix1=numIntegrate(tmp,x);
      tmp1[1]=phi1xc(xv)[0]*xix1*tmp1[3]; // 1 4 sl
      
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[2]=numIntegrate(tmp,x)*tmp1[3]; // 1 4 ss
      temp=temp+(tmp1[0]+tmp1[1]+tmp1[2])*theta[4]/w2[4]; // 1 4 
      
      tmp1[0]=v[0];
      tmp1[1]=xv[0];
      for(int k=0; k < m; ++k){  
        tmp1[2]=x[k];
        tmp1[3]=xi[k];
        tmp[k]=R12c(tmp1);
        tmp2[k]=R21c(tmp1);
        tmp3[k]=R22c(tmp1);
      }
      tmp1[0]=numIntegrate(tmp,x)*y1x;
      tmp1[1]=numIntegrate(tmp2,x)*y1x;
      tmp1[2]=numIntegrate(tmp3,x)*y1x;
      temp=temp+(tmp1[0]+tmp1[1]+tmp1[2])*theta[5]/w2[5]; // 2  4
      
      tmp=phi1xc(xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=tmp[k]*y1[k];
      }
      tmp1[0]=numIntegrate(tmp,x)*phi1xc(xu)[0];
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[0]=numIntegrate(tmp,x)*tmp1[0];  // 3  4 ls
      
      x1=R1c(xu,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      tmp1[1]=numIntegrate(tmp,x)*xix1*phi1xc(xv)[0]; // 3 4 sl
      
      x1=R1c(xu,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      tmp1[2]=numIntegrate(tmp,x);
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[2]=numIntegrate(tmp,x)*tmp1[2];  // 3 4 ss 
      temp=temp+(tmp1[0]+tmp1[1]+tmp1[2])*theta[6]/w2[6]; // 3 4
      
      
      xij=xij+temp*cmat[(i*ny)+j];
    }
  }
  total=total+xij;
  return total;
}


// L1 prediciton for simulation
// [[Rcpp::export]]
NumericVector Funpred124l1(NumericVector x1test, NumericVector x2test, NumericVector x1int, NumericVector x2int, NumericVector x, NumericMatrix Xt, NumericMatrix Psi, NumericVector theta, NumericVector dmat, NumericVector cmat, NumericVector w2) {
  int x1n = x1test.size(), x1m = x1int.size();
  NumericVector xtest(4), tmp(x1m), res(x1n);
  
  for(int i = 0; i < x1n; ++i) {
    for(int j = 0; j < x1m; ++j) {
      
      xtest[0]=x1test[i];
      xtest[1]=x1int[j];
      xtest[2]=x2test[i];
      xtest[3]=x2int[j];
      tmp[j]=gestimation124cl1(xtest,x,Xt,Psi,theta,dmat,cmat,w2);
    }
    
    res[i]=numIntegrate(tmp,x1int);
    
  }
  
  return res;
}


// L2 estimation for simulaiton
// [[Rcpp::export]]
double gestimation124cl2(NumericVector xtest,  NumericVector x, NumericMatrix Xt, NumericMatrix Psi, NumericVector theta, NumericVector dmat, NumericVector cmat) {
  int m = Xt.nrow(), n = Xt.ncol(), ny = Psi.ncol();
  
  double total, temp, xij, y1x, xy1x, xix1;
  NumericVector u(1), v(1), xu(1), xv(1), y1(m), xi(m);
  NumericVector tmp0(2), tmp(m), tmp1(4), tmp2(m), tmp3(m);
  NumericMatrix x1(1,m);
  
  u[0]=xtest[0];
  v[0]=xtest[1];
  xu[0]=xtest[2];
  xv[0]=xtest[3];
  
  total=0;
  total=total+dmat[0];
  total=total+phi1xc(u)[0]*dmat[1];
  total=total+phi1xc(xu)[0]*dmat[2];
  total=total+phi1xc(xv)[0]*dmat[3];
  tmp0[0]=u[0];
  tmp0[1]=xu[0];
  total=total+phi2xc(tmp0)[0]*dmat[4];
  tmp0[0]=u[0];
  tmp0[1]=xv[0];
  total=total+phi2xc(tmp0)[0]*dmat[5];
  tmp0[0]=xu[0];
  tmp0[1]=xv[0];
  total=total+phi2xc(tmp0)[0]*dmat[6];
  tmp0[0]=v[0];
  tmp0[1]=xv[0];
  total=total+phi2xc(tmp0)[0]*dmat[7];
  
  xij=0;
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < ny; ++j) {
      for(int k=0; k < m; ++k){
        y1[k]=Psi(k,j);
        xi[k]=Xt(k,i);
      }
      x1=R1c(u,x);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      temp=numIntegrate(tmp,x)*theta[0]; // 1s
      
      x1=R1c(xu,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      temp=temp+numIntegrate(tmp,x)*theta[1];  //3s
      
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      y1x=numIntegrate(y1,x);
      temp=temp+numIntegrate(tmp,x)*y1x*theta[2]; //4s
      
      tmp1[0]=u[0];
      tmp1[1]=xu[0];
      for(int k=0; k < m; ++k){  
        tmp1[2]=x[k];
        tmp1[3]=xi[k];
        tmp[k]=R12c(tmp1)*y1[k];
        tmp2[k]=R21c(tmp1)*y1[k];
        tmp3[k]=R22c(tmp1)*y1[k];
      }
      temp=temp+(numIntegrate(tmp,x)+numIntegrate(tmp2,x)+numIntegrate(tmp3,x))*theta[3]; //1 3
      
      tmp=phi1xc(x);
      for(int k=0; k < m; ++k){  
        tmp[k]=tmp[k]*y1[k];
      }
      xy1x=numIntegrate(tmp,x);
      tmp1[0]=phi1xc(u)[0]*xy1x;
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[0]=numIntegrate(tmp,x)*tmp1[0]; // 1 4 ls
      
      x1=R1c(u,x);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      tmp1[3]=numIntegrate(tmp,x);
      tmp=phi1xc(xi);
      xix1=numIntegrate(tmp,x);
      tmp1[1]=phi1xc(xv)[0]*xix1*tmp1[3]; // 1 4 sl
      
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[2]=numIntegrate(tmp,x)*tmp1[3]; // 1 4 ss
      temp=temp+(tmp1[0]+tmp1[1]+tmp1[2])*theta[4]; // 1 4 
      
      tmp1[0]=v[0];
      tmp1[1]=xv[0];
      for(int k=0; k < m; ++k){  
        tmp1[2]=x[k];
        tmp1[3]=xi[k];
        tmp[k]=R12c(tmp1);
        tmp2[k]=R21c(tmp1);
        tmp3[k]=R22c(tmp1);
      }
      tmp1[0]=numIntegrate(tmp,x)*y1x;
      tmp1[1]=numIntegrate(tmp2,x)*y1x;
      tmp1[2]=numIntegrate(tmp3,x)*y1x;
      temp=temp+(tmp1[0]+tmp1[1]+tmp1[2])*theta[5]; // 2  4
      
      tmp=phi1xc(xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=tmp[k]*y1[k];
      }
      tmp1[0]=numIntegrate(tmp,x)*phi1xc(xu)[0];
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[0]=numIntegrate(tmp,x)*tmp1[0];  // 3  4 ls
      
      x1=R1c(xu,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      tmp1[1]=numIntegrate(tmp,x)*xix1*phi1xc(xv)[0]; // 3 4 sl
      
      x1=R1c(xu,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      tmp1[2]=numIntegrate(tmp,x);
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[2]=numIntegrate(tmp,x)*tmp1[2];  // 3 4 ss 
      temp=temp+(tmp1[0]+tmp1[1]+tmp1[2])*theta[6]; // 3 4
      
      
      xij=xij+temp*cmat[(i*ny)+j];
    }
  }
  total=total+xij;
  return total;
}

// L2 prediction for simulation
// [[Rcpp::export]]
NumericVector Funpred124l2(NumericVector x1test, NumericVector x2test, NumericVector x1int, NumericVector x2int, NumericVector x, NumericMatrix Xt, NumericMatrix Psi, NumericVector theta, NumericVector dmat, NumericVector cmat) {
  int x1n = x1test.size(), x1m = x1int.size();
  NumericVector xtest(4), tmp(x1m), res(x1n);
  
  for(int i = 0; i < x1n; ++i) {
    for(int j = 0; j < x1m; ++j) {
      
      xtest[0]=x1test[i];
      xtest[1]=x1int[j];
      xtest[2]=x2test[i];
      xtest[3]=x2int[j];
      tmp[j]=gestimation124cl2(xtest,x,Xt,Psi,theta,dmat,cmat);
    }
    
    res[i]=numIntegrate(tmp,x1int);
    
  }
  
  return res;
}

// L2 EFLM  estimation for simulation
// [[Rcpp::export]]
double gestimationEFLM124cl2(NumericVector xtest,  NumericVector x, NumericMatrix Xt, NumericMatrix Psi, NumericVector theta, NumericVector dmat, NumericVector cmat) {
  int m = Xt.nrow(), n = Xt.ncol(), ny = Psi.ncol();
  
  double total, temp, xij, y1x, xy1x, xix1;
  NumericVector u(1), v(1), xv(1), y1(m), xi(m);
  NumericVector tmp0(2), tmp(m), tmp1(4), tmp2(m), tmp3(m);
  NumericMatrix x1(1,m);
  
  u[0]=xtest[0];
  v[0]=xtest[1];
  xv[0]=xtest[2];
  
  total=0;
  total=total+dmat[0];
  total=total+phi1xc(u)[0]*dmat[1];
  total=total+phi1xc(xv)[0]*dmat[2];
  tmp0[0]=u[0];
  tmp0[1]=xv[0];
  total=total+phi2xc(tmp0)[0]*dmat[3];
  tmp0[0]=v[0];
  tmp0[1]=xv[0];
  total=total+phi2xc(tmp0)[0]*dmat[4];
  
  xij=0;
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < ny; ++j) {
      for(int k=0; k < m; ++k){
        y1[k]=Psi(k,j);
        xi[k]=Xt(k,i);
      }
      x1=R1c(u,x);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      temp=numIntegrate(tmp,x)*theta[0]; // 1 s
      
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      y1x=numIntegrate(y1,x);
      temp=temp+numIntegrate(tmp,x)*y1x*theta[1];// 4 s
      
      tmp=phi1xc(x);
      for(int k=0; k < m; ++k){  
        tmp[k]=tmp[k]*y1[k];
      }
      xy1x=numIntegrate(tmp,x);
      tmp1[0]=phi1xc(u)[0]*xy1x;
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[0]=numIntegrate(tmp,x)*tmp1[0];
      //temp=temp+numIntegrate(tmp,x)*tmp1[0]*theta[2];// 1 4 ls
      
      x1=R1c(u,x);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k)*y1[k];
      }
      tmp1[3]=numIntegrate(tmp,x);
      tmp=phi1xc(xi);
      xix1=numIntegrate(tmp,x);
      tmp1[1]=phi1xc(xv)[0]*xix1*tmp1[3];
      // temp=temp+phi1xc(xv)[0]*xix1*tmp1[3]*theta[3];// 1 4 sl
      
      x1=R1c(xv,xi);
      for(int k=0; k < m; ++k){  
        tmp[k]=x1(0,k);
      }
      tmp1[2]=numIntegrate(tmp,x)*tmp1[3]; // 1 4 ss
      temp=temp+tmp1[0]*theta[2]+tmp1[1]*theta[2]+tmp1[2]*theta[2];
      // temp=temp+tmp1[1]*theta[2]+(tmp1[0]+tmp1[2])*theta[3];
      
      tmp1[0]=v[0];
      tmp1[1]=xv[0];
      for(int k=0; k < m; ++k){  
        tmp1[2]=x[k];
        tmp1[3]=xi[k];
        tmp[k]=R12c(tmp1);
        tmp2[k]=R21c(tmp1);
        tmp3[k]=R22c(tmp1);
      }
      tmp1[0]=numIntegrate(tmp,x)*y1x;
      tmp1[1]=numIntegrate(tmp2,x)*y1x;
      tmp1[2]=numIntegrate(tmp3,x)*y1x;
      temp=temp+tmp1[0]*theta[3]+tmp1[1]*theta[3]+tmp1[2]*theta[3];
      // temp=temp+(tmp1[0]+tmp1[1]+tmp1[2])*theta[4];
      
      xij=xij+temp*cmat[(i*ny)+j];
    }
  }
  total=total+xij;
  return total;
}

// L2 EFLM prediction for simulation
// [[Rcpp::export]]
NumericVector FunpredEFLM124l2(NumericVector x1test, NumericVector x1int, NumericVector x2int, NumericVector x, NumericMatrix Xt, NumericMatrix Psi, NumericVector theta, NumericVector dmat, NumericVector cmat) {
  int x1n = x1test.size(), x1m = x1int.size();
  NumericVector xtest(3), tmp(x1m), res(x1n);
  
  for(int i = 0; i < x1n; ++i) {
    for(int j = 0; j < x1m; ++j) {
      
      xtest[0]=x1test[i];
      xtest[1]=x1int[j];
      xtest[2]=x2int[j];
      tmp[j]=gestimationEFLM124cl2(xtest,x,Xt,Psi,theta,dmat,cmat);
    }
    
    res[i]=numIntegrate(tmp,x1int);
    
  }
  
  return res;
}

// L2 NCM estimation for simulation
// [[Rcpp::export]]
double gestimationNCM124cl2(NumericVector xtest,  NumericVector x, NumericMatrix Xt, NumericMatrix Psi, NumericVector theta, NumericVector dmat, NumericVector cmat) {
  int m = Xt.nrow(), n = Xt.ncol(), ny = Psi.ncol();
  
  double total, temp, xij;
  NumericVector u(1), xu(1), y1(m), xi(m);
  NumericVector tmp0(2), tmp(m), tmp1(4), tmp2(m), tmp3(m);
  NumericMatrix x1(1,m);
  
  u[0]=xtest[0];
  xu[0]=xtest[1];
  
  total=0;
  total=total+dmat[0];
  total=total+phi1xc(u)[0]*dmat[1];
  total=total+phi1xc(xu)[0]*dmat[2];
  
  tmp0[0]=u[0];
  tmp0[1]=xu[0];
  total=total+phi2xc(tmp0)[0]*dmat[3];
  
  xij=0;
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < ny; ++j) {
      for(int k=0; k < m; ++k){
        y1[k]=Psi(k,j);
        xi[k]=Xt(k,i);
      }
      x1=R1c(u,x);
      for(int k=0; k < m; ++k){
        tmp[k]=x1(0,k)*y1[k];
      }
      temp=numIntegrate(tmp,x)*theta[0]; // 1s
      
      x1=R1c(xu,xi);
      for(int k=0; k < m; ++k){
        tmp[k]=x1(0,k)*y1[k];
      }
      temp=temp+numIntegrate(tmp,x)*theta[1];  //3s
      
      tmp1[0]=u[0];
      tmp1[1]=xu[0];
      for(int k=0; k < m; ++k){
        tmp1[2]=x[k];
        tmp1[3]=xi[k];
        tmp[k]=R12c(tmp1)*y1[k];
        tmp2[k]=R21c(tmp1)*y1[k];
        tmp3[k]=R22c(tmp1)*y1[k];
      }
      temp=temp+numIntegrate(tmp,x)*theta[2];
      temp=temp+numIntegrate(tmp2,x)*theta[2]; 
      temp=temp+numIntegrate(tmp3,x)*theta[2]; //1 3
      
      xij=xij+temp*cmat[(i*ny)+j];
    }
  }
  total=total+xij;
  return total;
}



// L2 NCM prediciton for simulation
// [[Rcpp::export]]
NumericVector FunpredNCM124l2(NumericVector x1test, NumericVector x2test, NumericVector x, NumericMatrix Xt, NumericMatrix Psi, NumericVector theta, NumericVector dmat, NumericVector cmat) {
  int x1n = x1test.size();
  NumericVector xtest(2), res(x1n);
  
  for(int i = 0; i < x1n; ++i) {
    
    xtest[0]=x1test[i];
    xtest[1]=x2test[i];
    
    res[i]=gestimationNCM124cl2(xtest,x,Xt,Psi,theta,dmat,cmat);
    
  }
  
  return res;
}
