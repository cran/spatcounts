#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

void betaindbisection(double* mu, int* b_n, double* beta1, double* phi, int* y_n, double* gammavec, double* xb, double* y, double* ex, double* P0, int* acceptb, double* beta_r)
{
  /* printf("%f\n",xb[y_n+b_n* *b_n]);*/
int i, j, z;
double *blo = (double*) R_alloc(*b_n,sizeof(double));
double *bup = (double*) R_alloc(*b_n,sizeof(double));
double *bm = (double*) R_alloc(*b_n,sizeof(double));
double *xblo = (double*) R_alloc(*y_n,sizeof(double));
double *xbup = (double*) R_alloc(*y_n,sizeof(double));
double *xbm = (double*) R_alloc(*y_n,sizeof(double));
double *mulo = (double*) R_alloc(*y_n,sizeof(double));
double *muup = (double*) R_alloc(*y_n,sizeof(double));
double *mum = (double*) R_alloc(*y_n,sizeof(double));
double *munew = (double*) R_alloc(*y_n,sizeof(double));
double dllo;
double dlup;
double dlm;
double propcov;
double u;
double *betanew = (double*) R_alloc(*b_n,sizeof(double));
double temp;
double lbeta1;

for (j=0; j < *b_n; j++){
   blo[j] = beta1[j];
   bup[j] = beta1[j];
   bm[j] = beta1[j];

}
for (j=0; j < *b_n; j++){
 /*Starting values for bisection*/
  if (j<3){
     blo[j]=-100;
     bup[j]=100;
  } else {
    blo[j]=beta1[j]-50;
    bup[j]=beta1[j]+50;
  }
   bm[j]=(blo[j]+bup[j])/2;

for (i=0; i < *y_n; i++){
       xblo[i]=0;
       xbup[i]=0;
       xbm[i]=0;
}
for (i=0; i < *y_n; i++){
for (z=0; z < *b_n; z++){
       xblo[i]+=xb[i+z* *y_n]*blo[z];
       xbup[i]+=xb[i+z* *y_n]*bup[z];
       xbm[i]+=xb[i+z* *y_n]*bm[z];  
}
       mulo[i]=ex[i]*exp(xblo[i]+gammavec[i]);
       muup[i]=ex[i]*exp(xbup[i]+gammavec[i]);
       mum[i]=ex[i]*exp(xbm[i]+gammavec[i]);
}
dllo=0;
dlup=0;
dlm=0;
for (i=0; i < *y_n; i++){
       dllo+=xb[i+j* *y_n]*(1+(y[i]-1)*mulo[i]/(mulo[i]+(*phi-1) * y[i])-mulo[i] / *phi);
       dlup+=xb[i+j* *y_n]*(1+(y[i]-1)*muup[i]/(muup[i]+(*phi-1) * y[i])-muup[i] / *phi);
       dlm+=xb[i+j* *y_n]*(1+(y[i]-1)*mum[i]/(mum[i]+(*phi-1) * y[i])-mum[i] / *phi);
}
dllo=dllo-blo[j] * P0[j+j* *b_n];
dlup=dlup-bup[j] * P0[j+j* *b_n];
dlm=dlm-bm[j] * P0[j+j* *b_n];

/*if ((dllo<0) || (dlup>0)){
     printf(" dllo<0 or dlup>0\n");
     }*/

   while(fabs(bup[j] - blo[j]) > 0.01*fabs(bup[j])) {

   if (dlm<0){
      dlup=dlm;
      bup[j]=bm[j];
      bm[j]=(blo[j]+bup[j])/2;
   } else {
      dllo=dlm;
      blo[j]=bm[j];
      bm[j]=(blo[j]+bup[j])/2;
   }

for (i=0; i < *y_n; i++){
       xbm[i]=0;
}

for (i=0; i < *y_n; i++){
for (z=0; z < *b_n; z++){
       xbm[i]+=xb[i+z* *y_n]*bm[z];
}
mum[i]=ex[i]*exp(xbm[i]+gammavec[i]);
}
dlm=0;
for (i=0; i < *y_n; i++){
      dlm+=xb[i+j* *y_n]*(1+(y[i]-1)*mum[i]/(mum[i]+(*phi-1)*y[i])-mum[i] / *phi);
}
dlm=dlm-bm[j] * P0[j+j* *b_n];
}

propcov=0;
for (i=0; i < *y_n; i++){
    propcov+=(-(pow(xb[i+j* *y_n],2)*mum[i]*((y[i]-1)*(*phi-1)*(y[i]/pow(mum[i]+(*phi-1)*y[i],2))-1/ *phi) ));
}
propcov=1/(propcov+P0[j+j* *b_n]);

GetRNGstate();
u=runif(0,1);
PutRNGstate();

for (z=0; z < *b_n; z++){
betanew[z] = beta1[z];
}

GetRNGstate();
betanew[j]=(bm[j]+pow(21*propcov/20, 0.5)*rt(20));
PutRNGstate();
for (i=0; i < *y_n; i++){
       xbm[i]=0;
}
for (i=0; i < *y_n; i++){
for (z=0; z < *b_n; z++){
       xbm[i]+=xb[i+z* *y_n]*betanew[z];
}
munew[i]=ex[i]*exp(xbm[i]+gammavec[i]);
}
temp=0;
for (i=0; i< *y_n; i++){
   temp+=(log(munew[i]/mu[i])+(y[i]-1)*log((munew[i]+y[i]*(*phi-1))/(mu[i]+y[i]*(*phi-1)))-1 / *phi*(munew[i]-mu[i]));
}

/*N(0,P0^-1) Prior for beta*/
 temp = temp - 0.5*((pow(betanew[j],2) - pow(beta1[j],2))*P0[j+j* *b_n]) + 10.5*log(1+pow(betanew[j]-bm[j],2)*(1/propcov)/20)-log(1+pow(beta1[j]-bm[j],2)*(1/propcov)/20);

lbeta1=(temp<0)*temp;

if ((log(u) < lbeta1) | (lbeta1 >= 0)){
for (z=0; z < *b_n; z++){
  beta1[z] = betanew[z];
}
for (i=0; i < *y_n; i++){
  mu[i] = munew[i];
}
  acceptb[j] = acceptb[j]+1;
} else {
for (z=0; z < *b_n; z++){
  beta1[z] = beta1[z];
}
for (i=0; i < *y_n; i++){
  mu[i] = mu[i];
}
}
}

for (z=0; z < *b_n; z++){
beta_r[z] = beta1[z];
beta_r[*b_n + z] = acceptb[z];
}
for (i=0; i < *y_n; i++){
beta_r[2 * *b_n + i] = mu[i];
}
}
