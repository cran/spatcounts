#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

void betanbind(double* mu, int* nb, double* beta1, double* gammavec, double* r, double* xb, int* ny, double* y, double* ex, int* acceptb, double* P0, double* beta_r)
{

int i, j, z;
double *blo = (double*) R_alloc(*nb,sizeof(double));
double *bup = (double*) R_alloc(*nb,sizeof(double));
double *bm = (double*) R_alloc(*nb,sizeof(double));
double *xblo = (double*) R_alloc(*ny,sizeof(double));
double *xbup = (double*) R_alloc(*ny,sizeof(double));
double *xbm = (double*) R_alloc(*ny,sizeof(double));
double *mulo = (double*) R_alloc(*ny,sizeof(double));
double *muup = (double*) R_alloc(*ny,sizeof(double));
double *mum = (double*) R_alloc(*ny,sizeof(double));
double *munew = (double*) R_alloc(*ny,sizeof(double));
double dllo;
double dlup;
double dlm;
double propcov;
double u;
double *betanew = (double*) R_alloc(*nb,sizeof(double));
double temp;
double lbeta1;

for (j=0; j < *nb; j++){
   blo[j] = beta1[j];
   bup[j] = beta1[j];
   bm[j] = beta1[j];

}
for (j=0; j < *nb; j++){
 /*Starting values for bisection*/
  if (j<3){
     blo[j]=-100;
     bup[j]=100;
  } else {
    blo[j]=beta1[j]-50;
    bup[j]=beta1[j]+50;
  }
   bm[j]=(blo[j]+bup[j])/2;

for (i=0; i < *ny; i++){
       xblo[i]=0;
       xbup[i]=0;
       xbm[i]=0;
}
for (i=0; i < *ny; i++){
for (z=0; z < *nb; z++){
       xblo[i]+=xb[i+z* *ny]*blo[z];
       xbup[i]+=xb[i+z* *ny]*bup[z];
       xbm[i]+=xb[i+z* *ny]*bm[z];
}
       mulo[i]=ex[i]*exp(xblo[i]+gammavec[i]);
       muup[i]=ex[i]*exp(xbup[i]+gammavec[i]);
       mum[i]=ex[i]*exp(xbm[i]+gammavec[i]);
}
dllo=0;
dlup=0;
dlm=0;
for (i=0; i < *ny; i++){
       dllo+=xb[i+j* *ny]*(y[i]-(*r+y[i])*mulo[i]/(mulo[i]+*r));
       dlup+=xb[i+j* *ny]*(y[i]-(*r+y[i])*muup[i]/(muup[i]+*r));
       dlm+=xb[i+j* *ny]*(y[i]-(*r+y[i])*mum[i]/(mum[i]+*r));
}
dllo=dllo-blo[j] * P0[j+j* *nb];
dlup=dlup-bup[j] * P0[j+j* *nb];
dlm=dlm-bm[j] * P0[j+j* *nb];

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

for (i=0; i < *ny; i++){
       xbm[i]=0;
}

for (i=0; i < *ny; i++){
for (z=0; z < *nb; z++){
       xbm[i]+=xb[i+z* *ny]*bm[z];
}
mum[i]=ex[i]*exp(xbm[i]+gammavec[i]);
}
dlm=0;
for (i=0; i < *ny; i++){
      dlm+=xb[i+j* *ny]*(y[i]-(*r+y[i])*mum[i]/(mum[i]+*r));
}
dlm=dlm-bm[j] * P0[j+j* *nb];
}

propcov=0;
for (i=0; i < *ny; i++){
    propcov+=( pow(xb[i+j* *ny],2)*(*r+y[i])* *r *mum[i]/pow(mum[i]+*r,2) );
}
propcov=1/(propcov+P0[j+j* *nb]);

GetRNGstate();
u=runif(0,1);
PutRNGstate();

for (z=0; z < *nb; z++){
betanew[z] = beta1[z];
}

GetRNGstate();
betanew[j]=(bm[j]+pow(21*propcov/20, 0.5)*rt(20));
PutRNGstate();
for (i=0; i < *ny; i++){
       xbm[i]=0;
}
for (i=0; i < *ny; i++){
for (z=0; z < *nb; z++){
       xbm[i]+=xb[i+z* *ny]*betanew[z];
}
munew[i]=ex[i]*exp(xbm[i]+gammavec[i]);
}
temp=0;
for (i=0; i< *ny; i++){
   temp+=(log((mu[i]+ *r)/(munew[i]+ *r))*( *r+y[i])+y[i]*(log(munew[i])-log(mu[i])));
}

/*N(0,P0^-1) Prior for beta*/
 temp = temp - 0.5*((pow(betanew[j],2) - pow(beta1[j],2))*P0[j+j* *nb]) + 10.5*log(1+pow(betanew[j]-bm[j],2)*(1/propcov)/20)-log(1+pow(beta1[j]-bm[j],2)*(1/propcov)/20);

lbeta1=(temp<0)*temp;

if ((log(u) < lbeta1) | (lbeta1 >= 0)){
for (z=0; z < *nb; z++){
  beta1[z] = betanew[z];
}
for (i=0; i < *ny; i++){
  mu[i] = munew[i];
}
  acceptb[j] = acceptb[j]+1;
} else {
for (z=0; z < *nb; z++){
  beta1[z] = beta1[z];
}
for (i=0; i < *ny; i++){
  mu[i] = mu[i];
}
}
}

for (z=0; z < *nb; z++){
beta_r[z] = beta1[z];
beta_r[*nb + z] = acceptb[z];
}
for (i=0; i < *ny; i++){
beta_r[2 * *nb + i] = mu[i];
}
}
