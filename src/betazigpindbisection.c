#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

void betazigpindbisection(double* mu, int* nb, double* beta1, double* phi, int* ny, double* gammavec, double* x, int* zfull, double* y, double* ex, double* P0, int* acceptb, double* beta_r)
{

int i, j, z, ny0;
double *y0 = (double*) R_alloc(*ny,sizeof(double));
double *ex0 = (double*) R_alloc(*ny,sizeof(double));
double *mu0 = (double*) R_alloc(*ny,sizeof(double));
double *g0 = (double*) R_alloc(*ny,sizeof(double)); 
double *x0b = (double*) R_alloc(*ny*(*nb),sizeof(double));
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

ny0=0;
for(z=0; z < *ny; z++){
if (zfull[z]==0){
y0[ny0]=y[z];
ex0[ny0]=ex[z];
mu0[ny0]=mu[z];
g0[ny0]=gammavec[z];
for(i=0; i < *nb; i++){
x0b[i+ny0*(*nb)]=x[*ny + (z + (i * *ny))];
}

ny0 += 1;
}
}

for (j=0; j < *nb; j++){
   blo[j] = beta1[j];
   bup[j] = beta1[j];
   bm[j] = beta1[j];

}
for (j=1; j < *nb; j++){
 /*Starting values for bisection*/
  if (j<3){
     blo[j]=-100;
     bup[j]=100;
  } else {
    blo[j]=beta1[j]-50;
    bup[j]=beta1[j]+50;
  }
   bm[j]=(blo[j]+bup[j])/2;
   
for (i=0; i < ny0; i++){
       xblo[i]=0;
       xbup[i]=0;
       xbm[i]=0;
}
for (i=0; i < ny0; i++){
for (z=0; z < *nb; z++){
       xblo[i]+=x0b[z+(i * *nb)] * blo[z];
       xbup[i]+=x0b[z+(i * *nb)]*bup[z];
       xbm[i]+=x0b[z+(i * *nb)]*bm[z];
}
       mulo[i]=ex0[i]*exp(xblo[i]+g0[i]);
       muup[i]=ex0[i]*exp(xbup[i]+g0[i]);
       mum[i]=ex0[i]*exp(xbm[i]+g0[i]);
}
dllo=0;
dlup=0;
dlm=0;
for (i=0; i < ny0; i++){
       dllo+=x0b[j+(i * *nb)]*(1+(y0[i]-1)*mulo[i]/(mulo[i]+(*phi-1) * y0[i])-mulo[i] / *phi);
       dlup+=x0b[j+(i * *nb)]*(1+(y0[i]-1)*muup[i]/(muup[i]+(*phi-1) * y0[i])-muup[i] / *phi);
       dlm+=x0b[j+(i * *nb)]*(1+(y0[i]-1)*mum[i]/(mum[i]+(*phi-1) * y0[i])-mum[i] / *phi);
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

for (i=0; i < ny0; i++){
       xbm[i]=0;
}
for (i=0; i < ny0; i++){
for (z=0; z < *nb; z++){
       xbm[i]+=x0b[z+(i * *nb)]*bm[z];
}
       mum[i]=ex0[i]*exp(xbm[i]+g0[i]);
}
dlm=0;
for (i=0; i < ny0; i++){
       dlm+=x0b[j+(i * *nb)]*(1+(y0[i]-1)*mum[i]/(mum[i]+(*phi-1) * y0[i])-mum[i] / *phi);
}
dlm=dlm-bm[j] * P0[j+j* *nb];
}

propcov=0;
for (i=0; i < ny0; i++){
    propcov+=(-(pow(x0b[j+i* *nb],2)*mum[i]*((y0[i]-1)*(*phi-1)*y0[i]/(pow(mum[i]+(*phi-1)*y[i],2))-1/ *phi) ));
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

for (i=0; i < ny0; i++){
       xbm[i]=0;
}
for (i=0; i < ny0; i++){
for (z=0; z < *nb; z++){
       xbm[i]+=x0b[(z+(i * *nb))]*betanew[z];
}
munew[i]=ex0[i]*exp(xbm[i]+g0[i]);
}

temp=0;
for (i=0; i< ny0; i++){
   temp+=(log(munew[i]/mu0[i])+(y0[i]-1)*log((munew[i]+y0[i]*(*phi-1))/(mu0[i]+y0[i]*(*phi-1)))-1 / *phi*(munew[i]-mu0[i]));
}

/*N(0,P0^-1) Prior for beta*/
 temp = temp - 0.5*((pow(betanew[j],2) - pow(beta1[j],2))*P0[j+j* *nb]) + 10.5*log(1+pow(betanew[j]-bm[j],2)*(1/propcov)/20)-log(1+pow(beta1[j]-bm[j],2)*(1/propcov)/20);

lbeta1=(temp<0)*temp;

if ((log(u) < lbeta1) | (lbeta1 >= 0)){
for (z=0; z < *nb; z++){
  beta1[z] = betanew[z];
}
for (i=0; i < *ny; i++){
  mu0[i] = munew[i];
}
  acceptb[j] = acceptb[j]+1;
} else {
for (z=0; z < *nb; z++){
  beta1[z] = beta1[z];
}
for (i=0; i < *ny; i++){
  mu0[i] = mu0[i];
}
}
}

for (i=0; i < *ny; i++){
       xbm[i]=0;
}
for (i=0; i < *ny; i++){
for (z=0; z < *nb; z++){
       xbm[i]+=x[*ny + (i + (z * *ny))]*beta1[z];
}
}

for (z=0; z < *nb; z++){
beta_r[z] = beta1[z];
beta_r[*nb + z] = acceptb[z];
}
for (i=0; i < *ny; i++){
beta_r[2 * *nb + i] = ex[i]*exp(xbm[i]+gammavec[i]);
}
}
