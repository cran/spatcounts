#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

void betaintercept(double* mu, int* nb, double* beta1, double* phi, int* ny, double* gammavec, double* omega, double* xb, double* y, double* ex, double* P0, int* acceptbint, double* beta_r)
{
  /* printf("%f\n",xb[ny+nb* *nb]);*/
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

 /*Starting values for bisection*/
   blo[0]=beta1[0]-50;
   bup[0]=beta1[0]+50;
   bm[0]=(blo[0]+bup[0])/2;

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
   if (y[i]==0){ 
       dllo+=(-(mulo[i]*(1-*omega)*1 / *phi*exp(-mulo[i] / *phi))/(*omega+(1-*omega)*exp(-mulo[i] / *phi)));
       dlup+=(-(muup[i]*(1-*omega)*1 / *phi*exp(-muup[i] / *phi))/(*omega+(1-*omega)*exp(-muup[i] / *phi)));
       dlm+=(-(mum[i]*(1-*omega)*1 / *phi*exp(-mum[i] / *phi))/(*omega+(1-*omega)*exp(-mum[i] / *phi)));
   } else {
       dllo+=(1+(y[i]-1)*mulo[i]/(mulo[i]+(*phi-1)*y[i])-mulo[i] / *phi);
       dlup+=(1+(y[i]-1)*muup[i]/(muup[i]+(*phi-1)*y[i])-muup[i] / *phi);
       dlm+=(1+(y[i]-1)*mum[i]/(mum[i]+(*phi-1)*y[i])-mum[i] / *phi);
}
}
dllo=dllo-blo[0] * P0[1+1* *nb];
dlup=dlup-bup[0] * P0[1+1* *nb];
dlm=dlm-bm[0] * P0[1+1* *nb];

/*if ((dllo<0) || (dlup>0)){
      printf(" dllo<0 or dlup>0\n");
      }*/

   while(fabs(bup[0] - blo[0]) > 0.01*fabs(bup[0])) {

   if (dlm<0){
      dlup=dlm;
      bup[0]=bm[0];
      bm[0]=(blo[0]+bup[0])/2;
   } else {
      dllo=dlm;
      blo[0]=bm[0];
      bm[0]=(blo[0]+bup[0])/2;
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
   if (y[i]==0){ 
       dlm+=(-(mum[i]*(1-*omega)*1 / *phi*exp(-mum[i] / *phi))/(*omega+(1-*omega)*exp(-mum[i] / *phi)));
   } else {
       dlm+=(1+(y[i]-1)*mum[i]/(mum[i]+(*phi-1)*y[i])-mum[i] / *phi);
}
}
       dlm=dlm-bm[0] * P0[1+1* *nb];
}
propcov=0;

for (i=0; i < *ny; i++){
   if (y[i]==0){
        propcov += (-mum[i]*(1-*omega) / *phi*exp(-mum[i] / *phi)*(*omega*(mum[i] / *phi-1)+exp(-mum[i] / *phi)*(*omega-1))/pow(*omega+(1-*omega)*exp(-mum[i] / *phi),2));
    } else {
        propcov += (-(y[i]-1)*(*phi-1)*y[i]*mum[i]/pow(mum[i]+(*phi-1)*y[i],2)+mum[i] / *phi);
    }
}

propcov=1/(propcov+P0[1+1* *nb]);

GetRNGstate();
u=runif(0,1);
PutRNGstate();

for (z=0; z < *nb; z++){
betanew[z] = beta1[z];
}

GetRNGstate();
betanew[0]=(bm[0]+pow(21*propcov/20, 0.5)*rt(20));
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
for (i=0; i < *ny; i++){
   if (y[i]==0){
      temp += (log((*omega+(1-*omega)*exp(-1 / *phi*munew[i]))/(*omega+(1-*omega)*exp(-1 / *phi*mu[i]))));
     } else {
      temp += (betanew[0]-beta1[0]+(y[i]-1)*log((munew[i]+y[i]*(*phi-1))/(mu[i]+y[i]*(*phi-1)))-1 / *phi*(munew[i]-mu[i]));
     }
}

/*N(0,P0^-1) Prior for beta*/
 temp = temp - 0.5*((pow(betanew[0],2) - pow(beta1[0],2))*P0[1+1* *nb]) + 10.5*log(1+pow(betanew[0]-bm[0],2)*(1/propcov)/20)-log(1+pow(beta1[0]-bm[0],2)*(1/propcov)/20);

lbeta1=(temp<0)*temp;

if ((log(u) < lbeta1) | (lbeta1 >= 0)){
for (z=0; z < *nb; z++){
  beta1[z] = betanew[z];
}
for (i=0; i < *ny; i++){
  mu[i] = munew[i];
}
  *acceptbint = *acceptbint+1;
} else {
for (z=0; z < *nb; z++){
  beta1[z] = beta1[z];
}
for (i=0; i < *ny; i++){
  mu[i] = mu[i];
}
}
for (z=0; z < *nb; z++){
beta_r[z] = beta1[z];
}
beta_r[*nb] = *acceptbint;
for (i=0; i < *ny; i++){
beta_r[2 * *nb + i] = mu[i];
}
}
