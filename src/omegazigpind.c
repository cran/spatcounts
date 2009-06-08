#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

void omegazigpind(double* omega, double* mu, double* phi, double* x, double* y, double* ex, int* acceptomega, int* ny, int* n1, double* omega_r)
{

int i;
double omegalo;
double omegaup;
double omegam;
double dllo;
double dlup;
double dlm;
double propcov;
double u;
double omeganew;
double temp;
double lomega;

/*Posterior Mode and inverse curvature at the mode for omega-update, using bisection*/
omegalo=0;
omegaup=0.9999;
omegam=(omegalo+omegaup)/2;

dllo=0;
dlup=0;
dlm=0;

for (i=0; i < *ny; i++){
   if (y[i]==0){
    dllo += (1-exp(-mu[i] / *phi))/(omegalo+(1-omegalo)*exp(-mu[i] / *phi));
    dlup += (1-exp(-mu[i] / *phi))/(omegaup+(1-omegaup)*exp(-mu[i] / *phi));
    dlm += (1-exp(-mu[i] / *phi))/(omegam+(1-omegam)*exp(-mu[i] / *phi));
  }
}

     dllo = dllo-*n1/(1-omegalo);
     dlup = dlup-*n1/(1-omegaup);
     dlm = dlm-*n1/(1-omegam);


if ((dllo<0) || (dlup>0)){
  /*printf(" dllo<0 or dlup>0\n");*/
omegam=0;
} else {
 while (fabs(omegaup - omegalo) > 0.001*fabs(omegaup)){
   if (dlm < 0){
      dlup = dlm;
      omegaup = omegam;
      omegam = (omegalo+omegaup)/2;
   } else {
      dllo = dlm;
      omegalo = omegam;
      omegam = (omegalo+omegaup)/2;
  }
      dlm = 0;
 for (i=0; i < *ny; i++){
   if (y[i]==0){
    dlm += (1-exp(-mu[i] / *phi))/(omegam+(1-omegam)*exp(-mu[i] / *phi));
  }
}
     dlm = dlm-*n1/(1-omegam);
}
}

propcov=0;
 for (i=0; i < *ny; i++){
   if (y[i]==0){
        propcov += pow(1-exp(-mu[i] / *phi),2)/pow(omegam+(1-omegam)*exp(-mu[i] / *phi),2);
    }
}

propcov = propcov+*n1/pow(1-omegam,2);
propcov = 1/propcov;

GetRNGstate();
u=runif(0,1);
PutRNGstate();

 /*t-Proposal for omega with v degrees of freedom*/
GetRNGstate();
  omeganew = (omegam + pow(21*propcov/20,0.5)*rt(20));
PutRNGstate();

GetRNGstate();
  while (omeganew >= 1 | omeganew < 0){
     omeganew = (omegam + pow(21*propcov/20,0.5)*rt(20));
  }
PutRNGstate();

/*Calculation of acceptance probability*/
 temp=0;
 for (i=0; i < *ny; i++){
   if (y[i]==0){
     temp += log((omeganew+(1-omeganew)*exp(-mu[i] / *phi))/(*omega+(1-*omega)*exp(-mu[i] / *phi)));
     } else {
     temp += log((1-omeganew)/(1-*omega));
     }
 }

/*Proposal in acceptance probability*/
  temp = temp + 10.5*log(1+pow(omeganew-omegam,2)*(1/propcov)/20)-log(1+pow(*omega-omegam,2)*(1/propcov)/20);

 lomega = (temp<0)*temp;
 if ((log(u) < lomega) | (lomega >= 0)){
 *omega = omeganew;
 *acceptomega = *acceptomega+1;
 } else {
 *omega = *omega;
 }
 omega_r[0] = *omega;
 omega_r[1] = *acceptomega;
}
