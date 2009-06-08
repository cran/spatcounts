#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

void phigpind(double* phi, double* mu, double* y, int* acceptphi, int* ny, double* phi_r)
{

int i;
double philo;
double phiup;
double phim;
double dllo;
double dlup;
double dlm;
double propcov;
double u;
double phinew;
double temp;
double lphi;

/*Posterior Mode and inverse curvature at the mode for phi-update, using bisection*/
philo=1;
phiup=10000;
phim=(philo+phiup)/2;

dllo=0;
dlup=0;
dlm=0;

for (i=0; i < *ny; i++){
  dllo+=((y[i]-1)*y[i]/(mu[i]+(philo-1)*y[i])-(y[i]/philo)-(1/pow(philo,2))*(y[i]-mu[i]));
    dlup+=((y[i]-1)*y[i]/(mu[i]+(phiup-1)*y[i])-(y[i]/phiup)-(1/pow(phiup,2))*(y[i]-mu[i]));
      dlm+=((y[i]-1)*y[i]/(mu[i]+(phim-1)*y[i])-(y[i]/phim)-(1/pow(phim,2))*(y[i]-mu[i]));
}
  dllo=dllo-(2/philo);
  dlup=dllo-(2/phiup);
  dlm=dlm-(2/phim);
if (dllo<0 & dlup<0){
phim=1;
} else {
 while (fabs(phiup-philo) > 0.001*fabs(phiup)){
   if (dlm < 0){
      dlup = dlm;
      phiup = phim;
      phim = (philo + phiup)/2;
   } else {
      dllo = dlm;
      philo = phim;
      phim = (philo + phiup)/2;
   }
   dlm=0;
      for (i=0; i < *ny; i++){
	dlm+=((y[i]-1)*y[i]/(mu[i]+(phim-1)*y[i])-(y[i]/phim)-(1/pow(phim,2))*(y[i]-mu[i]));
      }
	dlm=dlm-(2/phim);
  }
}

propcov=0;
for (i=0; i < *ny; i++){
  propcov+=(pow(y[i],2)*(y[i]-1)/pow(mu[i]+(phim-1)*y[i],2)-y[i]/pow(phim,2)-2/pow(phim,3)*(y[i]-mu[i]));
}
  propcov = 1/(propcov-2/pow(phim,2));

GetRNGstate();
u=runif(0,1);
PutRNGstate();

 /*t-Proposal for phi with v degrees of freedom*/
GetRNGstate();
  phinew = phim + pow(21*propcov/20,0.5)*rt(20);
PutRNGstate();

GetRNGstate();
  while (phinew < 1){
     phinew = phim + pow(21*propcov/20,0.5)*rt(20);
  }
PutRNGstate();

/*Calculation of acceptance probability*/
 temp=0;
 for (i=0; i < *ny; i++){
    temp+=((y[i]-1)*log(mu[i]+(phinew-1)*y[i]/(mu[i]+(*phi-1)*y[i]))-y[i]*log(phinew / *phi)+(y[i]-mu[i])*(1/phinew-1 / *phi));
 }

   /*Prior for phi*/
   temp=temp-2*log(phinew / *phi);

/*Proposal in acceptance probability*/
  temp = temp + 10.5*log(1+pow(phinew-phim,2)*(1/propcov)/20)-log(1+pow(*phi-phim,2)*(1/propcov)/20);

 lphi = (temp<0)*temp;
 if ((log(u) < lphi) | (lphi >= 0)){
 *phi = phinew;
 *acceptphi = *acceptphi+1;
 } else {
 *phi = *phi;
 }
 phi_r[0] = *phi;
 phi_r[1] = *acceptphi;
}
