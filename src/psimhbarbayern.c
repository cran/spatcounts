#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

void psimhbarbayern(double* psi, int* J, double* ga, double* eigenv, double* invsigmasq, double* gKerng, int* acceptpsi, double* psivar, double* alpha, double* psi_r)
{
/*Prior of psi:1/(1+psi)^alpha*/
/*Proposal for zeta <- psi/(1+psi)*/
/*eigenv <- eig(-KERN)*/

int i;

double u;
double psibarnew;
double psinew;
double temp;
double temp1;
double temp2;
double lpsi;

GetRNGstate();
u=runif(0,1);
PutRNGstate();

/*generate proposal for psibar*/
GetRNGstate();
psibarnew = rnorm(*psi/(1+*psi),*psivar);
PutRNGstate();

while (psibarnew >= 1 | psibarnew < 0){
GetRNGstate();
psibarnew = rnorm(*psi/(1+*psi),*psivar);
PutRNGstate();
}
psinew = psibarnew/(1-psibarnew);

/*Calculation acceptance probability*/
for (i=0; i < *J; i++){
temp1 += (log(1-psinew*eigenv[i])-log(1-*psi*eigenv[i]));
}

temp2 = *invsigmasq * *gKerng*(psinew-*psi);

/*Prior of psi*/
/*temp = (temp1-temp2)/2+alpha*log((1+psi)/(1+psinew))*/
temp = (temp1 - temp2)/2;

lpsi = (temp < 0)*temp;

if ((log(u) < lpsi) | (lpsi >= 0)){
  *psi = psinew;
  *acceptpsi = *acceptpsi+1;
} else {
*psi = *psi;
}
psi_r[0] = *psi;
psi_r[1] = *acceptpsi;
}

