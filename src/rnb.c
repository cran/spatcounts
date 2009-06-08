#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

void rnb(double* mu, double* r, double* x, int* ny, double* y, double* ex, int* acceptr, double* rvar, double* a, double* b, double* r_r)
{

int i;

double u;
double rnew;
double temp;
double lr;

/*a,b: Parameter of Gamma(a,b)-Prior of r*/

GetRNGstate();
u=runif(0,1);
PutRNGstate();

/*Proposal forr r: truncated normal*/
  rnew = rnorm(*r,*rvar);
    while (rnew < 0){  /*| rnew > 100){*/
      rnew = rnorm(*r,*rvar);
    }

/*Calculation of acceptance probability*/
 temp=0;
 for (i=0; i < *ny; i++){
    temp+=((lgammafn(y[i]+rnew)+lgammafn(*r))-(lgammafn(y[i]+ *r)+lgammafn(rnew))+rnew*log(rnew/(mu[i]+rnew))-(*r)*log(*r/(mu[i]+*r))+y[i]*log((mu[i]+*r)/(mu[i]+rnew)));
}

/*Prior for r*/
temp = temp + (*a-1)*log((rnew)/(*r)) - *b * ((rnew)-(*r));

/*Proposal Ratio for gamma proposal for r*/
/*temp=temp+((*r)-(rnew))/ *rvar*(log(*rvar)-1)+log(gammafn((*r)/ *rvar)/gammafn((rnew)/ *rvar))-((*r)/ *rvar-1)*log((rnew))+((rnew)/ *rvar-1)*log((*r));*/

lr = (temp<0)*temp;

 if ((log(u) < lr) | (lr >= 0)){
 *r = rnew;
 *acceptr = *acceptr+1;
 } else {
 *r = *r;
 }
 r_r[0] = *r;
 r_r[1] = *acceptr;
}
