#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

void gammazigpindbisection(int* nb, double* beta1, int* J, double* ga,  int* ny, double* mu,  double* phi, double* psi, double* invsigmasq, double* x, double* y, int* zfull, double* ex, int* acceptg, int* nmat, int* maxindex, double* gamma_r)
{

int i, j, z, w, n, nj, ny0;

double *y0 = (double*) R_alloc(*ny,sizeof(double));
double *ex0 = (double*) R_alloc(*ny,sizeof(double));
double *mu0 = (double*) R_alloc(*ny,sizeof(double));
double *g0 = (double*) R_alloc(*ny,sizeof(double)); 
double *x0b = (double*) R_alloc(*ny*(*nb),sizeof(double));
double *x0 = (double*) R_alloc(*ny*(*nb+1),sizeof(double));
double *yj = (double*) R_alloc(*ny,sizeof(double));
double *xjb = (double*) R_alloc(*ny*(*nb),sizeof(double));
double *exj = (double*) R_alloc(*ny,sizeof(double));
double *mujold = (double*) R_alloc(*ny,sizeof(double));
double *ganew = (double*) R_alloc(*J,sizeof(double));
double pnj;
double P0g;
double gamsum;
double galo;
double gaup;
double gam1;
double *xjbeta = (double*) R_alloc(*ny,sizeof(double));
double *mulo = (double*) R_alloc(*ny,sizeof(double));
double *muup = (double*) R_alloc(*ny,sizeof(double));
double *mum = (double*) R_alloc(*ny,sizeof(double));
double dllo;
double dlup;
double dlm;
double propcov;
double u;
double *mucan = (double*) R_alloc(*ny,sizeof(double));
double temp;
double lgamma1;

ny0=0;
for(z=0; z < *ny; z++){
if (zfull[z]==0){
y0[ny0]=y[z];
ex0[ny0]=ex[z];
mu0[ny0]=mu[z];
for(i=0; i < (*nb+1); i++){
x0[i+ny0*(*nb+1)]=x[(z + (i * *ny))];
}

ny0 += 1;
}
}

for(j=0; j < *J; j++){
/*only data in region j*/
n=0;
for(z=0; z < ny0; z++){
if (x0[z*(*nb+1)]==(j+1)){
yj[n]=y0[z];
exj[n]=ex0[z];
mujold[n]=mu0[z];
for(i=0; i < *nb; i++){
xjb[i+n*(*nb)]=x0[(i+1)+z*(*nb+1)];
}
n+=1;
}
}

for(z=0; z < *J; z++){
   ganew[z] = ga[z];
}
   nj = nmat[*J * (*maxindex-1) + j];
   pnj = 1 + *psi*nj;
   P0g = *invsigmasq*pnj;

   /*sum of all gammas which are neighbours of j*/
   gamsum=0;

   for (i=0; i < nj; i++){
    gamsum+=ga[(nmat[j + (*J * i)]-1)];
   }

   /*galo, gaup: lower and upper starting values for bisection algorithm*/
   /*galo <- ga[j] - 3*/
   /*gaup <- ga[j] + 3*/
   galo = ga[j] - 15;
   gaup = ga[j] + 15;
   gam1 = (galo + gaup)/2;

for (i=0; i < *ny; i++){
       xjbeta[i]=0;
}

for(i=0; i < *nb; i++){
for(z=0; z < n; z++){
   xjbeta[z]+=xjb[i+(z * *nb)] * beta1[i];
}
}

dllo=0;
dlup=0;
dlm=0;
for(z=0; z < n; z++){
   mulo[z] = exj[z]*exp(xjbeta[z]+galo);
   muup[z] = exj[z]*exp(xjbeta[z]+gaup);
   mum[z] = exj[z]*exp(xjbeta[z]+gam1);

    dllo+=1 + (yj[z]-1)*mulo[z]/(mulo[z]+(*phi-1)*yj[z])-mulo[z] / *phi;
    dlup+= 1 + (yj[z]-1)*muup[z]/(muup[z]+(*phi-1)*yj[z])-muup[z] / *phi;
    dlm+=1 + (yj[z]-1)*mum[z]/(mum[z]+(*phi-1)*yj[z])-mum[z] / *phi;
}

dllo = dllo - P0g*(galo-*psi/pnj*gamsum);
dlup = dlup - P0g*(gaup-*psi/pnj*gamsum);
dlm = dlm - P0g*(gam1-*psi/pnj*gamsum);

/*if ((dllo<0) || (dlup>0)){
      printf(" dllo<0 or dlup>0\n");
      }*/

while(fabs(gaup - galo) > 0.55*fabs(gaup)) {
   if (dlm < 0){
      dlup = dlm;
      gaup = gam1;
      gam1 = (galo + gaup)/2;
   } else {
      dllo = dlm;
      galo = gam1;
      gam1 = (galo + gaup)/2;
   }


   for (i=0; i < *ny; i++){
       xjbeta[i]=0;
   }


   for (i=0; i < *nb; i++){
   for (z=0; z < n; z++){
       xjbeta[z]+=xjb[i+(z * *nb)] * beta1[i];
   }
   }

   dlm=0;
   for (z=0; z < n; z++){
       mum[z] = exj[z]*exp(xjbeta[z]+gam1);

       dlm += 1 + (yj[z]-1)*mum[z]/(mum[z]+(*phi-1)*yj[z])-mum[z] / *phi;
    }

    dlm = dlm - P0g*(gam1-*psi/pnj*gamsum);
    }

propcov=0; 
for (z=0; z < n; z++){
    propcov+=(- mum[z]*(((yj[z]-1)*(*phi-1)*yj[z])/pow(mum[z]+(*phi-1)*yj[z],2)-1 / *phi));
}
propcov=1/(propcov+P0g);

GetRNGstate();
u = runif(0,1);
PutRNGstate();

/*Proposal for gamma:*/
 /*ganew[j] <- mvrnorm(n = 1, propmean,propcov)*/
 /*t-Proposal for gamma*/
GetRNGstate();
 ganew[j] = gam1+pow(21*propcov/20,0.5)*rt(20);
PutRNGstate();

/*Calculation of new mu*/
  for (i=0; i < *ny; i++){
       xjbeta[i]=0;
   }
   for (i=0; i < *nb; i++){
   for (z=0; z < n; z++){
       xjbeta[z]+=xjb[i+(z * *nb)] * beta1[i];
   }
   }

   dlm=0;
   for (z=0; z < n; z++){
       mucan[z] = exj[z]*exp(xjbeta[z]+ganew[j]);
   }

/*Calculation of acceptance probability*/
temp=0;
for (z=0; z< n; z++){
 temp+=(ganew[j]-ga[j]+(yj[z]-1)*log((mucan[z]+yj[z]*(*phi-1))/(mujold[z]+yj[z]*(*phi-1)))-1 / *phi*(mucan[z]-mujold[z]));
}

/*CAR-Prior for gammas*/
 temp = temp - 0.5 * *invsigmasq * ((ganew[j]-ga[j])*((ganew[j]+ga[j])*pnj-2 * *psi*gamsum));

/*Proposal in acceptance probability*/
 /*temp = temp + log(mvnpdf(ga[j],propmean,propcov))-log(mvnpdf(ganew[j],propmean,propcov))*/
 temp = temp + 10.5*(log(1+(pow(ganew[j]-gam1,2))*(1/propcov)/20)-log(1+(pow(ga[j]-gam1,2))*(1/propcov)/20));

lgamma1=(temp<0)*temp;

if ((log(u)<lgamma1) | (lgamma1>=0)){
for (z=0; z < *J; z++){
   ga[z] = ganew[z];
}
   acceptg[j] = acceptg[j]+1;
} else {
for (z=0; z < *J; z++){
   ga[z]= ga[z];
  }
 }

gamma_r[j] = ga[j];
gamma_r[*J + j] = acceptg[j];

}
}
