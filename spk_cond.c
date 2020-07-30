#include "mex.h"
#include "fastexponent.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

/* Constants for ran2 random number generator. From Numerical Recipes. */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* declare functions used other than mexFunction*/
double ran2(long *idum);    /* Ran2 uniform random deviate generator from Numerical Recipes */
double gasdev(long *idum);
int myComparamsisonFunction(const void *x, const void *y);

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{

// declare variables
    double *params,*seed,*spk_pad,*xi_c;
    double taueff,sigmaeff,Veff,tref,Vth,Vr,Delta,VT;
    double C,gx,tau_x,Vx,Vxh,Dx;
    double c,xi,cc,ci,tstop,trans,dt,t,sqrtdt;
   
    long idum;
    int countM,R,r,Ni,Ne,N,n,ne,ni,i,Tloop,count;
    
    // Get pointers to input arrays
    params=mxGetPr(prhs[0]);
    seed=mxGetPr(prhs[1]);  
            
    idum=seed[0]<0?(long)seed[0]:-(long)seed[0]; /* Ran2 needs a negative value for a seed. If the seed is not negative, make it negative. */
    
    /* Get paramsameters */
    taueff = params[0];
    Veff = params[1];
    Vth = params[2];
    Vr = params[3];
    VT = params[4];
    Delta = params[5];
    C = params[6];
    tref = params[7];
    gx = params[8];
    Vx = params[9];
    Vxh = params[10];
    tau_x = params[11];
    Dx = params[12];
    sigmaeff = params[13];
    c = params[14];
    tstop = params[15];
    trans = params[16];
    dt = params[17]; 
    Ne = params[18];

    Tloop=(int)(tstop/dt);
    N=Ne;
    countM=100*Ne*2*tstop/1000;
    

    /* Create an mxArray for the output data */
    plhs[0]=mxCreateDoubleMatrix(2,countM,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(1,Tloop,mxREAL);
    
    /* Create a pointer to the output data */
    spk_pad=mxGetPr(plhs[0]);
    xi_c=mxGetPr(plhs[1]);
    
    double V[N];
    double V0[N];
    double tlast[N];
    double x[N];
    double x0[N];
  
    for(n=0;n<N;n++){ //ICs
        V0[n]=Vr;
        V[n]=Vr;
        tlast[n]=-tstop;
        x0[n]=0.;   
        x[n]=0.;     
    }

    t=0;
    count=0;
    cc=sqrt(c);
    ci=sqrt(1-c);
    sqrtdt = sqrt(dt);

    // pi=3.14159265;

    // mexPrintf("%f\n",tref);
    
    for(i=0;i<Tloop;i++) {
        
        t=t+dt;
        xi_c[i]=cc*gasdev(&idum);
        
        for(ne=0;ne<N;ne++){ //fusiform cells
            
            xi=ci*gasdev(&idum);

            x[ne] = x0[ne] + (dt/tau_x)*(1/(1+EXP(-(V0[ne]-Vxh)/Dx))-x0[ne]);
            V[ne] = V0[ne] + (dt/taueff)*(Veff-V0[ne]+Delta*EXP((V0[ne]-VT)/Delta)) + dt/C*gx*x0[ne]*(Vx-V0[ne]) + sqrtdt/taueff*sqrt(2*taueff)*sigmaeff*(xi_c[i]+xi); 

            if(t-tlast[ne]<=tref){
                V[ne]=Vth-(t-tlast[ne])/tref*(Vth-Vr); //box spike
            }
            
            if(V[ne]>Vth){
                V[ne]=Vth;
                tlast[ne]=t;
                // mexPrintf("%f\n",t);

                if(count<=countM && t>=trans){
                    spk_pad[count]=t;
                    spk_pad[count+1]=ne+1;
                    count+=2;
                }
            }
             
            V0[ne]=V[ne];
            x0[ne]=x[ne];

        } //end neuron loop
        
    } //end time loop
} //end mex function
      
int myComparamsisonFunction(const void *x, const void *y) {

    // x and y are pointers to doubles.

    // Returns -1 if x < y
    //          0 if x == y
    //         +1 if x > y

    double dx, dy;

    dx = *(double *)x;
    dy = *(double *)y;

    if (dx < dy) {
        return -1;
    } 
    else if (dx > dy) {
        return +1;
    }
    return 0;
}


double ran2(long *idum)
{
 int j;
 long k;
 static long idum2=123456789;
 static long iy=0;
 static long iv[NTAB];
 float temp;

 if (*idum <= 0) {
  if (-(*idum) < 1) *idum=1;
  else *idum = -(*idum);
  idum2=(*idum);
  for (j=NTAB+7;j>=0;j--) {
   k=(*idum)/IQ1;
   *idum=IA1*(*idum-k*IQ1)-k*IR1;
   if (*idum < 0) *idum += IM1;
   if (j < NTAB) iv[j] = *idum;
  }
  iy=iv[0];
 }
 k=(*idum)/IQ1;
 *idum=IA1*(*idum-k*IQ1)-k*IR1;
 if (*idum < 0) *idum += IM1;
 k=idum2/IQ2;
 idum2=IA2*(idum2-k*IQ2)-k*IR2;
 if (idum2 < 0) idum2 += IM2;
 j=iy/NDIV;
 iy=iv[j]-idum2;
 iv[j] = *idum;
 if (iy < 1) iy += IMM1;
 if ((temp=AM*iy) > RNMX) return RNMX;
 else return temp;
}

double gasdev(long *idum)
{
    double ran2(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;

    if (*idum < 0) iset=0;
    if  (iset == 0) {
        do {
            v1=2.0*ran2(idum)-1.0;
            v2=2.0*ran2(idum)-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    } else {
        iset=0;
        return gset;
    }
}