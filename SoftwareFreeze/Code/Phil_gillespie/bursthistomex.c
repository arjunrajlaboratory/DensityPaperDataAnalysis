#include <math.h>


//Here's the stuff for the random number generator
static unsigned long jz,jsr=123456789;

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3*.2328306e-9)
#define IUNI SHR3

static long hz;
static unsigned long iz, kn[128], ke[256];
static float wn[128],fn[128], we[256],fe[256];

#define RNOR (hz=SHR3, iz=hz&127, (fabs(hz)<kn[iz])? hz*wn[iz] : nfix())
#define REXP (jz=SHR3, iz=jz&255, (    jz <ke[iz])? jz*we[iz] : efix())


float nfix(void);
float efix(void);
void zigset(unsigned long jsrseed);



//END random number generator

#ifndef __STDIO_H__
#include <stdio.h>
#endif

#ifndef mex_h
#include "mex.h"
#endif

#define TIMES_OUT     plhs[0]
#define SPECIES_OUT   plhs[1]

#define CURRT_IN          prhs[0]
#define SPECIES_IN        prhs[1]
#define RATES_IN          prhs[2]
#define PROPENSITIES_IN   prhs[3]
#define SEED_IN           prhs[4]
#define M_IN              prhs[5]
#define TOTALT_IN         prhs[6]

//INCLUDE NSPECIES HERE
#define NSPECIES 4

//INCLUDE NUMRXNS HERE
#define NUMRXNS 7


//void gillespie(long m, double *times_out,double *species_out,double currT,double *species,double *rates,double *propensities,double *rand1,double *rand2)

void gillespie(long m, double *times_out,double *species_out,double currT,double *species,double *rates,double *propensities,long seed, double totalt)
{
  double alpha, deltaT, p;
  double cumpropensities[NUMRXNS];

  double deltaTsave;
  long savecount;
  long totaliterations;

//INSERT ALL VARIABLE DECLARATIONS HERE
  long P,Pactive,Pinactive,mRNA;
  double gam ,lam,k2 ,mu ,del ,mu_p ,del_p ;

  long i, j, k;
//UNPACK ALL SPECIES HERE
  P = (long)species[0];
  Pactive = (long)species[1];
  Pinactive = (long)species[2];
  mRNA = (long)species[3];
//UNPACK ALL RATES HERE
  gam  = rates[0];
  lam = rates[1];
  k2  = rates[2];
  mu  = rates[3];
  del  = rates[4];
  mu_p  = rates[5];
  del_p  = rates[6];

  zigset(seed);
  totaliterations = 0;

  savecount = 0;
  deltaTsave = totalt/(m-1);


  //  for (i=0; i<m; i++) {
  while (currT<totalt) {
    totaliterations++;
    cumpropensities[0] = propensities[0];
    for (j=1; j<NUMRXNS; j++)
      cumpropensities[j] = cumpropensities[j-1]+propensities[j];
    
    alpha = cumpropensities[NUMRXNS-1];
    //deltaT = -1.0/alpha*log(rand1[i]);
    deltaT = -1.0/alpha*log(UNI);
    currT += deltaT;
    //p = rand2[i]*alpha;
    p = UNI*alpha;
    
    
        
    while (currT > savecount*deltaTsave && savecount < m) {
//INSERT SAVE HERE
species_out[savecount*NSPECIES+0] = P;
species_out[savecount*NSPECIES+1] = Pactive;
species_out[savecount*NSPECIES+2] = Pinactive;
species_out[savecount*NSPECIES+3] = mRNA;
		times_out[savecount] = savecount*deltaTsave;
		savecount++;
    }
        
//INSERT IF STATEMENT HERE
if (p<cumpropensities[0]) {
  // rxn: Pactive=Pinactive
  Pactive=Pactive + -1;
  Pinactive=Pinactive + 1;

  //update propensity for Pactive=Pinactive
  propensities[0] = gam *Pactive;
  //update propensity for Pinactive=Pactive
  propensities[1] = lam*Pinactive;
  //update propensity for Pinactive=Pinactive+mRNA
  propensities[2] = k2 *Pinactive;
  //update propensity for Pactive=Pactive+mRNA
  propensities[3] = mu *Pactive;
} else if (p<cumpropensities[1]) {
  // rxn: Pinactive=Pactive
  Pactive=Pactive + 1;
  Pinactive=Pinactive + -1;

  //update propensity for Pactive=Pinactive
  propensities[0] = gam *Pactive;
  //update propensity for Pinactive=Pactive
  propensities[1] = lam*Pinactive;
  //update propensity for Pinactive=Pinactive+mRNA
  propensities[2] = k2 *Pinactive;
  //update propensity for Pactive=Pactive+mRNA
  propensities[3] = mu *Pactive;
} else if (p<cumpropensities[2]) {
  // rxn: Pinactive=Pinactive+mRNA
  mRNA=mRNA + 1;

  //update propensity for mRNA=
  propensities[4] = del *mRNA;
  //update propensity for mRNA=mRNA+P
  propensities[5] = mu_p *mRNA;
} else if (p<cumpropensities[3]) {
  // rxn: Pactive=Pactive+mRNA
  mRNA=mRNA + 1;

  //update propensity for mRNA=
  propensities[4] = del *mRNA;
  //update propensity for mRNA=mRNA+P
  propensities[5] = mu_p *mRNA;
} else if (p<cumpropensities[4]) {
  // rxn: mRNA=
  mRNA=mRNA + -1;

  //update propensity for mRNA=
  propensities[4] = del *mRNA;
  //update propensity for mRNA=mRNA+P
  propensities[5] = mu_p *mRNA;
} else if (p<cumpropensities[5]) {
  // rxn: mRNA=mRNA+P
  P=P + 1;

  //update propensity for P=
  propensities[6] = del_p *P;
} else if (p<cumpropensities[6]) {
  // rxn: P=
  P=P + -1;

  //update propensity for P=
  propensities[6] = del_p *P;
}


  }
  printf("Total iterations = %d\n",totaliterations);


}
/*-------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{ 
  long i;

  double *rand1, *rand2, *propensities_in, *species, *rates;
  double propensities[NUMRXNS];
  double *currTpr;
  double *seedpr, *mpr, *totaltpr;
  double seed, m, totalt;

  double *times_out, *species_out;

  double currT;

  //Load scalars first...
  currTpr = mxGetPr(CURRT_IN);
  currT = currTpr[0];
  
  seedpr = mxGetPr(SEED_IN);
  seed = seedpr[0];
  mpr = mxGetPr(M_IN);
  m = mpr[0];

  totaltpr = mxGetPr(TOTALT_IN);
  totalt = totaltpr[0];

  //Now load matrices/vectors...
  species = mxGetPr(SPECIES_IN);
  propensities_in = mxGetPr(PROPENSITIES_IN);
  rates = mxGetPr(RATES_IN);

  //Now let's create the output matrices
  TIMES_OUT = mxCreateDoubleMatrix(1,(long)m,mxREAL);
  times_out = mxGetPr(TIMES_OUT);
  SPECIES_OUT = mxCreateDoubleMatrix(NSPECIES,(long)m,mxREAL);
  species_out = mxGetPr(SPECIES_OUT);

  // Copy propensities_in into propensities
  for (i = 0; i<NUMRXNS; i++)
    propensities[i] = propensities_in[i];

  
  gillespie((long)m,times_out,species_out,currT,species,rates,propensities,(long)seed, totalt);

} /* end function mexFunction */

/*-------------------------------------------------------------------*/
// Random number stuff below


/* nfix() generates variates from the residue when rejection in RNOR occurs. */

float nfix(void)
{
const float r = 3.442620f;     /* The start of the right tail */
static float x, y;
 for(;;)
  {  x=hz*wn[iz];      /* iz==0, handles the base strip */
     if(iz==0)
       { do{ x=-log(UNI)*0.2904764; y=-log(UNI);}	/* .2904764 is 1/r */
        while(y+y<x*x);
        return (hz>0)? r+x : -r-x;
       }
                         /* iz>0, handle the wedges of other strips */
      if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;

     /* initiate, try to exit for(;;) for loop*/
      hz=SHR3;
      iz=hz&127;
      if(fabs(hz)<kn[iz]) return (hz*wn[iz]);
  }

}

/* efix() generates variates from the residue when rejection in REXP occurs. */
float efix(void)
{ float x;
 for(;;)
  {  if(iz==0) return (7.69711-log(UNI));          /* iz==0 */
     x=jz*we[iz]; if( fe[iz]+UNI*(fe[iz-1]-fe[iz]) < exp(-x) ) return (x);

      /* initiate, try to exit for(;;) loop */
   jz=SHR3;
   iz=(jz&255);
   if(jz<ke[iz]) return (jz*we[iz]);
  }
}
/*--------This procedure sets the seed and creates the tables------*/

void zigset(unsigned long jsrseed)
{  const double m1 = 2147483648.0, m2 = 4294967296.;
   double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
   double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
   int i;
   jsr^=jsrseed;

/* Set up tables for RNOR */
   q=vn/exp(-.5*dn*dn);
   kn[0]=(dn/q)*m1;
   kn[1]=0;

   wn[0]=q/m1;
   wn[127]=dn/m1;

   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);

    for(i=126;i>=1;i--)
    {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
     kn[i+1]=(dn/tn)*m1;
     tn=dn;
     fn[i]=exp(-.5*dn*dn);
     wn[i]=dn/m1;
    }

/* Set up tables for REXP */
    q = ve/exp(-de);
    ke[0]=(de/q)*m2;
    ke[1]=0;

    we[0]=q/m2;
    we[255]=de/m2;

    fe[0]=1.;
    fe[255]=exp(-de);

   for(i=254;i>=1;i--)
  {de=-log(ve/de+exp(-de));
   ke[i+1]= (de/te)*m2;
   te=de;
   fe[i]=exp(-de);
   we[i]=de/m2;
  }
}

