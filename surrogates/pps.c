/*=================================================================
 *
 * PPS.C	.MEX file corresponding to PPS.M
 *              Computes a pseudo periodic surrogate (noisy local
 *              linear model) from data y. Use embedding dimension 
 *              de and embedding lag tau. Generate an n point 
 *              surrogate with noise radius rad.
 *
 * The calling syntax is:
 *
 *		yp = pps(y,de,tau,n,rad,y1), or
 *              yp = pps(y,emb,[],n,rad,y1)
 *
 *
 * This is a MEX-file for MATLAB.  
 *=================================================================
 *
 * Original code, copyright 2001 
 * Michael Small (ensmall@polyu.edu.hk)
 * Hong Kong Polytechnic University
 * Hong Kong SAR, PR China 
 * For comments/problems/redistribution, email the author.
 *
 *=================================================================*/
/* $Revision: 1.5 $ */
#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <math.h>

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

float ran2(long *idum);
int surrogate(double *ys, double *yi, double *y, int lengthy, int de, int* emb, int n, double rho, int xinit);

float ran2(long *idum)
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

int surrogate(double *ys,
	      double *yi,
	      double *y,
	      int lengthy, 
	      int de, 
	      int* emb,
	      int n, 
	      double rho,
	      int xinit)
     /* compute the pls surrogate */
{
  const int x0=*(emb+de-1);
  const int ny=lengthy-x0+1;
  const int ny_1=ny-2;
  long idum;
  
  int i,j,k,xi;
  double prob,sum,temp,ri;
  double *pdf;

  idum=-abs(rand()); /*random number seed*/

  pdf = (double *)malloc(sizeof(double)*ny);

  if (xinit==0)
    xi= (int)(rand()*ny)/(RAND_MAX+1.0)+x0; /*random integer < ny*/
else
	xi=xinit;

  *ys=*(y+xi);
  *yi=xi;

  for (i=1; i<n; i++){
    /* repeat n times, to get n point surrogate */

    /* find the current image */
    sum=0;
    for (j=0; j<ny_1; j++){         /*compute the distribution of interpoint distances*/
      temp = 0;
      for (k=0; k<de; k++) {      /* L2-norm^2 between 2 points */ 
	prob = ( *(y+xi-*(emb+k)) - *(y+x0+j-*(emb+k)) );
	temp += prob*prob;
      }
      prob = exp( -0.5*sqrt(temp)/rho );/*pow(exp(0.5/rho),-temp);  */
      if (1 /*!isnan(prob)*/) {
	*(pdf+j) = prob;
	sum += prob;                /* cumsum */
      } else {
	*(pdf+j) = 0;
      }

    }
    ri = rand()/(RAND_MAX+1.0);  /*random number in ]0,1[ */
    ri=ri*sum;                   /* random number in ]0,sum[ */
    xi=0;
    temp=0;
    while (temp<ri) {
      temp += *(pdf+xi); 
      xi++;
    }
    xi += x0;

    /* and add it to the list */
    *(ys+i) = *(y+xi);
    *(yi+i) = xi;

  }
  free((double *)pdf);
  return 0;
}



  
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     /* the MATLAB mex wrapper function */     
{ 
  const int dtau=1,dde=5;
  const double drho=0.01;
  double *y,*yi,*ys,*yout,*temp;
  int i,nrows,ncols,lengthy;
  int *emb;
  double tau,de,n,xinit,rho;
    
    /* Check for proper number of arguments */

    
    if (nrhs > 6) { 
	mexErrMsgTxt("Too many input arguments."); 
    } 
    if (nrhs < 1) {
        mexErrMsgTxt("Insufficient input arguments.");
    } 
    if (nlhs > 2) {
	mexErrMsgTxt("Too many output arguments."); 
    } 

    
    /*Assign a pointer to the input vector y*/
    y = mxGetPr(prhs[0]);
    /* Check the dimensions of y */     
    nrows = mxGetM(prhs[0]); 
    ncols = mxGetN(prhs[0]);
    /* check that y is a vector */
    if ((nrows!=1) && (ncols!=1)) {
      mexWarnMsgTxt("First input should be a vector.");
    }
    lengthy = nrows*ncols;

    /*Get the parameters tau, de, n and rho*/
    if (nrhs>1) {
      /* get tau */
      de=mxGetScalar(prhs[1]);
      nrows = mxGetM(prhs[1]);
      ncols = mxGetN(prhs[1]);
      /* check de is a scalar */
      if ((nrows!=1) || (ncols!=1)) {
	/*	mexWarnMsgTxt("Second input should be scalar.");*/
      }
      /* check de is an integer */
      if (floor(de)!=de) {
	tau=(int)floor(de);
	mexWarnMsgTxt("Second input rounded to integer.");
      }
      if (nrhs>2) {
	/* get tau */
	tau = mxGetScalar(prhs[2]);
	nrows = mxGetM(prhs[2]);
	ncols = mxGetN(prhs[2]);
	/* check de is a scalar */
	if ((nrows>1) || (ncols>1)) {
	  mexWarnMsgTxt("Third input should be scalar.");
	}
	/* check de is an integer */
	if (floor(tau)!=tau) {
	  tau=(int)floor(tau);
	  mexWarnMsgTxt("Third input rounded to integer.");
	}
	if (nrhs>3) {
	  /* get n */
	  n=mxGetScalar(prhs[3]);
	  nrows = mxGetM(prhs[3]);
	  ncols = mxGetN(prhs[3]);
	  /* check n is a scalar */
	  if ((nrows!=1) || (ncols!=1)) {
	    mexWarnMsgTxt("Fourth input should be scalar.");
	  }
	  /* check n in an integer */
	  if (floor(n)!=n) {
	  n = (int)floor(n);
	  mexWarnMsgTxt("Fourth input rounded to integer.");
	  }
	  if (nrhs>4) {
	    /* get rho */
	    rho = mxGetScalar(prhs[4]);
	    nrows = mxGetM(prhs[4]);
	    ncols = mxGetN(prhs[4]);
	    /* check rho is a scalar */
	    if ((nrows!=1) || (ncols!=1)) {
	      mexWarnMsgTxt("Fifth input should be scalar.");
	    }
            if (nrhs>5) {
              /* get the initial state, if any */
	      xinit=mxGetScalar(prhs[5]);
	      nrows = mxGetM(prhs[4]);
	      ncols = mxGetN(prhs[4]);
	      /* check rho is a scalar */
	      if ((nrows!=1) || (ncols!=1)) {
	        mexWarnMsgTxt("Sixth input should be scalar.");
	      }	  
              /* check xinit in an integer */
	      if (floor(xinit)!=xinit) {
	        xinit = (int)floor(xinit);
	        mexWarnMsgTxt("Sixth input rounded to integer.");
              }
            } else {
              xinit=0;
            }

	  } else {
	    /* default rho */
	    rho=drho;
	  }
	} else {
	  /* default rho */
	  rho=drho;
	  /*default n */
	  n=lengthy;
	}
      } else {
	/* default rho */
	rho=drho;
	/*default n */
	n=lengthy;
	/*default de */
	de=dde;
      }
    } else {
      /* default rho */
      rho=drho;
      /*default n */
      n=lengthy;
      /*default de */
      de=dde;
      /* default tau */
      tau=dtau;
    }
     
    /* unless second arg is a vector, use [0 tau 2tau ... (de-1)tau] otherwise, use whats given*/
    /* check de is a scalar */
    nrows = mxGetM(prhs[1]);
    ncols = mxGetN(prhs[1]);
    if ((nrows!=1) || (ncols!=1)) {
      mexWarnMsgTxt("Second arg is a vector, ignoring third arg.");
      /*Assign a pointer to the input vector y*/
      temp = mxGetPr(prhs[1]);
      de = nrows*ncols;
      emb = (int *)malloc(sizeof(int)*de);
      for (i=0; i<de; i++) {
	*(emb+i)=(int)*(temp+i);
      }
    } else {
      emb = (int *)malloc(sizeof(int)*de);
      for (i=0; i<de; i++) {
	*(emb+i)=(int)(i*tau);
      }
    }
    /*now, de is the embedding dimension (length of emb), and emb if the vector of embedding lags */
    

    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1, n, mxREAL);


    /* Assign pointer to the ouput matrix */ 
    ys = mxGetPr(plhs[0]);
    yi = mxGetPr(plhs[1]);

    /* Do the actual computations in a subroutine */
    surrogate(ys,yi,y,lengthy,(int)de,emb,(int)n,rho,(int)xinit);
    
    
    /* mexAddFlops(1); */ 
    /* I'm sure I can figure that out properly, later.*/
    free((int *)emb);

    return;
    
}


#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* Random number code (C) Copr. 1986-92 Numerical Recipes Software D*V. */
