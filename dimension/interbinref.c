/*=================================================================
 *
 * INTERBINREF.C	
 *              .MEX file corresponding to INTERBINREF.M
 *              returns the number of interpoint distance of the 
 *              embedded time series in each bin
 *
 * The calling syntax is:
 *
 *		np = entropy(x,bins,nref,nt)
 *
 *
 * This is a MEX-file for MATLAB.  
 *=================================================================*/
/* $Revision: 1.5 $ */
#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <time.h>


void ipb(double	*data, /* the embedded data (d-by-n) */
	 int d, /* d-dimensional embedding  */
	 int n, /* n-points */
	 double *bins, /* the bins */
	 int nbins, /* the number of bins */
	 int nref, /* number of reference vectors to choose */
	 int nt, /* the Thieler band (forbidden points) */
	 double *binpop) /* the population in each bin (to be computed) */
     /*actually calculate the interpoint binning*/
{
  int xi,i,j,k;
  double dist,val;
  bool toosmall;

  /*initialise random number generator */
  srand((unsigned)time(NULL));

  for (i=0; i<nref; i++) {
    /* choose xi randomly between 0 and n */
    xi=(int)(rand()/(RAND_MAX+1.0)*n)+1; 
    /* I know that the C/C++ rand functions are crap, but it'll do */
    for (j=0; j<n; j++) {
      
      /* proceed if the points aren't too close  */
      if ( abs(xi-j)>nt ) {

	/* compute the square of the distance for these points  */
	dist=0;
	for (k=0; k<d; k++) {
	  val = *(data+d*xi+k)-*(data+d*j+k);
	  dist += val*val;
	}
	
	/* find the bin it fits in  */
	k=0;
	toosmall=1;
	while (toosmall && k<nbins) {
	  if ( (*(bins+k))>dist) {
	    toosmall=0;
	    (*(binpop+k))++;
	  }
	  k++;
	}

      }/* if abs(...) */

    }/* for j */
  }/* for i */

}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     /* the MATLAB mex wrapper function */     
{ 
    double *data,*bins,*thearg3,*thearg4,*localbins,*binpop;
    int mrows,ncols,i,j;
    int d,n,nbins,nref,nt;
    bool warnThem;
    
    /* Check for proper number of arguments */
    
    if (nrhs > 4) { 
	mexErrMsgTxt("Too many input arguments."); 
    }
    if (nrhs < 2) {
        mexErrMsgTxt("Insufficient input arguments.");
    }
    if (nrhs < 3) {
      /* set nref=0 */
      nref=1000;
    } 
    if (nrhs < 4) {
      /* set nt=0 */
      nt=0;
    } 
    if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /*Assign a pointer to the input matrix*/
    data = mxGetPr(prhs[0]);

    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */     
    mrows = mxGetM(prhs[0]); 
    ncols = mxGetN(prhs[0]);
    d=mrows; /* dimension of embedding */
    n=ncols; /* number of data points */
    
    /* Get the size of the "forbidden zone" --- the second input arg. */
    bins=mxGetPr(prhs[1]);
    mrows = mxGetM(prhs[1]);
    ncols = mxGetN(prhs[1]);
    /* check that x is a vector */
    if ((mrows!=1) && (ncols!=1)) {
      mexWarnMsgTxt("Second input should be a vector");
    }
    nbins = mrows*ncols;
    /* check that all bins ascending */
    warnThem=0;
    for (i=0; i<(nbins-1); i++) 
      if ( *(bins+i+1)<*(bins+i)) 
	warnThem=1;
    if (warnThem) {
      mexErrMsgTxt("Error in interbin: Second input not ascending");
    }

    /*alloc memory and make a copy of bins (squaring as we go)*/
    localbins = (double *) calloc(nbins,sizeof(double));
    for (i=0; i<nbins; i++) 
      *(localbins+i) = (*(bins+i))*(*(bins+i));

    /* */
    if (nrhs>2){
      thearg3=mxGetPr(prhs[2]);
      mrows = mxGetM(prhs[2]);
      ncols = mxGetN(prhs[2]);
      /* check that x is a scalar */
      if ((mrows!=1) || (ncols!=1)) {
	mexWarnMsgTxt("Third input should be a scalar");
      }
      /* check that it is an integer */
      warnThem=0;
      if (floor(*(thearg3))!=*(thearg3)) {
	mexWarnMsgTxt("Third input being rounded to an integer value");
	}
      nref=(int)(*(thearg3));
    }

    /* */
    if (nrhs==4){
      thearg4=mxGetPr(prhs[3]);
      mrows = mxGetM(prhs[3]);
      ncols = mxGetN(prhs[3]);
      /* check that x is a scalar */
      if ((mrows!=1) || (ncols!=1)) {
	mexWarnMsgTxt("Fourth input should be a scalar");
      }
      /* check that it is an integer */
      warnThem=0;
      if (floor(*(thearg4))!=*(thearg4)) {
	mexWarnMsgTxt("Fourth input being rounded to an integer value");
	}
      nt=(int)(*(thearg4));
    }    

    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(nbins+1,1, mxREAL); 

    /* Assign pointer to the ouput matrix */ 
    binpop = mxGetPr(plhs[0]);
        
    /* Do the actual computations in a subroutine */
    ipb(data,d,n,localbins,nbins,nref,nt,binpop);
    
    /*free allocated memory (if any) */
    /*    if (nrhs>0)
      free(data);
    if (nrhs>1)
      free(bins);
    if (nrhs>2)
      free(thearg3);
    if (nrhs>3)
      free(thearg4);*/
    free(localbins);

    return;
    
}



