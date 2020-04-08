/*=================================================================
 *
 * INTERPOINT.C	
 *              .MEX file corresponding to INTERPOINT.M
 *              returns the number of interpoint distance of the 
 *              embedded time series in each bin
 *
 * The calling syntax is:
 *
 *		np = interpoint(y,de,tau,bins,nref,nt)
 *
 *
 * This is a MEX-file for MATLAB.  
 *=================================================================*/
/* $Revision: 1.5 $ */
#include <math.h>
#include <mex.h>
#include <stdlib.h>
#include <time.h>


void ipb(double	*data, /* the embedded data (d-by-n) */
	 int ndata, /* n-points */
	 int *de, /* embedding dimensions  */
	 int nde, /* number of embedding dimensions */
	 int tau, /* embedding lag */
	 double *bins, /* the bins */
	 int nbins, /* the number of bins */
	 long int nref, /* number of reference vectors to choose */
	 int nt, /* the Thieler band (forbidden points) */
	 double *binpop) /* the population in each bin (to be computed) */
     /*actually calculate the interpoint binning*/
{
  long int i,idum;
  int xi,yi,k,p,q; /* overdoes on indices */
  int nextde,maxde,maxn;
  double dist,val;
  bool toosmall,useall;
   char str[80];
  /*seed the RNG using the standard one*/
  idum=rand();

  /* determnie the maximum embedding window */
  maxde=*(de+nde-1);
  maxn=ndata-(maxde-1)*tau;
  
  /*initialise random number generator */
  srand((unsigned)time(NULL));

  /*use all the data or random reference vectors? */
  if (nref==0) {
    useall=1;
    if (nt>0) {
      nref=(maxn-2*nt+2);
      nref=nref*(maxn-2*nt+1);;
      nref += 2*(nt-1)*(maxn-nt+1)-(nt-1)*nt;
    } else {
      nref=maxn;
      nref=nref*nref;
    }
    xi=0;
    yi=0+nt-1;
    /*  nref=maxn*(maxn-nt); */
  } else {
    useall=0;
  }

  for (i=0; i<nref; i++) {

    /* get xi and yi */
    if (useall) {
      /* use all data and therefore just increment xi and yi */
      yi++;
      if (yi==maxn) {
	xi++;
	yi=0;
      }
      /* ensure that the points aren't too close  */     
      while (abs(xi-yi)<nt) {
	yi++;
	if (yi==maxn) {
	  yi=0;
	  xi++;
	}
      }
    } else {
      /* choose xi and yi randomly between 0 and maxn */
		xi=rand()%maxn+1;
		yi=rand()%maxn+1;	
      /* check that the points aren't too close  */     
      while (abs(xi-yi)<nt) {
		yi=rand()%maxn+1;
	}

    } /* if (useall) */
        

 /*  mexPrintf("xi = %i and yi = %i.\n", xi, yi);*/
   
    /* compute the square of the distance for these points  */
    dist=0;
    nextde=*(de);
    p=0;
    q=0;
    k=0;
    for (k=0; k<maxde; k++) {
      val = *(data+xi+k*tau)-*(data+yi+k*tau);
      dist += val*val;
      if (k==nextde-1) { /*is this one of the de we should check */
	
	/* find the bin this distance fits in  */
	if ( (*(bins+p))>dist) {
	  (*(binpop+p+q*(nbins+1)))++; /*increment the relevant bin */
	} else {
	  toosmall=1;
	  while (toosmall && p<nbins) {
	    if ( (*(bins+p))>dist) {
	      toosmall=0;
	      (*(binpop+p+q*(nbins+1)))++; /*increment the relevant bin */
	    }
	    p++;
	  } /* while */
	  p--;
	}
	/*increment nextde*/
	q++;
	nextde=*(de+q);
	
      } /* if (k==nextde) */
	
    } /* for (k= ...) */
    
  }/* for i */
  
}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     /* the MATLAB mex wrapper function */     
{ 
    double *data,*bins,*thearg,*binpop;
    int mrows,ncols,i,j;
    int d,n,nbins,tau,nt,nde;
    int *de;
    long int nref;
    bool warnThem;
    
    /* Check for proper number of arguments */
    
    if (nrhs > 6) { 
	mexErrMsgTxt("Too many input arguments."); 
    }
    if (nrhs < 4) {
        mexErrMsgTxt("Insufficient input arguments.");
    }
    if (nrhs < 5) {
      /* set nref=0 */
      nref=0;
    } 
    if (nrhs < 6) {
      /* set nt=0 */
      nt=0;
    } 
    if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /*Assign a pointer to the input matrix*/
    data = mxGetPr(prhs[0]);

    /* The first input argument (y) */     
    mrows = mxGetM(prhs[0]); 
    ncols = mxGetN(prhs[0]);
    /* check that y is a vector */
    if ((mrows!=1) && (ncols!=1)) {
      mexWarnMsgTxt("First input should be a vector");
    }
    n=mrows*ncols; /* number of data points */
    
    /* The second input argument (de) */
    thearg =mxGetPr(prhs[1]);
    mrows = mxGetM(prhs[1]);
    ncols = mxGetN(prhs[1]);
    /* check that all bins ascending */
    warnThem=0;
    nde=mrows*ncols;
    for (i=0; i<(nde-1); i++) 
      if ( *(thearg+i+1)<*(thearg+i)) 
	warnThem=1;
    if (warnThem) {
      mexErrMsgTxt("Error in interpoint: Second input not ascending");
    }
    /* alloc memory and make a copy of de (squaring as we go)*/
    /* and,  ... check that it is an integer */
    warnThem=0;
    de = (int *) calloc(nde,sizeof(int));
    for (i=0; i<nde; i++) {
      *(de+i) = (int)(*(thearg+i));
      if (floor(*(thearg+i))!=*(thearg+i)) {
	warnThem=1;
	}
    }
    if (warnThem) {
      mexWarnMsgTxt("Second input being rounded to integer values");
    }

    /* The third input argument (tau) */
    thearg = mxGetPr(prhs[2]);
    mrows = mxGetM(prhs[2]);
    ncols = mxGetN(prhs[2]);
    /* check that x is a scalar */
    if ((mrows!=1) || (ncols!=1)) {
      mexWarnMsgTxt("Third input should be a scalar");
    }
    /* check that it is an integer */
    if (floor(*(thearg))!=*(thearg)) {
      mexWarnMsgTxt("Third input being rounded to an integer value");
    }
    tau = (int)(*(thearg));
    
    /* The fourth input arg (bins) */
    thearg = mxGetPr(prhs[3]);
    mrows = mxGetM(prhs[3]);
    ncols = mxGetN(prhs[3]);
    /* check that x is a vector */
    if ((mrows!=1) && (ncols!=1)) {
      mexWarnMsgTxt("Fourth input should be a vector");
    }
    nbins = mrows*ncols;
    /* check that all bins ascending */
    warnThem=0;
    for (i=0; i<(nbins-1); i++) 
      if ( *(thearg+i+1)<*(thearg+i)) 
	warnThem=1;
    if (warnThem) {
      mexErrMsgTxt("Error in interpoint: Second input not ascending");
    }
    /*alloc memory and make a copy of bins (squaring as we go)*/
    bins = (double *) calloc(nbins,sizeof(double));
    for (i=0; i<nbins; i++) 
      *(bins+i) = (*(thearg+i))*(*(thearg+i));

    /* The fifth input argument (nref) */
    if (nrhs>4){
      thearg =mxGetPr(prhs[4]);
      mrows = mxGetM(prhs[4]);
      ncols = mxGetN(prhs[4]);
      /* check that x is a scalar */
      if ((mrows!=1) || (ncols!=1)) {
	mexWarnMsgTxt("Fifth input should be a scalar");
      }
      /* check that it is an integer */
      warnThem=0;
      if (floor(*(thearg))!=*(thearg)) {
	mexWarnMsgTxt("Fifth input being rounded to an integer value");
	}
      nref=(long int)(*(thearg));
    }

    /* The sixth input argument (nt) */
    if (nrhs==6){
      thearg =mxGetPr(prhs[5]);
      mrows = mxGetM(prhs[5]);
      ncols = mxGetN(prhs[5]);
      /* check that x is a scalar */
      if ((mrows!=1) || (ncols!=1)) {
	mexWarnMsgTxt("Sixth input should be a scalar");
      }
      /* check that it is an integer */
      warnThem=0;
      if (floor(*(thearg))!=*(thearg)) {
	mexWarnMsgTxt("Fourth input being rounded to an integer value");
	}
      nt=(int)(*(thearg));
    }    

    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(nbins+1,nde, mxREAL); 

    /* Assign pointer to the ouput matrix */ 
    binpop = mxGetPr(plhs[0]);
        
    /* Do the actual computations in a subroutine */
    ipb(data,n,de,nde,tau,bins,nbins,nref,nt,binpop);

    return;
    
}



