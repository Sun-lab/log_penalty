/*
 * This c program is orginally downloaded from http://www.netlib.org/c/
 * The original code is modified to fit in other parameters of the 
 * function of interest
 */


/**********************************************************************
 *
 *  function to calculate ||b||
 *  equation (10) in supplementary materials
 **********************************************************************/

double fb(double b_norm_j, int h, double kappa_j, double *sigma2, 
		  double X2_j, double *bj_bar, double lambda)
{
	double val = 1/(4.0*kappa_j*kappa_j);
	double tmp = 2.0*kappa_j*b_norm_j;
	
	int s;
	double tmp1;
	
	for(s=0; s<h; s++){
		tmp1 = 1.0/(tmp + lambda*sigma2[s]/X2_j);
		val -= tmp1*tmp1*bj_bar[s]*bj_bar[s];
	}
	
	return(val);
}

/*
 ************************************************************************
 *	    		    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,tol)
 *	double ax; 			Root will be seeked for within
 *	double bx;  			a range [ax,bx]
 *	double (*f)(double x);		Name of the function whose zero
 *					will be seeked for
 *	double tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bissection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bissection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is 
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bissection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */

#include "math.h"

double zeroin(double ax, double bx, double (*f)(), double tol, int h, 
 double kappa_j, double *sigma2, double X2_j, double *bj_bar, double lambda)
/***
 * ax: Left border of the range where the root is seeked
 * bx: Left border of the range 
 * *f: Function under investigation
 * tol: Acceptable tolerance
 ***/
{
  double a,b,c;			/* Abscissae, descr. see above	*/
  double fa;				/* f(a)				*/
  double fb;				/* f(b)				*/
  double fc;				/* f(c)				*/
  double EPSILON = 1e-10;
  
  a = ax;  b = bx;
  fa = (*f)(a, h, kappa_j, sigma2, X2_j, bj_bar, lambda);
  fb = (*f)(b, h, kappa_j, sigma2, X2_j, bj_bar, lambda);
  c = a;   fc = fa;

  for(;;)		/* Main iteration loop	*/
  {
    /* Distance from the last but one to the last approximation */
    double prev_step = b-a;
    /* Actual tolerance		*/
    double tol_act;			
    /* Interpolation step is calculated in the form p/q; 
    division operations is delayed until the last moment */
    double p; 			
    double q;
    /* Step at this iteration */		    
    double new_step;
   
    if( fabs(fc) < fabs(fb) )
    { /* Swap data for b to be the best approximation	*/
      a = b;  b = c;  c = a;
      fa=fb;  fb=fc;  fc=fa;
    }
    tol_act = 2*EPSILON*fabs(b) + tol/2;
    new_step = (c-b)/2;

    if( fabs(new_step) <= tol_act || fb == (double)0 )
      return b;				/* Acceptable approx. is found	*/

    /**
     * Decide if the interpolation can be tried	
     * If prev_step was large enough, and was in true direction, 
     * Interpolatiom may be tried
     */
    if( fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb) )
    {
      register double t1,cb,t2;
      cb = c-b;
      if( a==c )	/* If we have only two distinct	*/
      {				    /* points linear interpolation can only be applied	*/
        t1 = fb/fa;
        p = cb*t1;
        q = 1.0 - t1;
      }
      else				/* Quadric inverse interpolation*/
      {
        q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
        p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
        q = (q-1.0) * (t1-1.0) * (t2-1.0);
      }
      if( p > (double)0 )
        /* p was calculated with the opposite sign; make p positive*/
        q = -q;
      else	/* and assign possible minus to	q */
        p = -p;
      
      /* If b+p/q falls in [b,c] and isn't too large*/
      if( p < (0.75*cb*q-fabs(tol_act*q)/2) && p < fabs(prev_step*q/2) ){
        new_step = p/q;
        /* If p/q is too large then the bissection procedure
         *  can reduce [b,c] range to more extent */
      }
    }

    if( fabs(new_step) < tol_act ){
      /* Adjust the step to be not less than tolerance */
      if( new_step > (double)0 )
      	new_step = tol_act;
      else
      	new_step = -tol_act;
    }
    
    a = b;  fa = fb;			/* Save the previous approx.	*/
    b += new_step;
    /* Do step to a new approxim.	*/
    fb = (*f)(b, h, kappa_j, sigma2, X2_j, bj_bar, lambda);	
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
    {  /* Adjust c for it to have a sign opposite to that of b	*/
      c = a;  fc = fa;
    }
  }
}

