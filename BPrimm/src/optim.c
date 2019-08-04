/*
 * =======================================================================
 *  Use Nelder-Mead method to find the MLE of b_j
 *  well, this is only used to verify the results based on zeroin
 *
 *  We could use thie general optimization algorithm to sovle b_j, but
 *  (1) it has a memory leaking problem in R
 *  (2) it is slow.. probably due to (1)
 * 
 *  any way, just for comparison and verify zeroin, this should be enough
 *
 * =======================================================================
 */

#include "stdio.h"
#include "math.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

/*
void nmmin(int n, double *xin, double *x, double *Fmin, optimfn fn,
           int *fail, double abstol, double intol, void *ex,
           double alpha, double beta, double gamma, 
           int trace, int *fncount, int maxit);
*/

double fn(int n, double *x, void *ex)
{
  int s;
  double tmp, b_norm=0.0, ff=0.0;
  double *Para = (double *) ex;
  double *bj_bar = Para;
  double *w_j    = &Para[n];
  double lambda  = Para[2*n];
  
  for(s=0; s<n; s++){
  	tmp = x[s] - bj_bar[s];
  	ff += tmp*tmp/w_j[s];
  	b_norm += x[s]*x[s];
  }
  
  ff += lambda*sqrt(b_norm);
  return(ff);
}

void solve_bj(int* dim, double *xin, double *x, double *Fmin, 
	int *fail, double* control, int *fncount, double *ex)
{
  int n, maxit, trace;
  double abstol, intol, alpha, beta, gamma;
  
  n = dim[0];
  maxit = dim[1];
  trace = dim[2];
  
  abstol = control[0];
  intol  = control[1];
  alpha  = control[2];
  beta   = control[3];
  gamma  = control[4];
  
	nmmin(n, xin, x, Fmin, fn, fail, abstol, intol, (void *)ex, 
	  alpha, beta, gamma, trace, fncount, maxit);
	
}
