#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

static double *vector(int n)
{
  double *v;
  v = Calloc(n, double);
  return v;
}

static void free_vector(double *v)
{
  Free(v);
}

static int checkConvergence(double *beta, double *beta_old,
                            double eps, int l, int J)
{
  int j;
  int converged = 1;
  for (j=0; j < J; j++)
    {
      if (beta[l*J+j]!=0 & beta_old[j]!=0)
	{
	  if (fabs((beta[l*J+j]-beta_old[j])/beta_old[j]) > eps)
	    {
	      converged = 0;
	      break;
	    }
	}
      else if (beta[l*J+j]==0 & beta_old[j]!=0)
	{
	  converged = 0;
	  break;
	}
      else if (beta[l*J+j]!=0 & beta_old[j]==0)
	{
	  converged = 0;
	  break;
	}
    }
  return(converged);
}

static double S(double z, double l)
{
  if (z > l)  return(z - l);
  if (z < -l) return(z + l);
  return(0);
}

static double MCP(double z, double l1, double gamma, double v, int Index)
{
  double s=0;
  
  if (Index==1) {

    if (z > 0) s = 1;
    else if (z < 0) s = -1;
    if (fabs(z) <= l1) return(0);
    else if (fabs(z) <= gamma*l1) return(s*(fabs(z)-l1)/(v*(1-1/gamma)));
    else return(z/(v));
  }else{
      return(z/(v));
  }
}




static double SCAD(double z, double l1, double gamma, double v, int Index)
{
  double s=0;
  if (Index==1) {

    if (z > 0) s = 1;
    else if (z < 0) s = -1;
    if (fabs(z) <= l1) return(0);
    else if (fabs(z) <= (l1+l1))
      return(s*(fabs(z)-l1)/(v));
    else if (fabs(z) <= gamma*l1)
      return(s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1))));
    else return(z/(v));
    
  }else{
    return(z/(v));
  }
}



  //* LOG2 and SICA2 is for lm;
static double LOGlm(double z, double lambda, double tau, double bhat,
                    double v, int Index)
{
  double s=0;
  double threshold;
  
  if (Index==1) {
    threshold = ((lambda)/(fabs(bhat) + tau));
    if (z > 0) s = 1;
    else if (z < 0) s = -1;
    
    if (fabs(z) <= threshold/v) return(0);
    else return(z - s*threshold/v );
  }else{
    return(z);
  }
}



static double SICAlm(double z, double lambda, double tau, double bhat,
                     double v, int Index)
{
  double s=0;
  double threshold;
  if (Index==1) {
    
    threshold = lambda*tau*(1 + tau)/((fabs(bhat) + tau)*(fabs(bhat) + tau));
    
    if (z > 0) s = 1;
    else if (z < 0) s = -1;
    if (fabs(z) <= threshold/v) return(0);
    else return(z - s*threshold/v );
    
  }else{
    return(z);
  }
}




static double LOGbinom(double z, double lambda, double tau, double bhat,
                       double v, double n, int Index)
{
  double s=0;
  double threshold;
  
  if (Index==1) {
    threshold = ((lambda)/(fabs(bhat) + tau));
    if (z > 0) s = 1;
    else if (z < 0) s = -1;
    if (fabs(z)/v <= threshold/(n*v)) return(0);
    else return(z/v - s*threshold/(n*v) );
    
  } else {
    return(z/v);
  }
}

static double SICAbinom(double z, double lambda, double tau, double bhat,
                        double v, double n, int Index)
{
  double s=0;
  double threshold;
  
  if (Index==1) {

    threshold = lambda*tau*(1 + tau)/((fabs(bhat) + tau)*(fabs(bhat) + tau));
    if (z > 0) s = 1;
    else if (z < 0) s = -1;
    if (fabs(z)/v <= threshold/(n*v)) return(0);
    else return(z/v - s*threshold/(n*v) );
    
  }else{
    return(z/v);
  }
}




static void gaussian_model(double *beta, int *iter, double *x, double *y,
                           int *n_, int *p_, char **penalty_, double *lambda,
                           double *tau,int *L_, double *eps_, int *max_iter_,
                           double *gamma_,  int *dfmax_, int *user_,
                           int *ChooseXindex_, int *restriction_)
{
  /* Declarations */
  int L=L_[0]; int p=p_[0]; int n=n_[0]; int max_iter=max_iter_[0];
  double eps=eps_[0]; double gamma=gamma_[0];
  char *penalty=penalty_[0];
  int dfmax=dfmax_[0]; int user=user_[0];
  int restriction = restriction_[0];
  
  int converged, active, lstart;
  double *r, *beta_old;
  
  r = vector(n); for (int i=0;i<n;i++) r[i] = y[i];
  beta_old = vector(p);
  if (user) lstart = 0;
  else lstart = 1;

  /* Path */
  for (int l=lstart;l<L;l++){
    if (l != 0) for (int j=0;j<p;j++) beta_old[j] = beta[(l-1)*p+j];
    
    while (iter[l] < max_iter){
      converged = 0;
      iter[l] = iter[l] + 1;

      /* Check dfmax */
      active = 0;
      for (int j=0;j<p;j++) if (beta[l*p+j]!=0) active++;
      
      if (active > dfmax){
        for (int ll=l;ll<L;ll++){
          for (int j=0;j<p;j++) beta[ll*p+j] = R_NaReal;
        }
        free_vector(beta_old);
        free_vector(r);
        return;
      }

      /* Covariates */
      for (int j=0;j<p;j++){
      
        /* Calculate z */
        double z = 0;
        for (int i=0;i<n;i++) z = z + x[j*n+i]*r[i];
        z = z/n + beta_old[j];
        
        if (restriction==1 && ChooseXindex_[j]==1 && z < 0) {
          beta[l*p+j] = 0.0;
          
        }else{
          /* Update beta_j */
          if (strcmp(penalty,"MCP")==0){
            beta[l*p+j] = MCP(z,lambda[l],gamma,1,ChooseXindex_[j]);
          }
          
          if (strcmp(penalty,"SCAD")==0){
            beta[l*p+j] = SCAD(z,lambda[l],gamma,1,ChooseXindex_[j]);
          }
          
          if (strcmp(penalty,"LOG")==0){
            beta[l*p+j] = LOGlm(z,lambda[l],tau[l],beta_old[j],n,
                                ChooseXindex_[j]);
          }
          
          if (strcmp(penalty,"SICA")==0){
            beta[l*p+j] = SICAlm(z,lambda[l],tau[l],beta_old[j],n,
                                 ChooseXindex_[j]);
          }
        }
        
        /* Update r */
        if (beta[l*p+j] != beta_old[j]){
          for (int i=0;i<n;i++)
            r[i] = r[i] - (beta[l*p+j] - beta_old[j])*x[j*n+i];
        }
      }

      /* Check for convergence */
      if (checkConvergence(beta,beta_old,eps,l,p)){
        converged  = 1;
        break;
      }
      
      for (int j=0;j<p;j++) beta_old[j] = beta[l*p+j];
      
    }
      /*if (converged==0) warning("Failed to converge");*/
  }

  free_vector(beta_old);
  free_vector(r);
}

static void binomial_model(double *beta0, double *beta, int *iter, double *x,
                           double *y, int *n_, int *p_, char **penalty_,
                           double *lambda, double *tau, int *L_, double *eps_,
                           int *max_iter_, double *gamma_, int *dfmax_,
                           int *user_, int *ChooseXindex_, int *restriction_)
{
  /* Declarations */
  int L=L_[0];int p=p_[0];int n=n_[0];int max_iter=max_iter_[0];
  double eps=eps_[0]; double gamma=gamma_[0];
  char *penalty=penalty_[0];
  int dfmax=dfmax_[0]; int user=user_[0];
  int restriction = restriction_[0];

  int converged, active, lstart;
  double beta0_old;
  double *r, *w, *beta_old;
  r = vector(n);
  w = vector(n);
  beta_old = vector(p);
  if (user) lstart = 0;
  else lstart = 1;

  /* Initialization */
  double ybar=0;
  for (int i=0;i<n;i++) ybar = ybar + y[i];
  ybar = ybar/n;
  if (lstart==0) beta0_old = log(ybar/(1-ybar));
  else beta0[0] = log(ybar/(1-ybar));

  /* Path */
  double xwr, xwx, eta, pi, yp, yy, z, v;
  
  for (int l=lstart;l<L;l++){
    
    if (l != 0){
      beta0_old = beta0[l-1];
      for (int j=0;j<p;j++) beta_old[j] = beta[(l-1)*p+j];
    }

    while (iter[l] < max_iter){
      converged = 0;
      iter[l] = iter[l] + 1;

      /* Check dfmax */
      active = 0;
      for (int j=0;j<p;j++) if (beta[l*p+j]!=0) active++;
      
      if (active > dfmax){
        for (int ll=l;ll<L;ll++){
          beta0[ll] = R_NaReal;
          for (int j=0;j<p;j++) beta[ll*p+j] = R_NaReal;
        }
        free_vector(beta_old);
        free_vector(w);
        free_vector(r);
        return;
      }

      /* Approximate L */
      yp=yy=0;
      
      for (int i=0;i<n;i++){
        eta = beta0_old;
        for (int j=0;j<p;j++) eta = eta + x[j*n+i]*beta_old[j];
        pi = exp(eta)/(1+exp(eta));
        if (pi > .9999){
          pi = 1;
          w[i] = .0001;
        }else if (pi < .0001){
          pi = 0;
          w[i] = .0001;
        }else{
          w[i] = pi*(1-pi);
          r[i] = (y[i] - pi)/w[i];
          yp = yp + pow(y[i]-pi,2);
          yy = yy + pow(y[i]-ybar,2);
        }
      }
      
      if (yp/yy < .01){
        /*warning("Model saturated; exiting...");*/
        for (int ll=l;ll<L;ll++){
          beta0[ll] = R_NaReal;
          for (int j=0;j<p;j++) beta[ll*p+j] = R_NaReal;
        }
        free_vector(beta_old);
        free_vector(w);
        free_vector(r);
        return;
      }

      /* Intercept */
      xwr = xwx = 0;
      
      for (int i=0;i<n;i++){
        xwr = xwr + w[i]*r[i];
        xwx = xwx + w[i];
      }
      
      beta0[l] = xwr/xwx + beta0_old;
      for (int i=0; i<n; i++) r[i] = r[i] - (beta0[l] - beta0_old);

      /* Covariates */
      for (int j=0;j<p;j++){
        /* Calculate z */
        xwr=0;
        xwx=0;
        
        for (int i=0;i<n;i++){
          xwr = xwr + x[j*n+i]*w[i]*r[i];
          xwx = xwx + x[j*n+i]*w[i]*x[j*n+i];
        }
        
        z = xwr/n + (xwx/n)*beta_old[j];
        v = xwx/n;

        /* Update beta_j */
        if (restriction==1 && ChooseXindex_[j]==1 && z < 0) {
          beta[l*p+j] = 0.0;
          
        }else{

          if (strcmp(penalty,"MCP")==0){
            beta[l*p+j] = MCP(z,lambda[l], gamma,v,ChooseXindex_[j]);
          }
          
          if (strcmp(penalty,"SCAD")==0){
            beta[l*p+j] = SCAD(z,lambda[l], gamma,v,ChooseXindex_[j]);
          }
          
          if (strcmp(penalty,"LOG")==0){
            beta[l*p+j] = LOGbinom(z,lambda[l],tau[l],beta_old[j],v,n,ChooseXindex_[j]);
          }
          
          if (strcmp(penalty,"SICA")==0){
            beta[l*p+j] = SICAbinom(z,lambda[l],tau[l],beta_old[j],v,n,ChooseXindex_[j]);
          }
        }

        /* Update r */
        if (beta[l*p+j] != beta_old[j]){
          for (int i=0;i<n;i++)
            r[i] = r[i] - (beta[l*p+j] - beta_old[j])*x[j*n+i];
        }
      }

      /* Check for convergence */
      if (checkConvergence(beta,beta_old,eps,l,p)){
        converged = 1;
        break;
      }
      
      beta0_old = beta0[l];
      for (int j=0;j<p;j++) beta_old[j] = beta[l*p+j];
    }
  }

  free_vector(beta_old);
  free_vector(w);
  free_vector(r);
}

static const R_CMethodDef cMethods[] = {
  {"gaussian_model", (DL_FUNC) &gaussian_model, 17},
  {"binomial_model", (DL_FUNC) &binomial_model, 18},
  NULL
};

void R_init_PEN(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

