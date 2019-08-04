/*
 *  mGAP.c
 *
 *  Created by Wei Sun on Sun Mar 8th 2009.
 * 
 *  Last Edit: Nov 28th, 2009
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "optim.h"
#include "mGAP.h"

/*********************************************************************
 *
 * mGAPr
 *
 * The Iterative Adaptive Lasso, with X input from R
 *
 *********************************************************************/

void mGAPr(double* Ry, double* RX, double* RB,
  double* Rlambda, double* Rtau, int* dims, double* Repsilon, int* nIter, 
  int* b_update_order, double* Rscore, int* RscoNA, double* score2use, 
  double* lambda2use, double* tau2use,
  double* RB2, double* Rscore2, int* RscoNA2, double* score2use2, 
  double* lambda2use2, double* tau2use2,
  double* RB3, double* Rscore3, int* RscoNA3, double* score2use3, 
  double* lambda2use3, double* tau2use3, double* extbicgamma)
{
  double **y, **X, **B, **B2, **B3; 
  int n = dims[0];
  int p = dims[1];
  int h = dims[2];

  /* 
   * reorganize vector into matrix 
   * NOTE, each row of y is one response variable
   * and each column of X is a covariate
   */
  reorg(Ry, &y, h, n);
  reorg(RX, &X, p, n);
  
  reorg(RB, &B, p, h);
  reorg(RB2, &B2, p, h);
  reorg(RB3, &B3, p, h);

  mGAP(y, X, B, Rlambda, Rtau, dims, Repsilon, nIter, b_update_order, 
        Rscore, RscoNA, score2use, lambda2use, tau2use, 
        B2, Rscore2, RscoNA2, score2use2, lambda2use2, tau2use2, 
        B3, Rscore3, RscoNA3, score2use3, lambda2use3, tau2use3, 
        extbicgamma);
}


/*********************************************************************
 *
 * mGAPc
 *
 * The Iterative Adaptive Lasso, with X input from text file
 *
 *********************************************************************/

void mGAPc(double* Ry, char** fname, double* RB, 
  double* Rlambda, double* Rtau, int* dims, double* Repsilon, int* nIter, 
  int* b_update_order, double* Rscore, int* RscoNA, double* score2use, 
  double* lambda2use, double* tau2use, 
  double* RB2, double* Rscore2, int* RscoNA2, double* score2use2, 
  double* lambda2use2, double* tau2use2,
  double* RB3, double* Rscore3, int* RscoNA3, double* score2use3, 
  double* lambda2use3, double* tau2use3, double* extbicgamma)
{
  double **y, **X, **B, **B2, **B3; 
  int i;
  int n = dims[0];
  int p = dims[1];
  int h = dims[2];
  int offsetrow  = dims[5];
  int offsetcol  = dims[6];
  int transposeX = dims[7];

  /* 
   * reorganize vector into matrix 
   * NOTE, each row of y is one response variable
   * and each column of X is a covariate
   */
  reorg(Ry, &y, h, n);
  reorg(RB, &B, p, h);
  reorg(RB2, &B2, p, h);
  reorg(RB3, &B3, p, h);

  /* 
   * allocate memory 
   * X need to be a matrix of size p*n
   */
  X = (double **) Calloc(p, double*);
  X[0] = (double *) Calloc(n*p, double);
  if(X[0] == NULL){ error("fail to allocate memory of size %d\n", n*p); }
  for(i=1; i<p; i++){
    X[i] = X[0] + i*n;
  }
  
  /* read in data into matrix */
  // Rprintf("fname=%s, n=%d, p=%d, offsetrow=%d, offsetcol=%d, 
  // transposeX=%d\n", fname[0], n, p, offsetrow, offsetcol, transposeX);
  readtext(X, fname[0], p, n, offsetrow, offsetcol, transposeX);

  mGAP(y, X, B, Rlambda, Rtau, dims, Repsilon, nIter, b_update_order, 
  Rscore, RscoNA, score2use, lambda2use, tau2use, 
  B2, Rscore2, RscoNA2, score2use2, lambda2use2, tau2use2, 
  B3, Rscore3, RscoNA3, score2use3, lambda2use3, tau2use3, 
  extbicgamma);
  
  Free(X[0]);
  Free(X);
}

/*********************************************************************
 *
 * mGAP
 *
 * Iterative Multivariate Adaptive Lasso 
 *
 *********************************************************************/

void mGAP(double** y, double** X, double** RB, double* Rlambda, 
          double* Rtau, int* dims, double* Repsilon, int* nIter, 
          int* b_update_order, double* Rscore, int* RscoNA, 
          double* score2use, double* lambda2use, double* tau2use,
          double** RB2, double* Rscore2, int* RscoNA2, 
          double* score2use2, double* lambda2use2, double* tau2use2,
          double** RB3, double* Rscore3, int* RscoNA3, 
          double* score2use3, double* lambda2use3, double* tau2use3, 
          double* extbicgamma)
{
  int i, j, s, s1, s2, n, p, h, L;
  int Nmax, n_lambda, n_tau, p_max, *Jset, **scoNA, recursive, bicvec;
  int k1, k2, k, w, jk, found, n_b2use, *w_b2use;
  double XjAve, *X2, **B, **B1, **B_bar, B_js, **res, **resCov, b0;
  double *b_norm, *kappa, *bj_bar, bj_bar_norm, *sigma2, tmp;
  double **score, epsilon, lambda1, tau1, BIC, dfBIC, tmp1;
  double *b_bar_norm, X2_j, kappa_j, b_norm_j, tmp2, eta_norm;
  double sigma2_max, sigma2_min, bj_norm_lower, bj_norm_upper;

  double *B_j_init, *B_j_solu, Fmin, *ex, flower, fupper, r_js;
  double control[5];

  double **score2, **score3, extBIC, extBICgg, ebg;
  int **scoNA2, **scoNA3;
  
  
  int isfail, fncount;
  int dim[3];
  
  time_t timer;

  // a vector indicating whether b is used or not
  int* b2use;
  
  n = dims[0];
  p = dims[1];
  h = dims[2];
  L = dims[3];
  Nmax = dims[4];
  n_lambda = dims[8];
  n_tau    = dims[9];
  p_max    = dims[10];
  recursive = dims[11];
  bicvec = dims[12];
  
  epsilon  = *Repsilon;
  ebg = *extbicgamma;
  
  int DEBUG = 0;
    
  if(DEBUG){
    Rprintf("n=%d p=%d h=%d L=%d Nmax=%d \n", n, p, h, L, Nmax);
    Rprintf("n_lambda=%d n_tau=%d p_max=%d\n", n_lambda, n_tau, p_max);
  }

  // BIC score
  *score2use = DBL_MAX; // 1.79769e+308
  *score2use2 = DBL_MAX; 
  *score2use3 = DBL_MAX; 

  /* reorganize vector into matrix */
  reorg(Rscore, &score, n_lambda, n_tau);
  reorg_int(RscoNA, &scoNA, n_lambda, n_tau);
  
  reorg(Rscore2, &score2, n_lambda, n_tau);
  reorg_int(RscoNA2, &scoNA2, n_lambda, n_tau);
  
  reorg(Rscore3, &score3, n_lambda, n_tau);
  reorg_int(RscoNA3, &scoNA3, n_lambda, n_tau);

  /* allocate memory */
  sigma2 = (double *)Calloc(h, double);
  kappa  = (double *)Calloc(p, double);
  X2     = (double *)Calloc(p, double);
  b_norm = (double *)Calloc(p, double);
  b_bar_norm = (double *)Calloc(p, double);
  ex         = (double *)Calloc(2*h+1, double);
  B_j_init   = (double *)Calloc(h, double);
  B_j_solu   = (double *)Calloc(h, double);
  
  Jset    = (int *)Calloc(p, int);
  b2use   = (int *)Calloc(p, int);
  w_b2use = (int *)Calloc(p_max, int);

  res = (double**) Calloc(h, double*);
  res[0] = (double*) Calloc(n*h, double);
  for(s=1; s<h; s++){
    res[s] = res[0] + s*n;
  }

  B = (double**) Calloc(p, double*);
  B[0] = (double*) Calloc(p*h, double);
  for(j=1; j<p; j++){
    B[j] = B[0] + j*h;
  }

  B1 = (double**) Calloc(p, double*);
  B1[0] = (double*) Calloc(p*h, double);
  for(j=1; j<p; j++){
    B1[j] = B1[0] + j*h;
  }

  B_bar = (double**) Calloc(p, double*);
  B_bar[0] = (double*) Calloc(p*h, double);
  for(j=1; j<p; j++){
    B_bar[j] = B_bar[0] + j*h;
  }
  
  /* allocate memory */
  resCov = (double **)Calloc(h, double*);
  resCov[0] = (double*) Calloc(h*h, double);
  for(j=1; j<h; j++){
    resCov[j] = resCov[0] + j*h;
  }

  // Rprintf("memory allocation is done now\n");

  /* dim is the input parameter for function solve_bj */
  dim[0] = h;    // number of responses
  dim[1] = 1000; // maximum iterations
  dim[2] = 0;    // whether to trace the numerical optimization
  
  /* control is the control parameters for function solve_bj */
  control[0] = -1.0/0.0; // abstol = -Inf
  control[1] = 1e-8; // reltol
  control[2] = 1.0;  // alpha
  control[3] = 0.5;  // beta
  control[4] = 2.0;  // gamma

  /* initialize the order to update coefficients */
  for(j=0; j<p; j++){
    Jset[j] = j;
  }

  /**
   * Remove mean values of Xj and Calculate X2
   * sum square for each marker, i.e., each column of X
   */
  for(j=0; j<p; j++){
    XjAve = 0.0;
    for(i=0; i<n; i++){
      XjAve += X[j][i];
    }
    XjAve /= n;
    
    tmp = 0.0;
    for(i=0; i<n; i++){
      X[j][i] -= XjAve;
      tmp += X[j][i]*X[j][i];
    }
    X2[j] = tmp;
  }
  
  GetRNGstate();
  
  if(DEBUG){
    Rprintf("finish initializing data \n");
  }
  
  /**********************************************************
   * The EM algorithm
   **********************************************************/
  /* initialize residuals */
  for(s=0; s<h; s++){
    for(i=0; i<n; i++){
      res[s][i] = y[s][i];
    }
  }
  
  for(k1=0; k1 < n_lambda; k1++){
    lambda1 = Rlambda[k1];

    for(k2=0; k2 < n_tau; k2++){
      tau1 = Rtau[k2];
      
      if(DEBUG){
      	Rprintf("\n------------------------------------\n");
        Rprintf("lambda1 = %f, tau1 = %f\n", lambda1, tau1);
      }
      
      /**
       * step 1. Initialization
       */
        
      timer=time(NULL);
      // Rprintf("%s\n", asctime(localtime(&timer)));
      
      if(recursive){ 
        // use the coefficients from the precious combination of lambda and tau
        for(j=0; j<p; j++){
          for(s=0; s<h; s++){
            B1[j][s] = RB[j][s];
            B[j][s]  = RB[j][s];
          }
        }
        
        for(s=0; s<h; s++){
          sigma2[s] = var(res[s], n);
          // Rprintf("s=%d, sigma2[s]=%f\n", s, sigma2[s]);
        }
        
      }else{
        // assign the initail values to 0
        for(j=0; j<p; j++){
          for(s=0; s<h; s++){
            B1[j][s] = 0.0;
            B[j][s]  = 0.0;
          }
        }
        
        /* initialize residuals */
        for(s=0; s<h; s++){
          for(i=0; i<n; i++){
            res[s][i] = y[s][i];
          }
        }
        
        for(s=0; s<h; s++){
          sigma2[s] = var(res[s], n);
          // Rprintf("s=%d, sigma2[s]=%f\n", s, sigma2[s]);
        }
        
      }
      
      
      for(j=0; j<p; j++){
        
        kappa[j] = tau1;
        
        /* include all the covariates in the beginning 
         * this is not neccesary if RB=0, but just in case 
         */
        b2use[j] = 1;
        
      }

      /* k is used to check the convergence */
      k = 0;
    

      for(w=1; w<=Nmax; w++){
        /**
         * step 3.1 choose the order of updating b[j]
         */

        if(*b_update_order == 1){
          /* need to do nothing */
        }else if(*b_update_order == 2){
          rsample(Jset, p);
        }else{
          error("invalid b_update_order\n");
        }

        sigma2_max  = sigma2_min = sigma2[0];
        
        for(s=0; s<h; s++){
          if (sigma2_max < sigma2[s]) {
            sigma2_max = sigma2[s];
          }
          
          if (sigma2_min > sigma2[s]) {
            sigma2_min = sigma2[s];
          }        
        }
        
        /**
         * step 3.2 Update B[j][i]
         */

        for(jk=0; jk<p; jk++){
          j       = Jset[jk];
          bj_bar  = B_bar[j];
          X2_j    = X2[j];
          kappa_j = kappa[j];

          if (b2use[j]) {
            /* remove the effect of Xj from the residual*/
            for(s=0; s<h; s++){
              B_js  = B[j][s];
              for(i=0; i<n; i++){
                res[s][i] += X[j][i]*B_js;
              }
            }
          }
          
          bj_bar_norm = 0.0;
          eta_norm    = 0.0;

          for(s=0; s<h; s++){
            tmp1 = 0.0;
            for(i=0; i<n; i++){
              tmp1 += X[j][i]*res[s][i];
            }
            bj_bar[s] = tmp1/X2_j;
            bj_bar_norm += bj_bar[s]*bj_bar[s];
            
            tmp1 = 2.0*X2_j*bj_bar[s]*kappa_j/sigma2[s];
            eta_norm += tmp1*tmp1;
          }
          
          bj_bar_norm   = sqrt(bj_bar_norm);
          b_bar_norm[j] = bj_bar_norm;
          eta_norm      = sqrt(eta_norm);
          
          if (eta_norm <= lambda1) {
            b_norm_j = 0.0;
          }else{
            tmp = lambda1/(2.0*kappa_j*X2_j);
            bj_norm_lower = bj_bar_norm - tmp*sigma2_max;
            bj_norm_upper = bj_bar_norm - tmp*sigma2_min;
            
            if (bj_norm_upper < 1e-10) {
              if(bj_norm_upper < -1e-10){
                Rprintf("Warning: bj_norm_upper=%2e\n", bj_norm_upper);
              }
              b_norm_j = 0.0;
            }else {
              if (bj_norm_lower < 1e-10) { bj_norm_lower = 1e-10; }
              
              if(bj_norm_upper - bj_norm_lower < 1e-10){
                b_norm_j = bj_norm_upper;
              }else{
              
                flower = fb(bj_norm_lower, h, kappa_j, sigma2, X2_j, bj_bar, lambda1);
                fupper = fb(bj_norm_upper, h, kappa_j, sigma2, X2_j, bj_bar, lambda1);

                if (fabs(flower) < 1e-15) {
                  b_norm_j = bj_norm_lower;
                }else if (fabs(fupper) < 1e-15) {
                  b_norm_j = bj_norm_upper;
                }else if (flower*fupper < 0.0) {
                  b_norm_j = zeroin(bj_norm_lower, bj_norm_upper, fb, 0.0, h,
                                    kappa_j, sigma2, X2_j, bj_bar, lambda1);
                  if (b_norm_j < 0) { error("b_norm_j=%.2e\n", b_norm_j); }
                }else{
                  b_norm_j = -1.0;
                }
                
              }
            }
          }
          
          if (b_norm_j <= 1e-10) {
            b_norm[j] = 0.0;
            b2use[j]  = 0;
            for(s=0; s<h; s++){ B[j][s] = 0.0; }            
          }else{
            b_norm[j] = b_norm_j;
            b2use[j]  = 1;
            for(s=0; s<h; s++){
              tmp  = 2.0*kappa_j*X2_j*b_norm_j;
              B_js = bj_bar[s]*tmp/(tmp + lambda1*sigma2[s]);
              B[j][s] = B_js;
              for(i=0; i<n; i++){ res[s][i] -= X[j][i]*B_js; }
            }
          }
          
        } // end of updating B[j]
          
        /**
         * step 2. Update b0
         */
          
        for(s=0; s<h; s++){
          b0 = mean(res[s], n);
          if(fabs(b0) > 1e-5){
            Rprintf("Warning: s=%d, intercept=%f, and it is not 0.0\n", s, b0);
            for(i=0; i<n; i++){ res[s][i] -= b0; }
          }
        }

        /**
         * step 4 Update sigma^2
         */
        // Rprintf("Update sigma^2\n");

        for(s=0; s<h; s++){
          tmp = 0.0;
          for(i=0; i<n; i++){
            tmp += res[s][i]*res[s][i];
          }
          sigma2[s] = tmp/n;
        }

        /**
         * step 5 Update kappa
         */
        // Rprintf("Update kappa\n");

        for(j=0; j<p; j++){
          kappa[j] = b_norm[j] + tau1;
        }
    
        /**
         * check convergence
         */
        // Rprintf("check convergence\n");
        
        found = 0;
        for(j=0; j<p; j++){
          for(s=0; s<h; s++){
            if(fabs(B1[j][s] - B[j][s]) > epsilon){
              found = 1; k = 0; break;
            }
          }
          if (found) { break; }
        }
        
        for(j=0; j<p; j++){
          for(s=0; s<h; s++){
          	B1[j][s] = B[j][s];
          }
        }
    
        if(found == 0){
          k += 1;
          if(k>=L){ break; } // converged :)
        }
      } // end of loop for one set of lambda and tau 

      n_b2use = 0;
      for(j=0; j < p; j++){
        if(b2use[j]){
          w_b2use[n_b2use] = j;
          n_b2use++;
          if(n_b2use >= p_max){
            break;
          }
        }
      }

      if(DEBUG){
        if(n_b2use >= p_max){
          n_b2use = 0;
          for(j=0; j < p; j++){ if(b2use[j]) n_b2use++; }
        }
        
        Rprintf("n_b2use = %d, coverged at %d iterations, k=%d \n", n_b2use, w, k);
      }

      /**
       * ignore this combinition of lambda and tau 
       * if too many covariates are chosen 
       */
      // if(n_b2use == 0 || n_b2use >= p_max){
      //   continue;
      // }
      
      if(n_b2use >= p_max){
        continue;
      }
      
      /**
       * calculate the residual covariance matrix
       */
      resCov[0][0] = sigma2[0];
      
      for(s1=1; s1<h; s1++){
        resCov[s1][s1] = sigma2[s1];
        for(s2=0; s2<s1; s2++){
          tmp = 0.0;
          for(i=0; i<n; i++){
            tmp += res[s1][i]*res[s2][i];
          }
          tmp /= n;
          resCov[s1][s2] = tmp;
          resCov[s2][s1] = tmp;
        }
      }
      
      dfBIC = 0.0;
      
      if(bicvec){
        for(j=0; j<n_b2use; j++){
          jk   = w_b2use[j];
          tmp  = lambda1/(2.0*(b_norm[jk] + tau1)*b_norm[jk]*X2[jk]);
          tmp1 = b_norm[jk]*b_norm[jk];
          for (s=0; s<h; s++) {
            r_js   = 1 + tmp*sigma2[s];
            dfBIC += r_js/(r_js*r_js - (r_js - 1)*B[jk][s]*B_bar[jk][s]/tmp1);
          }
        }
        // note resCov will be changed in the function determinant
        determinant(resCov, h, &BIC);
        BIC = n*log(BIC) + dfBIC*log(n);
        extBIC = BIC + 2.0*ebg*lchoose(p, n_b2use); //
        extBICgg = BIC + 4.0*ebg*n_b2use*log(p);
      }else{
        for(j=0; j<p; j++){
          for(s=0; s<h; s++){
            if(B[j][s] > 1e-10){
             	dfBIC++;
             }
          }
        }
        // note resCov will be changed in the function determinant
        determinant(resCov, h, &BIC);
        BIC = n*log(BIC) + dfBIC*log(n);
        extBIC = BIC + 2.0*ebg*lchoose(p*h, dfBIC); //
        extBICgg = BIC + 4.0*ebg*dfBIC*log(p*h);
      }      
      
      if(BIC < *score2use){
        *score2use  = BIC;
        *lambda2use = lambda1;
        *tau2use    = tau1;
        *nIter = w;
        for(j=0; j<p; j++){
          for(s=0; s<h; s++){
            RB[j][s] = B[j][s];
          }
        }
      }
      
      if(extBIC < *score2use2){
        *score2use2  = extBIC;
        *lambda2use2 = lambda1;
        *tau2use2    = tau1;
        for(j=0; j<p; j++){
          for(s=0; s<h; s++){
            RB2[j][s] = B[j][s];
          }
        }
      }
      
      if(extBICgg < *score2use3){
        *score2use3  = extBICgg;
        *lambda2use3 = lambda1;
        *tau2use3    = tau1;
        for(j=0; j<p; j++){
          for(s=0; s<h; s++){
            RB3[j][s] = B[j][s];
          }
        }
      }
      
      if(DEBUG){  
        Rprintf("n_b2use=%d, dfBIC = %f, BIC=%f\n", n_b2use, dfBIC, BIC);
        Rprintf("lambda = %f, tau = %f\n", lambda1, tau1);
        Rprintf("=================================================\n\n");
      }

      score[k1][k2] = BIC;
      scoNA[k1][k2] = 0;
      
      score2[k1][k2] = extBIC;
      scoNA2[k1][k2] = 0;
      
      score3[k1][k2] = extBICgg;
      scoNA3[k1][k2] = 0;
      
    }
  }

  Free(sigma2);
  Free(kappa);
  Free(X2);
  Free(b_norm);
  Free(b_bar_norm);
  Free(ex);
  Free(B_j_init);
  Free(B_j_solu);
  
  Free(Jset);
  Free(b2use);
  Free(w_b2use);

  Free(res[0]);
  Free(res);

  Free(resCov[0]);
  Free(resCov);

  Free(B[0]);
  Free(B);

  Free(B1[0]);
  Free(B1);
  
  Free(B_bar[0]);
  Free(B_bar);
}
