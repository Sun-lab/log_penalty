/*
 *  bayes_shrinkage.c
 *
 *  Created by Wei Sun on Aug 15 2007.
 *
 *  Last updated by Wei Sun on Nov 18 2009.
 */
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "bayes_shrinkage.h"
#include "arms.h"

/*********************************************************************
 *
 * get_b_order
 *
 * the order of b for updating
 *
 *********************************************************************/
void get_b_order(int b_update_order, int* Jset, int w, int n, int p,
  double* y, double** X, double* b){

  int i, j;
  double ave_y, ave_Xj, Syx, Sxx, tmp, *btmp;

  btmp = (double *)calloc(p, sizeof(double));

  if(b_update_order == 1){
    for(j=0; j<p; j++){
      Jset[j] = j;
    }
  }else if(b_update_order == 2){
    rsample(Jset, p);
  }else if(b_update_order == 3){
    if(w==1){
      for(j=0; j<p; j++){
        Jset[j] = j;
        /* l = lm(y~X[,i]) */
        ave_y  = mean(y, n);
        ave_Xj = mean_j(X, j, n);
        Syx = Sxx = 0.0;
        for(i=0; i<n; i++){
          tmp = X[i][j] - ave_Xj;
          Syx += (y[i] - ave_y)*tmp;
          Sxx += tmp*tmp;
        }
        btmp[j] = fabs(Syx/Sxx);
      }
      revsort(btmp, Jset, p);
    }else{
      for(j=0; j<p; j++){
        Jset[j] = j;
        btmp[j] = fabs(b[j]);
      }
      revsort(btmp, Jset, p);
    }
  }else if(b_update_order == 4){
    for(j=0; j<p; j++){
      Jset[j] = p-j-1;
    }
  }else{
    error("invalid value for b_update_order!\n");
  }

  free(btmp);
}


/*********************************************************************
 *
 * IALr
 *
 * The Iterative Adaptive Lasso, with X input from R
 *
 *********************************************************************/

void IALr(double* y, double* RX, double* bInit, double* RbSample, double* Rb0,
  double* Rdelta, double* Rtau, int* dims, double* Repsilon, int* nIter, 
  int* b_update_order, double* Rscore, int* RscoNA, double* RBIC_gamma,
  int* Ridc, char** Rcriterion, double* score2use, double* delta2use, 
  double* tau2use, double* Rratio)
{
  double **X;
  int **IDC;

  int n = dims[0];
  int p = dims[1];
  int Times_cv = dims[12];

  /* reorganize vector into matrix */
  reorg_int(Ridc, &IDC, n, Times_cv);
  reorg(RX, &X, n, p);

  IAL(y, X, bInit, RbSample, Rb0, Rdelta, Rtau, dims, Repsilon, nIter, 
    b_update_order, Rscore, RscoNA, RBIC_gamma, IDC, Rcriterion, 
    score2use, delta2use, tau2use, Rratio);
}

/*********************************************************************
 *
 * IALc
 *
 * The Iterative Adaptive Lasso, with X input from text file
 *
 *********************************************************************/

void IALc(double* y, char** fname, double* bInit, double* RbSample, double* Rb0, 
  double* Rdelta, double* Rtau, int* dims, double* Repsilon, int* nIter, 
  int* b_update_order, double* Rscore, int* RscoNA, double* RBIC_gamma, 
  int* Ridc, char** Rcriterion, double* score2use, double* delta2use, 
  double* tau2use, double* Rratio)
{
  double **X;
  int ** IDC;
  int i;
  int n = dims[0];
  int p = dims[1];
  int offsetrow  = dims[5];
  int offsetcol  = dims[6];
  int transposeX = dims[7];
  int Times_cv   = dims[12];

  /* reorganize vector into matrix */
  reorg_int(Ridc, &IDC, n, Times_cv);

  /* allocate memory */
  X = (double **) malloc(n * sizeof(double*));
  X[0] = (double *) calloc(n*p, sizeof(double));
  if(X[0] == NULL){ error("fail to allocate memory of size %d\n", n*p); }
  for(i=1; i<n; i++){
    X[i] = X[0] + i*p;
  }
  
  /* read in data into matrix */
  // Rprintf("fname=%s, n=%d, p=%d, offsetrow=%d, offsetcol=%d, 
  // transposeX=%d\n", fname[0], n, p, offsetrow, offsetcol, transposeX);
  readtext(X, fname[0], n, p, offsetrow, offsetcol, transposeX);

  IAL(y, X, bInit, RbSample, Rb0, Rdelta, Rtau, dims, Repsilon, nIter, 
    b_update_order, Rscore, RscoNA, RBIC_gamma, IDC, Rcriterion,
    score2use, delta2use, tau2use, Rratio);
  
  free(X[0]);
  free(X);
}

/*********************************************************************
 *
 * IAL
 *
 * The Iterative Adaptive Lasso
 *
 *********************************************************************/

void IAL(double* y, double** X, double* bInit, double* RbSample, double* Rb0,
  double* Rdelta, double* Rtau, int* dims, double* Repsilon, int* nIter, 
  int* b_update_order, double* Rscore, int* RscoNA, double* RBIC_gamma, 
  int** IDC, char** Rcriterion, double* score2use, double* delta2use, 
  double* tau2use, double* Rratio)
{
  int i, j, jk, k, c1, t1, w, n, p, L, init_b, *Jset, n_b2use, *w_b2use;
  int s1=0, Times_cv, total;
  int Nmax, found=0, n_delta, n_tau, p_max, k1, k2, ncv, **scoNA;
  int *i_test, nTrain, nTest, checkThreshold;
  double xij, b0, bj_bar, v0, sj, tmp, epsilon, tau1, delta1, delta2;
  double diff=0.0, tsumrescv;
  double *b, *b1, *v, *res, *Xj2, **score, *rssCV, rss1, *rescv;
  double *y_test, **X_train, *y_train, *Threshold;
  double mse1, BIC_gamma = *RBIC_gamma, sigma, ratio=*Rratio;
  char* criterion = Rcriterion[0];
  
  

  // Rprintf("criterion = %s\n", criterion);
  
  n = dims[0];
  p = dims[1];
  L = dims[2];
  Nmax = dims[3];
  init_b   = dims[4];
  n_delta  = dims[8];
  n_tau    = dims[9];
  p_max    = dims[10];
  ncv      = dims[11];
  Times_cv = dims[12];
  epsilon  = *Repsilon;
  total    = ncv*Times_cv;

  // socre can be eitehr BIC or rss/n for CV
  *score2use = DBL_MAX;
  
  /* reorganize vector into matrix */
  reorg(Rscore, &score, n_delta, n_tau);
  reorg_int(RscoNA, &scoNA, n_delta, n_tau);

  /* allocate memory */
  w_b2use = (int *)calloc(p_max, sizeof(int));
  rssCV   = (double *)calloc(ncv*Times_cv, sizeof(double));

  if(strcmp(criterion, "cv")==0){
    y_train = (double *)calloc(n, sizeof(double));
    y_test  = (double *)calloc(n, sizeof(double));
    X_train = (double **)calloc(n, sizeof(double*));
    i_test  = (int *)calloc(n, sizeof(int));
  }
  
  Threshold = (double *)calloc(p, sizeof(double));

  b     = (double *)calloc(p, sizeof(double));
  b1    = (double *)calloc(p, sizeof(double));
  v     = (double *)calloc(p, sizeof(double));
  res   = (double *)calloc(n, sizeof(double));
  Xj2   = (double *)calloc(p, sizeof(double));
  Jset  = (int *)calloc(p, sizeof(int));
  rescv = (double *)calloc(n, sizeof(double));

  for(j=0; j<p; j++){
    Jset[j] = j;
  }

  // Rprintf("start IAL, n_delta=%d, n_tau=%d\n", n_delta, n_tau);
  
  /**
   * Calculate Xj2
   * sum square for each marker, i.e., each column of X
   * note this quantity is re-calculated if we use cv to 
   * select delta and tau
   */
  for(j=0; j<p; j++){
    Xj2[j] = 0.0;
    for(i=0; i<n; i++){
      xij = X[i][j];
      Xj2[j] += xij*xij;
    }
  }

  GetRNGstate();
  
  /* IAL calculation */
  for(k1=0; k1 < n_delta; k1++){
    delta1 = Rdelta[k1];
    delta2 = delta1 + 1.0;
    
    for(k2=0; k2 < n_tau; k2++){
      tau1 = Rtau[k2];
      
      /**********************************************************
       * BIC
       **********************************************************/
      if(strcmp(criterion, "BIC")==0){
  
        /**
         * step 1. Initialization
         */

        // Rprintf("k1=%d, k2=%d, delta1=%f, tau1=%f\n", k1, k2, delta1, tau1);

        b0 = 0.0;

        if(init_b){
          for(j=0; j<p; j++){
            b[j] = b1[j] = bInit[j];
            v[j] = (fabs(b[j]) + tau1)/delta2;
          }    
        }else{
          for(j=0; j<p; j++){
            b[j] = b1[j] = 0.0;
            v[j] = tau1/delta2;
          }
        }
        
        /* k is used to check the convergence */
        k = 0;
        
        /* initialize residuals */
        for(i=0; i<n; i++){
          res[i] = y[i];
          for(j=0; j<p; j++){
            res[i] -= X[i][j]*b[j];
          }
        }
        
        b0 = mean(res, n);
        
        for(i=0; i<n; i++){
          res[i] -= b0;
        }
        sigma =var(res, n); 
        v0 = var(res, n);

        for(w=1; w<=Nmax; w++){
                    
          /**
           * step 3.1 choose the order of updating b[j]
           */
          if(*b_update_order > 1){
            get_b_order(*b_update_order, Jset, w, n, p, y, X, b);
          }
      
          /**
           * step 3.2 Update b[j]
           */
          for(jk=0; jk<p; jk++){
            j = Jset[jk];
      
            /* remove the effect of Xj from the residual*/
            for(i=0; i<n; i++){
              res[i] = res[i] + X[i][j]*b[j];
            }
      
            bj_bar = 0.0;
            for(i=0; i<n; i++){
              bj_bar += X[i][j]*res[i];
            }
            bj_bar /= (Xj2[j]);
            sj   = (ratio*sigma + v0)/(Xj2[j]);
                        
            tmp = sj/v[j];
            if(bj_bar > tmp){
              b[j] = bj_bar - tmp;
            }else if(bj_bar < -tmp){
              b[j] = bj_bar + tmp;
            }else{
              b[j] = 0.0;
            }
      
            /* add the effect of Xj back into the residual */
            for(i=0; i<n; i++){
              res[i] = res[i] - X[i][j]*b[j];
            }
          }
      
          /**
           * step 2. Update b0 and res
           */
          
          diff = mean(res, n);
          
          for(i=0; i<n; i++){
            res[i] = res[i] - diff;
          }
          
          b0 = b0 + diff;
          
          
          /**
           * step 4 Update v0, i.e. sigma^2
           */
          v0 = 0.0;
          for(i=0; i<n; i++){
            v0 += res[i]*res[i];
          }
          v0 /= n;
          
          /**
           * step 5 Update v[j], i.e. kappa_j
           */
          // Rprintf("w=%d, tau=%e, delta=%e\n", w, tau, delta);
          for(j=0; j<p; j++){
            v[j] = (fabs(b[j]) + tau1)/delta2;
          }          
          
          /**
           * check convergence
           */          
          found = 0;
          for(j=0; j<p; j++){
            if(fabs(b1[j] - b[j]) > epsilon){
              found = 1;
              k = 0;
              break;
            }
          }
          
          for(j=0; j<p; j++){
            b1[j] = b[j];
          }
      
          if(found == 0){
            k += 1;
            if(k>=L){
              break; // converged :)
            }
          }
        }
        
        /*if(v0!=sigma){
          Rprintf("v0 = %f, sigma=%f\n", v0, sigma);
          checkThreshold = 1;}*/
          
        
                
        n_b2use = 0;
        for(j=0; j < p; j++){
          if(b[j]>1e-10){
            b[j]=b[j]+(ratio*sigma + v0)/(Xj2[j]*v[j]) - v0/(Xj2[j]*v[j]);
          }else { if(b[j]<-1e-10){
            b[j]=b[j]-(ratio*sigma + v0)/(Xj2[j]*v[j])  + v0/(Xj2[j]*v[j]) ;}
          }

          if(fabs(b[j]) > 1e-10){
            // w_b2use[n_b2use] = j;
            n_b2use++;
            if(n_b2use >= p_max){
              break;
            }
          }
        }
        for(i=0; i<n; i++){
          tmp = 0;
          for(j=0; j<p; j++){
            tmp = tmp + X[i][j]*b[j];
          }
          res[i] = y[i] - tmp - b0;
        }
        v0 = var(res, n);
        //*if(checkThreshold==1){Rprintf("v0 = %f\n", v0);}
        
        
        
        if(n_b2use == 0 || n_b2use >= p_max){
          continue;
        }
        
        tmp = 2*BIC_gamma*lchoose((double)p, (double)n_b2use);
        tmp = n*log(v0) + n_b2use*log(n) + tmp;
        
        if(tmp < *score2use){
          *score2use = tmp;
          *delta2use = delta1;
          *tau2use   = tau1;
          *nIter = w;
          *Rb0   = b0;
          for(j=0; j<p; j++){
            RbSample[j] = b[j];
          }
        }
        // Rprintf("delta = %f, tau = %f, BIC=%f\n", delta1, tau1, tmp);
        score[k1][k2] = tmp;
        scoNA[k1][k2] = 0;

      /**********************************************************
       * cross-validation
       **********************************************************/
      }else if(strcmp(criterion, "cv")==0){
      
        s1=0;
        
        for (t1=0; t1 < Times_cv; t1++){

          for(c1=0; c1 < ncv; c1++){
        
            /**
             * step 1.1 Initialization
             */
            b0 = 0.0;
            
            if(init_b){
              for(j=0; j<p; j++){
                b[j] = b1[j] = bInit[j];
                v[j] = (fabs(b[j]) + tau1)/delta2;
              }    
            }else{
              for(j=0; j<p; j++){
                b[j] = b1[j] = 0.0;
                v[j] = tau1/delta2;
              }
            }
            
            /**
             * step 1.2 Split the data
             */

            k = 0;
            nTrain = 0;
            nTest  = 0;
              
            for(i=0; i<n; i++){
              if(IDC[i][t1]!=c1){
                X_train[nTrain] = X[i];
                y_train[nTrain] = y[i];
                nTrain++;
              }else{
                i_test[nTest] = i;
                y_test[nTest] = y[i];
                nTest++;
              }
            }
              
            if(nTest + nTrain != n){
              error("nTest + nTrain != n: %d + %d != %d\n", nTest, nTrain, n);
            }
            
            for(j=0; j<p; j++){
              Xj2[j] = 0.0;
              for(i=0; i<nTrain; i++){
                xij = X_train[i][j];
                Xj2[j] += xij*xij;
              }
            }

            /* initialize residuals */
            b0 = mean(y_train, nTrain);

            for(i=0; i<nTrain; i++){
              rescv[i] = y_train[i] - b0;
            }

            v0 = var(y_train, nTrain);

            for(w=1; w<=Nmax; w++){
          
              /**
               * step 3.1 choose the order of updating b[j]
               */
              if(*b_update_order > 1)
              get_b_order(*b_update_order, Jset, w, nTrain, p, y_train, X_train, b);
          
              /**
               * step 3.2 Update b[j]
               */
              for(jk=0; jk<p; jk++){
                j = Jset[jk];
          
                /* remove the effect of Xj from the residual*/
                for(i=0; i<nTrain; i++){
                  rescv[i] = rescv[i] + X_train[i][j]*b[j];
                }
          
                bj_bar = 0.0;
                for(i=0; i<nTrain; i++){
                  bj_bar += X_train[i][j]*rescv[i];
                }
                
                bj_bar /= (Xj2[j]);
                sj      = v0/(Xj2[j]);
                tmp     = sj/v[j];
                
                if(bj_bar > tmp){
                  b[j] = bj_bar - tmp;
                }else if(bj_bar < -tmp){
                  b[j] = bj_bar + tmp;
                }else{
                  b[j] = 0.0;
                }
          
                /* add the effect of Xj back into the residual */
                for(i=0; i<nTrain; i++){
                  rescv[i] = rescv[i] - X_train[i][j]*b[j];
                }
              }
          
              tsumrescv = 0;
              for (i=0; i<nTrain; i++){
                tsumrescv = tsumrescv + rescv[i];
              }
              diff = tsumrescv / nTrain;

              for(i=0; i<nTrain; i++){
                rescv[i] = rescv[i] - diff;
              }
              
              b0 = b0 + diff;   
              
              /**
               * step 4 Update v0, i.e. sigma^2
               */
              v0 = 0.0;
              for(i=0; i<nTrain; i++){
                v0 += rescv[i]*rescv[i];
              }
              v0 /= nTrain;
          
              /**
               * step 5 Update v[j], i.e. kappa_j
               */
          
              // Rprintf("w=%d, tau=%e, delta=%e\n", w, tau, delta);
              for(j=0; j<p; j++){
                v[j] = (fabs(b[j]) + tau1)/delta2;
              }
          
              /**
               * check convergence
               */
              
              found = 0;
              for(j=0; j<p; j++){
                if(fabs(b1[j] - b[j]) > epsilon){
                  found = 1;
                  k = 0;
                  break;
                }
              }
              
              for(j=0; j<p; j++){
                b1[j] = b[j];
              }
          
              if(found == 0){
                k += 1;
                if(k>=L){
                  break;
                }
              }
            }
            
            n_b2use = 0;
            for(j=0; j < p; j++){
              if(fabs(b1[j]) > 1e-10){
                w_b2use[n_b2use] = j;
                n_b2use++;
                if(n_b2use >= p_max){
                  break;
                }
              }
            }
            
            /***
             * if no covarite or too many covariates are selected, 
             * just record the rss as the rss of y_test
             */
            if(n_b2use == 0 || n_b2use >= p_max){
              rss1 = 0.0;
              for(i=0; i < nTest; i++){
                tmp = mean(y_train,nTrain);              
                tmp  = y_test[i] -  tmp;
                rss1 += tmp*tmp;
              }
              
              rssCV[s1] = rss1/nTest;
              s1 = s1 +1;
            }else{
              rss1 = 0.0;
              for(i=0; i < nTest; i++){
                tmp = 0.0;
                for(j=0; j<n_b2use; j++){
                  tmp += X[i_test[i]][w_b2use[j]]*b[w_b2use[j]];
                }
                tmp  = y_test[i] - b0 - tmp;
                rss1 += tmp*tmp;
              }
              rssCV[s1] = rss1/nTest;
              s1 += 1;
            }
          }
        }

        if(s1 > total){
          error("s1 > total: %d > %d\n", s1, total);
        }
          
        tmp = mean(rssCV, total);
        
        if(tmp < *score2use){
          *score2use = tmp;
          *delta2use = delta1;
          *tau2use   = tau1;
        }
        /***
         * hm, so score would not be NA
         * since we take the variance of y_test as rss
         * even if too many or no varaible is selected.
         */
        score[k1][k2] = tmp;
        scoNA[k1][k2] = 0;
      }else{
        error("wrong criterion: %s\n", criterion);
      }
      // Rprintf("delta=%f, tau=%f, score=%f, b0=%f\n", delta1, tau1, tmp, b0);
    }
  }
  
  /**********************************************************
   * cross-validation
   *
   * if the criterion is cross-validation, we only selectedd 
   * delta and tau by cv in previous code, now we need to 
   * compute b0 and b using all the data again
   **********************************************************/
  if(strcmp(criterion, "cv")==0){
    delta1 = *delta2use;
    tau1   = *tau2use;
    delta2 = delta1 + 1;
    
    // Rprintf("\n\ndelta=%f, tau=%f\n", delta1, tau1);
    
    /**
     * Calculate Xj2
     * sum square for each marker, i.e., each column of X
     * note this quantity is re-calculated if we use cv to 
     * select delta and tau
     */
    for(j=0; j<p; j++){
      Xj2[j] = 0.0;
      for(i=0; i<n; i++){
        xij = X[i][j];
        Xj2[j] += xij*xij;
      }
    }
    
    /**
     * step 1. Initialization
     */
    b0 = 0.0;
    v0 = var(y, n);
    
    if(init_b){
      for(j=0; j<p; j++){
        b[j] = b1[j] = bInit[j];
        v[j] = (fabs(b[j]) + tau1)/delta2;
      }    
    }else{
      for(j=0; j<p; j++){
        b[j] = b1[j] = 0.0;
        v[j] = tau1/delta2;
      }
    }
    
    /* k is used to check the convergence */
    k = 0;

    /* initialize residuals */
    for(i=0; i<n; i++){
      res[i] = y[i];
    }

    for(w=1; w<=Nmax; w++){
      /**
       * step 2. Update b0
       */
      b0 = mean(res, n);
  
      for(i=0; i<n; i++){
        res[i] -= b0;
      }
  
      /**
       * step 3.1 choose the order of updating b[j]
       */
      if(*b_update_order > 1){
        get_b_order(*b_update_order, Jset, w, n, p, y, X, b);
      }
  
      /**
       * step 3.2 Update b[j]
       */
      for(jk=0; jk<p; jk++){
        j = Jset[jk];
  
        /* remove the effect of Xj from the residual*/
        for(i=0; i<n; i++){
          res[i] = res[i] + X[i][j]*b[j];
        }
  
        bj_bar = 0.0;
        for(i=0; i<n; i++){
          bj_bar += X[i][j]*res[i];
        }
        bj_bar /= (Xj2[j]);
        sj   = v0/(Xj2[j]);
        tmp  = sj/v[j];
        
        if(bj_bar > tmp){
          b[j] = bj_bar - tmp;
        }else if(bj_bar < -tmp){
          b[j] = bj_bar + tmp;
        }else{
          b[j] = 0.0;
        }
  
        /* add the effect of Xj back into the residual */
        for(i=0; i<n; i++){
          res[i] = res[i] - X[i][j]*b[j];
        }
      }
  
      /**
       * step 4 Update v0, i.e. sigma^2
       */
      v0 = 0.0;
      for(i=0; i<n; i++){
        v0 += res[i]*res[i];
      }
      v0 /= n;
  
      /**
       * step 5 Update v[j], i.e. kappa_j
       */
      // Rprintf("w=%d, tau=%e, delta=%e\n", w, tau, delta);
      for(j=0; j<p; j++){
        v[j] = (fabs(b[j]) + tau1)/delta2;
      }
  
      /**
       * check convergence
       */          
      found = 0;
      for(j=0; j<p; j++){
        if(fabs(b1[j] - b[j]) > epsilon){
          found = 1;
          k = 0;
          break;
        }
      }
      
      for(j=0; j<p; j++){
        b1[j] = b[j];
      }
  
      if(found == 0){
        k += 1;
        if(k>=L){
          break; // converged :)
        }
      }
    }

    *nIter = w;
    *Rb0   = b0;
    for(j=0; j<p; j++){
      RbSample[j] = b[j];
    }
    // Rprintf("nIter = %d, b0=%f\n", 2, b0);
  }
  
  PutRNGstate();
  
  free(b);
  free(b1);
  free(v);
  free(res);
  free(Xj2);
  free(Jset);
  free(w_b2use);
  free(rssCV);
  
  if(strcmp(criterion, "cv")==0){
    free(y_train);
    free(y_test);
    free(X_train);
    free(i_test);
  }

}

/*********************************************************************
 * end of IAL
 *********************************************************************/



/*********************************************************************
 *
 * ld_bj
 *
 * log density function for b_j used in Bayesian Adaptive Lasso
 *
 *********************************************************************/
struct bj_parm {
  double bj_bar;
  double sigma_j2;
  double kappa_j;
};

double ld_bj(double x, void *bj_parm)
{
	struct bj_parm *d;
	double y, tmp;

	d   = bj_parm;
	tmp = x - d->bj_bar;
	y   = -tmp*tmp/(2*d->sigma_j2) - fabs(x)/d->kappa_j;
	return y;
}

void bj_sampler(int* nsample, double* bj_bar,
  double* sigma_j2, double* kappa_j, double* bj, int* neval)
{
  /* variables for Adpative Rejection Metroplois Sampling (arms) */
  int err, ninit = 4, npoint = 50, ncent = 0;
  double xinit[10], xl = -100.0, xr = 100.0;
  xinit[0] = *bj_bar - 0.5;
  xinit[1] = *bj_bar - 0.05;
  xinit[2] = *bj_bar + 0.05;
  xinit[3] = *bj_bar + 0.5;
  double xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.0;
  int dometrop  = 0;   // Whether metropolis step is required
  double xprev  = 0.0; // Previous value from the Markov chain
  struct bj_parm d;

	d.bj_bar   = *bj_bar;
	d.sigma_j2 = *sigma_j2;
	d.kappa_j  = *kappa_j;

  GetRNGstate();
  err = arms(xinit,ninit,&xl,&xr,ld_bj,&d,&convex,npoint,
    dometrop,&xprev,bj,*nsample,qcent,xcent,ncent,neval);
  PutRNGstate();

  if(err>0){
    error("error in arms, error code = %d\n",err);
  }
}

/*********************************************************************
 *
 * ld_delta
 *
 * log density function for delta used in Bayesian AL or Bayesan t
 *
 *********************************************************************/
struct delta_parm {
  double tau;
  double sum_log_v; // sum_{j=0}^{p-1} log(sigma_j^2)
  int p;
};

double ld_delta(double x, void *delta_parm)
{
	struct delta_parm *d;
	double y;

	d = delta_parm;
	y = x*(d->p*log(d->tau) - d->sum_log_v) - d->p*lgammafn(x);
	return y;
}

/*********************************************************************
 *
 * bayes_AL
 *
 * The Bayesian version of adaptive Lasso
 *
 *********************************************************************/
void bayes_AL(double* y, double* RX, double* bInit, double* RbSample, 
  double* deltaSam, double* tauSam, int* dims, 
  double* Rdelta, double* Rtau, int* b_update_order)
{
  int i, j, jk, k, w, n, p, Nburn, Nthin, Nsample, Nall, init_b, *Jset;
  int update_delta, update_tau;
  double xij, b0, b0_bar, bj_bar, v0, s0, sj, delta, tau;
  double *b, bj, *v, *res, *Xj2, **X, **bSample, bj_t, sj_t, rate, tmp;

  /* variables for Adpative Rejection Metroplois Sampling (arms) */
  int err, neval=0, ninit = 4, npoint = 50, nsamp=1, ncent = 4;
  double binit[10];
  double bl = -100.0, br = 100.0;
  double bprev = 0.0;
  double xinit[10]={0.01, 2.0, 10.0, 20.0}, xl = 0.0, xr = 500.0;
  double xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.0;
  int dometrop  = 0;
  double xprev  = 0.01;
  
  struct delta_parm delta_data;
  struct bj_parm d;

  n = dims[0];
  p = dims[1];
  Nburn = dims[2];
  Nthin = dims[3];
  Nsample = dims[4];
  update_delta = dims[5];
  update_tau   = dims[6];
  init_b  = dims[7];

  delta_data.p = p;

  /**
   * if update_delta *Rdelta is the initial value of delta
   * otherwise it is the fixed value of detla
   * similar for tau
   */
  delta = *Rdelta;
  tau = *Rtau;

  // Rprintf("delta = %f, tau = %e\n", delta, tau);

  /* allocate memory */
  b    = (double *)calloc(p, sizeof(double));
  v    = (double *)calloc(p, sizeof(double));
  res  = (double *)calloc(n, sizeof(double));
  Xj2  = (double *)calloc(p, sizeof(double));
  Jset = (int *)calloc(p, sizeof(int));

  /* reoganize vector into matrix */
  reorg(RX, &X, n, p);
  reorg(RbSample, &bSample, p, Nsample);

  /**
   * Calculate Xj2
   * sum square for each marker, i.e., each column of X
   */
  for(j=0; j<p; j++){
    Xj2[j] = 0.0;
    for(i=0; i<n; i++){
      xij = X[i][j];
      Xj2[j] += xij*xij;
    }
  }

  /**
   * step 1. Initialization
   */
  b0 = 0.0;
  v0 = 0.01;
  if(init_b){
    for(j=0; j<p; j++){
      b[j] = bInit[j];
      v[j] = 0.01;
    }    
  }else{
    for(j=0; j<p; j++){
      b[j] = 0.0;
      v[j] = 0.01;
    }
  }

  Nall = Nburn + Nthin*Nsample;
  k = -1;

  for(j=0; j<p; j++){
    Jset[j] = j;
  }

  GetRNGstate();
  for(w=1; w<=Nall; w++){    
    /**
     * step 2. Update b0
     */
    for(i=0; i<n; i++){
      res[i] = y[i];
      for(j=0; j<p; j++){
        res[i] -= X[i][j]*b[j];
      }
    }
    b0_bar = mean(res, n);
    s0 = v0/n;
    b0 = rnorm(b0_bar, sqrt(s0));

    for(i=0; i<n; i++){
      res[i] -= b0;
    }

    /**
     * step 3.1 choose the order of updating b[j]
     */

    if(*b_update_order > 1){
      get_b_order(*b_update_order, Jset, w, n, p, y, X, b);
    }

    /**
     * step 3.2 Update b[j]
     */

    for(jk=0; jk<p; jk++){
      j = Jset[jk];

      /* remove the effect of Xj from the residual*/
      for(i=0; i<n; i++){
        res[i] = res[i] + X[i][j]*b[j];
      }

      bj_bar = 0.0;
      for(i=0; i<n; i++){
        bj_bar += X[i][j]*res[i];
      }
      bj_bar /= Xj2[j];
      sj   = v0/Xj2[j];

    	d.bj_bar   = bj_bar;
    	d.sigma_j2 = sj;
    	d.kappa_j  = v[j];
    	xprev      = b[j];
    	
      if(bj_bar > sj/v[j]){
        bj_t = bj_bar - sj/v[j];
      }else if(bj_bar < -sj/v[j]){
        bj_t = bj_bar + sj/v[j];
      }else{
      	bj_t = 0.0;
      }
      sj_t = sqrt(sj);
      
      binit[0] = bj_t - 2*sj_t;
      binit[1] = bj_t - 0.5*sj_t;
      binit[2] = bj_t + 0.5*sj_t;
      binit[3] = bj_t + 2*sj_t;
      
      err = arms(binit,ninit,&bl,&br,ld_bj,&d,&convex,npoint,
        dometrop,&xprev,&bj,nsamp,qcent,xcent,ncent,&neval);
      if(neval > 1000){
        warning("In total %d evaluations is used in ARMS\n", neval);
      }
      if(err>0){
        warning("error in arms, error code = %d, b[j] is not updated\n", err);
        bj = b[j];
      }else{
        b[j] = bj;
      }
      /* add the effect of Xj back into the residual */
      for(i=0; i<n; i++){
        res[i] = res[i] - X[i][j]*bj;
      }
    }

    /**
     * step 4 Update v0, i.e. sigma_0^2
     */
    v0 = 0.0;
    for(i=0; i<n; i++){
      v0 += res[i]*res[i];
    }
    v0 /= rchisq(n);

    /**
     * step 5 Update v[j], i.e. kappa_j
     */
    for(j=0; j<p; j++){
      v[j] = 2*(fabs(b[j]) + tau)/rchisq(2+2*delta);
    }

    /**
    * step 6 Update tau, if needed
    */
    if(update_tau){
      rate = 0.0;
      for(j=0; j<p; j++){
      	rate += 1/v[j];
      }
      tau = rgamma(p*delta, 1/rate);
      if(tau == 0){
      	error("tau=0, rate=%e\n", rate);
      }
    }

    /**
    * step 7 Update delta, if needed
    */
    if(update_delta){
    	delta_data.tau = tau;
    	tmp = 0.0;
    	for(j=0; j<p; j++){
    		tmp += log(v[j]);
    	}
    	delta_data.sum_log_v = tmp;
    	xprev = delta;
      err = arms(xinit,ninit,&xl,&xr,ld_delta,&delta_data,&convex,
        npoint,dometrop,&xprev,&delta,nsamp,qcent,xcent,ncent,&neval);
      if(err>0){
        error("error in arms, error code = %d\n",err);
      }
      for(j=0; j<ninit; j++){
      	xinit[j] = xcent[j];
      }
    }

    /**
     * last step Record b when neccesary
     */
    if(w > Nburn){
      if((w - Nburn) % Nthin == 0){
        k = k + 1;
        for(j=0; j<p; j++){
          bSample[j][k] = b[j];
        }
        deltaSam[k] = delta;
        tauSam[k] = tau;
      }
    }

  }

  PutRNGstate();

  free(b);
  free(v);
  free(res);
  free(Xj2);
  free(Jset);
}

/*********************************************************************
 *
 * bayes_t
 *
 * The Bayesian shrinkage mehtod
 *
 *********************************************************************/
void bayes_t(double* y, double* RX, double* bInit, double* RbSample,
  double* deltaSam, double* tauSam, int* dims,
  double* Rdelta, double* Rtau, int* b_update_order)
{
  int i, j, jk, k, w, n, p, Nburn, Nthin, Nsample, Nall, init_b, *Jset;
  int update_delta, update_tau;
  double xij, b0, b0_bar, bj_bar, v0, s0, sj, delta, tau;
  double *b, *v, *res, *Xj2, **X, **bSample, rate, tmp;

  /* variables for Adpative Rejection Metroplois Sampling (arms) */
  int err, neval = 0, ninit = 4, npoint = 50, nsamp = 1, ncent = 4;
  double xinit[10]={0.01, 2.0, 10.0, 20.0}, xl = 0.0, xr = 500.0;
  double xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.0;
  int dometrop  = 0;
  double xprev  = 0.01;

  struct delta_parm delta_data;

  n = dims[0];
  p = dims[1];
  Nburn = dims[2];
  Nthin = dims[3];
  Nsample = dims[4];
  update_delta = dims[5];
  update_tau   = dims[6];
  init_b  = dims[7];

  delta_data.p = p;
  
  /**
   * if update_delta *Rdelta is the initial value of delta
   * otherwise it is the fixed value of detla
   * similar for tau
   */
  delta = *Rdelta;
  tau = *Rtau;

  /* allocate memory */
  b    = (double *)calloc(p, sizeof(double));
  v    = (double *)calloc(p, sizeof(double));
  res  = (double *)calloc(n, sizeof(double));
  Xj2  = (double *)calloc(p, sizeof(double));
  Jset = (int *)calloc(p, sizeof(int));

  /* reoganize vector into matrix */
  reorg(RX, &X, n, p);
  reorg(RbSample, &bSample, p, Nsample);

  /**
   * Calculate Xj2
   * sum square for each marker, i.e., each column of X
   */
  for(j=0; j<p; j++){
    Xj2[j] = 0.0;
    for(i=0; i<n; i++){
      xij = X[i][j];
      Xj2[j] += xij*xij;
    }
  }

  /**
   * step 1. Initialization
   */
  b0 = 0.0;
  v0 = 0.01;
  if(init_b){
    for(j=0; j<p; j++){
      b[j] = bInit[j];
      v[j] = 0.01;
    }    
  }else{
    for(j=0; j<p; j++){
      b[j] = 0.0;
      v[j] = 0.01;
    }
  }

  Nall = Nburn + Nthin*Nsample;
  k = -1;

  for(j=0; j<p; j++){
    Jset[j] = j;
  }

  GetRNGstate();
  for(w=1; w<=Nall; w++){
    /**
     * step 2. Update b0
     */
    for(i=0; i<n; i++){
      res[i] = y[i];
      for(j=0; j<p; j++){
        res[i] -= X[i][j]*b[j];
      }
    }
    b0_bar = mean(res, n);
    s0 = v0/n;
    b0 = rnorm(b0_bar, sqrt(s0));

    for(i=0; i<n; i++){
      res[i] -= b0;
    }

    /**
     * step 3.1 choose the order of updating b[j]
     */
    if(*b_update_order > 1){
      get_b_order(*b_update_order, Jset, w, n, p, y, X, b);
    }

    /**
     * step 3.2 Update b[j]
     */

    for(jk=0; jk<p; jk++){
      j = Jset[jk];

      /* remove the effect of Xj from the residual*/
      for(i=0; i<n; i++){
        res[i] = res[i] + X[i][j]*b[j];
      }

      bj_bar = 0.0;
      for(i=0; i<n; i++){
        bj_bar += X[i][j]*res[i];
      }
      bj_bar /= (Xj2[j] + v0/v[j]);
      sj   = v0/(Xj2[j] + v0/v[j]);
      b[j] = rnorm(bj_bar, sqrt(sj));

      /* add the effect of Xj back into the residual*/
      for(i=0; i<n; i++){
        res[i] = res[i] - X[i][j]*b[j];
      }
    }

    /**
     * step 4 Update v0, i.e. sigma_0^2
     */
    v0 = 0.0;
    for(i=0; i<n; i++){
      v0 += res[i]*res[i];
    }
    v0 /= rchisq(n);

    /**
     * step 5 Update v[j], i.e. sigma_j^2
     */
    for(j=0; j<p; j++){
      v[j] = (b[j]*b[j] + 2*tau)/rchisq(1+2*delta);
    }

    /**
    * step 6 Update tau, if needed
    */
    if(update_tau){
      rate = 0.0;
      for(j=0; j<p; j++){
      	rate += 1/v[j];
      }
      tau = rgamma(p*delta, 1/rate);
      if(tau == 0){
      	error("tau=0, rate=%e\n", rate);
      }
    }

    /**
    * step 7 Update delta, if needed
    */
    if(update_delta){
    	delta_data.tau = tau;
    	tmp = 0.0;
    	for(j=0; j<p; j++){
    		tmp += log(v[j]);
    	}
    	delta_data.sum_log_v = tmp;
    	xprev = delta;
      err = arms(xinit,ninit,&xl,&xr,ld_delta,&delta_data,&convex,
        npoint,dometrop,&xprev,&delta,nsamp,qcent,xcent,ncent,&neval);
      if(err>0){
        error("error in arms, error code = %d\n",err);
      }
      for(j=0; j<ninit; j++){
      	xinit[j] = xcent[j];
      }
    }

    /**
     * step 8 Record b when neccesary
     */
    if(w > Nburn){
      if((w - Nburn) % Nthin == 0){
        k = k + 1;
        for(j=0; j<p; j++){
          bSample[j][k] = b[j];
        }
        deltaSam[k] = delta;
        tauSam[k] = tau;
      }
    }

  }

  PutRNGstate();

  free(b);
  free(v);
  free(res);
  free(Xj2);
  free(Jset);
}

/*********************************************************************
 *
 * bayes_Lasso
 *
 * The Bayesian version of Lasso
 *
 *********************************************************************/

void bayes_Lasso(double* y, double* RX, double* bInit, double* RbSample, 
 double* kappaSam, int* dims, double* Rs, double* Rr, int* b_update_order)
{
  int i, j, jk, k, w, n, p, Nburn, Nthin, Nsample, Nall, init_b, *Jset, method;
  double xij, b0, b0_bar, bj_bar, v0, s0, sj, s, r, rate, kappa2;
  double *b, *v, *res, *Xj2, **X, **bSample, tmp;

  n = dims[0];
  p = dims[1];
  Nburn = dims[2];
  Nthin = dims[3];
  Nsample = dims[4];
  method  = dims[5];
  init_b  = dims[6];
  s = *Rs;
  r = *Rr;

  // Rprintf("s = %f, r = %f, method=%d\n", s, r, method);

  /* allocate memory */
  b    = (double *)calloc(p, sizeof(double));
  v    = (double *)calloc(p, sizeof(double));
  res  = (double *)calloc(n, sizeof(double));
  Xj2  = (double *)calloc(p, sizeof(double));
  Jset = (int *)calloc(p, sizeof(int));

  /* reoganize vector into matrix */
  reorg(RX, &X, n, p);
  reorg(RbSample, &bSample, p, Nsample);

  /**
   * Calculate Xj2
   * sum square for each marker, i.e., each column of X
   */
  for(j=0; j<p; j++){
    Xj2[j] = 0.0;
    for(i=0; i<n; i++){
      xij = X[i][j];
      Xj2[j] += xij*xij;
    }
  }

  /**
   * step 1. Initialization
   */
  b0 = 0.0;
  v0 = 0.01;
  
  if(init_b){
    for(j=0; j<p; j++){
      b[j] = bInit[j];
      v[j] = 0.01;
    }    
  }else{
    for(j=0; j<p; j++){
      b[j] = 0.0;
      v[j] = 0.01;
    }
  }

  kappa2 = 3.0;
  //Rprintf("kappa2 = %e\n", kappa2);

  Nall = Nburn + Nthin*Nsample;
  k = -1;

  for(j=0; j<p; j++){
    Jset[j] = j;
  }

  GetRNGstate();
  for(w=1; w<=Nall; w++){
  	//Rprintf("w=%d\n", w);

    /**
     * step 2. Update b0
     */
    for(i=0; i<n; i++){
      res[i] = y[i];
      for(j=0; j<p; j++){
        res[i] -= X[i][j]*b[j];
      }
    }
    b0_bar = mean(res, n);
    s0 = v0/n;
    b0 = rnorm(b0_bar, sqrt(s0));

    for(i=0; i<n; i++){
      res[i] -= b0;
    }

    /**
     * step 3.1 choose the order of updating b[j]
     */

    if(*b_update_order > 1){
      get_b_order(*b_update_order, Jset, w, n, p, y, X, b);
    }

    /**
     * step 3.2 Update b[j]
     */

    for(jk=0; jk<p; jk++){
      j = Jset[jk];

      /* remove the effect of Xj from the residual*/
      for(i=0; i<n; i++){
        res[i] = res[i] + X[i][j]*b[j];
      }

      bj_bar = 0.0;
      for(i=0; i<n; i++){
        bj_bar += X[i][j]*res[i];
      }

      if(method==1){
      	tmp = Xj2[j] + v0/v[j];
        bj_bar /= tmp;
        sj  = v0/tmp;
      }else{
      	tmp = Xj2[j] + 1/v[j];
        bj_bar /= tmp;
        sj  = v0/tmp;
		  }

      b[j] = rnorm(bj_bar, sqrt(sj));

      /* add the effect of Xj back into the residual */
      for(i=0; i<n; i++){
        res[i] = res[i] - X[i][j]*b[j];
      }
    }

    /**
     * step 4 Update v0, i.e. sigma_0^2
     */

    v0 = 0.0;
    for(i=0; i<n; i++){
      v0 += res[i]*res[i];
    }

    if(method == 1){
      v0 /= rchisq(n);
    }else if(method==2){
    	for(j=0; j<p; j++){
    		v0 += b[j]*b[j]/v[j];
    	}
    	v0 /= rchisq(n+p);
		}else{
			error("invalid value of method\n");
		}

    /**
     * step 5 Update v[j], i.e. simga_j^2
     */
    if(method == 1){
      for(j=0; j<p; j++){
        tmp = rinvGauss1(sqrt(kappa2)/fabs(b[j]), kappa2);
        if(tmp < 1e-10){
          warning("random number from inverse Gaussion < 1e-10, replaced it by 1e-10\n");
          tmp = 1e-10;
        }
        v[j] = 1/tmp;
      }
    }else if(method == 2){
      for(j=0; j<p; j++){
        tmp = rinvGauss1(sqrt(kappa2*v0)/fabs(b[j]), kappa2);
        if(tmp < 1e-10){
          warning("random number from inverse Gaussion < 1e-10, replaced it by 1e-10\n");
          tmp = 1e-10;
        }
        v[j] = 1/tmp;
      }
    }

    /**
     * step 6 Update kappa
     */
    rate = 0.0;
    for(j=0; j<p; j++){
    	rate += v[j];
    }
    rate += r;

    kappa2 = 2*rgamma(p+s, 1/rate);

    // Rprintf("w=%d, rate=%e, kappa2=%e\n", w, rate, kappa2);

    /**
     * step 7 Record b when neccesary
     */
    if(w > Nburn){
      if((w - Nburn) % Nthin == 0){
        k = k + 1;
        for(j=0; j<p; j++){
          bSample[j][k] = b[j];
        }
        kappaSam[k] = sqrt(kappa2);
      }
    }

  }

  PutRNGstate();

  free(b);
  free(v);
  free(res);
  free(Xj2);
  free(Jset);
}

