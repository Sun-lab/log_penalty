
double fb(double b_norm_j, int h, double kappa_j, double *sigma2, 
		  double Xj2_j, double *bj_bar, double lambda);

double zeroin(double ax, double bx, double (*f)(), double tol, int H, 
 double kappa_j, double *sigma2, double Xj2_j, double *bj_bar, double lambda);

void mGAPr(double* Ry, double* RX, double* RB,
  double* Rlambda, double* Rtau, int* dims, double* Repsilon, int* nIter, 
  int* b_update_order, double* Rscore, int* RscoNA, double* score2use, 
  double* lambda2use, double* tau2use,
  double* RB2, double* Rscore2, int* RscoNA2, double* score2use2, 
  double* lambda2use2, double* tau2use2,
  double* RB3, double* Rscore3, int* RscoNA3, double* score2use3, 
  double* lambda2use3, double* tau2use3, double* extbicgamma);
  
void mGAPc(double* Ry, char** fname, double* RB, 
  double* Rlambda, double* Rtau, int* dims, double* Repsilon, int* nIter, 
  int* b_update_order, double* Rscore, int* RscoNA, double* score2use, 
  double* lambda2use, double* tau2use, 
  double* RB2, double* Rscore2, int* RscoNA2, double* score2use2, 
  double* lambda2use2, double* tau2use2,
  double* RB3, double* Rscore3, int* RscoNA3, double* score2use3, 
  double* lambda2use3, double* tau2use3, double* extbicgamma);

void mGAP(double** y, double** X, double** RB, double* Rlambda, 
          double* Rtau, int* dims, double* Repsilon, int* nIter, 
          int* b_update_order, double* Rscore, int* RscoNA, 
          double* score2use, double* lambda2use, double* tau2use,
          double** RB2, double* Rscore2, int* RscoNA2, 
          double* score2use2, double* lambda2use2, double* tau2use2,
          double** RB3, double* Rscore3, int* RscoNA3, 
          double* score2use3, double* lambda2use3, double* tau2use3, 
          double* extbicgamma);
