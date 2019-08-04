
void get_b_order(int b_update_order, int* Jset, int w, int n, int p, 
  double* y, double** X, double* b);

void IAL(double* y, double** X, double* bInit, double* RbSample, double* Rb0,
  double* Rdelta, double* Rtau, int* dims, double* Repsilon, int* nIter, 
  int* b_update_order, double* Rscore, int* RscoNA, double* RBIC_gamma, 
  int** IDC, char** Rcriterion, double* score2use, double* delta2use, 
  double* tau2use, double* Rratio);
      
void bayes_AL(double* y, double* RX, double* bInit, double* RbSample, 
  double* deltaSam, double* tauSam, int* dims, 
  double* Rdelta, double* Rtau, int* b_update_order);
  
void bayes_t(double* y, double* RX, double* bInit, double* RbSample,
  double* deltaSam, double* tauSam, int* dims,
  double* Rdelta, double* Rtau, int* b_update_order);
  
void bayes_Lasso(double* y, double* RX, double* bInit, double* RbSample, 
  double* kappaSam, int* dims, double* Rs, double* Rr, int* b_update_order);
