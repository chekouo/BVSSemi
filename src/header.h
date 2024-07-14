#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define PI 3.1415926535897932384
#include <stdbool.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

void nrerror(char error_text[]);
double norm(int n, double * x);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int facto(int n);
void MaximizeBeta2(double *Xty,double nu,double s2,double *Ainv,int maxiter,double stop,double *betaInit,int k,int r,_Bool Posbeta);
double * Ain(int k,int K,int ng,int N,double h1, double h0,double hg,double **Ps);
void LogLikNonLocal(int *k1, int *K1, int *ng1, int *N1, double *alpha1, double *psi1, double *y, double *Ps1, 
double *betaMax, int *r1, double *h0v, double *hg1, int *maxiter1, double *stop1, int *Posbeta1, double *LogLik1);
