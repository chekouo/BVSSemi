// R  CMD SHLIB  LogLikNonLocal.c -lgsl -lgslcblas -o NonLocal.so
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>
#include "header.h"

int facto(int n)
{
    int f = 1;
    int i;
    for (i = 1; i <= n; i++)
    {
        f *= i;
    }
    return f;
}
double *Ain(int k, int K, int ng, int N, double h1, double h0, double hg, double **Ps)
{
    // K=number of prognosic factors(covariates)
    double *Ainv = malloc(k * k * sizeof(double));
    int i, j, l;
    // Ainv+(2*r*nu*s2/(nu-2)).*diag(1./(betaMax.^2))
    // printf("N===%d\n",N);
    for (i = 0; i < k; i++)
    {
        for (j = 0; j <= i; j++)
        {
            double a = 0;
            if ((i == 0) && (j == 0))
            {
                a = N + (1 / h0);
            }
            else if (j == 0)
            {
                for (l = 0; l < N; l++)
                    a += Ps[l][i - 1];
            }
            else if ((i != 0) && (j != 0))
            {
                for (l = 0; l < N; l++)
                    a += Ps[l][i - 1] * Ps[l][j - 1];
                if (i == j)
                {
                    if (i < K)
                        a += (1.0 / h1);
                    else
                        a += (1.0 / hg);
                }
            }
            Ainv[i * k + j] = Ainv[j * k + i] = a;
            // printf("%f %d %d \n",Ainv[i*k+j],i,j);
        }
    }
    return Ainv;
}
void LogLikNonLocal(int *k1, int *K1, int *ng1, int *N1, double *alpha1, double *psi1, double *y, double *Ps1, 
double *betaMax, int *r1, double *h0v, double *hg1, int *maxiter1, double *stop1, int *Posbeta1, double *LogLik1)
{
    // Posbeta=1 if the regression coefficients are positifs
    // Posbeta=0 if the regression coefficients are positifs or negatifs
    // K=number of covariates (prognostic factors)
    // k=nbr total of selected covariates
    // ng=nbr of selected genes
    // Ps matrix of significant covariates
    // K1 is the number of features that we are not selecting
    // ng is the number of feature we want to select
    // k=1+K1+ng
    // alpha and psi are hyperparameters for sigma2
    // y is the dependent var
    // Ps1 is the matrix of selected covariates including the constant
    // betaMax is the estimated beta
    // r1 is the r paramater in the pMOM
    // h1v ann hg1 should be equal, hov should be large (hov), it's the tau for the constant as the constant should be noninformative
    // stop is the error term to stop with the algorithm and can be set as stop1=10^-3
    // maxiter1 is the number of iterations to maximize beta
    //

    int k = k1[0];
    int K = K1[0];
    int ng = ng1[0];
    int N = N1[0];
    double alpha = alpha1[0];
    double psi = psi1[0];
    _Bool Posbeta = Posbeta1[0];
    int r = r1[0];
    double h0 = h0v[0];
    double h1 = h0;
    double hg = hg1[0];
    int maxiter = maxiter1[0];
    double stop = stop1[0];
    double **Ps = dmatrix(0, N - 1, 0, K + ng - 1);
    int ll = 0;
    int i, j, l;
    double a;

    for (l = 0; l < N; l++)
    {
        for (i = 0; i < K + ng; i++)
        {
            Ps[l][i] = Ps1[ll];
            ll = ll + 1;
        }
    }

    double *Ainv = Ain(k, K, ng, N, h1, h0, hg, Ps);
    double *Ainv1 = malloc(k * k * sizeof(double));
    for (i = 0; i < k; i++)
        for (j = 0; j <= i; j++)
            Ainv1[i * k + j] = Ainv1[j * k + i] = Ainv[i * k + j];
    gsl_matrix_view m = gsl_matrix_view_array(Ainv1, k, k);
    gsl_linalg_cholesky_decomp(&m.matrix);

    double *Ainvp = malloc(k * k * sizeof(double));
    double *Xty = malloc(k * sizeof(double));
    for (i = 0; i < k; i++)
    {
        a = 0;
        for (l = 0; l < N; l++)
        {
            if (Posbeta == 0)
            {
                if (i == 0)
                    a += y[l];
                else
                    a += Ps[l][i - 1] * y[l];
            }
            else
            { // Posbeta==1 and without constant
                a += Ps[l][i] * y[l];
            }
        }
        Xty[i] = a;
    }

    gsl_vector_view b = gsl_vector_view_array(Xty, k);

    gsl_vector *x = gsl_vector_alloc(k);

    gsl_linalg_cholesky_solve(&m.matrix, &b.vector, x);

    double nu = N + 2 * r * k + 2 * alpha;
    double s2 = 0;
    double betaHat[k];

    for (i = 0; i < k; i++)
    {
        betaHat[i] = gsl_vector_get(x, i);
        s2 += Xty[i] * betaHat[i];
    }
    gsl_vector_free(x);
    double yy = 0;
    for (i = 0; i < N; i++)
    {
        yy += pow(y[i], 2);
    }
    s2 = (2 * psi + yy - s2) / nu;
    if (Posbeta == 0)
    {
        for (i = 0; i < k; i++)
        {
            betaMax[i] = betaHat[i];
        }
    }
    else
    {
        for (i = 0; i < k; i++)
        {
            betaMax[i] = betaHat[i] * (betaHat[i] > 0);
        }
    }
    MaximizeBeta2(Xty, nu, s2, Ainv, maxiter, stop, betaMax, k, r, Posbeta);
    for (i = 0; i < k; i++)
    {
        for (j = 0; j <= i; j++)
        {
            Ainvp[i * k + j] = Ainvp[j * k + i] = Ainv[i * k + j];
            if (i == j)
                Ainvp[i * k + i] += (2 * r * nu * s2 / (nu - 2)) * (1 / pow(betaMax[i], 2));
        }
    }
    double L1 = gsl_sf_lngamma(nu / 2) + (alpha * log(psi)) + (nu / 2) * log(2);
    if (Posbeta == 1)
        L1 += k * log(2);
    double betaAibeta = 0;
    double sumlogbeta = 0;
    double difbetaAibeta = 0;
    for (i = 0; i < k; i++)
    {
        sumlogbeta += log(pow(betaMax[i], 2));
        for (j = 0; j < i; j++)
        {
            betaAibeta += 2 * Ainv[i * k + j] * betaHat[i] * betaHat[j];
            difbetaAibeta += 2 * Ainv[i * k + j] * (betaHat[i] - betaMax[i]) * (betaHat[j] - betaMax[j]);
        }
        betaAibeta += pow(betaHat[i], 2) * Ainv[i * k + i];
        difbetaAibeta += pow(betaHat[i] - betaMax[i], 2) * Ainv[i * k + i];
    }

    double L2 = -(nu / 2) * log(2 * psi + yy - betaAibeta) + r * sumlogbeta;
    double L3 = -((nu - 2) / (2 * nu * s2)) * difbetaAibeta;

    gsl_matrix_view Aip = gsl_matrix_view_array(Ainvp, k, k);
    gsl_linalg_cholesky_decomp(&Aip.matrix);

    double logdet = 0;
    for (i = 0; i < k; i++)
    {
        logdet += 2 * log(gsl_matrix_get(&Aip.matrix, i, i));
    }
    double L4 = -0.5 * logdet;
    // double L4=-0.5*gsl_linalg_CD_lndet(&Aip.matrix);

    double doublefact = facto(2 * r - 1) / ((1 << (r - 1)) * facto(r - 1));
    double L5 = -gsl_sf_lngamma(alpha) - k * log(doublefact) - (N / 2.0) * log(2 * PI)  - (K / 2.0 + r * K) * log(h1) - (ng / 2.0 + r * ng) * log(hg) - (0.5 + r) * log(h0);
    double LogLik = L1 + L2 + L3 + L4 + L5;
    free(Xty);
    free(Ainvp);
    free_dmatrix(Ps, 0, N - 1, 0, k - 1);
    free(Ainv1);
    // printf("L1 L2 L3 L4 L5 =%f %f %f %f %f\n",L1,L2,L3,L4,L5);
    LogLik1[0] = LogLik;
}
void MaximizeBeta2(double *Xty, double nu, double s2, double *Ainv, int maxiter, double stop, double *betaInit, int k, int r, _Bool Posbeta)
{
    int i, m, m1;
    double a = (nu * s2) / (nu - 2);
    int converged = 0;
    double *beta = malloc(k * sizeof(double));
    for (i = 0; i < k; i++)
    {
        beta[i] = betaInit[i];
        // printf("am==%f \n",Ainv[i*k+i]);
    }
    i = 0;
    while ((i < maxiter) && (converged == 0))
    {
        for (m = 0; m < k; m++)
        {
            double am = Ainv[m * k + m];
            // printf("am==%f \n",am);
            double AinvBetaMinus_m = 0;
            if (k > 0)
            {
                for (m1 = 0; m1 < k; m1++)
                    if (m1 != m)
                        AinvBetaMinus_m += Ainv[m * k + m1] * beta[m1];
                // AinvBetaMinus_m=Ainv(m,:)*beta-am*beta(m);
            }
            double delta = pow(AinvBetaMinus_m - Xty[m], 2) + (8 * r * a * am);
            double f1 = (-AinvBetaMinus_m + Xty[m] + sqrt(delta)) / (2 * am);
            // diffValues is the difference between the values of the objective function for the two roots f1 and f2;
            beta[m] = f1;
            if (Posbeta == 0)
            {
                double f2 = -2 * r * a / (am * f1);
                double diffValues = betaInit[m] * ((f1 - f2) / (f1 * f2)) - log(pow(f1, 2)) + log(pow(f2, 2));
                if (diffValues > 0)
                    beta[m] = f2;
            }
        } // end of m
        double diffBeta[k];
        for (m1 = 0; m1 < k; m1++)
        {
            diffBeta[m1] = betaInit[m1] - beta[m1];
            // printf("diiVal==%f \n",beta[m1]);
        }

        if ((norm(k, diffBeta) / (norm(k, beta) + norm(k, betaInit))) < stop)
            converged = 1;
        for (m1 = 0; m1 < k; m1++)
        {
            betaInit[m1] = beta[m1];
            // printf("Beta=%f ",beta[m1]);
        }
        i++;
    }
    free(beta);
    // return beta;
}
double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;

    m = (double **)malloc((unsigned)(nrh - nrl + 1) * sizeof(double *));
    if (!m)
        nrerror("allocation failure 1 in dmatrix()");
    m -= nrl;

    for (i = nrl; i <= nrh; i++)
    {
        m[i] = (double *)malloc((unsigned)(nch - ncl + 1) * sizeof(double));
        if (!m[i])
            nrerror("allocation failure 2 in dmatrix()");
        m[i] -= ncl;
    }
    return m;
}
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
    int i;

    for (i = nrh; i >= nrl; i--)
        free((char *)(m[i] + ncl));

    free((char *)(m + nrl));
}
double norm(int n, double *x)
{
    double normx = 0;
    int i;
    for (i = 0; i < n; i++)
    {
        normx += pow(x[i], 2);
    }
    return sqrt(normx);
}

void nrerror(char error_text[])
{
    printf("Utils run-time error...\n");
    printf("%s\n", error_text);
    printf("...now exiting to system...\n");
    exit(1);
}
