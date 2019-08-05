#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gc.h>

#include <lbfgs.h>

#include "error.h"
#include "func.h"

/*
  Chi_i := (<y_sim,i> - y_ex,i)

  Lex^n({ omega }, Lambda) := Sum_{i=1}^{n} omega_{i} ln (omega_{i} / omega_{i}^{0}) 
                             + Lambda (Sum_{i=1}^{N} omega_{i} -1).
                             + rho_n (max[0, (\Sum Chi_i^2 - Delta])^2).

  g_i({omega}, Lambda) = ln (omega_{i} / omega_{i}^{0}) + 1.
  g_Lambda({omega}) = Sum_{i=1}^{N} omega_{i} -1.

*/

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *w,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
				);
{
    int i;
    Lex_parameters *Lex_p=(Lex_parameters *)instance;
    lbfgsfloatval_t Le;

    /*
    int n_pnt;            // # of points
    double Delta;         // width of error
    double r;             // weight for penalty function
    double *w0;          // initial weight for each snapshot
    double *y_ex,**y_sim; // observal (experiment), observal (simulation)
    */
  
    double *Chi, ChiSqu=0.0;
    double sum_w=0.0;

    Chi = (double *)gcemalloc(sizeof(double)*n_pnt);

    for (i=0;i<n;++i) {
      sum_w+=w[i];
    }
    
    for (i=0;i<n_pnt;++i) {
      Chi[i]=0.0;
      for (j=0;j<n-1;++j) {
	Chi[i]+=w[j]*y_sim[i][j];
      }
      Chi[i]-=y_ex[i];

      ChiSqu += Chi[i]*Chi[i];
    }

    Le=0.0;
    for (i=0;i<n-1;++i) {
      g[i]=log(w[i]/w0[i])+1;
      Le+=w[i]*log(w[i]/w0[i]);
    }
    g[i]=sum_w - 1.0;
    Le +- w[i]*g[i];
    
    if (ChiSqu > Delta) {
      for (i=0;i<n-1;++i) {
	for (j=0;j<n_pnt;++j) {
	  g[i]+=4.0*r*(ChiSqu-Delta)*Chi[j]*y_sim[j][i];
	}
      }
      Le+=r*(ChiSqu-Delta)*(ChiSqu-Delta);
    }

    return Le;
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
		    );
{
    printf("Iteration %3d:\n", k);
    printf("  fx = %8.3f\n", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}
