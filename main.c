/* 
 Ensemble refinment with experimental data.

 INPUT;
 { Y } : experimental data. M dimension.
 { y } : simulation data projected to the observal space. M x N dimension.

 OUTPUT;
 omega: weights for each snapshots.

 Where,
     rho_sim : distribution derived from simulation.
     rho_sim : optimal distribution.

     Chi := \int \rho y dx - Y.

 { omega^{opt.} } <- arg min_{ omega } D(rho || rho_sim),
 subject to Chi^2 < Delta, and Sigma (omega) = 1.
*/

#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <lbfgs.h>

#include "error.h"
//#include "func.h"

#define M 5

struct Lex_parameters {
  int n_pnt;            // # of points
  double Delta;         // width of error
  double r;             // weight for penalty function
  double *wi0;          // initial weight for each snapshot
  double *y_ex,**y_sim; // observal (experiment), observal (simulation)
};

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step);

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
    int ls);

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,ret = 0;
  double dummy,f;

  int n_snp, n_pnt=1, N;   // # of snapshots, # of data points

  double rk[4];            // weight for penalty function
  
  double Del=10e-2;        // width of error
  double r;                // weight for penalty function
  double *wi0;             // initial weight for each snapshot
  double *y_ex,**y_sim;    // observal (experiment), observal (simulation)

  struct Lex_parameters Lex_p;    // parametrs for extended Lagrangian
  
  lbfgsfloatval_t Lex;     // extended Lagrangian
  lbfgsfloatval_t *wopt;   // (optimal) weights for snapshots
  lbfgs_parameter_t param; // parametrs for LBFGS

  char *inputfileExpDataname;
  char *inputfileSimDataname;
  char *outputfilename;

  FILE *inputfileExpData;
  FILE *inputfileSimData;
  FILE *outputfile;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;

  rk[0]=1.0e-1; rk[1]=1.0e0; rk[2]=1.0e1; rk[3]=1.0e2;
    
  struct option long_opt[] = {
    {"Del",1,NULL,'d'},
    {"r1",1,NULL,'1'},
    {"r2",1,NULL,'2'},
    {"r3",1,NULL,'3'},
    {"r4",1,NULL,'4'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hd:1:2:3:4:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'd':
      Del=atof(optarg);  break;
    case '1':
      rk[0]=atof(optarg);  break;
    case '2':
      rk[1]=atof(optarg);  break;
    case '3':
      rk[2]=atof(optarg);  break;
    case '4':
      rk[3]=atof(optarg);  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfileExpDataname = *argv;
  inputfileSimDataname = *++argv;
  outputfilename       = *++argv;

  printf("inputfile of Experimental Data=%s\n",inputfileExpDataname);
  printf("inputfile of Simulation   Data=%s\n",inputfileSimDataname);
  printf("outputfile      =%s\n",outputfilename);
  
  y_ex=(double *)gcemalloc(sizeof(double)*1);

  inputfileExpData=efopen(inputfileExpDataname,"r");
  n_pnt=0;
  while ( (d=fscanf(inputfileExpData,"%lf",&f)) != EOF )  {
    y_ex[n_pnt]=f;
    ++n_pnt;
    y_ex=(double *)gcerealloc(y_ex,sizeof(double)*(n_pnt));
  }
  fclose(inputfileExpData);

  printf("# of data points : %4d\n",n_pnt);

  /* debug printf */
  for (i=0;i<n_pnt;++i) {
    printf("%8.3lf\n",y_ex[i]);
  }
  printf("\n");

  n_snp = 1;
  y_sim   =(double **)gcemalloc(sizeof(double *)*n_snp);
  y_sim[0]=(double  *)gcemalloc(sizeof(double)*n_pnt);
  
  inputfileSimData=efopen(inputfileSimDataname,"r");
  i=0;
  d = 1;
  while ( (d=fscanf(inputfileSimData,"%lf",&f)) != EOF )  {
    if (i<n_pnt) {
      y_sim[n_snp-1][i] = f;

      ++i;
    }
    else {
      ++n_snp;
      i=0;

      y_sim=(double **)gcerealloc(y_sim,sizeof(double *)*(n_snp));
      y_sim[n_snp-1]=(double *)gcemalloc(sizeof(double)*n_pnt);
      y_sim[n_snp-1][i] = f;

      ++i;
    }
  }
  fclose(inputfileSimData);

  printf("# of snap shots  : %4d\n",n_snp);
  
  /* debug printf */
  for (i=0;i<n_snp;++i) {
    for (j=0;j<n_pnt;++j) {
      printf("%8.3lf\n",y_sim[i][j]);
    }
    printf("\n");
  }

  wi0 = (double *)gcemalloc(sizeof(double)*n_snp);

  for (i=0;i<n_snp;++i) {
    wi0[i] = 1.0/n_snp;
  }

  wopt = lbfgs_malloc(n_snp+1);

  for (i=0;i<n_snp;++i) {
    wopt[i] = wi0[i];
  }
  wopt[i] = 0.0;

  /* Initialize the parameters for the L-BFGS optimization. */
  lbfgs_parameter_init(&param);

  Lex_p.Delta=Del;
  Lex_p.wi0=wi0;
  Lex_p.y_ex=y_ex;
  Lex_p.y_sim=y_sim;

  N = n_snp + 1;

  for (i=0;i<M;++i) {
    
    Lex_p.r=rk[i];
  
    /*Start the L-BFGS optimization.*/
    ret = lbfgs(N, wopt, &Lex, evaluate, progress, &(Lex_p), &param);
    printf("n=%3d  Ln = %8.3f\n", i+1,Lex);
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<n_snp;++i) {
    fprintf(outputfile,"%4d %10.8lf\n",i+1,wopt[i]);
  }
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("[-Del <Delta Value> ] \n");
  printf("[-r1  <rho1> ] \n");
  printf("[-r2  <rho2> ] \n");
  printf("[-r3  <rho3> ] \n");
  printf("[-r4  <rhor> ] \n");
  printf("%s inputname1(exp data) inputname2(simu data) outfilename\n",progname);
}

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
				) {
  int i, j;
  struct Lex_parameters *Lex_p=(struct Lex_parameters *)instance;
  lbfgsfloatval_t Le;

  double *Chi, ChiSqu=0.0;
  double sum_w=0.0;
  
  Chi = (double *)gcemalloc(sizeof(double)*Lex_p->n_pnt);

  for (i=0;i<n;++i) {
    sum_w+=w[i];
  }
    
  for (i=0;i<Lex_p->n_pnt;++i) {
    Chi[i]=0.0;
    for (j=0;j<n-1;++j) {
      Chi[i]+=w[j]*Lex_p->y_sim[i][j];
    }
    Chi[i]-=Lex_p->y_ex[i];
    
    ChiSqu += Chi[i]*Chi[i];
  }

  Le=0.0;
  for (i=0;i<n-1;++i) {
    g[i]=log(w[i]/Lex_p->wi0[i])+1;
    Le+=w[i]*log(w[i]/Lex_p->wi0[i]);
  }
  g[i]=sum_w - 1.0;
  Le += w[i]*g[i];
  
  if (ChiSqu > Lex_p->Delta) {
    for (i=0;i<n-1;++i) {
      for (j=0;j<Lex_p->n_pnt;++j) {
	g[i]+=4.0*Lex_p->r*(ChiSqu-Lex_p->Delta)*Chi[j]*Lex_p->y_sim[j][i];
      }
    }
    Le+=Lex_p->r*(ChiSqu-Lex_p->Delta)*(ChiSqu-Lex_p->Delta);
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
		    )
{
  printf("Iteration %3d:\n", k);
  printf("  fx = %8.3f\n", fx);
  printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
  printf("\n");
  return 0;
}
